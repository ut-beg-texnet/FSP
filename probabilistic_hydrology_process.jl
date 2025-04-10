using CSV
using DataFrames
using PrettyTables
using Distributions
using Statistics
using Dates
using LinearAlgebra




include("core/hydrology_calculations.jl")
include("core/utilities.jl")
include("TexNetWebToolLauncherHelperJulia.jl")
include("core/bill_pfront.jl")
include("graphs/julia_fsp_graphs.jl")

using .HydroCalculations
using .Utilities
using .TexNetWebToolLauncherHelperJulia
using .BillPFront
using .JuliaFSPGraphs


const ARGS_FILE_NAME = "args.json"
const RESULTS_FILE_NAME = "results.json"


struct HydrologyParams
    aquifer_thickness::Float64
    porosity::Float64
    permeability::Float64
    fluid_density::Float64
    dynamic_viscosity::Float64
    fluid_compressibility::Float64
    rock_compressibility::Float64
    plus_minus::Dict{String, Float64}
    n_iterations::Int64
end



function run_monte_carlo_hydrology(helper::TexNetWebToolLaunchHelperJulia, 
                                  params::HydrologyParams, 
                                  distribution_type::String="uniform",
                                  extrapolate_injection_rates::Bool=false)

    if distribution_type == "uniform"
        # check if we are missing the plus_minus value for a parameter, set it to 0.0
        for (key, value) in params.plus_minus
            if isnan(value)
                params.plus_minus[key] = 0.0
            end
        end

        # if for some reason we are missing the number of iterations, set it to 750
        # that shouldn't happen since the portal sets a default of 750 if it's missing
        if params.n_iterations == 0 || isnothing(params.n_iterations)
            params.n_iterations = 750
        end

        # check that the plus_minus values are not greater than the base values
        for (key, value) in params.plus_minus
            if key == "aquifer_thickness" && value > params.aquifer_thickness ||
                key == "porosity" && value > params.porosity ||
                key == "permeability" && value > params.permeability ||
                key == "fluid_density" && value > params.fluid_density ||
                key == "dynamic_viscosity" && value > params.dynamic_viscosity ||
                key == "fluid_compressibility" && value > params.fluid_compressibility ||
                key == "rock_compressibility" && value > params.rock_compressibility
                throw(ArgumentError("plus_minus value for $key cannot be greater than the base value."))
            end
        end

        # create the Uniform distributions for each parameter
        distributions = Dict(
            "aquifer_thickness" => create_uniform_distribution(params.aquifer_thickness, params.plus_minus["aquifer_thickness"]),
            "porosity" => create_uniform_distribution(params.porosity, params.plus_minus["porosity"]),
            "permeability" => create_uniform_distribution(params.permeability, params.plus_minus["permeability"]),
            "fluid_density" => create_uniform_distribution(params.fluid_density, params.plus_minus["fluid_density"]),
            "dynamic_viscosity" => create_uniform_distribution(params.dynamic_viscosity, params.plus_minus["dynamic_viscosity"]),
            "fluid_compressibility" => create_uniform_distribution(params.fluid_compressibility, params.plus_minus["fluid_compressibility"]),
            "rock_compressibility" => create_uniform_distribution(params.rock_compressibility, params.plus_minus["rock_compressibility"])
        )
    elseif distribution_type == "gaussian" # ADITTIONAL FEATURE, NOT USED YET AND WE MIGHT NEED TO ADD A DIFFERENT DISTRIBUTION TYPE
        # create the Gaussian distributions for each parameter
        distributions = Dict(
            "aquifer_thickness" => create_gaussian_distribution(params.aquifer_thickness, params.plus_minus["aquifer_thickness"]),
            "porosity" => create_gaussian_distribution(params.porosity, params.plus_minus["porosity"]),
            "permeability" => create_gaussian_distribution(params.permeability, params.plus_minus["permeability"]),
            "fluid_density" => create_gaussian_distribution(params.fluid_density, params.plus_minus["fluid_density"]),
            "dynamic_viscosity" => create_gaussian_distribution(params.dynamic_viscosity, params.plus_minus["dynamic_viscosity"]),
            "fluid_compressibility" => create_gaussian_distribution(params.fluid_compressibility, params.plus_minus["fluid_compressibility"]),
            "rock_compressibility" => create_gaussian_distribution(params.rock_compressibility, params.plus_minus["rock_compressibility"])
        )
    else
        throw(ArgumentError("Invalid distribution type."))
    end

    # Get fault data
    fault_data_path = get_dataset_file_path(helper, 5, "faults_model_inputs_output")
    if fault_data_path === nothing
        error("Required fault dataset not found or accessible.")
    end
    fault_df = CSV.read(fault_data_path, DataFrame)
    num_faults = nrow(fault_df)
    println("Read faults_model_inputs_output.csv: $fault_df")

    # create the matrix to store the results of the Monte Carlo simulations
    ppOnFaultMC = zeros(params.n_iterations, num_faults) # rows are iterations, columns are faults

    # get the year of interest
    year_of_interest = get_parameter_value(helper, 5, "year_of_interest")
    #year_of_interest = 2022
    
    if year_of_interest === nothing
        year_of_interest = Dates.year(today())  # default to current year
        add_message_with_step_index!(helper, 5, "Year of interest was not provided, using the current year ($year_of_interest) as the default value", 0)
    else
        # Check if year_of_interest is already an Int before parsing
        if isa(year_of_interest, String)
            year_of_interest = parse(Int, year_of_interest)
        elseif !isa(year_of_interest, Int)
            # Handle unexpected types if necessary, e.g., throw an error or log a warning
            error("Unexpected type for year_of_interest: $(typeof(year_of_interest))")
        end
        # If it's already an Int, we don't need to do anything
        println("Year of interest: $year_of_interest")
    end

    # Get injection well data
    injection_wells_filepath, injection_data_type = get_injection_dataset_path(helper, 5)
    if injection_wells_filepath === nothing
        error("No injection well dataset found. Please provide injection well data.")
    end
    injection_wells_df = CSV.read(injection_wells_filepath, DataFrame)
    
    # Get unique well IDs based on data format - moved outside the loops for efficiency
    well_id_col = injection_data_type == "injection_tool_data" ? "API Number" : "APINumber"
    well_ids = unique(injection_wells_df[!, well_id_col])
    
    # Pre-process well data outside the Monte Carlo loop
    well_info = Dict{String, Dict{String, Any}}()
    for well_id in well_ids
        if injection_data_type == "injection_tool_data"
            well_data = injection_wells_df[injection_wells_df[!, well_id_col] .== well_id, :]
            if isempty(well_data)
                # Try UIC Number
                well_data = injection_wells_df[injection_wells_df[!, "UIC Number"] .== well_id, :]
                if isempty(well_data)
                    continue
                end
            end
            
            # Get well coordinates - use lat/lon directly
            well_lat = first(well_data[!, "Surface Latitude"])
            well_lon = first(well_data[!, "Surface Longitude"])
            
            # Determine injection period
            dates = Date[]
            try
                dates = Date.(well_data[!, "Date of Injection"], dateformat"y-m-d")
            catch
                try
                    dates = Date.(well_data[!, "Date of Injection"], dateformat"m/d/y")
                catch
                    try
                        dates = Date.(well_data[!, "Date of Injection"], dateformat"m/d/yyyy")
                    catch e
                        @warn "Could not parse dates for well $well_id: $e"
                        continue
                    end
                end
            end
            
            if isempty(dates)
                continue
            end
            
            inj_start_year = year(minimum(dates))
            inj_end_year = year(maximum(dates))
        else
            # Annual or monthly format
            well_data = injection_wells_df[injection_wells_df[!, well_id_col] .== well_id, :]
            if isempty(well_data)
                continue
            end
            
            # Get well coordinates directly
            well_lat = first(well_data[!, "Latitude(WGS84)"])
            well_lon = first(well_data[!, "Longitude(WGS84)"])
            
            # Get injection period
            if "StartYear" in names(well_data)
                inj_start_year = first(well_data[!, "StartYear"])
                inj_end_year = first(well_data[!, "EndYear"])
            elseif "Year" in names(well_data)
                inj_start_year = minimum(well_data[!, "Year"])
                inj_end_year = maximum(well_data[!, "Year"])
            else
                continue
            end
        end
        
        # Check if the well's injection period overlaps with year of interest
        if year_of_interest < inj_start_year
            continue
        end
        actual_end = min(inj_end_year, year_of_interest)
        if actual_end <= inj_start_year
            continue
        end
        
        # Store well information for use in Monte Carlo loop
        well_info[well_id] = Dict{String, Any}(
            "lat" => well_lat,
            "lon" => well_lon,
            "inj_start_year" => inj_start_year,
            "inj_end_year" => inj_end_year,
            "actual_end" => actual_end
        )
    end
    
    println("Processed $(length(well_info)) wells that are active in year $year_of_interest")

    # PRE-PROCESS WELL DATA WITH EXTRAPOLATION (for efficiency)
    println("Pre-processing well injection data with extrapolation = $(extrapolate_injection_rates)...")
    prepared_well_data = Dict{String, Dict{String, Any}}()
    for (well_id, info) in well_info
        # Prepare injection data with extrapolation once
        days, rates = prepare_well_data_for_pressure_scenario(
            injection_wells_df,
            well_id,
            info["inj_start_year"],
            info["actual_end"],
            injection_data_type,
            year_of_interest,
            extrapolate_injection_rates # Use extrapolation parameter
        )
        
        if !isempty(days) && !isempty(rates)
            # Store the prepared data
            prepared_well_data[well_id] = Dict{String, Any}(
                "days" => days,
                "rates" => rates,
                "lat" => info["lat"],
                "lon" => info["lon"]
            )
        end
    end
    

    println("Running Monte Carlo simulations for hydrology...")
    # run the Monte Carlo simulations
    for i in 1:params.n_iterations
        # sample the parameters from the distributions
        sampled_params = Dict(
            "aquifer_thickness" => rand(distributions["aquifer_thickness"]),
            "porosity" => rand(distributions["porosity"]),
            "permeability" => rand(distributions["permeability"]),
            "fluid_density" => rand(distributions["fluid_density"]),
            "dynamic_viscosity" => rand(distributions["dynamic_viscosity"]),
            "fluid_compressibility" => rand(distributions["fluid_compressibility"]),
            "rock_compressibility" => rand(distributions["rock_compressibility"])
        )

        # calculate storativity and transmissivity
        S, T, rho = calcST(
            sampled_params["aquifer_thickness"], 
            sampled_params["porosity"], 
            sampled_params["permeability"], 
            sampled_params["fluid_density"], 
            sampled_params["dynamic_viscosity"], 
            9.81, 
            sampled_params["fluid_compressibility"], 
            sampled_params["rock_compressibility"]
        )

        STRho = (S, T, rho) 

        
        for f in 1:num_faults
            # Get fault coordinates
            x_fault_km = fault_df[f, "Longitude(WGS84)"]
            y_fault_km = fault_df[f, "Latitude(WGS84)"]
            
            # initialize pp
            ppOnFault = 0.0
            
            # Loop over wells using PRE-PROCESSED DATA
            for (well_id, data) in prepared_well_data
                # Use pre-processed injection data
                pressure_contribution = pfieldcalc_all_rates(
                    x_fault_km,
                    y_fault_km,
                    STRho,
                    data["days"],
                    data["rates"],
                    data["lon"],
                    data["lat"],
                    Date(year_of_interest, 1, 1)
                )
                
                # Add to total pressure for this fault
                ppOnFault += pressure_contribution
            end
            
            # Store the result for this fault and iteration
            ppOnFaultMC[i, f] = ppOnFault
        end
    end

    # Convert results to a DataFrame
    result_rows = []
    for i in 1:params.n_iterations
        for f in 1:num_faults
            push!(result_rows, (
                IterationID = i,
                FaultID = string(f),
                Pressure = ppOnFaultMC[i, f]
            ))
        end
    end
    
    results_df = DataFrame(result_rows)
    
    # Calculate statistics
    stats_df = DataFrame(
        FaultID = String[],
        Mean = Float64[],
        StdDev = Float64[],
        Median = Float64[],
        Min = Float64[],
        Max = Float64[]
    )
    
    for f in 1:num_faults
        fault_data = results_df[results_df.FaultID .== string(f), :Pressure]
        push!(stats_df, (
            string(f),
            mean(fault_data),
            std(fault_data),
            median(fault_data),
            minimum(fault_data),
            maximum(fault_data)
        ))
    end
    
    return results_df, stats_df
end

"""
Get the injection dataset path based on available data types
"""
function get_injection_dataset_path(helper::TexNetWebToolLaunchHelperJulia, step_index::Int)
    for param_name in ["injection_wells_annual", "injection_wells_monthly", "injection_tool_data"]
        filepath = get_dataset_file_path(helper, step_index, param_name)
        if filepath !== nothing
            if param_name == "injection_wells_annual"
                injection_data_type = "annual_fsp"
                return filepath, injection_data_type
            elseif param_name == "injection_wells_monthly"
                injection_data_type = "monthly_fsp"
                return filepath, injection_data_type
            elseif param_name == "injection_tool_data"
                injection_data_type = "injection_tool_data"
                return filepath, injection_data_type
            end
        end
    end
    
    return nothing, nothing
end

"""
Calculate slip potential by combining probabilistic geomechanics CDF with deterministic hydrology


Using Deterministic Hydrology to Calculate Slip Potential:
Once we have calculated the pore pressure added to a given fault (from deterministic hydrology) and the 
probability of fault slip as a function of pore pressure increase (from probabilistic geomechanics), we can calculate the 
cumulative probability of fault slip. We do this simply by slicing the probability shown 
by the CDF curve in probabilistic geomechanics for the appropriate pore pressure increase

"""
function calculate_deterministic_slip_potential(prob_geo_cdf::DataFrame, det_hydro_pressures::DataFrame, year_of_interest::Union{Int, Nothing}=nothing)
    # results df
    slip_potential_df = DataFrame(
        FaultID = String[],
        Year = Int[],
        Date = Date[],
        slip_pressure = Float64[],
        probability = Float64[]
    )
    
    # Get unique fault IDs
    fault_ids = unique(prob_geo_cdf.FaultID)
    
    # Filter deterministic hydrology data by year if specified
    filtered_hydro_pressures = det_hydro_pressures
    if !isnothing(year_of_interest) && "Date" in names(det_hydro_pressures)
        filtered_hydro_pressures = det_hydro_pressures[year.(det_hydro_pressures.Date) .== year_of_interest, :]
        println("Filtered deterministic hydrology data to year $year_of_interest")
    end
    
    for fault_id in fault_ids
        # Get fault's geomechanics CDF
        fault_cdf = prob_geo_cdf[prob_geo_cdf.FaultID .== fault_id, :]
        
        # Get all deterministic pressure entries for this fault
        fault_pressures = filtered_hydro_pressures[filtered_hydro_pressures.FaultID .== fault_id, :]
        
        if isempty(fault_pressures)
            continue
        end

        if "Date" in names(det_hydro_pressures)
            if !all(x -> x isa Date, det_hydro_pressures.Date)
                error("det_hydro_pressures.Date contains non-Date values")
            end
        else
            error("det_hydro_pressures is missing Date column")
        end
        
        # Process each pressure entry (could be multiple years)
        for row in eachrow(fault_pressures)
            pressure = row.slip_pressure
            evaluation_date = row.Date
            evaluation_year = year(evaluation_date)
            
            # Sort the CDF by pore pressure
            # Make sure we use the correct column names from prob_geo_cdf
            slipPressureCol = "slip_pressure"
            cumulativeProbCol = "probability"
            sort!(fault_cdf, slipPressureCol)
            
            # Find slip potential by interpolating the CDF
            slip_potential = 0.0
            
            # If pressure is less than minimum in CDF, slip potential is 0
            if pressure < minimum(fault_cdf[!, slipPressureCol])
                slip_potential = 0.0
            # If pressure is greater than maximum in CDF, slip potential is 100%
            elseif pressure > maximum(fault_cdf[!, slipPressureCol])
                slip_potential = 100.0
            else
                # Interpolate to find slip potential
                # CDF uses normalized cumulative probability (0-100%)
                slip_potential = interpolate_cdf(fault_cdf[!, slipPressureCol], fault_cdf[!, cumulativeProbCol], pressure)
            end
            
            # Convert fault_id to String before adding to DataFrame
            push!(slip_potential_df, (string(fault_id), evaluation_year, evaluation_date, pressure, slip_potential))
        end
    end

    println("Slip potential df:")
    pretty_table(slip_potential_df)
    
    
    return slip_potential_df
end

"""
Calculate slip potential by combining probabilistic geomechanics CDF with probabilistic hydrology
"""
function calculate_probabilistic_slip_potential(prob_geo_cdf::DataFrame, prob_hydro_results::DataFrame)
    # Create result dataframe
    slip_potential_df = DataFrame(
        FaultID = String[],
        MeanPorePressure = Float64[],
        MeanSlipPotential = Float64[],
        StdDevSlipPotential = Float64[],
        MinSlipPotential = Float64[],
        MaxSlipPotential = Float64[]
    )
    
    # Get unique fault IDs
    fault_ids = unique(prob_geo_cdf.FaultID)
    
    for fault_id in fault_ids
        # Get fault's geomechanics CDF
        fault_cdf = prob_geo_cdf[prob_geo_cdf.FaultID .== fault_id, :]
        sort!(fault_cdf, :slip_pressure)
        
        # Get all probabilistic hydrology pressures for this fault
        fault_pressures = prob_hydro_results[prob_hydro_results.FaultID .== fault_id, :Pressure]
        
        if isempty(fault_pressures)
            continue
        end
        
        # Calculate slip potential for each pressure value
        slip_potentials = Float64[]
        
        for pressure in fault_pressures
            # If pressure is less than minimum in CDF, slip potential is 0
            if pressure < minimum(fault_cdf.slip_pressure)
                push!(slip_potentials, 0.0)
            # If pressure is greater than maximum in CDF, slip potential is 100%
            elseif pressure > maximum(fault_cdf.slip_pressure)
                push!(slip_potentials, 100.0)
            else
                # Interpolate to find slip potential
                # CDF uses normalized cumulative probability (0-100%)
                slip_potential = interpolate_cdf(fault_cdf.slip_pressure, fault_cdf.cumulative_probability, pressure)
                push!(slip_potentials, slip_potential)
            end
        end
        
        # Calculate statistics
        mean_pressure = mean(fault_pressures)
        mean_slip = mean(slip_potentials)
        std_slip = std(slip_potentials)
        min_slip = minimum(slip_potentials)
        max_slip = maximum(slip_potentials)
        
        push!(slip_potential_df, (fault_id, mean_pressure, mean_slip, std_slip, min_slip, max_slip))
    end
    
    return slip_potential_df
end

"""
Interpolate a value on a CDF
"""
function interpolate_cdf(x_values::Vector{Float64}, y_values::Vector{Float64}, x::Float64)
    # Handle case where vectors are empty
    if isempty(x_values) || isempty(y_values)
        @warn "Empty vectors provided to interpolate_cdf"
        return 0.0
    end
    
    # Find indices where x falls between x_values
    i = findfirst(v -> v >= x, x_values)
    
    # Handle edge cases
    if i === nothing
        return last(y_values)
    elseif i == 1
        return first(y_values)
    end
    
    # Linear interpolation
    """
    If x falls between two points (x_values[i-1] and x_values[i]), these lines define the coordinates of the interval:
        x1 and y1 are the left endpoint.
        x2 and y2 are the right endpoint.

    """
    x1 = x_values[i-1]
    x2 = x_values[i]
    y1 = y_values[i-1]
    y2 = y_values[i]
    
    # Linear interpolation formula
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
end

# Add an overload to handle DataFrameColumns
function interpolate_cdf(x_values::Any, y_values::Any, x::Float64)
    # Convert to Vector{Float64} if needed
    x_vec = convert(Vector{Float64}, x_values)
    y_vec = convert(Vector{Float64}, y_values)
    
    return interpolate_cdf(x_vec, y_vec, x)
end

function main()

    println("\n=== Starting Probabilistic Hydrology Process ===")

    scratchPath = ARGS[1]

    helper = TexNetWebToolLaunchHelperJulia(scratchPath)

    
    hydro_model_type = get_parameter_value(helper, 5, "hydro_model_type")

    println("Hydro model type from the portal: $hydro_model_type")

    extrapolate_injection_rates = get_parameter_value(helper, 5, "extrapolate_injection_rates")

    # if it's not provided, we use the default value
    if extrapolate_injection_rates === nothing
        extrapolate_injection_rates = false
    else
        extrapolate_injection_rates = parse(Bool, extrapolate_injection_rates)
    end

    
    # Default to probabilistic if not specified
    if hydro_model_type === nothing
        hydro_model_type = "probabilistic"  # Default to probabilistic
        println("Hydrology model type not specified, defaulting to probabilistic")
    end

    # Get year of interest
    year_of_interest = get_parameter_value(helper, 5, "year_of_interest")
    if year_of_interest === nothing
        year_of_interest = year(today())  
        println("Year of interest not specified, defaulting to $year_of_interest")
    else
        # Check if year_of_interest is already an Int before parsing
        if isa(year_of_interest, String)
            year_of_interest = parse(Int, year_of_interest)
        elseif !isa(year_of_interest, Int)
            # Handle unexpected types if necessary, e.g., throw an error or log a warning
            error("Unexpected type for year_of_interest: $(typeof(year_of_interest))")
        end
        # If it's already an Int, we don't need to do anything
        println("Year of interest: $year_of_interest")
    end

    # Read probabilistic geomechanics CDF data
    # CONTINUE FROM HERE: configure 'prob_geomechanics_cdf_graph_data' in the portal
    prob_geo_results = get_dataset_file_path(helper, 5, "prob_geomechanics_cdf_graph_data")
    prob_geo_cdf = CSV.read(prob_geo_results, DataFrame)
    println("Loaded probabilistic geomechanics CDF data with $(nrow(prob_geo_cdf)) rows")

    if hydro_model_type == "deterministic"
        # deterministic hydrology
        println("Running deterministic hydrology model...")
        
        # Check if we have deterministic hydrology results from step 4
        # TO DO: verify that this properly parses that data from the previous step
        det_hydro_results = get_dataset_file_path(helper, 5, "deterministic_hydrology_results")
        det_hydro_df = CSV.read(det_hydro_results, DataFrame)
        println("Using existing deterministic hydrology results")
        
        
        # Calculate slip potential by combining with probabilistic geomechanics CDF
        #=
        Using Deterministic Hydrology to Calculate Slip Potential
        Once we have calculated the pore pressure added to a given fault and the 
        probability of fault slip as a function of pore pressure increase, we can calculate the 
        cumulative probability of fault slip. We do this simply by slicing the probability shown 
        by the CDF curve (from Figure 9) for the appropriate pore pressure increase
        =#
        
        # Calculate slip potential for all years in the dataset
        all_years_slip_potential = calculate_deterministic_slip_potential(prob_geo_cdf, det_hydro_df, nothing)
        
        # Save all years slip potential results
        save_dataframe_as_parameter!(helper, 5, "hydro_slip_potential_results_all_years", all_years_slip_potential)
        
        # Filter for the specific year of interest
        year_specific_slip_potential = all_years_slip_potential[all_years_slip_potential.Year .== year_of_interest, :]
        
        # Save the year-specific slip potential results
        #save_dataframe_as_parameter!(helper, 5, "hydro_slip_potential_results", year_specific_slip_potential)
        
        # Save the slip potential results
        #save_dataframe_as_parameter!(helper, 5, "slip_potential_results", slip_potential)
        
        # Pretty print the results
        println("\nSlip Potential Results (Deterministic Hydrology) for year $year_of_interest:")
        #pretty_table(year_specific_slip_potential)
        
        println("\nTotal Slip Potential Results (Deterministic Hydrology) - All Years:")
        println("Years: $(minimum(all_years_slip_potential.Year)) to $(maximum(all_years_slip_potential.Year))")
        println("Number of data points: $(nrow(all_years_slip_potential))")
        
        
        
        # Plot results (could call a function from julia_fsp_graphs.jl)
        # plot_det_hydro_combined_cdf(prob_geo_cdf, det_hydro_df)
        
    elseif hydro_model_type == "probabilistic"
        # probabilistic hydrology
        println("Running probabilistic hydrology model...")
        
        # Get hydrology parameters from the portal
        aquifer_thickness = get_parameter_value(helper, 5, "aquifer_thickness_ft")
        porosity = get_parameter_value(helper, 5, "porosity")
        permeability = get_parameter_value(helper, 5, "permeability_md")
        fluid_density = get_parameter_value(helper, 5, "fluid_density")
        dynamic_viscosity = get_parameter_value(helper, 5, "dynamic_viscosity")
        fluid_compressibility = get_parameter_value(helper, 5, "fluid_compressibility")
        rock_compressibility = get_parameter_value(helper, 5, "rock_compressibility")
        
        # Get uncertainty parameters
        aquifer_thickness_uncertainty = get_parameter_value(helper, 5, "aquifer_thickness_uncertainty")
        porosity_uncertainty = get_parameter_value(helper, 5, "porosity_uncertainty")
        permeability_uncertainty = get_parameter_value(helper, 5, "permeability_uncertainty")
        fluid_density_uncertainty = get_parameter_value(helper, 5, "fluid_density_uncertainty")
        dynamic_viscosity_uncertainty = get_parameter_value(helper, 5, "dynamic_viscosity_uncertainty")
        fluid_compressibility_uncertainty = get_parameter_value(helper, 5, "fluid_compressibility_uncertainty")
        rock_compressibility_uncertainty = get_parameter_value(helper, 5, "rock_compressibility_uncertainty")
        
        # Get number of iterations
        n_iterations = get_parameter_value(helper, 5, "hydro_mc_iterations")
        if n_iterations === nothing
            n_iterations = 750  # default
        end
        
        # Create the parameters structure
        params = HydrologyParams(
            aquifer_thickness,
            porosity,
            permeability,
            fluid_density,
            dynamic_viscosity,
            fluid_compressibility,
            rock_compressibility,
            Dict(
                "aquifer_thickness" => aquifer_thickness_uncertainty,
                "porosity" => porosity_uncertainty,
                "permeability" => permeability_uncertainty,
                "fluid_density" => fluid_density_uncertainty,
                "dynamic_viscosity" => dynamic_viscosity_uncertainty,
                "fluid_compressibility" => fluid_compressibility_uncertainty,
                "rock_compressibility" => rock_compressibility_uncertainty
            ),
            n_iterations
        )
        
        # Run Monte Carlo simulation
        prob_hydro_results, prob_hydro_stats = run_monte_carlo_hydrology(helper, params, "uniform", extrapolate_injection_rates)
        
        # Save the results
        #save_dataframe_as_parameter!(helper, 5, "prob_hydrology_results", prob_hydro_results)
        #save_dataframe_as_parameter!(helper, 5, "prob_hydrology_stats", prob_hydro_stats)

        # Generate probabilistic hydrology CDF data
        prob_hydro_cdf_data = prob_hydrology_cdf(prob_hydro_results)
        save_dataframe_as_parameter!(helper, 5, "prob_hydrology_cdf_graph_data", prob_hydro_cdf_data)
        
        # Calculate slip potential by combining with probabilistic geomechanics CDF
        slip_potential = calculate_probabilistic_slip_potential(prob_geo_cdf, prob_hydro_results)
        
        # Save the slip potential results
        save_dataframe_as_parameter!(helper, 5, "hydro_slip_potential_results", slip_potential)
        
        # Pretty print the statistics
        println("\nProbabilistic Hydrology Statistics:")
        pretty_table(prob_hydro_stats)
        
        println("\nSlip Potential Results (Probabilistic Hydrology):")
        pretty_table(slip_potential)
        
        # Plot results (could call a function from julia_fsp_graphs.jl)
        # plot_prob_hydro_combined_cdf(prob_geo_cdf, prob_hydro_results)
    end
    
    # Write the final args file
    write_final_args_file(helper, joinpath(helper.scratch_path, ARGS_FILE_NAME))

    # explicitly set this step's success state to true
    set_success_for_step_index!(helper, 5, true)
    
    # Write the results file
    write_results_file(helper)
    
    println("=== Probabilistic Hydrology Process Completed ===\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end



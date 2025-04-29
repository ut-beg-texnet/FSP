using CSV
using DataFrames
using PrettyTables
using Distributions
using Statistics
using Dates
using LinearAlgebra
using Interpolations




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
    aquifer_thickness::Union{Float64, Int64}
    porosity::Union{Float64, Int64}
    permeability::Union{Float64, Int64}
    fluid_density::Union{Float64, Int64}
    dynamic_viscosity::Union{Float64, Int64}
    fluid_compressibility::Union{Float64, Int64, Nothing}
    rock_compressibility::Union{Float64, Int64, Nothing}
    plus_minus::Dict{String, Union{Float64, Int64, Nothing}}
    n_iterations::Union{Int64, Nothing, Float64}
end



function run_monte_carlo_hydrology(helper::TexNetWebToolLaunchHelperJulia, 
                                  params::HydrologyParams, 
                                  distribution_type::String="uniform",
                                  extrapolate_injection_rates::Bool=false,
                                  year_of_interest_date::Date=Date(year_of_interest-1, 12, 31),
                                  year_of_interest::Int64=Dates.year(today()))

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
    #fault_data_path = get_dataset_file_path(helper, 2, "det_geomechanics_results")
    if fault_data_path === nothing
        error("Required fault dataset not found or accessible.")
    end

    fault_df = CSV.read(fault_data_path, DataFrame)

    
    
    

    
    num_faults = nrow(fault_df)
    println("Read faults_model_inputs_output.csv: $fault_df")
    
    # Extract actual fault IDs from the input file
    fault_id_col = "FaultID" in names(fault_df) ? "FaultID" : "ID"
    if !(fault_id_col in names(fault_df))
        # If neither FaultID nor ID exists, use a fallback column or create numeric IDs
        if "Name" in names(fault_df)
            fault_id_col = "Name"
        else
            println("Warning: No FaultID, ID, or Name column found in fault data. Using row indices as IDs.")
            fault_df[!, :AutoID] = string.(1:num_faults)
            fault_id_col = "AutoID"
        end
    end
    
    # Extract the actual fault IDs
    fault_ids = string.(fault_df[!, fault_id_col])
    println("Using fault IDs: $(join(fault_ids, ", "))")

    # create the matrix to store the results of the Monte Carlo simulations
    ppOnFaultMC = zeros(params.n_iterations, num_faults) # rows are iterations, columns are faults

    # Create storage for parameter samples across all Monte Carlo iterations
    hydro_samples = Dict{String, Vector{Float64}}()
    hydro_samples["aquifer_thickness"] = Float64[]
    hydro_samples["porosity"] = Float64[]
    hydro_samples["permeability"] = Float64[]
    hydro_samples["fluid_density"] = Float64[]
    hydro_samples["dynamic_viscosity"] = Float64[]
    hydro_samples["fluid_compressibility"] = Float64[]
    hydro_samples["rock_compressibility"] = Float64[]
    

    # Get injection well data
    injection_wells_filepath, injection_data_type = get_injection_dataset_path(helper, 5)
    if injection_wells_filepath === nothing
        error("No injection well dataset found. Please provide injection well data.")
    end

    
    
    # Explicitly read API Number as String to preserve leading zeros
    if injection_data_type == "injection_tool_data"
        injection_wells_df = CSV.read(injection_wells_filepath, DataFrame, types = Dict(
            "API Number" => String,
            #"APINumber" => String,
            "UIC Number" => String
        ), validate = false)
    else
        injection_wells_df = CSV.read(injection_wells_filepath, DataFrame)
    end

    
    
    # Get unique well IDs based on data format - moved outside the loops for efficiency
    well_id_col = injection_data_type == "injection_tool_data" ? "API Number" : "WellID"
    well_ids = unique(injection_wells_df[!, well_id_col])
    
    # Pre-process well data outside the Monte Carlo loop
    println("Starting to process $(length(well_ids)) wells for pre-processing")
    well_info = Dict{String, Dict{String, Any}}()
    for well_id in well_ids
        println("Processing well ID: $well_id")
        
        if injection_data_type == "injection_tool_data"
            well_data = injection_wells_df[string.(injection_wells_df[!, well_id_col]) .== string(well_id), :]
            if isempty(well_data)
                # Try UIC Number
                well_data = injection_wells_df[string.(injection_wells_df[!, "UIC Number"]) .== string(well_id), :]
                if isempty(well_data)
                    continue
                end
            end
            
            # Get well coordinates 
            well_lat = first(well_data[!, "Surface Latitude"])
            well_lon = first(well_data[!, "Surface Longitude"])
            
            # find injection period
            dates = Date[]
            
            # Check the type of the date values first
            if eltype(well_data[!, "Date of Injection"]) <: Date
                dates = well_data[!, "Date of Injection"]
            else
                # Need to parse from strings
                try
                    # Try different date formats
                    dates = Date.(well_data[!, "Date of Injection"], dateformat"y-m-d")
                catch
                    try
                        dates = Date.(well_data[!, "Date of Injection"], dateformat"m/d/y")
                    catch
                        try
                            dates = Date.(well_data[!, "Date of Injection"], dateformat"m/d/yyyy")
                        catch e
                            continue
                        end
                    end
                end
            end
            
            if isempty(dates)
                continue
            end
            
            println("Successfully parsed dates for well $well_id")
            println("Minimum date: $(minimum(dates)), Maximum date: $(maximum(dates))")
            
            # Convert years to integers properly
            inj_start_year = Dates.year(minimum(dates))
            inj_end_year = Dates.year(maximum(dates))
            inj_start_date = minimum(dates)
            inj_end_date = min(maximum(dates), year_of_interest_date)
            
            println("Well $well_id: start_year=$inj_start_year ($(typeof(inj_start_year))), end_year=$inj_end_year ($(typeof(inj_end_year)))")
        else
            # Annual or monthly format
            well_data = injection_wells_df[string.(injection_wells_df[!, well_id_col]) .== string(well_id), :]
            if isempty(well_data)
                continue
            end
            
            # Get well coordinates directly
            well_lat = first(well_data[!, "Latitude(WGS84)"])
            well_lon = first(well_data[!, "Longitude(WGS84)"])
            
            # Get injection period
            # FSP Annual format has StartYear and EndYear columns
            # FSP Monthly format has Year column
            if "StartYear" in names(well_data)
                inj_start_year = first(well_data[!, "StartYear"])
                inj_end_year = first(well_data[!, "EndYear"])
                inj_start_date = Date(inj_start_year, 1, 1)
                inj_end_date = Date(inj_end_year-1, 12, 31)
            elseif "Year" in names(well_data)
                inj_start_year = minimum(well_data[!, "Year"])
                inj_start_month = minimum(well_data[well_data[!, "Year"] .== inj_start_year, "Month"])
                inj_start_date = Date(inj_start_year, inj_start_month, 1)
                inj_end_year = maximum(well_data[!, "Year"])
                inj_end_month = maximum(well_data[well_data[!, "Year"] .== inj_end_year, "Month"])
                inj_end_date = Date(inj_end_year, inj_end_month, 1)
                inj_end_date = lastdayofmonth(inj_end_date)
            else
                continue
            end
        end
        
        # Check if the well's injection period started after the year of interest
        if inj_start_date > year_of_interest_date
            continue
        end
        actual_end_date = min(inj_end_date, year_of_interest_date)
        if actual_end_date < inj_start_date
            continue
        end
        
        # Store well information for use in Monte Carlo loop
        well_info[string(well_id)] = Dict{String, Any}(
            "lat" => well_lat,
            "lon" => well_lon,
            "inj_start_year" => inj_start_year,
            "inj_start_date" => inj_start_date,
            "inj_end_year" => inj_end_year,
            "inj_end_date" => inj_end_date,
            "actual_end_date" => actual_end_date
        )
    end
    
    

    # PRE-PROCESS WELL DATA WITH EXTRAPOLATION (for efficiency)
    println("Pre-processing well injection data with extrapolation = $(extrapolate_injection_rates)...")
    prepared_well_data = Dict{String, Dict{String, Any}}()
    for (well_id, info) in well_info
        #println("Calling prepare_well_data_for_pressure_scenario for well $well_id")
        #println("Types: well_id=$(typeof(well_id)), inj_start_year=$(typeof(info["inj_start_year"])), inj_start_date=$(typeof(info["inj_start_date"])), inj_end_year=$(typeof(info["inj_end_year"])), actual_end_date=$(typeof(info["actual_end_date"])), injection_data_type=$(typeof(injection_data_type)), year_of_interest=$(typeof(year_of_interest)), extrapolate_injection_rates=$(typeof(extrapolate_injection_rates)), year_of_interest_date=$(typeof(year_of_interest_date))")
        # Prepare injection data with extrapolation once
        days, rates = prepare_well_data_for_pressure_scenario(
            injection_wells_df,
            well_id,
            info["inj_start_year"],
            info["inj_start_date"],
            info["inj_end_year"],
            info["actual_end_date"],
            injection_data_type,
            year_of_interest,
            extrapolate_injection_rates,
            year_of_interest_date
        )
        
        #println("Returned from prepare_well_data_for_pressure_scenario for well $well_id: got $(length(days)) days and $(length(rates)) rates")
        
        if !isempty(days) && !isempty(rates)
            # Store the data for pressure scenario 
            prepared_well_data[string(well_id)] = Dict{String, Any}(
                "days" => days,
                "rates" => rates,
                "lat" => info["lat"],
                "lon" => info["lon"]
            )
        end
    end
    

    println("Running Monte Carlo simulations for hydrology with $(length(prepared_well_data)) prepared wells...")
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
        
        # Store the sampled parameters for histogram visualization
        push!(hydro_samples["aquifer_thickness"], sampled_params["aquifer_thickness"])
        push!(hydro_samples["porosity"], sampled_params["porosity"])
        push!(hydro_samples["permeability"], sampled_params["permeability"])
        push!(hydro_samples["fluid_density"], sampled_params["fluid_density"])
        push!(hydro_samples["dynamic_viscosity"], sampled_params["dynamic_viscosity"])
        push!(hydro_samples["fluid_compressibility"], sampled_params["fluid_compressibility"])
        push!(hydro_samples["rock_compressibility"], sampled_params["rock_compressibility"])

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
            
            # Get fault ID for logging
            fault_id = fault_ids[f]
            
            # initialize pp
            
            ppOnFault = 0.0
            
            #println("  * Processing fault ID: $fault_id")
            
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
                )

                
                
                # Add to total pressure for this fault
                ppOnFault += pressure_contribution
            end

            # Ensure pressure is not negative
            ppOnFault = max(0.0, ppOnFault)
            
            
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
                ID = fault_ids[f],  # Use actual fault ID instead of numeric index
                Pressure = ppOnFaultMC[i, f]
            ))
        end
    end
    
    results_df = DataFrame(result_rows)
    
    # Return both the results dataframe and the parameter samples
    return results_df, hydro_samples
end

"""
Get the injection dataset path based on available data types
"""
function get_injection_dataset_path(helper::TexNetWebToolLaunchHelperJulia, step_index::Int)
    for param_name in ["injection_wells_annual_prob_hydro", "injection_wells_monthly_prob_hydro", "injection_tool_data_prob_hydro"]
        filepath = get_dataset_file_path(helper, step_index, param_name)
        if filepath !== nothing
            if param_name == "injection_wells_annual_prob_hydro"
                injection_data_type = "annual_fsp"
                return filepath, injection_data_type
            elseif param_name == "injection_wells_monthly_prob_hydro"
                injection_data_type = "monthly_fsp"
                return filepath, injection_data_type
            elseif param_name == "injection_tool_data_prob_hydro"
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
        ID = String[],
        Year = Int[],
        Date = Date[],
        slip_pressure = Float64[],
        probability = Float64[]
    )

    
    # Get unique fault IDs
    fault_ids = unique(prob_geo_cdf.ID)

    
    # Filter deterministic hydrology data by year if specified
    filtered_hydro_pressures = det_hydro_pressures
    if !isnothing(year_of_interest) && "Date" in names(det_hydro_pressures)
        filtered_hydro_pressures = det_hydro_pressures[year.(det_hydro_pressures.Date) .== year_of_interest, :]
        println("Filtered deterministic hydrology data to year $year_of_interest")
    end
    
    for fault_id in fault_ids
        # Get fault's geomechanics CDF
        fault_cdf = prob_geo_cdf[prob_geo_cdf.ID .== fault_id, :]
        
        # Get all deterministic pressure entries for this fault
        fault_pressures = filtered_hydro_pressures[filtered_hydro_pressures.ID .== fault_id, :]
        
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
            
            # Ensure pressure is not negative
            if pressure < 0.0
                slip_potential = 0.0
            # If pressure is less than minimum in CDF, slip potential is 0
            elseif pressure < minimum(fault_cdf[!, slipPressureCol])
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

    
    
    return slip_potential_df
end

"""
Calculate slip potential by combining probabilistic geomechanics CDF with probabilistic hydrology exceedance curve.
This integration calculates the overall fault slip potential.
"""
function calculate_probabilistic_slip_potential(prob_geo_cdf::DataFrame, prob_hydro_results::DataFrame)
    # this dataframe will store the results
    slip_potential_df = DataFrame(
        ID = String[],
        MeanPorePressure = Float64[],
        SlipPotential = Float64[]
    )

    # Get unique fault IDs
    fault_ids = unique(prob_geo_cdf.ID)
    
    # prob_hydro_results are the raw monte carlo results, so we need to convert them to an exceedance curve
    prob_hydro_cdf_data = prob_hydrology_cdf(prob_hydro_results)

    #=
    println("inside calculate_probabilistic_slip_potential, prob_hydro_results:")
    pretty_table(prob_hydro_results)

    println("inside calculate_probabilistic_slip_potential, prob_hydro_cdf_data:")
    pretty_table(prob_hydro_cdf_data)

    println("inside calculate_probabilistic_slip_potential, prob_geo_cdf:")
    pretty_table(prob_geo_cdf)
    =#
    
    # verify column names
    geo_pressure_col = "slip_pressure" in names(prob_geo_cdf) ? "slip_pressure" : "pressure"
    geo_prob_col = "probability" in names(prob_geo_cdf) ? "probability" : "cumulative_probability"
    
    for fault_id in fault_ids
        # Get fault's geomechanics CDF
        fault_geo_cdf = prob_geo_cdf[prob_geo_cdf.ID .== fault_id, :]
        sort!(fault_geo_cdf, geo_pressure_col)
        
        # Get hydrology exceedance data for this fault
        fault_hydro_exceedance = prob_hydro_cdf_data[prob_hydro_cdf_data.ID .== fault_id, :]
        
        if isempty(fault_hydro_exceedance)
            println("No hydrology exceedance data for fault $fault_id")
            continue
        end
        
        # Sort by pressure for intersection finding
        sort!(fault_hydro_exceedance, :slip_pressure)
        
        # Get mean pore pressure
        fault_pressures = prob_hydro_results[prob_hydro_results.ID .== fault_id, :Pressure]
        mean_pressure = mean(fault_pressures)
        # round this to 2 decimal places
        mean_pressure = round(mean_pressure, digits=2)
        #println("fault_pressures max for fault $fault_id: $(maximum(fault_pressures))")
        #println("fault_pressures min for fault $fault_id: $(minimum(fault_pressures))")
        
        # Check for curve overlap - only calculate non-zero FSP if curves intersect
        hydro_max_pressure = maximum(fault_hydro_exceedance.slip_pressure)
        hydro_min_pressure = minimum(fault_hydro_exceedance.slip_pressure)
        geo_max_pressure = maximum(fault_geo_cdf[!, geo_pressure_col])
        geo_min_pressure = minimum(fault_geo_cdf[!, geo_pressure_col])
        
        # Check if hydrology is entirely to the left of geomechanics
        if hydro_max_pressure < geo_min_pressure
            #println("Fault $fault_id: Hydrology curve entirely to left of geomechanics curve - FSP = 0.0")
            push!(slip_potential_df, (fault_id, mean_pressure, 0.0))
            continue
        end
        
        # Check if hydrology is entirely to the right of geomechanics
        if hydro_min_pressure > geo_max_pressure
            #println("Fault $fault_id: Hydrology curve entirely to right of geomechanics curve - FSP = 1.0")
            push!(slip_potential_df, (fault_id, mean_pressure, 1.0))
            continue
        end
        
        # Find the intersection point between the curves
        # Combine all pressure points for evaluation
        all_pressures = unique(vcat(fault_geo_cdf[!, geo_pressure_col], fault_hydro_exceedance.slip_pressure))
        sort!(all_pressures)
        
        # Filter to pressures where both curves are defined
        valid_pressures = filter(p -> 
            p >= max(geo_min_pressure, hydro_min_pressure) && 
            p <= min(geo_max_pressure, hydro_max_pressure), 
            all_pressures)
        
        # Check each pressure point to find where curves cross
        intersection_found = false
        intersection_pressure = 0.0
        intersection_probability = 0.0
        
        for i in 1:(length(valid_pressures)-1)
            p1 = valid_pressures[i]
            p2 = valid_pressures[i+1]
            
            # Evaluate both curves at p1
            hydro_prob1 = interpolate_cdf_new(
                fault_hydro_exceedance.slip_pressure, 
                fault_hydro_exceedance.probability, 
                p1)
            
            geo_prob1 = interpolate_cdf_new(
                fault_geo_cdf[!, geo_pressure_col], 
                fault_geo_cdf[!, geo_prob_col], 
                p1)
            
            # Evaluate both curves at p2
            hydro_prob2 = interpolate_cdf_new(
                fault_hydro_exceedance.slip_pressure, 
                fault_hydro_exceedance.probability, 
                p2)
            
            geo_prob2 = interpolate_cdf_new(
                fault_geo_cdf[!, geo_pressure_col], 
                fault_geo_cdf[!, geo_prob_col], 
                p2)
            
            # Check if curves cross between p1 and p2
            if (hydro_prob1 - geo_prob1) * (hydro_prob2 - geo_prob2) <= 0
                # Found a crossing point
                # Use linear interpolation to find exact intersection
                if hydro_prob1 == geo_prob1
                    # Exact intersection at p1
                    intersection_pressure = p1
                    intersection_probability = hydro_prob1
                elseif hydro_prob2 == geo_prob2
                    # Exact intersection at p2
                    intersection_pressure = p2
                    intersection_probability = hydro_prob2
                else
                    # Interpolate to find crossing point
                    t = (geo_prob1 - hydro_prob1) / ((hydro_prob2 - hydro_prob1) - (geo_prob2 - geo_prob1))
                    intersection_pressure = p1 + t * (p2 - p1)
                    
                    # Calculate probability at intersection point - using geomechanics curve
                    intersection_probability = geo_prob1 + t * (geo_prob2 - geo_prob1)
                end
                
                intersection_found = true
                #println("Fault $fault_id: Intersection found at pressure = $intersection_pressure psi, probability = $intersection_probability%")
                break
            end
        end
        
        # If no intersection found, use the higher of the two curves where they're closest
        if !intersection_found
            # Find the point where the curves are closest
            min_diff = Inf
            closest_pressure = 0.0
            closest_probability = 0.0
            
            for p in valid_pressures
                hydro_prob = interpolate_cdf_new(
                    fault_hydro_exceedance.slip_pressure, 
                    fault_hydro_exceedance.probability, 
                    p)
                
                geo_prob = interpolate_cdf_new(
                    fault_geo_cdf[!, geo_pressure_col], 
                    fault_geo_cdf[!, geo_prob_col], 
                    p)
                
                diff = abs(hydro_prob - geo_prob)
                
                if diff < min_diff
                    min_diff = diff
                    closest_pressure = p
                    closest_probability = max(hydro_prob, geo_prob)
                end
            end
            
            intersection_pressure = closest_pressure
            intersection_probability = closest_probability


            
            # round this to 2 decimal places
            #intersection_probability = round(intersection_probability, digits=2)
            #println("Fault $fault_id: No exact intersection found. Using closest point at pressure = $intersection_pressure psi, probability = $intersection_probability%")
        end
        
        # Add to results - use the intersection probability as the FSP
        push!(slip_potential_df, (fault_id, mean_pressure, intersection_probability))
    end
    
    return slip_potential_df
end







function interpolate_cdf_new(x_values::Vector{Float64}, y_values::Vector{Float64}, x::Float64)
    # handle edge cases
    if isempty(x_values) || isempty(y_values)
        @warn "Empty vectors provided to interpolate_cdf"
        return 0.0
    end

    # sort x-values and reorder y-values accordingly
    p = sortperm(x_values)
    x_sorted = x_values[p]
    y_sorted = y_values[p]

    # Explicitly deduplicate knots to avoid warning
    x_deduplicated = copy(x_sorted)
    y_deduplicated = copy(y_sorted)
    indices = Interpolations.deduplicate_knots!(x_deduplicated)
    # Only keep corresponding y values
    if length(indices) < length(y_deduplicated)
        y_deduplicated = y_deduplicated[indices]
    end

    # create interpolation object
    # we use Flat() which returns the end
    itp = LinearInterpolation(x_deduplicated, y_deduplicated, extrapolation_bc=Flat())

    # return the interpolated value at x
    return itp(x)
end

function interpolate_cdf_new(x_values::Any, y_values::Any, x::Float64)
    # Convert inputs to Vector{Float64} and use the new function
    x_vec = convert(Vector{Float64}, x_values)
    y_vec = convert(Vector{Float64}, y_values)
    return interpolate_cdf_new(x_vec, y_vec, x)
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
    
    # We need to find the indices where x falls between x_values
    i = findfirst(v -> v >= x, x_values)
    
    # edge cases
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

    year_of_interest_date = Date(year_of_interest-1, 12, 31)

    # Read probabilistic geomechanics CDF data
    # CONTINUE FROM HERE: configure 'prob_geomechanics_cdf_graph_data' in the portal
    prob_geo_results = get_dataset_file_path(helper, 5, "prob_geomechanics_cdf_graph_data_prob_hydro")
    prob_geo_cdf = CSV.read(prob_geo_results, DataFrame)
    println("Loaded original probabilistic geomechanics CDF data: (last 10 rows)")
    pretty_table(prob_geo_cdf[end-10:end, :])

    if hydro_model_type == "deterministic"
        # deterministic hydrology
        println("Running deterministic hydrology model...")
        
        # Check if we have deterministic hydrology results from step 4
        # TO DO: verify that this properly parses that data from the previous step
        det_hydro_results = get_dataset_file_path(helper, 5, "deterministic_hydrology_results")
        det_hydro_df = CSV.read(det_hydro_results, DataFrame)
        println("Using existing deterministic hydrology results:")
        
        

        # rename the 'FaultID' column to 'ID'
        rename!(det_hydro_df, :FaultID => :ID)
        pretty_table(det_hydro_df)
        
        
        # Calculate slip potential by combining with probabilistic geomechanics CDF
        #=
        Using Deterministic Hydrology to Calculate Slip Potential
        Once we have calculated the pore pressure added to a given fault and the 
        probability of fault slip as a function of pore pressure increase, we can calculate the 
        cumulative probability of fault slip. We do this simply by slicing the prob geomechanics 
        CDF for the appropriate pore pressure increase
        =#
        
        # Calculate slip potential for all years in the dataset
        deterministic_slip_potential = calculate_deterministic_slip_potential(prob_geo_cdf, det_hydro_df, nothing)

        println("deterministic_slip_potential:")
        pretty_table(deterministic_slip_potential)
        
       
        
        # Filter for the specific year of interest
        #year_specific_slip_potential = deterministic_slip_potential[deterministic_slip_potential.Year .== year_of_interest, :]
        
        # Save the year-specific slip potential results
        save_dataframe_as_parameter!(helper, 5, "prob_hydrology_cdf_graph_data", det_hydro_df)
        


        
        
        
        
        # set the model run parameter (read by the summare step)
        set_parameter_value!(helper, 5, "model_run", 0)
       
        
        
    elseif hydro_model_type == "probabilistic"

        # set the model run parameter (read by the summary step)
        set_parameter_value!(helper, 5, "model_run", 1)

        # probabilistic hydrology
        println("Running probabilistic hydrology model...")
        
        # Get hydrology parameters from the portal
        aquifer_thickness = get_parameter_value(helper, 4, "aquifer_thickness_ft")
        porosity = get_parameter_value(helper, 4, "porosity")
        permeability = get_parameter_value(helper, 4, "permeability_md")
        fluid_density = get_parameter_value(helper, 4, "fluid_density")
        dynamic_viscosity = get_parameter_value(helper, 4, "dynamic_viscosity")
        fluid_compressibility = get_parameter_value(helper, 4, "fluid_compressibility")
        rock_compressibility = get_parameter_value(helper, 4, "rock_compressibility")
        
        # Get uncertainty parameters
        aquifer_thickness_uncertainty = get_parameter_value(helper, 5, "aquifer_thickness_uncertainty")
        porosity_uncertainty = get_parameter_value(helper, 5, "porosity_uncertainty")
        permeability_uncertainty = get_parameter_value(helper, 5, "permeability_uncertainty")
        fluid_density_uncertainty = get_parameter_value(helper, 5, "fluid_density_uncertainty")
        dynamic_viscosity_uncertainty = get_parameter_value(helper, 5, "dynamic_viscosity_uncertainty")
        fluid_compressibility_uncertainty = get_parameter_value(helper, 5, "fluid_compressibility_uncertainty")
        rock_compressibility_uncertainty = get_parameter_value(helper, 5, "rock_compressibility_uncertainty")

        println("input parameters:")
        println("aquifer_thickness: $aquifer_thickness")
        println("porosity: $porosity")
        println("permeability: $permeability")
        println("fluid_density: $fluid_density")
        println("dynamic_viscosity: $dynamic_viscosity")
        println("fluid_compressibility: $fluid_compressibility")
        println("rock_compressibility: $rock_compressibility")

        
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
        
        # Run Monte Carlo simulation - updated to only receive one return value
        prob_hydro_results, hydro_samples = run_monte_carlo_hydrology(helper, params, "uniform", extrapolate_injection_rates, year_of_interest_date, year_of_interest)
        

        # save the prob_hydro_results as a CSV in the current directory
        #CSV.write("prob_hydro_results.csv", prob_hydro_results)
        
        # Create histogram data for hydrology parameter distributions
        # Use the actual samples from the Monte Carlo simulation
        
        # Get fault IDs to use for the histogram data
        fault_dataset_path = get_dataset_file_path(helper, 5, "faults_model_inputs_output")
        faults_df = CSV.read(fault_dataset_path, DataFrame)
        fault_ids = string.(faults_df.FaultID)
        
        # Create a structure for each fault that references the global parameters
        fault_param_values = Dict{String, Dict{String, Vector{Float64}}}()
        for fault_id in fault_ids
            fault_param_values[fault_id] = Dict{String, Vector{Float64}}()
            # Include all the relevant distributions for this fault
            fault_param_values[fault_id]["aquifer_thickness"] = hydro_samples["aquifer_thickness"]
            fault_param_values[fault_id]["porosity"] = hydro_samples["porosity"]
            fault_param_values[fault_id]["permeability"] = hydro_samples["permeability"]
            fault_param_values[fault_id]["fluid_density"] = hydro_samples["fluid_density"]
            fault_param_values[fault_id]["dynamic_viscosity"] = hydro_samples["dynamic_viscosity"]
            fault_param_values[fault_id]["fluid_compressibility"] = hydro_samples["fluid_compressibility"]
            fault_param_values[fault_id]["rock_compressibility"] = hydro_samples["rock_compressibility"]
            
            # Add the slip pressure for each fault
            fault_pressures = prob_hydro_results[prob_hydro_results.ID .== fault_id, :Pressure]
            fault_param_values[fault_id]["slip_pressure"] = fault_pressures
        end
        
        # Generate histogram data
        histogram_d3_data = input_distribution_histograms_to_d3(fault_param_values, Dict{String, Vector{Float64}}(), "hydrology", nbins=25)
        
        # Save histogram data
        # TO DO: uncomment those in production and configrue the graph in the portal
        #CSV.write("hydro_histogram_sample_data.csv", histogram_d3_data)
        #save_dataframe_as_parameter!(helper, 5, "prob_hydrology_histogram_data", histogram_d3_data)
        
        # Generate probabilistic hydrology CDF data
        prob_hydro_cdf_data = prob_hydrology_cdf(prob_hydro_results)
        save_dataframe_as_parameter!(helper, 5, "prob_hydrology_cdf_graph_data", prob_hydro_cdf_data)
        println("prob_hydrology_cdf_graph_data (after prob_hydrology_cdf) (last 10 rows):")
        pretty_table(prob_hydro_cdf_data[end-10:end, :])

        # save the prob_hydro_cdf_data as a CSV in the current directory
        #CSV.write("prob_hydro_cdf_data.csv", prob_hydro_cdf_data)

        
        
        
        
        
        # Calculate slip potential by combining with probabilistic geomechanics CDF
        println("calculating slip potential by combining with probabilistic geomechanics CDF...")
        slip_potential = calculate_probabilistic_slip_potential(prob_geo_cdf, prob_hydro_results)

        # add the fsp to the faults dataframe
        fault_dataset_path = get_dataset_file_path(helper, 1, "faults_model_inputs_output")
        faults_df = CSV.read(fault_dataset_path, DataFrame)

        # Always recreate the column to ensure it's mutable (it's initialized in the model inputs process with missing values)
        faults_df[!, :prob_hydro_fsp] = zeros(nrow(faults_df))
        
        # update 'prob_hydro_fsp' column with the fsp values
        for row in eachrow(slip_potential)
            fault_id = row.ID
            fsp = row.SlipPotential
            
            # find the matching row in the faults_df
            idx = findfirst(id -> string(id) == string(fault_id), faults_df.FaultID)
            if !isnothing(idx)
                faults_df[idx, :prob_hydro_fsp] = fsp
            end

        end

        save_dataframe_as_parameter!(helper, 5, "faults_with_prob_hydro_fsp", faults_df)
        
        

        println("slip_potential:")
        pretty_table(slip_potential)
        
        # Save the slip potential results
        #save_dataframe_as_parameter!(helper, 5, "hydro_slip_potential_results", slip_potential)
        
        
        
        println("\nSlip Potential Results (Probabilistic Hydrology):")
        # from the 'SlipPotential' column, round to 2 decimal places
        slip_potential.SlipPotential = round.(slip_potential.SlipPotential, digits=2)
        pretty_table(slip_potential)
        println("size of slip_potential: $(size(slip_potential))")
        save_dataframe_as_parameter!(helper, 5, "slip_potential_results", slip_potential)
        
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



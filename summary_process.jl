include("core/utilities.jl")
include("core/hydrology_calculations.jl")
include("TexNetWebToolLauncherHelperJulia.jl")
include("core/bill_pfront.jl")
include("graphs/julia_fsp_graphs.jl")

using DataFrames
using CSV
using Dates
using PrettyTables
using Statistics
using Plots
using Interpolations
using ProgressMeter
using Random

using .Utilities
using .HydroCalculations
using .TexNetWebToolLauncherHelperJulia
using .BillPFront
using .JuliaFSPGraphs

const RESULTS_FILE_NAME = "results.json"
const ARGS_FILE_NAME = "args.json"

"""
HydrologyParams struct to hold all hydrology parameters for Monte Carlo simulations
"""
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

"""
Run Monte Carlo hydrology simulations for all years up to year_of_interest
Returns a DataFrame with fault ID, pressure, and year columns
"""
function run_mc_hydrology_time_series(
    params::HydrologyParams,
    fault_df::DataFrame,
    injection_wells_df::DataFrame,
    years_to_analyze::Vector{Int},
    year_of_interest::Int,
    injection_data_type::String,
    distribution_type::String="uniform"
)
    # Create distributions for MC sampling
    distributions = Dict{String, Distribution}()
    
    if distribution_type == "uniform"
        distributions = Dict(
            "aquifer_thickness" => Utilities.create_uniform_distribution(params.aquifer_thickness, params.plus_minus["aquifer_thickness"]),
            "porosity" => Utilities.create_uniform_distribution(params.porosity, params.plus_minus["porosity"]),
            "permeability" => Utilities.create_uniform_distribution(params.permeability, params.plus_minus["permeability"]),
            "fluid_density" => Utilities.create_uniform_distribution(params.fluid_density, params.plus_minus["fluid_density"]),
            "dynamic_viscosity" => Utilities.create_uniform_distribution(params.dynamic_viscosity, params.plus_minus["dynamic_viscosity"]),
            "fluid_compressibility" => Utilities.create_uniform_distribution(params.fluid_compressibility, params.plus_minus["fluid_compressibility"]),
            "rock_compressibility" => Utilities.create_uniform_distribution(params.rock_compressibility, params.plus_minus["rock_compressibility"])
        )
    elseif distribution_type == "gaussian"
        # Example for future extension
        @warn "Gaussian distribution not fully implemented, defaulting to uniform"
        return run_mc_hydrology_time_series(params, fault_df, injection_wells_df, years_to_analyze, year_of_interest, injection_data_type, "uniform")
    else
        error("Unsupported distribution type: $distribution_type")
    end
    
    # Get well IDs
    well_id_col = injection_data_type == "injection_tool_data" ? "API Number" : "WellID"
    if !(well_id_col in names(injection_wells_df))
        # Try alternates
        if "APINumber" in names(injection_wells_df)
            well_id_col = "APINumber"
        elseif "UIC Number" in names(injection_wells_df)
            well_id_col = "UIC Number"
        else
            error("Could not identify well ID column in injection data")
        end
    end
    well_ids = unique(injection_wells_df[!, well_id_col])
    
    # Get fault IDs
    num_faults = nrow(fault_df)
    fault_id_col = "FaultID" in names(fault_df) ? "FaultID" : "ID"
    if !(fault_id_col in names(fault_df))
        @warn "No FaultID column found, using sequential IDs"
        fault_df[!, :TempID] = string.(1:num_faults)
        fault_id_col = "TempID"
    end
    fault_ids = string.(fault_df[!, fault_id_col])
    
    # Pre-process well data to get date boundaries
    inj_start_date, inj_end_date = Utilities.get_date_bounds(injection_wells_df)
    
    # Container for results
    # Structure: year -> iteration -> fault -> pressure
    results = Dict{Int, Dict{Int, Dict{String, Float64}}}()
    
    # Create years to analyze
    if isempty(years_to_analyze)
        years_to_analyze = year(inj_start_date):year_of_interest
    end
    
    println("Running Monte Carlo simulations for years: $years_to_analyze")
    println("Number of iterations: $(params.n_iterations)")
    
    # Main Monte Carlo loop
    p = Progress(params.n_iterations, desc="MC Simulations: ", dt=1.0)
    for i in 1:params.n_iterations
        # Sample parameters from distributions
        sampled_params = Dict(
            "aquifer_thickness" => rand(distributions["aquifer_thickness"]),
            "porosity" => rand(distributions["porosity"]),
            "permeability" => rand(distributions["permeability"]),
            "fluid_density" => rand(distributions["fluid_density"]),
            "dynamic_viscosity" => rand(distributions["dynamic_viscosity"]),
            "fluid_compressibility" => rand(distributions["fluid_compressibility"]),
            "rock_compressibility" => rand(distributions["rock_compressibility"])
        )
        
        # Calculate storativity and transmissivity
        S, T, rho = HydroCalculations.calcST(
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
        
        # Process each year
        for analysis_year in years_to_analyze
            # Set up year cutoff date (Dec 31 of the analysis year)
            cutoff_date = Date(analysis_year, 12, 31)
            
            # Initialize year results if needed
            if !haskey(results, analysis_year)
                results[analysis_year] = Dict{Int, Dict{String, Float64}}()
            end
            
            # Initialize iteration results
            results[analysis_year][i] = Dict{String, Float64}()
            
            # Process each fault
            for f in 1:num_faults
                fault_id = fault_ids[f]
                fault_lat = fault_df[f, "Latitude(WGS84)"]
                fault_lon = fault_df[f, "Longitude(WGS84)"]
                
                # Initialize total pressure for this fault
                total_pressure = 0.0
                
                # Process each well's contribution
                for well_id in well_ids
                    well_data = injection_wells_df[string.(injection_wells_df[!, well_id_col]) .== well_id, :]
                    
                    if isempty(well_data)
                        continue
                    end
                    
                    # Get well coordinates
                    lat_col = injection_data_type == "injection_tool_data" ? "Surface Latitude" : "Latitude(WGS84)"
                    lon_col = injection_data_type == "injection_tool_data" ? "Surface Longitude" : "Longitude(WGS84)"
                    
                    well_lat = first(well_data[!, lat_col])
                    well_lon = first(well_data[!, lon_col])
                    
                    # Get injection period
                    if injection_data_type == "annual_fsp"
                        inj_start_year = first(well_data[!, "StartYear"])
                        inj_start_date = Date(inj_start_year, 1, 1)
                        inj_end_year = first(well_data[!, "EndYear"])
                        inj_end_date = Date(inj_end_year-1, 12, 31)
                    elseif injection_data_type == "monthly_fsp"
                        inj_start_year = minimum(well_data[!, "Year"])
                        inj_start_month = minimum(well_data[well_data[!, "Year"] .== inj_start_year, "Month"])
                        inj_start_date = Date(inj_start_year, inj_start_month, 1)
                        inj_end_year = maximum(well_data[!, "Year"])
                        inj_end_month = maximum(well_data[well_data[!, "Year"] .== inj_end_year, "Month"])
                        inj_end_date = Date(inj_end_year, inj_end_month, 1)
                        inj_end_date = lastdayofmonth(inj_end_date)
                    elseif injection_data_type == "injection_tool_data"
                        dates = try
                            Date.(well_data[!, "Date of Injection"], dateformat"y-m-d")
                        catch
                            try
                                Date.(well_data[!, "Date of Injection"], dateformat"m/d/y")
                            catch
                                try
                                    Date.(well_data[!, "Date of Injection"], dateformat"m/d/yyyy")
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
                        inj_start_date = minimum(dates)
                        inj_end_date = maximum(dates)
                    else
                        error("Unsupported injection data type: $injection_data_type")
                    end
                    
                    # Skip if the well hasn't started injecting by the analysis year
                    if year(inj_start_date) > analysis_year
                        continue
                    end
                    
                    # Limit end date to the analysis year cutoff
                    actual_end_date = min(inj_end_date, cutoff_date)
                    actual_end_year = year(actual_end_date)
                    
                    # Prepare injection data
                    days, rates = Utilities.prepare_well_data_for_pressure_scenario(
                        well_data,
                        well_id,
                        inj_start_year,
                        inj_start_date,
                        actual_end_year,
                        actual_end_date,
                        injection_data_type,
                        analysis_year,
                        false,  # Don't extrapolate
                        cutoff_date
                    )
                    
                    if isempty(days) || isempty(rates)
                        continue
                    end
                    
                    # Calculate pressure contribution from this well
                    pressure_contribution = HydroCalculations.pfieldcalc_all_rates(
                        fault_lon,
                        fault_lat,
                        STRho,
                        days,
                        rates,
                        well_lon,
                        well_lat
                    )
                    
                    # Add to total pressure for this fault
                    total_pressure += pressure_contribution
                end
                
                # Store result for this fault and iteration
                results[analysis_year][i][fault_id] = total_pressure
            end
        end
        
        # Update progress bar
        next!(p)
    end
    
    # Convert nested dictionary to DataFrame
    result_rows = []
    for year in sort(collect(keys(results)))
        for iter in 1:params.n_iterations
            for (fault_id, pressure) in results[year][iter]
                push!(result_rows, (
                    IterationID = iter,
                    ID = fault_id,
                    Pressure = pressure,
                    Year = year
                ))
            end
        end
    end
    
    # Create DataFrame
    results_df = DataFrame(result_rows)
    return results_df
end

"""
Calculate fault slip potential by combining probabilistic geomechanics CDF with hydrology results
"""
function calculate_fault_slip_potential(prob_geo_cdf::DataFrame, prob_hydro_df::DataFrame)
    # Group data by year and fault ID
    grouped_hydro = groupby(prob_hydro_df, [:Year, :ID])
    
    # Container for results
    fsp_results = DataFrame(
        ID = String[],
        Year = Int[],
        FSP = Float64[]
    )
    
    # Process each group (year-fault combination)
    for group in grouped_hydro
        if isempty(group)
            continue
        end
        
        year = first(group.Year)
        fault_id = first(group.ID)
        
        # Get fault's geomechanics CDF
        fault_cdf = prob_geo_cdf[string.(prob_geo_cdf.ID) .== fault_id, :]
        
        if isempty(fault_cdf)
            @warn "No geomechanics CDF data found for fault ID $fault_id"
            continue
        end
        
        # Sort by slip pressure
        sort!(fault_cdf, :slip_pressure)
        
        # Calculate slip probabilities for each hydrology pressure value
        slip_probs = Float64[]
        
        for pressure in group.Pressure
            # If pressure is less than minimum in CDF, slip potential is 0
            if pressure < minimum(fault_cdf.slip_pressure)
                push!(slip_probs, 0.0)
            # If pressure is greater than maximum in CDF, slip potential is 1.0 (100%)
            elseif pressure > maximum(fault_cdf.slip_pressure)
                push!(slip_probs, 1.0)
            else
                # Interpolate to find slip potential
                slip_prob = Utilities.interpolate_cdf(
                    fault_cdf.slip_pressure, 
                    fault_cdf.probability ./ 100.0,  # Convert from percentage to 0-1 scale
                    pressure
                )
                push!(slip_probs, slip_prob)
            end
        end
        
        # Calculate mean slip probability for this fault and year
        mean_slip_prob = mean(slip_probs)
        
        # Add to results
        push!(fsp_results, (fault_id, year, mean_slip_prob))
    end
    
    # Sort results by year and ID
    sort!(fsp_results, [:Year, :ID])
    
    return fsp_results
end

"""
Generate combined summary of results including both geomechanics and hydrology
"""
function generate_summary_report(fsp_results::DataFrame, prob_hydro_df::DataFrame, fault_df::DataFrame)
    # Calculate statistics for hydrology results by year and fault
    hydro_stats = combine(groupby(prob_hydro_df, [:Year, :ID]), 
        :Pressure => mean => :MeanPressure,
        :Pressure => std => :StdDevPressure,
        :Pressure => minimum => :MinPressure,
        :Pressure => maximum => :MaxPressure
    )
    
    # Join with FSP results
    summary_report = innerjoin(fsp_results, hydro_stats, on=[:Year, :ID])
    
    # Add fault metadata if available
    if !isempty(fault_df)
        fault_id_col = "FaultID" in names(fault_df) ? "FaultID" : "ID"
        if fault_id_col in names(fault_df)
            # Select relevant columns from fault data
            fault_metadata = select(fault_df, 
                fault_id_col => :ID, 
                ["Strike", "Dip", "FrictionCoefficient", "slip_pressure"], 
                skipmissing=true
            )
            
            # Rename slip_pressure column to avoid confusion
            if "slip_pressure" in names(fault_metadata)
                rename!(fault_metadata, :slip_pressure => :DeterministicSlipPressure)
            end
            
            # Join with summary
            summary_report = leftjoin(summary_report, fault_metadata, on=:ID)
        end
    end
    
    return summary_report
end

# Function that parses the injection well dataset filepath from the portal
# Use the helper function to get the file path for the given parameter name
# we want ot accept three possible formats
# 1) FSP (annual)
# 2) FSP (monthly)
# 3) Injection Tool Data
function get_injection_dataset_path(helper::TexNetWebToolLaunchHelperJulia, step_index::Int)
    println("heyyyyyyyyyyyyyyy")
    #println("DEBUG: get_injection_dataset_path called with step_index = $step_index")
    for param_name in ["injection_wells_annual_summary", "injection_wells_monthly_summary", "injection_tool_data_summary"]
        #println("DEBUG: Trying to get file path for param_name = $param_name")
        filepath = get_dataset_file_path(helper, step_index, param_name)
        #println("DEBUG: filepath = $filepath, type = $(typeof(filepath))")
        if filepath !== nothing
            if param_name == "injection_wells_annual_summary"
                injection_data_type = "annual_fsp"
                println("DEBUG: Returning filepath = $filepath, type = $(typeof(filepath)), injection_data_type = $injection_data_type")
                return filepath, injection_data_type
            elseif param_name == "injection_wells_monthly_summary"
                injection_data_type = "monthly_fsp"
                println("DEBUG: Returning filepath = $filepath, type = $(typeof(filepath)), injection_data_type = $injection_data_type")
                return filepath, injection_data_type
            elseif param_name == "injection_tool_data_summary"
                injection_data_type = "injection_tool_data"
                println("DEBUG: Returning filepath = $filepath, type = $(typeof(filepath)), injection_data_type = $injection_data_type")
                return filepath, injection_data_type
            end
        end
    end
    
    println("DEBUG: No injection dataset found, returning nothing")
    return nothing, nothing
end

function main()
    println("\n=== Starting FSP Summary Process ===")

    # 1) Get the inputs from the args.json file
    scratchPath = ARGS[1]
    helper = TexNetWebToolLaunchHelperJulia(scratchPath)

    # Get year of interest
    year_of_interest = get_parameter_value(helper, 6, "year_of_interest_summary")
    if year_of_interest === nothing
        year_of_interest = Dates.year(Dates.today())
        println("Year of interest not provided, using current year: $year_of_interest")
    elseif !isa(year_of_interest, Int)
        year_of_interest = parse(Int, year_of_interest)
    end
    
    year_of_interest_date = Date(year_of_interest - 1, 12, 31)
    println("Using year of interest: $year_of_interest, cutoff date: $year_of_interest_date")

    # 2) Read injection wells data
    injection_wells_csv_filepath, injection_data_type = get_injection_dataset_path(helper, 6)
    if injection_wells_csv_filepath === nothing
        error("No injection wells dataset provided.")
    end
    injection_wells_df = CSV.read(injection_wells_csv_filepath, DataFrame)
    println("Loaded injection wells data: $(nrow(injection_wells_df)) records")

    # 3) Read fault data
    fault_data_path = get_dataset_file_path(helper, 6, "faults")
    if fault_data_path === nothing
        error("Required fault dataset not found or accessible.")
    end
    fault_df = CSV.read(fault_data_path, DataFrame)
    println("Loaded fault data: $(nrow(fault_df)) faults")

    # 4) Read probabilistic geomechanics results
    prob_geo_cdf_path = get_dataset_file_path(helper, 6, "prob_geomechanics_cdf_graph_data_summary")
    if prob_geo_cdf_path === nothing
        error("Probabilistic geomechanics CDF data not found.")
    end
    prob_geo_cdf = CSV.read(prob_geo_cdf_path, DataFrame)
    println("Loaded probabilistic geomechanics data: $(nrow(prob_geo_cdf)) rows")

    # 5) Get hydrology parameters
    aquifer_thickness = get_parameter_value(helper, 6, "aquifer_thickness_ft_summary")
    porosity = get_parameter_value(helper, 6, "porosity_summary")
    permeability = get_parameter_value(helper, 6, "permeability_md_summary")
    fluid_density = get_parameter_value(helper, 6, "fluid_density_summary")
    dynamic_viscosity = get_parameter_value(helper, 6, "dynamic_viscosity_summary")
    fluid_compressibility = get_parameter_value(helper, 6, "fluid_compressibility_summary")
    rock_compressibility = get_parameter_value(helper, 6, "rock_compressibility_summary")
    
    # Get uncertainty parameters
    aquifer_thickness_uncertainty = get_parameter_value(helper, 6, "aquifer_thickness_uncertainty_summary")
    porosity_uncertainty = get_parameter_value(helper, 6, "porosity_uncertainty_summary")
    permeability_uncertainty = get_parameter_value(helper, 6, "permeability_uncertainty_summary")
    fluid_density_uncertainty = get_parameter_value(helper, 6, "fluid_density_uncertainty_summary")
    dynamic_viscosity_uncertainty = get_parameter_value(helper, 6, "dynamic_viscosity_uncertainty_summary")
    fluid_compressibility_uncertainty = get_parameter_value(helper, 6, "fluid_compressibility_uncertainty_summary")
    rock_compressibility_uncertainty = get_parameter_value(helper, 6, "rock_compressibility_uncertainty_summary")

    # Get number of iterations
    n_iterations = get_parameter_value(helper, 6, "hydro_mc_iterations_summary")
    if n_iterations === nothing
        n_iterations = 500  # default
    elseif !isa(n_iterations, Int)
        n_iterations = parse(Int, n_iterations)
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

    # 6) Get injection date bounds and determine years to analyze
    inj_start_date, inj_end_date = Utilities.get_date_bounds(injection_wells_df)
    start_year = year(inj_start_date)
    end_year = min(year(inj_end_date), year_of_interest)
    years_to_analyze = start_year:end_year
    
    println("Injection period: $inj_start_date to $inj_end_date")
    println("Years to analyze: $years_to_analyze")

    #=

    # Get the injection data type directly
    injection_wells_filepath, injection_data_type = Utilities.get_injection_dataset_path(helper, 6)
    if injection_wells_filepath === nothing
        error("No injection well dataset found. Please provide injection well data.")
    end
    =#

    # 7) Run Monte Carlo simulation for probabilistic hydrology for each year
    println("\nRunning probabilistic hydrology Monte Carlo simulation for multiple years...")
    prob_hydro_results = run_mc_hydrology_time_series(
        params, 
        fault_df, 
        injection_wells_df, 
        collect(years_to_analyze),
        year_of_interest,
        injection_data_type
    )
    
    # Save hydrology results
    save_dataframe_as_parameter!(helper, 6, "prob_hydrology_time_series", prob_hydro_results)
    println("Saved probabilistic hydrology results: $(nrow(prob_hydro_results)) records")

    # 8) Calculate fault slip potential by combining with geomechanics CDF
    println("\nCalculating fault slip potential...")
    fsp_results = calculate_fault_slip_potential(prob_geo_cdf, prob_hydro_results)
    
    # Save FSP results
    save_dataframe_as_parameter!(helper, 6, "fsp_summary", fsp_results)
    println("Saved fault slip potential summary: $(nrow(fsp_results)) records")

    # 9) Generate comprehensive summary report
    summary_report = generate_summary_report(fsp_results, prob_hydro_results, fault_df)
    save_dataframe_as_parameter!(helper, 6, "summary_report", summary_report)
    println("Saved comprehensive summary report: $(nrow(summary_report)) records")

    # 10) Print summary statistics
    println("\nSummary Statistics by Year:")
    yearly_stats = combine(groupby(fsp_results, :Year), 
        :FSP => mean => :MeanFSP,
        :FSP => maximum => :MaxFSP,
        nrow => :NumFaults
    )
    pretty_table(yearly_stats)

    # 11) Create visualizations
    year_color_map = Dict(zip(sort(collect(years_to_analyze)), 
                             range(colorant"blue", colorant"red", length(years_to_analyze))))
    
    # Plot FSP over time for each fault
    p1 = plot(title="Fault Slip Potential Over Time", 
         xlabel="Year", ylabel="Slip Potential", 
         legend=:outertopright, legendfontsize=8)
    
    for fault_id in unique(fsp_results.ID)
        fault_data = fsp_results[fsp_results.ID .== fault_id, :]
        plot!(p1, fault_data.Year, fault_data.FSP, 
              label="Fault $fault_id", marker=:circle, linewidth=2)
    end
    
    # Save the plots
    savefig(p1, "fsp_time_series.png")
    println("\nVisualization saved as fsp_time_series.png")

    # Finalize
    write_final_args_file(helper, joinpath(helper.scratch_path, ARGS_FILE_NAME))
    set_success_for_step_index!(helper, 6, true)
    write_results_file(helper)
    
    println("=== FSP Summary Process Completed ===\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
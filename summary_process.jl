include("core/utilities.jl")
include("core/hydrology_calculations.jl")
include("TexNetWebToolLauncherHelperJulia.jl")
include("core/bill_pfront.jl")
include("graphs/julia_fsp_graphs.jl")

using DataFrames
using CSV
using Dates
#using PrettyTables
using Statistics
using Interpolations
using Random
using Distributions
using Base.Threads  # Add Threads for parallelization

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
Pre-process well data to avoid expensive operations in Monte Carlo loops
Returns a dictionary with well_id as key and processed well info as value
"""
function preprocess_well_data(injection_wells_df::DataFrame, well_id_col::String, injection_data_type::String)
    #println("Pre-processing well data for optimization...")
    
    # Determine column names once
    lat_col = injection_data_type == "injection_tool_data" ? "Surface Latitude" : "Latitude(WGS84)"
    lon_col = injection_data_type == "injection_tool_data" ? "Surface Longitude" : "Longitude(WGS84)"
    
    # Group wells by ID and extract all needed info
    well_info = Dict{String, NamedTuple}()
    
    for well_id in unique(injection_wells_df[!, well_id_col])
        well_id_str = string(well_id)
        
        # Filter once per well (not millions of times)
        well_data = injection_wells_df[string.(injection_wells_df[!, well_id_col]) .== well_id_str, :]
        
        if isempty(well_data)
            continue
        end
        
        # Extract coordinates once
        well_lat = first(well_data[!, lat_col])
        well_lon = first(well_data[!, lon_col])
        
        # Process dates once based on injection type
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
            # Parse dates once
            dates = []
            if eltype(well_data[!, "Date of Injection"]) <: Date
                dates = well_data[!, "Date of Injection"]
            else
                # Try parsing with different formats
                try
                    dates = Date.(well_data[!, "Date of Injection"], dateformat"y-m-d")
                catch
                    try
                        dates = Date.(well_data[!, "Date of Injection"], dateformat"m/d/y")
                    catch
                        try
                            dates = Date.(well_data[!, "Date of Injection"], dateformat"m/d/yyyy")
                        catch e
                            @warn "Could not parse dates for well $well_id_str: $e"
                            continue
                        end
                    end
                end
            end
            
            if isempty(dates)
                continue
            end
            
            inj_start_date = minimum(dates)
            inj_end_date = maximum(dates)
        else
            error("Unsupported injection data type: $injection_data_type")
        end
        
        # Store all processed info
        well_info[well_id_str] = (
            data = well_data,
            latitude = well_lat,
            longitude = well_lon,
            start_date = inj_start_date,
            end_date = inj_end_date,
            start_year = year(inj_start_date),
            end_year = year(inj_end_date)
        )
    end
    
    #println("Pre-processed $(length(well_info)) wells successfully")
    return well_info
end

"""
Run Monte Carlo hydrology simulations for all years up to year_of_interest
Returns a DataFrame with fault ID, pressure, and year columns
Uses multi-threading with a ReentrantLock to avoid race conditions
"""
function run_mc_hydrology_time_series(
    params::HydrologyParams,
    fault_df::DataFrame,
    injection_wells_df::DataFrame,
    years_to_analyze::Vector{Int},
    year_of_interest::Union{Int, Nothing},
    injection_data_type::String,
    distribution_type::String="uniform"
)
    
    #println("Running Monte Carlo simulations with $(nthreads()) threads")
    
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
        elseif "Well ID" in names(injection_wells_df)
            well_id_col = "Well ID"
        elseif "UWI" in names(injection_wells_df)
            well_id_col = "UWI"
        else
            error("Could not identify well ID column in injection data")
        end
    end
    
    # PRE-PROCESS WELL DATA FOR OPTIMIZATION
    well_info = preprocess_well_data(injection_wells_df, well_id_col, injection_data_type)
    well_ids = collect(keys(well_info))  # Use pre-processed well IDs
    
    # For injection tool data format, verify date column
    if injection_data_type == "injection_tool_data"
        if !("Date of Injection" in names(injection_wells_df))
            error("'Date of Injection' column not found in injection well data")
        end
    end
    
    
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
    
    # Find the max end date of all injections
    max_injection_year = year(inj_end_date)
    
    
    # Structure: year -> iteration -> fault -> pressure
    results = Dict{Int, Dict{Int, Dict{String, Float64}}}()
    
    # Create range of years to analyze
    if isempty(years_to_analyze)
        years_to_analyze = year(inj_start_date):year_of_interest
    end
    
    # Initialize the results structure to avoid race conditions during parallel writing
    for analysis_year in years_to_analyze
        results[analysis_year] = Dict{Int, Dict{String, Float64}}()
        for i in 1:params.n_iterations
            results[analysis_year][i] = Dict{String, Float64}()
        end
    end
    
    # We'll use a mutex to prevent race conditions when updating the results
    # this allows only one thread/process to update the results at a time
    results_lock = ReentrantLock()
    
    # Prepare iteration indices for parallelization
    iterations = collect(1:params.n_iterations)
    
    # Process Monte Carlo iterations in parallel
    @threads for i in iterations
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
        
        # Process each year for this iteration
        for analysis_year in years_to_analyze
            # Set up year cutoff date (Dec 31 of the analysis year)
            cutoff_date = Date(analysis_year, 12, 31)
            
            # Local results for this iteration and year
            local_results = Dict{String, Float64}()
            
            # Process each fault
            for f in 1:num_faults
                fault_id = fault_ids[f]
                fault_lat = fault_df[f, "Latitude(WGS84)"]
                fault_lon = fault_df[f, "Longitude(WGS84)"]
                
                # Initialize total pressure for this fault
                total_pressure = 0.0
                
                # Process each well's contribution
                for well_id in well_ids
                    
                    # Use pre-processed well data (OPTIMIZED)
                    well = well_info[well_id]
                    
                    # Skip if the well hasn't started injecting by the analysis year
                    if well.start_year > analysis_year
                        continue
                    end
                    
                    # Use pre-processed coordinates and dates
                    well_lat = well.latitude
                    well_lon = well.longitude
                    inj_start_date = well.start_date
                    inj_end_date = well.end_date
                    
                    # Limit end date to the analysis year cutoff
                    actual_end_date = min(inj_end_date, cutoff_date)
                    actual_end_year = year(actual_end_date)
                    
                    # Prepare injection data using pre-processed well data
                    days, rates = Utilities.prepare_well_data_for_pressure_scenario(
                        well.data,  # Use pre-filtered DataFrame
                        String(well_id),
                        well.start_year,  # Use pre-processed start year
                        inj_start_date,
                        actual_end_year,
                        actual_end_date,
                        injection_data_type,
                        analysis_year,
                        false,  # Don't extrapolate
                        cutoff_date # December 31 of the analysis year
                    )
                    
                    if isempty(days) || isempty(rates)
                        continue
                    end
                    
                    # Calculate days from injection start to analysis date
                    evaluation_days_from_start = Float64((cutoff_date - inj_start_date).value + 1)
                    
                    # Calculate pressure contribution from this well
                    pressure_contribution = HydroCalculations.pfieldcalc_all_rates(
                        fault_lon, #fault longitude
                        fault_lat, #fault latitude
                        STRho, #storativity, transmissivity, fluid density
                        days, #days of injection
                        rates, #rates of injection
                        well_lon, #well longitude
                        well_lat, #well latitude
                        evaluation_days_from_start #days from injection start to evaluation date
                    )
                    
                    # Add to total pressure for this fault
                    total_pressure += pressure_contribution
                end
                
                # Ensure no negative pressure values
                total_pressure = max(0.0, total_pressure)
                
                # Store result for this fault and iteration
                local_results[fault_id] = total_pressure
            end
            
            # We update the results dictionary using a lock to avoid race conditions
            lock(results_lock) do
                for (fault_id, pressure) in local_results
                    results[analysis_year][i][fault_id] = pressure
                end
            end
        end
    end
    
    # Convert dictionary to df
    result_rows = []
    for year in sort(collect(keys(results)))
        for iter in 1:params.n_iterations
            for (fault_id, pressure) in results[year][iter]
                push!(result_rows, (
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
Calculate fault slip potential for each year using probabilistic hydrology results
"""
function calculate_fault_slip_potential(prob_geo_cdf::DataFrame, prob_hydro_df::DataFrame)
    # Container for the time series results
    fsp_through_time = DataFrame(
        ID = String[],
        Year = Int[],
        FSP = Float64[],
        epoch_time = Float64[]
    )
    
    # Group data by year
    years = unique(prob_hydro_df.Year)
    
    for year in years
        # Filter data for this year
        year_data = filter(row -> row.Year == year, prob_hydro_df)
        
        # Get unique fault IDs for this year
        fault_ids = unique(year_data.ID)
        
        for fault_id in fault_ids
        # Get fault's geomechanics CDF
            fault_geo_cdf = prob_geo_cdf[prob_geo_cdf.ID .== fault_id, :]

            if isempty(fault_geo_cdf)
                @warn "No geomechanics CDF data found for fault ID $fault_id"
                continue
            end
        
            # Get hydrology data for this fault and year
            fault_hydro_data = year_data[year_data.ID .== fault_id, :]
        
            if isempty(fault_hydro_data)
                @warn "No hydrology data found for fault ID $fault_id in year $year"
            continue
        end
        
            # Convert Monte Carlo results to an exceedance curve
            fault_hydro_exceedance = JuliaFSPGraphs.prob_hydrology_cdf(fault_hydro_data)
            
            if isempty(fault_hydro_exceedance)
                @warn "Failed to generate exceedance curve for fault ID $fault_id in year $year"
                continue
            end
            
            # Calculate mean pore pressure for this fault and year
            fault_pressures = fault_hydro_data.Pressure
            mean_pressure = mean(fault_pressures)
            
            # Verify column names
            geo_pressure_col = "slip_pressure" in names(fault_geo_cdf) ? "slip_pressure" : "pressure"
            geo_prob_col = "probability" in names(fault_geo_cdf) ? "probability" : "cumulative_probability"
            
            # Sort both datasets for intersection finding
            sort!(fault_geo_cdf, geo_pressure_col)
            sort!(fault_hydro_exceedance, :slip_pressure)
            
            # Check for curve overlap
            hydro_max_pressure = maximum(fault_hydro_exceedance.slip_pressure)
            hydro_min_pressure = minimum(fault_hydro_exceedance.slip_pressure)
            geo_max_pressure = maximum(fault_geo_cdf[!, geo_pressure_col])
            geo_min_pressure = minimum(fault_geo_cdf[!, geo_pressure_col])
            
            # Default FSP value
            slip_potential = 0.0
            
            # Check if hydrology is entirely to the left of geomechanics
            if hydro_max_pressure < geo_min_pressure
                slip_potential = 0.0
            # Check if hydrology is entirely to the right of geomechanics
            elseif hydro_min_pressure > geo_max_pressure
                slip_potential = 1.0
            else
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
                intersection_probability = 0.0
                
                for i in 1:(length(valid_pressures)-1)
                    p1 = valid_pressures[i]
                    p2 = valid_pressures[i+1]
                    
                    # Evaluate both curves at p1 and p2
                    # Use interpolate_cdf which now handles deduplication internally
                    hydro_prob1 = Utilities.interpolate_cdf(
                        fault_hydro_exceedance.slip_pressure, 
                        fault_hydro_exceedance.probability, 
                        p1)
                    
                    geo_prob1 = Utilities.interpolate_cdf(
                        fault_geo_cdf[!, geo_pressure_col], 
                        fault_geo_cdf[!, geo_prob_col], 
                        p1)
                    
                    hydro_prob2 = Utilities.interpolate_cdf(
                        fault_hydro_exceedance.slip_pressure, 
                        fault_hydro_exceedance.probability, 
                        p2)
                    
                    geo_prob2 = Utilities.interpolate_cdf(
                        fault_geo_cdf[!, geo_pressure_col], 
                        fault_geo_cdf[!, geo_prob_col], 
                        p2)
                    
                    # Check if curves cross between p1 and p2
                    if (hydro_prob1 - geo_prob1) * (hydro_prob2 - geo_prob2) <= 0
                        # Found a crossing point - use linear interpolation
                        if hydro_prob1 == geo_prob1
                            # Exact intersection at p1
                            intersection_probability = hydro_prob1
                        elseif hydro_prob2 == geo_prob2
                            # Exact intersection at p2
                            intersection_probability = hydro_prob2
                        else
                            # Interpolate to find crossing point
                            t = (geo_prob1 - hydro_prob1) / ((hydro_prob2 - hydro_prob1) - (geo_prob2 - geo_prob1))
                            # Calculate probability at intersection point
                            intersection_probability = geo_prob1 + t * (geo_prob2 - geo_prob1)
                        end
                        
                        intersection_found = true
                        break
                    end
                end
                
                # If no intersection found, use the higher of the two curves where they're closest
                if !intersection_found
                    # Find the point where the curves are closest
                    min_diff = Inf
                    closest_probability = 0.0
                    
                    for p in valid_pressures
                        # Use interpolate_cdf which now handles deduplication internally
                        hydro_prob = Utilities.interpolate_cdf(
                            fault_hydro_exceedance.slip_pressure, 
                            fault_hydro_exceedance.probability, 
                            p)
                        
                        geo_prob = Utilities.interpolate_cdf(
                            fault_geo_cdf[!, geo_pressure_col], 
                            fault_geo_cdf[!, geo_prob_col], 
                            p)
                        
                        diff = abs(hydro_prob - geo_prob)
                        
                        if diff < min_diff
                            min_diff = diff
                            closest_probability = max(hydro_prob, geo_prob)
                        end
                    end
                    
                    intersection_probability = closest_probability
                end
                
                slip_potential = intersection_probability
            end
        
        # Add to time series results with epoch timestamp
        epoch_timestamp = JuliaFSPGraphs.date_to_js_timestamp(Date(year, 1, 1))
            push!(fsp_through_time, (fault_id, year, slip_potential, epoch_timestamp))
        end
    end
    
    # Sort results by ID and Year
    sort!(fsp_through_time, [:ID, :Year])
    
    return fsp_through_time
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
                ["Strike", "Dip", "FrictionCoefficient", "slip_pressure"]
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



function get_injection_dataset_path_summary_step(helper::TexNetWebToolLaunchHelperJulia, step_index::Int)
    for param_name in ["injection_wells_annual_summary", "injection_wells_monthly_summary", "injection_tool_data_summary"]
        
        filepath = get_dataset_file_path(helper, step_index, param_name)
        if filepath !== nothing
            if param_name == "injection_wells_annual_summary"
                injection_data_type = "annual_fsp"
                return filepath, injection_data_type
            elseif param_name == "injection_wells_monthly_summary"
                
                injection_data_type = "monthly_fsp"
                return filepath, injection_data_type
            elseif param_name == "injection_tool_data_summary"
                injection_data_type = "injection_tool_data"
                return filepath, injection_data_type
            end
        else
            #println("$param_name not found")
        end
    end
            
            return nothing, nothing
        end

"""
Calculate fault slip potential with deterministic hydrology results
"""
function calculate_deterministic_fault_slip_potential(prob_geo_cdf::DataFrame, det_hydro_df::DataFrame, years_to_analyze::Vector{Int})
    # Container for the time series results
    fsp_through_time = DataFrame(
        ID = String[],
        Year = Int[],
        FSP = Float64[],
        epoch_time = Float64[]
    )
    
    # Get unique fault IDs
    fault_ids = unique(det_hydro_df.ID)
    
    for fault_id in fault_ids
        # Get fault's geomechanics CDF
        fault_cdf = prob_geo_cdf[string.(prob_geo_cdf.ID) .== fault_id, :]
        
        if isempty(fault_cdf)
            @warn "No geomechanics CDF data found for fault ID $fault_id"
            continue
        end
        
        # Sort by slip pressure
        sort!(fault_cdf, :slip_pressure)
        
        # Get hydrology results for this fault
        fault_pressures = det_hydro_df[det_hydro_df.ID .== fault_id, :]
        
        # Process each year's pressure for this fault
        for year in years_to_analyze
            # Find the row for this year, if it exists
            year_pressure_row = filter(row -> Dates.year(row.Date) == year, fault_pressures)
            
            if isempty(year_pressure_row)
                # No data for this year, skip or use 0
                push!(fsp_through_time, (
                    fault_id, 
                    year, 
                    0.0, 
                    JuliaFSPGraphs.date_to_js_timestamp(Date(year, 1, 1))
                ))
                continue
            end
            
            # Get pressure value for this year
            pressure = first(year_pressure_row.slip_pressure)
            
            # Calculate slip potential
            slip_prob = 0.0
            
            # If pressure is less than minimum in CDF, slip potential is 0
            if pressure < minimum(fault_cdf.slip_pressure)
                slip_prob = 0.0
            # If pressure is greater than maximum in CDF, slip potential is 1.0
            elseif pressure > maximum(fault_cdf.slip_pressure)
                slip_prob = 1.0
            else
                # Interpolate to find slip potential
                slip_prob = Utilities.interpolate_cdf(
                    fault_cdf.slip_pressure, 
                    fault_cdf.probability ./ 100.0,  # Convert from percentage to 0-1 scale
                    pressure
                )
            end
            
            # Add to time series results with epoch timestamp
            epoch_timestamp = JuliaFSPGraphs.date_to_js_timestamp(Date(year, 1, 1))
            push!(fsp_through_time, (fault_id, year, slip_prob, epoch_timestamp))
        end
    end
    
    # Sort results by ID and Year
    sort!(fsp_through_time, [:ID, :Year])
    
    return fsp_through_time
end

"""
Run deterministic hydrology calculations for all years up to year_of_interest
Returns a DataFrame with fault ID, pressure, and year columns
"""
function run_deterministic_hydrology_time_series(
    aquifer_thickness::Float64,
    porosity::Float64,
    permeability::Float64,
    fluid_density::Float64,
    dynamic_viscosity::Float64,
    fluid_compressibility::Float64,
    rock_compressibility::Float64,
    fault_df::DataFrame,
    injection_wells_df::DataFrame,
    years_to_analyze::Vector{Int},
    year_of_interest::Int,
    injection_data_type::String
)
    # Print key input parameter values for debugging
    #println("===== DEBUG: Input Parameters for Deterministic Hydrology =====")
    #println("aquifer_thickness: ", aquifer_thickness)
    #println("porosity: ", porosity)
    #println("permeability: ", permeability)
    #println("fluid_density: ", fluid_density)
    #println("dynamic_viscosity: ", dynamic_viscosity)
    #println("fluid_compressibility: ", fluid_compressibility)
    #println("rock_compressibility: ", rock_compressibility)
    #println("===================================")
    
    # Calculate storativity and transmissivity
    S, T, rho = HydroCalculations.calcST(
        aquifer_thickness,
        porosity,
        permeability,
        fluid_density,
        dynamic_viscosity,
        9.81,
        fluid_compressibility,
        rock_compressibility
    )
    
    STRho = (S, T, rho)
    #println("DEBUG: Calculated Storativity = $S, Transmissivity = $T")
    
    # Get well IDs
    well_id_col = injection_data_type == "injection_tool_data" ? "API Number" : "WellID"
    if !(well_id_col in names(injection_wells_df))
        # Try alternates
        if "APINumber" in names(injection_wells_df)
            well_id_col = "APINumber"
        elseif "UIC Number" in names(injection_wells_df)
            well_id_col = "UIC Number"
        elseif "Well ID" in names(injection_wells_df)
            well_id_col = "Well ID"
        elseif "UWI" in names(injection_wells_df)
            well_id_col = "UWI"
        else
            error("Could not identify well ID column in injection data")
        end
    end
    
    # PRE-PROCESS WELL DATA FOR OPTIMIZATION (deterministic version)
    well_info = preprocess_well_data(injection_wells_df, well_id_col, injection_data_type)
    well_ids = collect(keys(well_info))  # Use pre-processed well IDs
    
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
    #println("Injection rate time window for all wells: inj_start_date = $inj_start_date, inj_end_date = $inj_end_date")
    
    # Container for results
    # Structure: year -> fault -> pressure
    results = Dict{Int, Dict{String, Float64}}()
    
    # Create years to analyze
    if isempty(years_to_analyze)
        years_to_analyze = year(inj_start_date):year_of_interest
    end
    
    # Find the max end date of all injections - we'll set the evaluation date to this for pressure diffusion
    max_injection_year = year(inj_end_date)
    #println("Maximum injection end date: $inj_end_date (year $max_injection_year)")
    
    # Process each year
    for analysis_year in years_to_analyze
        # Set up year cutoff date (Dec 31 of the analysis year)
        cutoff_date = Date(analysis_year, 12, 31)
        
        # Initialize year results if needed
        if !haskey(results, analysis_year)
            results[analysis_year] = Dict{String, Float64}()
        end
        
        # Process each fault
        for f in 1:num_faults
            fault_id = fault_ids[f]
            fault_lat = fault_df[f, "Latitude(WGS84)"]
            fault_lon = fault_df[f, "Longitude(WGS84)"]
            
            # Initialize total pressure for this fault
            total_pressure = 0.0
            
            # Process each well's contribution
            for well_id in well_ids
                # Use pre-processed well data (OPTIMIZED - deterministic)
                well = well_info[well_id]
                
                # Skip if the well hasn't started injecting by the analysis year
                if well.start_year > analysis_year
                    continue
                end
                
                # Use pre-processed coordinates and dates
                well_lat = well.latitude
                well_lon = well.longitude
                inj_start_date = well.start_date
                inj_end_date = well.end_date
                
                # Limit end date to the analysis year cutoff
                actual_end_date = min(inj_end_date, cutoff_date)
                actual_end_year = year(actual_end_date)
                
                # Prepare injection data using pre-processed well data
                days, rates = Utilities.prepare_well_data_for_pressure_scenario(
                    well.data,  # Use pre-filtered DataFrame
                    well_id,
                    well.start_year,  # Use pre-processed start year
                    inj_start_date,
                    actual_end_year,
                    actual_end_date,
                    injection_data_type,
                    analysis_year,
                    false,  # Don't extrapolate
                    cutoff_date # December 31 of the analysis year
                )
                
                if isempty(days) || isempty(rates)
                    continue
                end
                
                # Calculate days from injection start to analysis date
                evaluation_days_from_start = Float64((cutoff_date - inj_start_date).value + 1)
                
                # For years after injection has stopped, evaluation date should use the current year
                # to properly model pressure diffusion over time
                if analysis_year > max_injection_year
                    #println("Calculating pressure diffusion for year $analysis_year (after injection end)")
                end
                
                # Calculate pressure contribution from this well
                pressure_contribution = HydroCalculations.pfieldcalc_all_rates(
                    fault_lon, #fault longitude
                    fault_lat, #fault latitude
                    STRho, #storativity, transmissivity, fluid density
                    days, #days of injection
                    rates, #rates of injection
                    well_lon, #well longitude
                    well_lat, #well latitude
                    evaluation_days_from_start #days from injection start to evaluation date
                )
                
                # Add to total pressure for this fault
                total_pressure += pressure_contribution
            end
            
            # Ensure total pressure is non-negative before storing
            total_pressure = max(0.0, total_pressure)
            
            # Store result for this fault and year
            results[analysis_year][fault_id] = total_pressure
        end
    end
    
    # Convert nested dictionary to DataFrame
    result_rows = []
    for year in sort(collect(keys(results)))
        for (fault_id, pressure) in results[year]
            push!(result_rows, (
                ID = fault_id,
                Pressure = pressure,
                Year = year
            ))
        end
    end
    
    # Create DataFrame
    results_df = DataFrame(result_rows)
    return results_df
end

function main()
    # Start timing the entire script execution
    script_start_time = time()
    
    #println("\n=== Starting FSP Summary Process ===")

    # 1) Get the inputs from the args.json file
    scratchPath = ARGS[1]
    helper = TexNetWebToolLaunchHelperJulia(scratchPath)

    # Get year of interest
    #year_of_interest = get_parameter_value(helper, 6, "year_of_interest_summary")
    #=
    if year_of_interest === nothing
        year_of_interest = Dates.year(Dates.today())
    elseif !isa(year_of_interest, Int)
        year_of_interest = parse(Int, year_of_interest)
    end
    =#

    # print the number of threads used 
    #add_message_with_step_index!(helper, 6, "Number of threads used: $(nthreads())", 0)
    
    #year_of_interest_date = Date(year_of_interest - 1, 12, 31)

    # 2) Read injection wells data
    injection_wells_csv_filepath, injection_data_type = get_injection_dataset_path_summary_step(helper, 6)
    if injection_wells_csv_filepath === nothing
        error("No injection wells dataset provided.")
    elseif "WellID" in names(CSV.read(injection_wells_csv_filepath, DataFrame))
        # check if we have a 'WellID' column
        injection_wells_df = CSV.read(injection_wells_csv_filepath, DataFrame, types=Dict("WellID" => String), pool=false)
    elseif "API Number" in names(CSV.read(injection_wells_csv_filepath, DataFrame))
        # check if we have a 'API Number' column
        injection_wells_df = CSV.read(injection_wells_csv_filepath, DataFrame, types=Dict("API Number" => String), pool=false)
    elseif "APINumber" in names(CSV.read(injection_wells_csv_filepath, DataFrame))
        # check if we have a 'APINumber' column
        injection_wells_df = CSV.read(injection_wells_csv_filepath, DataFrame, types=Dict("APINumber" => String), pool=false)
    else
        error("No valid well ID column found in the injection wells dataset.")
    end

    # check the type of the WellID column
    #println("DEBUG: WellID column type = $(typeof(injection_wells_df[!, "WellID"]))")

    # Get the model that was run from the previous step
    model_run = get_parameter_value(helper, 6, "model_run_summary")
    if model_run === nothing
        # Default to probabilistic if not specified
        model_run = 1   
        #println("Model run type not specified, defaulting to probabilistic")
    end
    
    #println("Using model type: $model_run")

    # 3) Read fault data
    fault_data_path = get_dataset_file_path(helper, 6, "faults")
    if fault_data_path === nothing
        error("Required fault dataset not found or accessible.")
    end
    fault_df = CSV.read(fault_data_path, DataFrame, types=Dict("FaultID" => String), pool=false)

    # check the type of the FaultID column
    #println("DEBUG: FaultID column type = $(typeof(fault_df[!, "FaultID"]))")

   
    
    


    # 4) Read probabilistic geomechanics results
    prob_geo_cdf_path = get_dataset_file_path(helper, 6, "prob_geomechanics_cdf_graph_data_summary")
    if prob_geo_cdf_path === nothing
        error("Probabilistic geomechanics CDF data not found.")
    end
    prob_geo_cdf = CSV.read(prob_geo_cdf_path, DataFrame, types=Dict("ID" => String), pool=false)
    #println("Loaded probabilistic geomechanics data (first 10 rows) out of $(nrow(prob_geo_cdf)) rows:")
    
    #pretty_table(prob_geo_cdf[1:10, :])

    # 5) Get hydrology parameters
    aquifer_thickness = get_parameter_value(helper, 4, "aquifer_thickness_ft")
    porosity = get_parameter_value(helper, 4, "porosity")
    permeability = get_parameter_value(helper, 4, "permeability_md")
    fluid_density = get_parameter_value(helper, 4, "fluid_density")
    dynamic_viscosity = get_parameter_value(helper, 4, "dynamic_viscosity")
    fluid_compressibility = get_parameter_value(helper, 4, "fluid_compressibility")
    rock_compressibility = get_parameter_value(helper, 4, "rock_compressibility")

    for param in [aquifer_thickness, porosity, permeability, fluid_density, dynamic_viscosity, fluid_compressibility, rock_compressibility]
        if param === nothing
            error("Parameter $param is not found.")
        end
    end

    
    
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
    end_year = (year(inj_end_date) + 1) #since we evaluate up to the end of the last year, we need to add 1 to include the last year
    # we are also adding 2 extra years to the end of the analysis to account for the future and show the diffusion of the pressure
    end_year = end_year + 2
    
    years_to_analyze = start_year:end_year
    
    #println("Injection period (for all wells): $inj_start_date to $inj_end_date")
    #println("Years to analyze: $years_to_analyze")

    # Depending on model_run, either run deterministic or probabilistic hydrology
    if model_run == 0
        #println("\nRunning deterministic hydrology calculation for all years...")
        pressure_through_time_results = run_deterministic_hydrology_time_series(
            aquifer_thickness,
            porosity,
            permeability,
            fluid_density,
            dynamic_viscosity,
            fluid_compressibility,
            rock_compressibility,
            fault_df,
            injection_wells_df,
            collect(years_to_analyze),
            nothing,
            injection_data_type
        )
        
        # Aggregate results - not needed for deterministic as we don't have multiple iterations
        pressure_through_time_results_aggregated = select(pressure_through_time_results, [:ID, :Year, :Pressure])
        
        # Add epoch timestamp and round pressure values
        pressure_through_time_results_aggregated[!, :epoch_time] = [
            Float64(JuliaFSPGraphs.date_to_js_timestamp(Date(year, 1, 1))) 
            for year in pressure_through_time_results_aggregated.Year
        ]
        pressure_through_time_results_aggregated[!, :Pressure] = round.(pressure_through_time_results_aggregated.Pressure, digits=2)
        
        # Calculate fault slip potential
        #println("\nCalculating fault slip potential...")
        
        fsp_through_time = calculate_deterministic_fault_slip_potential(prob_geo_cdf, pressure_through_time_results_aggregated, years_to_analyze)
    elseif model_run == 1
        # Probabilistic model
        #println("\nRunning probabilistic hydrology Monte Carlo simulation for all years...")
        
        # Add timing measurement (use only in performance testing)
        start_time = time()
        pressure_through_time_results = run_mc_hydrology_time_series(
            params, 
            fault_df, 
            injection_wells_df, 
            collect(years_to_analyze),
            nothing,
            injection_data_type
        )
        end_time = time()
        elapsed_seconds = end_time - start_time
        
        # Store timing information for later output
        mc_elapsed_seconds = elapsed_seconds
        mc_total_iterations = params.n_iterations * length(collect(years_to_analyze)) * nrow(fault_df)

        # print the pressure_through_time results, without omitting any rows
        #println("pressure_through_time_results: ")
        #pretty_table(pressure_through_time_results)

            # Aggregate results
        pressure_through_time_results_aggregated = combine(
        groupby(pressure_through_time_results, [:ID, :Year]),
        :Pressure => mean => :Pressure
        )

            # Add epoch timestamp and round pressure values
        pressure_through_time_results_aggregated[!, :epoch_time] = [
            JuliaFSPGraphs.date_to_js_timestamp(Date(year, 1, 1)) 
            for year in pressure_through_time_results_aggregated.Year
        ]
        pressure_through_time_results_aggregated[!, :Pressure] = round.(pressure_through_time_results_aggregated.Pressure, digits=2)

        # Calculate fault slip potential using probabilistic approach
        #println("\nCalculating fault slip potential...")
        fsp_through_time = calculate_fault_slip_potential(prob_geo_cdf, pressure_through_time_results)
        
        # Round FSP values to 2 decimal places
        fsp_through_time[!, :FSP] = round.(fsp_through_time.FSP, digits=2)
        
        
    end

    # try to parse the 'ID' column as an integer, if that fails, keep the original 'ID' column
    try
        fsp_through_time[!, :ID] = parse.(Int, fsp_through_time.ID)
    catch
        fsp_through_time[!, :ID] = fsp_through_time.ID
    end
    
    # sort both dataframes by the 'ID' column

    # try to parse the 'ID' column as an integer, if it doesn't fail, then sort by the 'ID' column numerically
    # if it fails, then sort by the 'ID' column treating it as a string
    try
        fsp_through_time[!, :ID] = parse.(Int, fsp_through_time.ID)
    catch
        fsp_through_time[!, :ID] = fsp_through_time.ID
    end
    


    sort!(pressure_through_time_results_aggregated, [:ID])
    sort!(fsp_through_time, [:ID])

    save_dataframe_as_parameter!(helper, 6, "pressure_through_time_results", pressure_through_time_results_aggregated)
    
    # Save FSP results
    save_dataframe_as_parameter!(helper, 6, "fsp_through_time_results", fsp_through_time)

    #=
    # create a df with columns 'ID', 'x', and 'y'
    year_of_interest_line = DataFrame(
        ID = String[],
        x = Float64[],
        y = Float64[]
    )

    # Get epoch timestamp for January 1st of year_of_interest
    year_of_interest_epoch = Float64(JuliaFSPGraphs.date_to_js_timestamp(Date(year_of_interest, 1, 1)))
    
    # this will draw a simple line segment at the year of interest
    push!(year_of_interest_line, ("year_of_interest_line", year_of_interest_epoch, 0.0))
    push!(year_of_interest_line, ("year_of_interest_line", year_of_interest_epoch, 1.0))

    # save the df
    save_dataframe_as_parameter!(helper, 6, "year_of_interest_line", year_of_interest_line)
    =#

    # Finalize
    write_final_args_file(helper, joinpath(helper.scratch_path, ARGS_FILE_NAME))
    set_success_for_step_index!(helper, 6, true)
    write_results_file(helper)
    
    # Print performance metrics if we ran the probabilistic model (USE ONLY IN PERFORMANCE TESTING)
    
    if model_run == 1
        println("\n=== Performance Metrics ===")
        println("Monte Carlo simulation completed in $(round(mc_elapsed_seconds, digits=2)) seconds ($(round(mc_elapsed_seconds/60, digits=2)) minutes)")
        println("Processed $(mc_total_iterations) simulation points at $(round(mc_total_iterations/mc_elapsed_seconds, digits=2)) points/second")
        println("Using $(nthreads()) threads")
        println("==========================")
    end
    


    # Calculate total script execution time
    script_end_time = time()
    total_execution_time = script_end_time - script_start_time
    
    # Print performance metrics (USE ONLY IN PERFORMANCE TESTING)
    
    println("\n=== Total Script Execution Summary ===")
    println("Total execution time: $(round(total_execution_time, digits=2)) seconds ($(round(total_execution_time/60, digits=2)) minutes)")
    println("Number of threads used: $(nthreads())")
    if model_run == 1
        println("Monte Carlo portion: $(round(mc_elapsed_seconds, digits=2)) seconds ($(round(100 * mc_elapsed_seconds / total_execution_time, digits=1))% of total time)")
    end
    println("======================================")
    
    #println("=== FSP Summary Process Completed ===\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
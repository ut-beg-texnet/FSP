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
using Random
using Distributions

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
    year_of_interest::Union{Int, Nothing},
    injection_data_type::String,
    distribution_type::String="uniform"
)
    # Print key input parameter values for debugging
    println("===== DEBUG: Input Parameters =====")
    println("aquifer_thickness: ", params.aquifer_thickness)
    println("porosity: ", params.porosity)
    println("permeability: ", params.permeability)
    println("fluid_density: ", params.fluid_density)
    println("dynamic_viscosity: ", params.dynamic_viscosity)
    println("fluid_compressibility: ", params.fluid_compressibility)
    println("rock_compressibility: ", params.rock_compressibility)
    println("===================================")
    
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
    
    # When we work with data from the injection reporting tool, always treat the well ID as a string
    # They are API Number so we need to preserve leading/trailing zeros
    well_ids = string.(unique(injection_wells_df[!, well_id_col]))
    
    # Debug - print well IDs and dataframe columns
    #=
    println("===== DEBUG: Well Data =====")
    println("Well ID column: ", well_id_col)
    println("Number of wells: ", length(well_ids))
    println("Well IDs: ", well_ids)
    println("Injection data columns: ", names(injection_wells_df))
    =#
    
    # For injection tool data format, verify date column
    # We should expect them to already be Date objects, but we'll check
    if injection_data_type == "injection_tool_data"
        if !("Date of Injection" in names(injection_wells_df))
            #println("WARNING: 'Date of Injection' column not found in injection data")
            error("'Date of Injection' column not found in injection well data")
        else
            date_col = "Date of Injection"
            #println("Date column type: ", eltype(injection_wells_df[!, date_col]))
            
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

    println("Injection rate time window for all wells: inj_start_date = $inj_start_date, inj_end_date = $inj_end_date")
    
    # Container for results
    # Structure: year -> iteration -> fault -> pressure
    results = Dict{Int, Dict{Int, Dict{String, Float64}}}()
    
    # Create years to analyze
    if isempty(years_to_analyze)
        years_to_analyze = year(inj_start_date):year_of_interest
    end
    
    
    # Main Monte Carlo loop
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
        
        # Debug print only for first iteration
        if i == 1
            println("DEBUG: Calculated Storativity = $S, Transmissivity = $T")
        end
        
        
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
                    
                    # filter the injection wells df for the well id
                    well_data = injection_wells_df[string.(injection_wells_df[!, well_id_col]) .== string(well_id), :]
                    
                    if isempty(well_data)
                        continue
                    end
                    
                    # Get well coordinates
                    # if it's injection tool data, we need to use the 'Surface Latitude' and 'Surface Longitude' columns
                    # otherwise, we can use the 'Latitude(WGS84)' and 'Longitude(WGS84)' columns
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
                        dates = []
                        # Check if dates are already Date objects
                        if eltype(well_data[!, "Date of Injection"]) <: Date
                            dates = well_data[!, "Date of Injection"]
                        else
                            # Need to parse from strings
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
                        cutoff_date # December 31 of the analysis year
                    )
                    
                    if isempty(days) || isempty(rates)
                        continue
                    end
                    
                    # Calculate pressure contribution from this well
                    pressure_contribution = HydroCalculations.pfieldcalc_all_rates(
                        fault_lon, #fault longitude
                        fault_lat, #fault latitude
                        STRho, #storativity, transmissivity, fluid density
                        days, #days of injection
                        rates, #rates of injection
                        well_lon, #well longitude
                        well_lat #well latitude
                    )

                    
                    
                    # Add to total pressure for this fault
                    total_pressure += pressure_contribution
                end
                
                # Ensure no negative pressure values
                total_pressure = max(0.0, total_pressure)
                
                # Store result for this fault and iteration
                results[analysis_year][i][fault_id] = total_pressure
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
Calculate fault slip potential by combining probabilistic geomechanics CDF with hydrology results
"""
function calculate_fault_slip_potential(prob_geo_cdf::DataFrame, prob_hydro_df::DataFrame)

    #println("inside calculate_fault_slip_potential, prob_geo_cdf (first 10 rows of fault with ID 'F4'):")
    #pretty_table(prob_geo_cdf[prob_geo_cdf.ID .== "F4", :])

    #println("inside calculate_fault_slip_potential, prob_hydro_df (first 10 rows of fault with ID 'F4'):")
    # sort the dataframe by ID and then by Year
    sort!(prob_hydro_df, [:ID, :Year])
    #pretty_table(prob_hydro_df[prob_hydro_df.ID .== "F4", :])

    # Group data by year and fault ID
    grouped_hydro = groupby(prob_hydro_df, [:Year, :ID])
    
    # dataframe that will be used by D3 to plot the time series
    fsp_through_time = DataFrame(
        ID = String[],
        Year = Int[],
        FSP = Float64[],
        epoch_time = Float64[]
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
            # check that pressure is not negative
            if pressure < 0.0
                push!(slip_probs, 0.0)
            # If pressure is less than minimum in CDF, slip potential is 0
            elseif pressure < minimum(fault_cdf.slip_pressure)
                push!(slip_probs, 0.0)
            # If pressure is greater than maximum in CDF, slip potential is 1.0 (100%)
            elseif pressure > maximum(fault_cdf.slip_pressure)
                push!(slip_probs, 1.0)
            else

                # deduplicate values
                slip_p = fault_cdf.slip_pressure
                Interpolations.deduplicate_knots!(slip_p; move_knots=false) #this will suppress the warning

                # Interpolate to find slip potential
                slip_prob = Utilities.interpolate_cdf(
                    fault_cdf.slip_pressure, 
                    fault_cdf.probability,  
                    pressure
                )
                push!(slip_probs, slip_prob)
            end
        end
        
        # Calculate mean slip probability for this fault and year
        mean_slip_prob = mean(slip_probs)
        
        # Add to time series results with epoch timestamp
        epoch_timestamp = JuliaFSPGraphs.date_to_js_timestamp(Date(year, 1, 1))
        push!(fsp_through_time, (fault_id, year, mean_slip_prob, epoch_timestamp))
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
                prtinln("hellooooooooo")
                injection_data_type = "monthly_fsp"
                return filepath, injection_data_type
            elseif param_name == "injection_tool_data_summary"
                injection_data_type = "injection_tool_data"
                return filepath, injection_data_type
            end
        else
            println("$param_name not found")
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
    println("===== DEBUG: Input Parameters for Deterministic Hydrology =====")
    println("aquifer_thickness: ", aquifer_thickness)
    println("porosity: ", porosity)
    println("permeability: ", permeability)
    println("fluid_density: ", fluid_density)
    println("dynamic_viscosity: ", dynamic_viscosity)
    println("fluid_compressibility: ", fluid_compressibility)
    println("rock_compressibility: ", rock_compressibility)
    println("===================================")
    
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
    println("DEBUG: Calculated Storativity = $S, Transmissivity = $T")
    
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
    
    well_ids = string.(unique(injection_wells_df[!, well_id_col]))
    
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
    println("Injection rate time window for all wells: inj_start_date = $inj_start_date, inj_end_date = $inj_end_date")
    
    # Container for results
    # Structure: year -> fault -> pressure
    results = Dict{Int, Dict{String, Float64}}()
    
    # Create years to analyze
    if isempty(years_to_analyze)
        years_to_analyze = year(inj_start_date):year_of_interest
    end
    
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
                # Filter the injection wells df for the well id
                well_data = injection_wells_df[string.(injection_wells_df[!, well_id_col]) .== string(well_id), :]
                
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
                    dates = []
                    # Check if dates are already Date objects
                    if eltype(well_data[!, "Date of Injection"]) <: Date
                        dates = well_data[!, "Date of Injection"]
                    else
                        # Need to parse from strings
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
                    cutoff_date # December 31 of the analysis year
                )
                
                if isempty(days) || isempty(rates)
                    continue
                end
                
                # Calculate pressure contribution from this well
                pressure_contribution = HydroCalculations.pfieldcalc_all_rates(
                    fault_lon, #fault longitude
                    fault_lat, #fault latitude
                    STRho, #storativity, transmissivity, fluid density
                    days, #days of injection
                    rates, #rates of injection
                    well_lon, #well longitude
                    well_lat #well latitude
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
    println("\n=== Starting FSP Summary Process ===")

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
    
    #year_of_interest_date = Date(year_of_interest - 1, 12, 31)

    # 2) Read injection wells data
    injection_wells_csv_filepath, injection_data_type = get_injection_dataset_path_summary_step(helper, 6)
    if injection_wells_csv_filepath === nothing
        error("No injection wells dataset provided.")
    else
                injection_wells_df = CSV.read(injection_wells_csv_filepath, DataFrame)
    end

    # Get the model that was run from the previous step
    model_run = get_parameter_value(helper, 6, "model_run_summary")
    if model_run === nothing
        # Default to probabilistic if not specified
        model_run = 1
        println("Model run type not specified, defaulting to probabilistic")
    end
    
    println("Using model type: $model_run")

    # 3) Read fault data
    fault_data_path = get_dataset_file_path(helper, 6, "faults")
    if fault_data_path === nothing
        error("Required fault dataset not found or accessible.")
    end
    fault_df = CSV.read(fault_data_path, DataFrame)
    

    # 4) Read probabilistic geomechanics results
    prob_geo_cdf_path = get_dataset_file_path(helper, 6, "prob_geomechanics_cdf_graph_data_summary")
    if prob_geo_cdf_path === nothing
        error("Probabilistic geomechanics CDF data not found.")
    end
    prob_geo_cdf = CSV.read(prob_geo_cdf_path, DataFrame)
    println("Loaded probabilistic geomechanics data (first 10 rows) out of $(nrow(prob_geo_cdf)) rows:")
    
    pretty_table(prob_geo_cdf[1:10, :])

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

    println("aquifer_thickness: $aquifer_thickness")
    
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
    end_year = year(inj_end_date)
    years_to_analyze = start_year:end_year
    
    println("Injection period (for all wells): $inj_start_date to $inj_end_date")
    println("Years to analyze: $years_to_analyze")

    # Depending on model_run, either run deterministic or probabilistic hydrology
    if model_run == 0
        println("\nRunning deterministic hydrology calculation for all years...")
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
        println("\nCalculating fault slip potential...")
        
        fsp_through_time = calculate_deterministic_fault_slip_potential(prob_geo_cdf, pressure_through_time_results_aggregated, years_to_analyze)
    elseif model_run == 1
        # Probabilistic model
        println("\nRunning probabilistic hydrology Monte Carlo simulation for all years...")
        pressure_through_time_results = run_mc_hydrology_time_series(
            params, 
            fault_df, 
            injection_wells_df, 
            collect(years_to_analyze),
            nothing,
            injection_data_type
        )

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
        println("\nCalculating fault slip potential...")
        fsp_through_time = calculate_fault_slip_potential(prob_geo_cdf, pressure_through_time_results)
        println("fsp_through_time:")
        pretty_table(fsp_through_time)
    end
    
    # Common code for both models
    println("pressure_through_time_results_aggregated (Pressure through time graph): ")
    pretty_table(pressure_through_time_results_aggregated)

    save_dataframe_as_parameter!(helper, 6, "pressure_through_time_results", pressure_through_time_results_aggregated)

    # Round FSP values to 2 decimal places
    fsp_through_time[!, :FSP] = round.(fsp_through_time.FSP, digits=2)
    
    #=
    
    # Get the FSP and pressure values for the year of interest
    year_of_interest_fsp = filter(row -> row.Year == year_of_interest, fsp_through_time)
    year_of_interest_pressure = filter(row -> row.Year == year_of_interest, pressure_through_time_results_aggregated)
    
    println("\nFSP values for year of interest ($year_of_interest):")
    pretty_table(year_of_interest_fsp)
    
    println("\nPressure values for year of interest ($year_of_interest):")
    pretty_table(year_of_interest_pressure)
    
    
    # Add FSP and pressure columns to fault_df
    # First, ensure fault_df has "ID" column (matches the column name in fsp_through_time)
    fault_id_col = "FaultID" in names(fault_df) ? "FaultID" : "ID"
    if fault_id_col != "ID"
        rename!(fault_df, fault_id_col => "ID")
    end
    
    # Initialize summary columns with default values
    fault_df[!, :summary_fsp] .= 0.0
    fault_df[!, :summary_pressure] .= 0.0
    
    # Update the values for each fault that has data
    for row in eachrow(year_of_interest_fsp)
        fault_idx = findfirst(x -> x == row.ID, fault_df.ID)
        if fault_idx !== nothing
            fault_df[fault_idx, :summary_fsp] = row.FSP
        end
    end
    
    for row in eachrow(year_of_interest_pressure)
        fault_idx = findfirst(x -> x == row.ID, fault_df.ID)
        if fault_idx !== nothing
            fault_df[fault_idx, :summary_pressure] = row.Pressure
        end
    end
    
    println("\nFault dataframe with summary values:")
    pretty_table(fault_df[!, [:ID, :summary_fsp, :summary_pressure]])
    
    # this is the table that the user can select faults and change the graph views
    save_dataframe_as_parameter!(helper, 6, "faults_fsp_summary", fault_df)
    =#
    
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
    
    println("=== FSP Summary Process Completed ===\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
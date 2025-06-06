using CSV 
using DataFrames 
using Dates 
using Statistics 
#using PrettyTables
#using Printf
using InlineStrings



include("TexNetWebToolLauncherHelperJulia.jl")
include("core/hydrology_calculations.jl")
include("core/utilities.jl")
include("core/bill_pfront.jl")
include("core/geomechanics_model.jl")
include("graphs/julia_fsp_graphs.jl")
include("deterministic_geomechanics_process.jl")

using .TexNetWebToolLauncherHelperJulia
using .HydroCalculations
using .Utilities
using .BillPFront
using .GeomechanicsModel
using .JuliaFSPGraphs
using .GeomechanicsDriver

const ARGS_FILE_NAME = "args.json"
const RESULTS_FILE_NAME = "results.json"



# Function that parses the injection well dataset filepath from the portal
# Use the helper function to get the file path for the given parameter name
# we want ot accept three possible formats
# 1) FSP (annual)
# 2) FSP (monthly)
# 3) Injection Tool Data
function get_injection_dataset_path(helper::TexNetWebToolLaunchHelperJulia, step_index::Int)
    #println("DEBUG: get_injection_dataset_path called with step_index = $step_index")
    for param_name in ["injection_wells_annual_hydrology", "injection_wells_monthly_hydrology", "injection_tool_data_hydrology"]
        #println("DEBUG: Trying to get file path for param_name = $param_name")
        filepath = get_dataset_file_path(helper, step_index, param_name)
        #println("DEBUG: filepath = $filepath, type = $(typeof(filepath))")
        if filepath !== nothing
            if param_name == "injection_wells_annual_hydrology"
                injection_data_type = "annual_fsp"
                #println("DEBUG: Returning filepath = $filepath, type = $(typeof(filepath)), injection_data_type = $injection_data_type")
                return filepath, injection_data_type
            elseif param_name == "injection_wells_monthly_hydrology"
                injection_data_type = "monthly_fsp"
                #println("DEBUG: Returning filepath = $filepath, type = $(typeof(filepath)), injection_data_type = $injection_data_type")
                return filepath, injection_data_type
            elseif param_name == "injection_tool_data_hydrology"
                injection_data_type = "injection_tool_data"
                #println("DEBUG: Returning filepath = $filepath, type = $(typeof(filepath)), injection_data_type = $injection_data_type")
                return filepath, injection_data_type
            end
        end
    end
    
    #println("DEBUG: No injection dataset found, returning nothing")
    return nothing, nothing
end

# function that normalizes column names
# remove spaces, and make lowercase
function normalize_name(name::Union{String, Vector{String}})
    if isa(name, String)
        return lowercase(replace(name, " " => ""))
    else
        # Handle vector of strings
        return [lowercase(replace(s, " " => "")) for s in name]
    end
end


function main()

    
    
    scratchPath = ARGS[1]
    helper = TexNetWebToolLaunchHelperJulia(scratchPath)

    
    #println("LOADING MODEL PARAMETERS FROM THE PORTAL...")
    
    # Fluid compressibility (1/psi)
    beta = get_parameter_value(helper, 4, "fluid_compressibility")
    if beta === nothing
        add_message_with_step_index!(helper, 4, "Fluid compressibility was not provided, using the default value of 3.6e-10 1/psi", 0)
        beta = 3.6e-10
    end
    
    # Rock compressibility (1/psi)
    alphav = get_parameter_value(helper, 4, "rock_compressibility")
    if alphav === nothing
        add_message_with_step_index!(helper, 4, "Rock compressibility was not provided, using the default value of 1.08e-09 1/psi", 0)
        alphav = 1.08e-09 
    end
    
    # Fluid density (kg/m³)
    rho = get_parameter_value(helper, 4, "fluid_density")
    if rho === nothing
        add_message_with_step_index!(helper, 4, "Fluid density was not provided, using the default value of 1000 kg/m³", 0)
        rho = 1000.0  
    end
    
    
    # Dynamic viscosity (Pa·s)
    mu = get_parameter_value(helper, 4, "dynamic_viscosity")
    if mu === nothing
        add_message_with_step_index!(helper, 4, "Dynamic viscosity was not provided, using the default value of 0.0008 Pa·s", 0)
        mu = 0.0008
    end
    
    
    

    
    year_of_interest_str = get_parameter_value(helper, 4, "year_of_interest")
    
    if year_of_interest_str === nothing
        year_of_interest = Dates.year(Dates.today())
        #println("- Year of interest: $year_of_interest (default value)")
        add_message_with_step_index!(helper, 4, "Year of interest was not provided, using the current year ($year_of_interest) as the default value", 0)
    else
        # Check if already an integer, if not, try to parse from string
        if typeof(year_of_interest_str) == Int
            year_of_interest = year_of_interest_str
        else
            year_of_interest = parse(Int, year_of_interest_str)
        end
        #println("- Year of interest: $year_of_interest")
    end

    # we don't include the year of interest in the analysis, we evaluate injection rates up to the end of the previous year
    # ex. if year_of_interest = 2025, we evaluate injection rates up to the end of 2024
    year_of_interest_date = Date(year_of_interest-1, 12, 31)

    
    # Aquifer thickness (ft)
    h_feet = get_parameter_value(helper, 4, "aquifer_thickness_ft")
    if h_feet === nothing
        add_message_with_step_index!(helper, 4, "Aquifer thickness was not provided, using the default value of 100 ft", 0)
        h_feet = 100.0
    end

    
    
    # Porosity
    phi = get_parameter_value(helper, 4, "porosity")
    if phi === nothing
        add_message_with_step_index!(helper, 4, "Porosity was not provided, using the default value of 0.1 (10%)", 0)
        phi = 0.1
    end
    
    # Permeability (mD)
    kap_md = get_parameter_value(helper, 4, "permeability_md")
    if kap_md === nothing
        add_message_with_step_index!(helper, 4, "Permeability was not provided, using the default value of 200 mD", 0)
        kap_md = 200.0  
    end

    #=
    # Summary of all hydrology parameters
    #println("\n------------------------------------------------------")
    #println("MODEL PARAMETERS SUMMARY")
    #println("------------------------------------------------------")
    #println("Aquifer thickness: $h_feet ft")
    #println("Porosity: $phi")
    #println("Permeability: $kap_md mD")
    #println("Fluid density: $rho kg/m³")
    #println("Dynamic viscosity: $mu Pa·s")
    #println("Fluid compressibility: $beta 1/psi")
    #println("Rock compressibility: $alphav 1/psi")
    #println("Year of interest: $year_of_interest")
    #println("------------------------------------------------------\n")
    =#
    #println("LOADING INJECTION WELL DATA...")

    # first we need to get the injection well data
    injection_wells_dataset_filepath, injection_data_type = get_injection_dataset_path(helper, 4)

    # Check if injection data was found
    if injection_wells_dataset_filepath === nothing
        error("No injection well dataset found. Please provide injection well data.")
    end

    # if we have injection tool data format, explicitly parse the 'API Number' column as a string
    if injection_data_type == "injection_tool_data"
        injection_wells_df = CSV.read(injection_wells_dataset_filepath, DataFrame, types = Dict(
            "API Number" => String,
            "UIC Number" => String
        ), validate = false)
        
    else
        # Parse injection well data
        injection_wells_df = CSV.read(injection_wells_dataset_filepath, DataFrame)
    end

    #println("- Loaded well dataset with $(nrow(injection_wells_df)) rows and $(ncol(injection_wells_df)) columns")

    
    #latitude = Float64[]
    #longitude = Float64[]
    well_ids = String[]

    # GET THE UNIQUE WELL IDS
    if injection_data_type == "annual_fsp" || injection_data_type == "monthly_fsp"
        

        # get all unique well ids
        well_ids = try
            # Check for "WellID" column first for FSP formats
            if "WellID" in names(injection_wells_df)
                unique(string.(injection_wells_df[!, "WellID"]))
            
            end
        catch e
            @warn "Error getting well IDs: $e, available columns: $(names(injection_wells_df))"
            String[]
        end
        #println("- Found $(length(well_ids)) unique wells in $(injection_data_type) format")
    elseif injection_data_type == "injection_tool_data"
        
        well_ids = try 
            # For injection tool data format, check for "API Number" first
            if "API Number" in names(injection_wells_df)
                unique(string.(injection_wells_df[!, "API Number"]))
            elseif "APINumber" in names(injection_wells_df)
                unique(string.(injection_wells_df[!, "APINumber"]))
            else
                unique(string.(injection_wells_df[!, "UIC Number"]))
            end
        catch e
            @warn "Error getting well IDs: $e, available columns: $(names(injection_wells_df))"
            String[]
        end
        #println("- Found $(length(well_ids)) unique wells in injection_tool_data format")
    else
        error("Unsupported injection data type: $injection_data_type")
    end

    # check if well_ids is String3, if so, convert to String
    if eltype(well_ids) == String3
        well_ids = String.(well_ids)
    end

    


    # Calculate storativity and transmissivity
    #println("\nCALCULATING AQUIFER PROPERTIES...")
    S, T, rho_return = calcST(h_feet, phi, kap_md, rho, mu, 9.81, beta, alphav)

    # tuple for the pressureScenario function
    STRho = (S, T, rho)
    #println("- Storativity (S): $S")
    #println("- Transmissivity (T): $T m²/s")

    #println("LOADING FAULT DATA...")
    # Load fault data
    fault_data_path = get_dataset_file_path(helper, 4, "faults")
    if fault_data_path === nothing
        error("Required fault dataset not found or accessible.")
    end

    fault_df = CSV.read(fault_data_path, DataFrame)

    # check if the 'Strike' and 'Dip' columns are floats, if not, convert them to floats
    
    if !(eltype(fault_df[!, "Strike"]) <: Float64)
        fault_df[!, "Strike"] = Float64.(fault_df[!, "Strike"])
    end
    if !(eltype(fault_df[!, "Dip"]) <: Float64)
        fault_df[!, "Dip"] = Float64.(fault_df[!, "Dip"])
    end
    

    

    

    # Get grid bounds from args.json if provided, otherwise calculate based on well locations
    #println("\nSETTING UP PRESSURE FIELD GRID BOUNDS...")
    
    # Try to get lat/lon grid bounds from args.json
    lat_min_param = get_parameter_value(helper, 4, "grid_lat_min")
    lat_max_param = get_parameter_value(helper, 4, "grid_lat_max")
    lon_min_param = get_parameter_value(helper, 4, "grid_lon_min")
    lon_max_param = get_parameter_value(helper, 4, "grid_lon_max")

    

    
    # Check if all lat/lon parameters were provided
    if lat_min_param !== nothing && lat_max_param !== nothing && 
       lon_min_param !== nothing && lon_max_param !== nothing
        # User has provided grid bounds in lat/lon
        local lat_min = Float64(lat_min_param)
        local lat_max = Float64(lat_max_param)
        local lon_min = Float64(lon_min_param)
        local lon_max = Float64(lon_max_param)
        
    else
        # Calculate grid bounds based on well lat/lon with a buffer
        #println("- Calculating lat/lon grid bounds based on well locations...")
        
        # Determine which lat/lon columns to use based on data type
        lat_column = injection_data_type == "injection_tool_data" ? "Surface Latitude" : "Latitude(WGS84)"
        lon_column = injection_data_type == "injection_tool_data" ? "Surface Longitude" : "Longitude(WGS84)"

        lat_column_faults = "Latitude(WGS84)"
        lon_column_faults = "Longitude(WGS84)"

        
        # Check if we have well data
        if size(injection_wells_df, 1) == 0
            error("No well data available. Please provide well data.")
        elseif size(fault_df, 1) == 0
            error("No fault data available. Please provide fault data.")
        else
            # for dynamic bounds, get the (lat_max - lat_min) and (lon_max - lon_min), and add 10% to each
            # First need to combine the arrays
            all_lats = vcat(injection_wells_df[!, lat_column], fault_df[!, lat_column_faults])
            all_lons = vcat(injection_wells_df[!, lon_column], fault_df[!, lon_column_faults])
            
            lat_range = maximum(all_lats) - minimum(all_lats)
            lon_range = maximum(all_lons) - minimum(all_lons)
            lat_range_buffer = lat_range * 0.3
            lon_range_buffer = lon_range * 0.3
            lat_min = minimum(all_lats) - lat_range_buffer
            lat_max = maximum(all_lats) + lat_range_buffer
            lon_min = minimum(all_lons) - lon_range_buffer
            lon_max = maximum(all_lons) + lon_range_buffer
            
            #println("- Grid bounds with lattitude $(lat_range_buffer)° and longitude $(lon_range_buffer)° buffer:")
            #println("  Latitude: [$lat_min, $lat_max]")
            #println("  Longitude: [$lon_min, $lon_max]")
        end
    end

    local number_points = 50
    
    #println("- Creating spatial grid with resolution = $number_points × $number_points points")
    
    # Create the grid using lat/lon coordinates
    LAT_grid, LON_grid, lat_range, lon_range = create_spatial_grid_latlon(
        lat_min, lat_max, lon_min, lon_max, number_points
    )
    
    #println("- Grid created with size: $(size(LAT_grid))")
    @assert size(LAT_grid) == (50, 50) "Grid size mismatch. Expected (50,50), got $(size(LAT_grid))"

    
    #println("- Latitude range: [$(minimum(LAT_grid)), $(maximum(LAT_grid))]")
    #println("- Longitude range: [$(minimum(LON_grid)), $(maximum(LON_grid))]")
    
    # 2D matrix representing the pressure field
    local total_pressure_2d = zeros(size(LAT_grid))

    # Initialize variables for radial curves
    local r_km = range(0.5, stop=20.0, length=50)
    local r_m = r_km .* 1000
    local radial_info = Vector{Tuple{String, Vector{Float64}, Vector{Float64}}}()

    # Store well locations for plotting/reference
    well_locations = Dict{String, Tuple{Float64, Float64}}()

    # Initialize min_injection_year to the year of interest
    min_injection_year = year_of_interest

    # Extrapolation option for missing injection data
    extrapolate_injection_rates = get_parameter_value(helper, 4, "extrapolate_injection_rates")
    if extrapolate_injection_rates === nothing
        extrapolate_injection_rates = false
    end

    #println("\n------------------------------------------------------")
    #println("CALCULATING PRESSURE FIELD AND RADIAL CURVES FOR YEAR OF INTEREST: $year_of_interest")
    #println("------------------------------------------------------")

    well_count = 0
    for well_id in well_ids
        well_count += 1
        #println("\nWell $well_count of $(length(well_ids)) - ID: $well_id")
        
        # Get well coordinates based on data type
        if injection_data_type == "annual_fsp" || injection_data_type == "monthly_fsp"
            # filter the dataframe for the well id
            well_data = injection_wells_df[string.(injection_wells_df[!, "WellID"]) .== well_id, :]
            if isempty(well_data)
                @warn "No data found for well $well_id"
                continue
            end
            
            # Get well lat/lon coordinates for this well
            well_lat = first(well_data[!, "Latitude(WGS84)"])
            well_lon = first(well_data[!, "Longitude(WGS84)"])
            
            # Get injection period
            if "StartYear" in names(well_data) # Annual format
                inj_start_year = first(well_data[!, "StartYear"])
                inj_start_date = Date(inj_start_year, 1, 1)
                inj_end_year = first(well_data[!, "EndYear"])
                inj_end_date = Date(inj_end_year-1, 12, 31)
            else
                # Monthly FSP format has a 'Year' and 'Month'column
                if "Year" in names(well_data) || "Month" in names(well_data)
                    inj_start_year = minimum(well_data[!, "Year"])
                    inj_start_month = minimum(well_data[well_data[!, "Year"] .== inj_start_year, "Month"])
                    inj_start_date = Date(inj_start_year, inj_start_month, 1)
                    inj_end_year = maximum(well_data[!, "Year"])
                    inj_end_month = maximum(well_data[well_data[!, "Year"] .== inj_end_year, "Month"])
                    inj_end_date = Date(inj_end_year, inj_end_month, 1)
                    inj_end_date = lastdayofmonth(inj_end_date)
                else
                    error("Injection start and end years could not be determined because the 'Year' column is missing.")
                end
            end
            
        elseif injection_data_type == "injection_tool_data"
            # Handle injection tool data
            well_data = injection_wells_df[string.(injection_wells_df[!, "API Number"]) .== well_id, :]
            if isempty(well_data)
                # Try UIC Number
                well_data = injection_wells_df[string.(injection_wells_df[!, "UIC Number"]) .== well_id, :]
                if isempty(well_data)
                    @warn "No data found for well $well_id"
                    continue
                end
            end
            
            # Get well lat/lon coordinates
            well_lat = first(well_data[!, "Surface Latitude"])
            well_lon = first(well_data[!, "Surface Longitude"])
            
            # Get injection period from 'Date of Injection' column
            if "Date of Injection" in names(well_data)
                # Parse dates to determine year range
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
                                error("Error: Could not parse dates from the injection tool dataset: $e")
                            end
                        end
                    end
                end
                
                if !isempty(dates)
                    # Extract the years using the Dates package
                    years = year.(dates)
                    inj_start_year = minimum(years)
                    inj_end_year = maximum(years)
                    
                    # Get the start and end dates
                    inj_start_date = minimum(dates)
                    inj_end_date = min(maximum(dates), year_of_interest_date)
                    
                    #println("  * Injection period determined from dates: $inj_start_year to $inj_end_year")
                else
                    error("Error: No valid dates found in the injection tool dataset")
                    continue
                end
            else
                error("Could not find 'Date of Injection' column for well $well_id")
                continue
            end
        else
            error("Unsupported injection well dataset format: $injection_data_type")
        end
        
        # Store well location
        well_locations[well_id] = (well_lat, well_lon)
        
        # Update min_injection_year for fault pressure calculations
        min_injection_year = min(min_injection_year, inj_start_year)
        
        # Check if well hasn't started injecting yet at the year of interest
        if inj_start_date > year_of_interest_date
            #println("- Status: NOT ACTIVE before year $year_of_interest (active $inj_start_year-$inj_end_year)")
            continue
        end
        
        # Calculate end year for pressure calculation
        actual_end_year = min(inj_end_year, year_of_interest)
        
        # Prepare injection data for pressure front calculation
        well_specific_data = well_data  # Already filtered above

        # convert well_id from String3 to String
        well_id = String(well_id)
        #println("type of well_id: $(typeof(well_id))")
        
        #println("- Processing well data for pressure calculation (start=$inj_start_year, end=$actual_end_year)...")
        days, rates = prepare_well_data_for_pressure_scenario(
            well_specific_data,
            well_id,
            inj_start_year,
            inj_start_date,
            actual_end_year,
            inj_end_date,
            injection_data_type,
            year_of_interest,
            extrapolate_injection_rates,
            year_of_interest_date
        )
        
        if isempty(days) || isempty(rates)
            @warn "No valid injection data for well $well_id"
            continue
        end
        
        #println("- Injection history: $(length(days)) step changes over $(maximum(days)) days")
        
        # 1. Calculate pressure field contribution for this well
        # Calculate days from injection start to analysis date
        evaluation_days_from_start = Float64((year_of_interest_date - inj_start_date).value + 1)
        
        pfield_this_well = pfieldcalc_all_rates(
            LON_grid, LAT_grid, STRho, days, rates,
            well_lon, well_lat, "latlon",
            evaluation_days_from_start
        )

        
        
        # Add to total pressure field (superposition)
        total_pressure_2d .+= pfield_this_well
        #println("- Maximum pressure contribution: $(maximum(pfield_this_well)) psi")
        
        # 2. Calculate radial curve data for this well (pressure vs distance)
        pressure_psi = pressureScenario_Rall(rates, days, collect(r_m), STRho)
        
        # Add to radial_info collection
        push!(radial_info, (well_id, collect(r_km), collect(pressure_psi)))
        #println("- Generated radial curve (max pressure: $(maximum(pressure_psi)) psi)")
    end

    #println("\n------------------------------------------------------")
    #println("PRESSURE FIELD RESULTS")
    #println("------------------------------------------------------")
    
    # 1. Calculate pressure field statistics
    min_pressure = minimum(total_pressure_2d)
    max_pressure = maximum(total_pressure_2d)
    mean_pressure = mean(total_pressure_2d)
    median_pressure = median(total_pressure_2d)
    std_pressure = std(total_pressure_2d)
    
    # 2. Find coordinates of maximum pressure point
    max_idx = argmax(total_pressure_2d)
    # Get lat/lon coordinates at the maximum pressure point
    # IMPORTANT: Due to the flipped meshgrid call in create_spatial_grid_latlon,
    # LAT_grid actually contains longitude values and LON_grid contains latitude values
    max_pressure_lat = LON_grid[max_idx]  # This is the actual latitude
    max_pressure_lon = LAT_grid[max_idx]  # This is the actual longitude
    
    # 3. Create a summary DataFrame for pressure statistics
    pressure_stats = DataFrame(
        Metric = ["Min Pressure (psi)", "Max Pressure (psi)", "Mean Pressure (psi)", 
                  "Median Pressure (psi)", "Standard Deviation (psi)",
                  "Max Pressure Latitude", "Max Pressure Longitude"],
        Value = [min_pressure, max_pressure, mean_pressure, 
                 median_pressure, std_pressure,
                 max_pressure_lat, max_pressure_lon]
    )
    
    

    

    

    # Load fault data to calculate pressure on faults
    #println("\n------------------------------------------------------")
    #println("CALCULATING PRESSURE ON FAULTS FOR YEAR $year_of_interest")
    #println("------------------------------------------------------")
    
    
    num_faults = nrow(fault_df)
    
    
    # Dataframe to store fault pressure
    fault_pressure_by_year = DataFrame(
        FaultID = String[],
        Date = Date[],
        slip_pressure = Float64[], # this is the pressure added to the fault, not the pore pressure required to slip
        probability = Float64[],
        Year = Int[]
    )
    
    # Create pressure array for faults at year_of_interest
    pressure_on_faults = zeros(num_faults)

    # Calculate pressure on each fault
    for f in 1:num_faults
        # Get fault coordinates
        fault_lat = fault_df[f, "Latitude(WGS84)"]
        fault_lon = fault_df[f, "Longitude(WGS84)"]
        
        # Get the actual fault ID from the fault dataset
        fault_id = if "FaultID" in names(fault_df)
            string(fault_df[f, "FaultID"])
        else
            # Fall back to using the index as the ID
            string(f)
        end
        
        #println("    - Fault $fault_id: lat=$(fault_lat)°, lon=$(fault_lon)°")
        
        # Process each well's contribution
        for well_id in well_ids
            # Skip wells with no location data
            if !haskey(well_locations, well_id)
                continue
            end

            
            
            # Get well coordinates
            well_lat, well_lon = well_locations[well_id]
            
            # Filter well-specific data
            local well_specific_data
            if injection_data_type == "annual_fsp" || injection_data_type == "monthly_fsp"
                well_specific_data = injection_wells_df[string.(injection_wells_df[!, "WellID"]) .== well_id, :]
            elseif injection_data_type == "injection_tool_data"
                if "API Number" in names(injection_wells_df)
                    well_specific_data = injection_wells_df[string.(injection_wells_df[!, "API Number"]) .== well_id, :]
                elseif "UIC Number" in names(injection_wells_df)
                    well_specific_data = injection_wells_df[string.(injection_wells_df[!, "UIC Number"]) .== well_id, :]
                else
                    continue
                end
            end
            
            if isempty(well_specific_data)
                continue
            end
            
            # Get injection period based on data type
            local inj_start_year, inj_end_year, inj_start_date, inj_end_date
            
            if injection_data_type == "annual_fsp"
                # Get injection period
                inj_start_year = minimum(well_specific_data[!, "StartYear"])
                inj_start_date = Date(inj_start_year, 1, 1)
                inj_end_year = maximum(well_specific_data[!, "EndYear"])
                inj_end_date = Date(inj_end_year-1, 12, 31)
            elseif injection_data_type == "monthly_fsp"
                # Get injection period
                inj_start_year = minimum(well_specific_data[!, "Year"])
                inj_start_month = minimum(well_specific_data[well_specific_data[!, "Year"] .== inj_start_year, "Month"])
                inj_start_date = Date(inj_start_year, inj_start_month, 1)
                inj_end_year = maximum(well_specific_data[!, "Year"])
                inj_end_month = maximum(well_specific_data[well_specific_data[!, "Year"] .== inj_end_year, "Month"])
                inj_end_date = Date(inj_end_year, inj_end_month, 1)
                inj_end_date = lastdayofmonth(inj_end_date)
            elseif injection_data_type == "injection_tool_data"
                # Get injection period from Date of Injection column
                if "Date of Injection" in names(well_specific_data)
                    # Parse dates to determine year range
                    dates = Date[]
                    
                    # Check the type of the date values first
                    if eltype(well_specific_data[!, "Date of Injection"]) <: Date
                        dates = well_specific_data[!, "Date of Injection"]
                    else
                        # Need to parse from strings
                        try
                            # Try different date formats
                            dates = Date.(well_specific_data[!, "Date of Injection"], dateformat"y-m-d")
                        catch
                            try
                                dates = Date.(well_specific_data[!, "Date of Injection"], dateformat"m/d/y")
                            catch
                                try
                                    dates = Date.(well_specific_data[!, "Date of Injection"], dateformat"m/d/yyyy")
                                catch e
                                    continue
                                end
                            end
                        end
                    end
                    
                    if !isempty(dates)
                        # Extract the years using the Dates package
                        years = year.(dates)
                        inj_start_year = minimum(years)
                        inj_end_year = maximum(years)
                        
                        # Get the start and end dates
                        inj_start_date = minimum(dates)
                        inj_end_date = min(maximum(dates), year_of_interest_date)
                    else
                        continue
                    end
                else
                    continue
                end
            else
                error("Unsupported data type: $injection_data_type")
            end
            
            # Check if well is active before the year of interest
            if inj_start_date > year_of_interest_date
                continue
            end
            
            # Calculate end year (limit to year_of_interest)
            actual_end_year = min(inj_end_year, year_of_interest)
            
            # Prepare injection data up to year_of_interest
            days, rates = prepare_well_data_for_pressure_scenario(
                well_specific_data,
                String(well_id),
                inj_start_year,
                inj_start_date,
                actual_end_year,
                inj_end_date,
                injection_data_type,
                year_of_interest,
                extrapolate_injection_rates,
                year_of_interest_date
            )
            
            if isempty(days) || isempty(rates)
                continue
            end
            
            # Calculate pressure contribution from this well at the fault location
            # Calculate days from injection start to analysis date
            evaluation_days_from_start = Float64((year_of_interest_date - inj_start_date).value + 1)
            
            pressure_contribution = pfieldcalc_all_rates(
                fault_lon,  # longitude is x
                fault_lat,  # latitude is y
                STRho,
                days,
                rates,
                well_lon,  # longitude is x
                well_lat,  # latitude is y
                evaluation_days_from_start
            )
            
            # Add to total pressure for this fault
            pressure_on_faults[f] += pressure_contribution
            
        end
        
        # Add results to the DataFrame
        push!(fault_pressure_by_year, (
            fault_id,
            year_of_interest_date,
            pressure_on_faults[f],
            0.0,
            year_of_interest
        ))

        # Add an extra row (so we can make it a line segment later)
        push!(fault_pressure_by_year, (
            fault_id,
            year_of_interest_date,
            pressure_on_faults[f],
            1.0,
            year_of_interest
        ))
        
       
    end
    
    
    
    
    # Save the fault pressure DataFrame as a parameter
    save_dataframe_as_parameter!(helper, 4, "deterministic_hydrology_results", fault_pressure_by_year)
    #println("- Saved fault pressure data as parameter 'deterministic_hydrology_results'")

    
    
    # 5. Create a grid statistics DataFrame
    grid_info = DataFrame(
        Parameter = ["Latitude Min", "Latitude Max", "Longitude Min", "Longitude Max",
                    "Grid Resolution", "Grid Points", "Year of Interest",
                    "Storativity", "Transmissivity"],
        Value = [lat_min, lat_max, lon_min, lon_max,
                number_points, number_points*number_points, year_of_interest,
                S, T]
    )
    
    
    
    
    
    
    
    # Create a DataFrame with lat, lon, and pressure values
    # this is the full grid that we reformat to the heatmap data
    grid_size = size(LAT_grid)
    n_points = prod(grid_size)
    
    full_grid_df = DataFrame(
        Latitude = vec(LON_grid),    # LON_grid contains latitude values
        Longitude = vec(LAT_grid),   # LAT_grid contains longitude values
        Pressure_psi = vec(total_pressure_2d)
    )
    
    # Save to CSV
    #=
    CSV.write(joinpath(results_dir, "full_grid.csv"), full_grid_df)
    #println("- Full grid saved to $(joinpath(results_dir, "full_grid.csv"))")
    =#
    


    # reformat the full grid to heatmap data
    #println("Reformatting full grid to heatmap data...")
    heatmap_data = reformat_pressure_grid_to_heatmap_data(full_grid_df, lat_range, lon_range)
    save_dataframe_as_parameter!(helper, 4, "hydrology_heatmap_data_arcgis", heatmap_data)
    #pretty_table(first(heatmap_data, 10))
    # print the dimensions of the heatmap data
    #println("Heatmap data dimensions: $(size(heatmap_data))")






    #println("\n------------------------------------------------------")
    #println("PROCESSING RADIAL CURVE DATA...")
    #println("------------------------------------------------------")

    if !isempty(radial_info)
        # Create a single DataFrame to hold all well data
        radial_df = DataFrame(
            Distance_km = Float64[],
            Pressure_psi = Float64[],
            ID = String[]
        )
        
        # Process each well's data
        for (well_id, distances, pressures) in radial_info
            # Add data for this well to the combined DataFrame
            for (dist, pres) in zip(distances, pressures)
                push!(radial_df, (dist, pres, string(well_id)))
            end
            
            #println("- Added radial data for well $well_id (max pressure: $(maximum(pressures)) psi)")
        end
        
       
        
        # Save to portal parameter
        save_dataframe_as_parameter!(helper, 4, "radial_curves_data", radial_df)
        #println("radial curves data:")
        #pretty_table(radial_df)
        #println("- Radial curve data saved to dataset")
    else
        #println("- No wells with valid data for radial curves")
    end


    
    #println("\n------------------------------------------------------")
    #println("PREPARING MOHR DIAGRAM DATA WITH UPDATED FAULT PRESSURES")
    #println("------------------------------------------------------")
    
    # Load geomechanics results with original stress state
    #=
    geo_results_path = get_dataset_file_path(helper, 4, "det_geomechanics_results")
    if geo_results_path === nothing
        error("Required deterministic geomechanics results not found or accessible.")
    end
    =#
    
    #geo_results_df = CSV.read(geo_results_path, DataFrame)
    #println("- Loaded geomechanics results with $(nrow(geo_results_df)) faults")
    
    # Extract stress state parameters from args.json or previous step
    #println("- Extracting stress state parameters...")
    stress_inputs = Dict(
        "reference_depth" => get_parameter_value(helper, 2, "reference_depth"),
        "vertical_stress" => get_parameter_value(helper, 2, "vertical_stress"),
        "min_horizontal_stress" => get_parameter_value(helper, 2, "min_horizontal_stress"),
        "max_horizontal_stress" => get_parameter_value(helper, 2, "max_horizontal_stress"),
        "pore_pressure" => get_parameter_value(helper, 2, "pore_pressure"),
        "max_stress_azimuth" => get_parameter_value(helper, 2, "max_stress_azimuth"),
        "model_type" => get_parameter_value(helper, 2, "stress_field_mode"),
        "aphi_value" => get_parameter_value(helper, 2, "aphi_value") === nothing ? nothing : get_parameter_value(helper, 2, "aphi_value"),
        "friction_coefficient" => get_parameter_value(helper, 2, "friction_coefficient")
    )

    #=
    # REMOVE THIS
    if stress_inputs["max_horizontal_stress"] === nothing
        add_message_with_step_index!(helper, 2, "Max Horizontal Stress Gradient is not provided, using default value of 1.22", 2)
        stress_inputs["max_horizontal_stress"] = 1.22
    end
    =#

    
    
    # Use friction coefficient from fault data
    # TO DO: we'll make the portal accepts this as a single scalar float value and not as part of the faults dataframe
    friction_coefficient = stress_inputs["friction_coefficient"]
    

    
    
    # Calculate stress state for the Mohr diagram
    
    
    # Calculate absolute stresses at reference depth
    stress_state, initial_pressure = GeomechanicsModel.calculate_absolute_stresses(
        stress_inputs, friction_coefficient, stress_inputs["model_type"]
    )
    
    # Get pressure changes for the year of interest from fault_pressure_by_year
    # we filter for probability = 0.0 because in the 'fault_pressure_by_year' dataframe,
    # for each fault and year we have two rows: one with probability = 0.0 and one with probability = 1.0
    # they both have the same data (pressure), but we only need one of them
    year_specific_data = fault_pressure_by_year[
        (fault_pressure_by_year.Date .== year_of_interest_date) .& 
        (fault_pressure_by_year.probability .== 0.0),
        :
    ]

    #println("pressure added to faults:")
    #pretty_table(year_specific_data)
    
    # Extract pore pressure changes for each fault
    pressure_changes_vec = zeros(nrow(fault_df))
    pressure_changes_df = DataFrame(FaultID = String[], Pressure = Float64[])

    

    # check if the fault_df 'FaultID' column is a string, if not, convert it to a string
    if !(eltype(fault_df[!, "FaultID"]) <: AbstractString)
        fault_df[!, "FaultID"] = string.(fault_df[!, "FaultID"])
    end

    


    for i in 1:nrow(fault_df)
        # Get the actual fault ID from the fault dataset
        fault_id = if "FaultID" in names(fault_df)
            fault_df[i, "FaultID"]
        else
            # Fall back to using the index as the ID
            string(i)
            #println("FALLBACK OPTION USED")
        end

        

        # TO DO: I removed the fallback option (make sure it's not needed)
        # also added the check if the fault_id is a string, if not, convert it to a string
        #=
        # check if the fault_id is a string, if not, convert it to a string
        if typeof(fault_id) != String
            fault_id = string(fault_id)
        end
        =#
        
        fault_rows = year_specific_data[year_specific_data.FaultID .== fault_id, :]
        if !isempty(fault_rows)
            pressure_value = first(fault_rows.slip_pressure)
            pressure_changes_vec[i] = pressure_value
            push!(pressure_changes_df, (FaultID = fault_id, Pressure = pressure_value))
        else
            push!(pressure_changes_df, (FaultID = fault_id, Pressure = 0.0))
        end
    end

    

    # We still need to keep pressure_changes_vec for compatibility with existing code
    pressure_changes = pressure_changes_vec

    # Process faults with updated pressure
    #println("- Processing faults with updated pressure...")
    faults_with_pressure = Vector{Dict{String, Any}}()
    for i in 1:nrow(fault_df)
        # Get the actual fault ID from the fault dataset
        fault_id = if "FaultID" in names(fault_df)
            string(fault_df[i, "FaultID"])
        else
            #println("No fault ID found for fault, using index as ID: $i")
            # Fall back to using the index as the ID
            string(i)
        end
        
        push!(faults_with_pressure, Dict{String, Any}(
            "strike" => fault_df[i, "Strike"],
            "dip" => fault_df[i, "Dip"],
            "friction_coefficient" => stress_inputs["friction_coefficient"],
            "fault_id" => fault_id
        ))
    end

    #println("faultdf inputs:")
    #pretty_table(fault_df)

    #println("- Pressure changes for year $year_of_interest:")
    #pretty_table(pressure_changes)

    #println("faults_with_pressure:")
    #pretty_table(faults_with_pressure)
    
    hydro_results = GeomechanicsDriver.process_faults(
        faults_with_pressure, 
        GeomechanicsDriver.GeomechanicsModel.StressState(stress_state.principal_stresses, stress_state.sH_azimuth), 
        initial_pressure,
        stress_inputs["friction_coefficient"]; 
        tab="det_hydro", 
        dp=pressure_changes
    )
    
    # Extract data for Mohr diagram
    tau_effective_faults = [result["shear_stress"] for result in hydro_results]
    sigma_n_faults = [result["normal_stress"] for result in hydro_results]
    slip_pressures = [result["slip_pressure"] for result in hydro_results]

    # Use actual fault IDs instead of indices
    fault_ids = [fault["fault_id"] for fault in faults_with_pressure]
    strikes = [fault["strike"] for fault in faults_with_pressure]
    
    # Determine stress regime
    # Normal: σV > σH > σh
    stress_regime = if stress_state.principal_stresses[1] >= stress_state.principal_stresses[3] &&
                      stress_state.principal_stresses[3] >= stress_state.principal_stresses[2]
        "Normal"
    # Reverse: σH > σh > σV
    elseif stress_state.principal_stresses[3] >= stress_state.principal_stresses[2] &&
           stress_state.principal_stresses[2] >= stress_state.principal_stresses[1]
        "Reverse"
    # Strike-Slip: σH > σV > σh
    elseif stress_state.principal_stresses[3] >= stress_state.principal_stresses[1] &&
           stress_state.principal_stresses[1] >= stress_state.principal_stresses[2]
        "Strike-Slip"
    else
        "Unknown"
    end
    

    # get arcs, slip, and faults dataframes from the geomechanics results
    geo_faults_df = get_dataset_file_path(helper, 4, "faultDF")
    geo_slip_df = get_dataset_file_path(helper, 4, "slipDF")
    geo_arcs_df = get_dataset_file_path(helper, 4, "arcsDF")

    geo_faults_df = CSV.read(geo_faults_df, DataFrame)
    geo_slip_df = CSV.read(geo_slip_df, DataFrame)
    geo_arcs_df = CSV.read(geo_arcs_df, DataFrame)

    


    fault_inputs_filepath = get_dataset_file_path(helper, 4, "faults")
    fault_inputs_df = CSV.read(fault_inputs_filepath, DataFrame)

    #println("fault_inputs_df: $(fault_inputs_df)")

    if typeof(fault_inputs_df.FaultID[1]) == String7
        fault_inputs_df.FaultID = String.(fault_inputs_df.FaultID)
    elseif typeof(fault_inputs_df.FaultID[1]) == InlineStrings.String1
        fault_inputs_df.FaultID = String.(fault_inputs_df.FaultID)
    elseif typeof(fault_inputs_df.FaultID[1]) == InlineStrings.String15
        fault_inputs_df.FaultID = String.(fault_inputs_df.FaultID)
    end



    # Get data for Mohr diagram with pressure changes
    arcsDF, slipDF, faultsDF = JuliaFSPGraphs.mohr_diagram_hydro_data_to_d3_portal(
        stress_state.principal_stresses[2], 
        stress_state.principal_stresses[3], 
        stress_state.principal_stresses[1], 
        tau_effective_faults, 
        sigma_n_faults, 
        initial_pressure, 
        1.0, 
        0.5, 
        pressure_changes, 
        strikes, 
        friction_coefficient, 
        stress_regime, 
        slip_pressures, 
        String.(fault_ids),
        geo_arcs_df,
        geo_faults_df,
        geo_slip_df,
        fault_inputs_df
    )

    # Add the slip pressures to fault_inputs_df
    fault_inputs_df[!, "pore_pressure_slip_det_hydro"] = zeros(nrow(fault_inputs_df))
    for (i, fault_id) in enumerate(fault_ids)
        idx = findfirst(id -> string(id) == fault_id, fault_inputs_df.FaultID)
        if !isnothing(idx)
            round(slip_pressures[i], digits=2)
            fault_inputs_df[idx, "pore_pressure_slip_det_hydro"] = slip_pressures[i]
        end
    end

    # Print the updated fault_inputs_df with new slip pressures
    #println("Updated fault_inputs_df with new slip pressures:")
    #pretty_table(fault_inputs_df, show_omitted_cell_summary=false, crop=:none)

    # Save the updated fault dataframe
    save_dataframe_as_parameter!(helper, 4, "faults_with_det_hydro_pp", fault_inputs_df)

    # Save datasets for visualization
    save_dataframe_as_parameter!(helper, 4, "arcsDF_hydro", arcsDF)
    save_dataframe_as_parameter!(helper, 4, "slipDF_hydro", slipDF)
    save_dataframe_as_parameter!(helper, 4, "faultsDF_hydro", faultsDF)
    
    

    # explicitly set this step's success state to true
    set_success_for_step_index!(helper, 4, true)

    # Save to results.json for the portal
    write_results_file(helper)
    
    #println("\nResults saved to CSV files in $(joinpath(@__DIR__, "output"))")
    #println("\n======================================================")
    #println("      DETERMINISTIC HYDROLOGY PROCESS COMPLETED        ")
    #println("======================================================\n")
end








if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end






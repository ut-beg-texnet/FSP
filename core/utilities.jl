module Utilities

include("../TexNetWebToolLauncherHelperJulia.jl")

using DataFrames
using Geodesy
using LinearAlgebra
using Dates
using MeshGrid
using Proj
using Distributions
using LibGEOS 
using CSV
using Interpolations
using .TexNetWebToolLauncherHelperJulia



import Proj: CRS # explicitly 

export latlon_to_wkt, convert_easting_northing_to_latlon, convert_latlon_to_easting_northing!
export prepare_well_data_for_pressure_scenario, create_spatial_grid_km, create_spatial_grid_latlon, create_uniform_distribution
export reformat_pressure_grid_to_heatmap_data, get_date_bounds, get_injection_dataset_path, interpolate_cdf, create_bounded_uniform_distribution



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


# Converts lat/lon to WKT format for polylines on the map
function latlon_to_wkt(faults_df::DataFrame,
    lat_column::String="Latitude(WGS84)",
    lon_column::String="Longitude(WGS84)",
    strike_column::String="Strike",
    length_column::String="LengthKm",
    id_column::String="FaultID"
)
    # Create a new column to add the WKT representation
    wkt_strings = String[]

    for row in eachrow(faults_df)
        lat = row[lat_column]
        lon = row[lon_column]
        strike = row[strike_column]
        length_km = row[length_column]
        id = row[id_column]

        # Create a LLA point (for WGS84) for the center of the fault
        center_point = LLA(lat, lon, 0.0)
        
        # Convert strike to radians (strike is measured clockwise from north)
        strike_rad = deg2rad(90 - strike) # Convert from azimuth to math angle
        
        # Calculate half length for extending in both directions
        half_length_km = length_km / 2.0
        
        # Create a local ENU coordinate system centered at the fault midpoint
        # This gives us a more accurate way to extend the fault in the proper direction
        enu_transform = ENUfromLLA(center_point, wgs84)
        lla_transform = LLAfromENU(center_point, wgs84)
        
        # Calculate the endpoints in the local ENU frame
        # In ENU, x is East, y is North, so we need to use the strike angle correctly
        dx = sin(strike_rad) * half_length_km * 1000.0  # Convert km to meters
        dy = cos(strike_rad) * half_length_km * 1000.0  # Convert km to meters
        
        # Start and end points in ENU coordinates (meters)
        start_point_enu = ENU(-dx, -dy, 0.0)  # Negative direction
        end_point_enu = ENU(dx, dy, 0.0)      # Positive direction
        
        # Convert back to LLA coordinates
        start_point_lla = lla_transform(start_point_enu)
        end_point_lla = lla_transform(end_point_enu)
        
        # Create standard WKT representation
        wkt_string = "LINESTRING ($(start_point_lla.lon) $(start_point_lla.lat),$(end_point_lla.lon) $(end_point_lla.lat))"
        push!(wkt_strings, wkt_string)
    end

    # Add the WKT column directly to the input DataFrame
    faults_df.wkt = wkt_strings
    return faults_df
end


function convert_easting_northing_to_latlon(easting::Float64, northing::Float64)
    println("lat: $(wgs84(easting, northing).lat), lon: $(wgs84(easting, northing).lon)")
    return wgs84(easting, northing)
end

# converts lat/lon in wgs-84 to easting northing for any dataframe with lat/lon columns
function convert_latlon_to_easting_northing!(df::DataFrame, lat_col::String, lon_col::String)
    try
        # Validate DataFrame
        if isnothing(df) || nrow(df) == 0
            throw(ArgumentError("Input DataFrame is empty or invalid"))
        end
        
        df_cols = try
            names(df)
        catch e
            throw(ArgumentError("Failed to get column names: $e"))
        end

        # Case-insensitive column name matching
        df_cols_lower = lowercase.(df_cols)
        
        # Find column indices using case-insensitive matching
        lat_column = try
            df_cols[findfirst(==(lowercase(lat_col)), df_cols_lower)]
        catch e
            throw(ArgumentError("Failed to find latitude column: $e"))
        end
        
        lon_column = try
            df_cols[findfirst(==(lowercase(lon_col)), df_cols_lower)]
        catch e
            throw(ArgumentError("Failed to find longitude column: $e"))
        end

        # Validate coordinate ranges
        try
            if any(lat -> !(-90 <= lat <= 90), df[:, lat_column])
                throw(ArgumentError("Latitude values must be in range [-90, 90] degrees"))
            end
            if any(lon -> !(-180 <= lon <= 180), df[:, lon_column])
                throw(ArgumentError("Longitude values must be in range [-180, 180] degrees"))
            end
        catch e
            throw(ArgumentError("Coordinate validation failed: $e"))
        end

        n_rows = nrow(df)
        eastings = Vector{Float64}(undef, n_rows)
        northings = Vector{Float64}(undef, n_rows)
        zones = Vector{String}(undef, n_rows)

        # Process each row
        for i in 1:n_rows
            try
                lat = df[i, lat_column]
                lon = df[i, lon_column]
                
                if isnothing(lat) || isnothing(lon) || ismissing(lat) || ismissing(lon)
                    throw(ArgumentError("Invalid coordinates in row $i"))
                end

                point_lla = try
                    LLA(lat, lon, 0.0)
                catch e
                    throw(ArgumentError("Failed to create LLA object: $e"))
                end

                zone = floor(Int, (lon + 180) / 6) + 1
                hemisphere = lat >= 0

                utm_transform = try
                    UTMfromLLA(zone, hemisphere, wgs84)
                catch e
                    throw(ArgumentError("Failed to create UTM transformation: $e"))
                end

                point_utm = try
                    utm_transform(point_lla)
                catch e
                    throw(ArgumentError("Failed to transform coordinates: $e"))
                end

                eastings[i] = point_utm.x / 1000.0  # Convert to km
                northings[i] = point_utm.y / 1000.0 # Convert to km
                zones[i] = "$(zone)$(hemisphere ? 'N' : 'S')"
                
            catch e
                println("\nError processing row $i:")
                println(sprint(showerror, e))
                throw(ArgumentError("Error converting coordinates for row $i: $e"))
            end
        end

        try
            df[!, "EastingKm"] = eastings
            df[!, "NorthingKm"] = northings
            df[!, "UTM_Zone"] = zones
        catch e
            throw(ArgumentError("Failed to add new columns to DataFrame: $e"))
        end
        
        return df
        
    catch e
        println("\nFatal error in coordinate conversion:")
        println(sprint(showerror, e))
        rethrow(e)
    end
end


# converts a dataframe with lat/lon columns to a dataframe with easting/northing columns
# uses FSP's assumption that we have cartesian coordinates in km
function convert_latlon_to_cartesian_km!(df::DataFrame, lat_col::String, lon_col::String)
    # set origin to be the first row of the dataframe (we can later make this a user input)
    origin_lat = df[1, lat_col]
    origin_lon = df[1, lon_col]

    # earth radius in km
    R = 6371.0008

    n_rows = nrow(df)

    eastings = Vector{Float64}(undef, n_rows)
    northings = Vector{Float64}(undef, n_rows)

    for i in 1:n_rows
        lat = df[i, lat_col]
        lon = df[i, lon_col]

        #get differences in degrees from origin and convert to radians
        dlat = deg2rad(lat - origin_lat)
        dlon = deg2rad(lon - origin_lon)

        #convert to cartesian coordinates (in km)
        x = R * dlon * cos(deg2rad(origin_lat))
        y = R * dlat

        eastings[i] = x
        northings[i] = y
    end

    # add the new columns to the dataframe
    df[!, "EastingKm"] = eastings
    df[!, "NorthingKm"] = northings

    return df
end
    
    

    
    


















"""
    prepare_well_data_for_pressure_scenario(df::DataFrame, well_id::String, start_year::Int, end_year::Int, data_type::String, year_of_interest::Int, extrapolate::Bool=false)

Prepares well data (days and injection rates) for pressure scenario calculations based on different dataframe types:
- annual_fsp: Annual injection data from FSP format
- monthly_fsp: Monthly injection data from FSP format
- injection_tool_data: Data from injection tool

Returns a tuple of (days, rates) where:
- days: Array of days relative to start_date (1-indexed)
- rates: Array of injection rates in barrels per day

Parameters:
- df: DataFrame containing the well data
- well_id: Identifier of the well
- start_year: Starting year for the injection period
- inj_start_date: Start date of the injection period for this well
- end_year: Ending year for the injection period
- inj_end_date: End date of the injection period for this well
- injection_data_type: Type of the data ("annual_fsp", "monthly_fsp", or "injection_tool_data")
- year_of_interest: The year for which to calculate pressure
- extrapolate: Whether to extrapolate missing data (default: false)
"""
function prepare_well_data_for_pressure_scenario(
    df::DataFrame, 
    well_id::String, 
    start_year::Int, 
    inj_start_date::Date,
    end_year::Int, 
    inj_end_date::Date,
    injection_data_type::String,
    year_of_interest::Int,
    extrapolate::Bool=false,
    year_of_interest_date::Date=Date(year_of_interest-1, 12, 31)
)
    # Validate inputs - start_year should not exceed end_year
    @assert start_year <= end_year "Start year must be <= end year"
    # Note: end_year is already capped at year_of_interest in the main script
    
    # Check if data is already filtered for the specific well
    # the way we currently call this function passes a dataframe that is already filtered for the specific well
    already_filtered = false
    if "APINumber" in names(df) && all(df.APINumber .== well_id)
        already_filtered = true
        #println("DEBUG: Data is already filtered for well $well_id")
        well_data = df
    elseif "API Number" in names(df) && all(df[!, "API Number"] .== well_id)
        already_filtered = true
        #println("DEBUG: Data is already filtered for well $well_id")
        well_data = df
    elseif "WellID" in names(df) && all(df.WellID .== well_id)
        already_filtered = true
        #println("DEBUG: Data is already filtered for well $well_id")
        well_data = df
    else
        # Filter data for the specific well - handle different possible well ID column names
        well_data = DataFrame()
        if "WellID" in names(df) && (injection_data_type == "annual_fsp" || injection_data_type == "monthly_fsp")
            well_data = df[string.(df[!, "WellID"]) .== well_id, :]
            #println("DEBUG: Filtered by WellID, got $(nrow(well_data)) rows")
        elseif "APINumber" in names(df)
            well_data = df[string.(df[!, "APINumber"]) .== well_id, :]
            #println("DEBUG: Filtered by APINumber, got $(nrow(well_data)) rows")
        elseif "API Number" in names(df)
            well_data = df[string.(df[!, "API Number"]) .== well_id, :]
            #println("DEBUG: Filtered by API Number, got $(nrow(well_data)) rows")
        elseif "Well ID" in names(df)
            well_data = df[string.(df[!, "Well ID"]) .== well_id, :]
            #println("DEBUG: Filtered by Well ID, got $(nrow(well_data)) rows")
        elseif "Well_ID" in names(df)
            well_data = df[string.(df[!, "Well_ID"]) .== well_id, :]
            #println("DEBUG: Filtered by Well_ID, got $(nrow(well_data)) rows")
        elseif "ID" in names(df)
            well_data = df[string.(df[!, "ID"]) .== well_id, :]
            #println("DEBUG: Filtered by ID, got $(nrow(well_data)) rows")
        else
            #println("DEBUG: Could not find well ID column in data. Available columns: $(join(names(df), ", "))")
            error("Could not find well ID column in data")
        end
    end
    
    if isempty(well_data)
        @warn "No data found for well ID: $well_id"
        return Float64[], Float64[]
    end
    
    # Select the appropriate preparation method based on data type
    if injection_data_type == "annual_fsp"
        #println("DEBUG: Calling prepare_annual_fsp_data")
        return prepare_annual_fsp_data(well_data, start_year, inj_start_date, end_year, inj_end_date, year_of_interest, year_of_interest_date)
    elseif injection_data_type == "monthly_fsp"
        #println("DEBUG: Calling prepare_monthly_fsp_data")
        return prepare_monthly_fsp_data(well_data, start_year, inj_start_date, end_year, inj_end_date, year_of_interest, extrapolate, year_of_interest_date)
    elseif injection_data_type == "injection_tool_data"
        #println("DEBUG: Calling prepare_injection_tool_data")
        return prepare_injection_tool_data(well_data, start_year, inj_start_date, end_year, inj_end_date, year_of_interest, extrapolate, year_of_interest_date)
    else
        error("Unsupported well dataset data type: $injection_data_type")
    end
end

"""
    prepare_annual_fsp_data(well_data::DataFrame, start_year::Int, end_year::Int, year_of_interest::Int)

Prepares injection data from annual FSP format by using a constant rate throughout the entire injection period.
This assumes the user provides a single injection rate for the entire period from start_year to end_year.
"""
function prepare_annual_fsp_data(
    well_data::DataFrame, 
    start_year::Int, 
    inj_start_date::Date,
    end_year::Int, 
    inj_end_date::Date,
    year_of_interest::Int,
    year_of_interest_date::Date
)
    # Use the provided end_year (already calculated as min(inj_end_year, year_of_interest))
    # If well starts after end year, return empty arrays
    if start_year > year(year_of_interest_date)
        println("DEBUG: Well starts after calculation end year - start_year ($start_year) > end_year ($end_year)")
        return Float64[], Float64[]
    end
    
    # Get the constant injection rate for the entire period
    if "InjectionRate(bbl/day)" in names(well_data)
        daily_rate_bbl = first(well_data[!, "InjectionRate(bbl/day)"])
    else
        error("Could not find injection rate in annual_fsp format. The 'InjectionRate(bbl/day)' column is required.")
    end
    
    @assert daily_rate_bbl >= 0 "Injection rate must be >= 0"
    
    # Define the time period
    
    global_start_date = inj_start_date
    global_end_date = min(inj_end_date, year_of_interest_date)
    
    
    
    # Calculate total days in the injection period
    days_total = (global_end_date - global_start_date).value + 1
    
    # Create arrays with just two points - start and end
    """
    1) step_times is an array containing only two elements:
        - 1.0: The first day of injection (day 1)
        - days_total: The total number of days in the injection period (converted to Float64)
    2) step_rates is a corresponding array that defines the injection rate at each time step:
        - daily_rate_bbl: The constant injection rate in barrels per day at the beginning
        - 0.0: The rate at the end (shutdown of injection)
    This is an optimization for constant-rate injection. Instead of creating an array with a value for every day of injection (which could be thousands of values), it uses just two points to represent:
        - A step up from 0 to daily_rate_bbl on day 1
        - A step down to 0 on day days_total
    """
    step_times = Float64[1.0, Float64(days_total)]
    step_rates = Float64[daily_rate_bbl, 0.0]
    
    println("  * Constant injection rate: $daily_rate_bbl bbl/day for $(days_total-1) days")
    println("  * Injection period: $(global_start_date) to $(global_end_date) ($(days_total) days)")
    
    return step_times, step_rates
end

"""
    prepare_monthly_fsp_data(well_data::DataFrame, start_year::Int, end_year::Int, year_of_interest::Int, extrapolate::Bool)

Prepares monthly injection data from FSP format with step changes at month boundaries.
"""
function prepare_monthly_fsp_data(
    well_data::DataFrame, 
    start_year::Int, 
    inj_start_date::Date,
    end_year::Int, 
    inj_end_date::Date,
    year_of_interest::Int, 
    extrapolate::Bool,
    year_of_interest_date::Date
)
    #println("\n===== DEBUG: Inside prepare_monthly_fsp_data =====")
    
    #println("DEBUG: start_year=$start_year, end_year=$end_year, year_of_interest=$year_of_interest")
    
    # If well starts after end_year, return empty arrays
    if start_year > year(year_of_interest_date)
        println("DEBUG: Well starts after calculation end year - start_year ($start_year) > end_year ($end_year)")
        return Float64[], Float64[]
    end
    
    # Define boundaries
    global_start_date = inj_start_date
    global_end_date = min(inj_end_date, year_of_interest_date)
    
    
    # Arrays for step changes
    step_times = Float64[]
    step_rates = Float64[]
    current_rate = 0.0
    
    # Helper function to get the first day of the next month
    function next_month(d::Date)
        y, m = year(d), month(d)
        if m < 12
            return Date(y, m+1, 1)
        else
            return Date(y+1, 1, 1)
        end
    end
    
    
    # Verify required columns exist
    required_cols = ["Year", "Month", "InjectionRate(bbl/month)"]
    missing_cols = filter(col -> !(col in names(well_data)), required_cols)
    
    if !isempty(missing_cols)
        println("DEBUG: Missing required columns: $(join(missing_cols, ", "))")
        
        # Try alternative column names for "InjectionRate(bbl/month)"
        if "InjectionRate(bbl/month)" in missing_cols
            alt_rate_cols = ["Injection Rate (bbl/month)", "MonthlyInjectionRate", "MonthlyVolume"]
            found_alt = false
            
            for alt_col in alt_rate_cols
                if alt_col in names(well_data)
                    println("DEBUG: Using alternative column '$alt_col' for injection rate")
                    rename!(well_data, alt_col => "InjectionRate(bbl/month)")
                    found_alt = true
                    missing_cols = filter(col -> !(col in names(well_data)), required_cols)
                    break
                end
            end
            
            if !found_alt
                println("DEBUG: No suitable monthly injection rate column found")
                return Float64[], Float64[]
            end
        end
        
        # If still missing required columns, return empty
        if !isempty(missing_cols)
            println("DEBUG: Still missing required columns after attempting alternatives: $(join(missing_cols, ", "))")
            return Float64[], Float64[]
        end
    end
    
    
    
    # Print the ranges we're looking for
    years_to_find = start_year:end_year-1
    months_to_find = 1:12
    
    
    # Print out all well data rows for debugging
    #=
    println("DEBUG: Well data rows (first 5 rows max):")
    for (i, row) in enumerate(eachrow(well_data))
        if i <= 5
            println("DEBUG: Row $i: Year=$(row.Year), Month=$(row.Month), Rate=$(row["InjectionRate(bbl/month)"])")
        else
            break
        end
    end
    =#
    
    current_month_date = global_start_date
    #println("DEBUG: Starting month loop with current_month_date = $current_month_date")
    
    # Process each month
    processed_months = 0
    while current_month_date <= global_end_date
        y, m = year(current_month_date), month(current_month_date)
        next_m = next_month(current_month_date)
        
        
        processed_months += 1
        
        # Calculate days in this month
        days_in_month = Dates.daysinmonth(current_month_date)
        
        
        # Try different approaches to find matching data
        month_data_numeric = well_data[(well_data.Year .== y) .& (well_data.Month .== m), :]
        month_data_string = well_data[(string.(well_data.Year) .== string(y)) .& (string.(well_data.Month) .== string(m)), :]
        
        
        
        # Determine which approach worked better
        month_data = nrow(month_data_numeric) > 0 ? month_data_numeric : month_data_string
        
        # Calculate daily rate for this month
        if !isempty(month_data)
            
            monthly_vol = first(month_data[!, "InjectionRate(bbl/month)"])
            daily_rate = monthly_vol / days_in_month
            
            
            
            # If rate changed, we add step point
            if abs(daily_rate - current_rate) > 1.0e-6 # check for a change of 0.000001
                days_since_start = (current_month_date - global_start_date).value + 1
                push!(step_times, Float64(days_since_start))
                push!(step_rates, daily_rate)
                current_rate = daily_rate
            else
                #println("DEBUG: No significant rate change (current=$current_rate, new=$daily_rate), skipping step point")
            end
        elseif extrapolate && !isnothing(current_rate) && current_rate > 0
            # Keep using current rate if extrapolating
            println("DEBUG: Extrapolating using current rate $current_rate for month $y-$m")
        else
            # No data for this month and not extrapolating
            # If current rate is non-zero, step down to zero
            if current_rate > 0
                days_since_start = (current_month_date - global_start_date).value + 1
                push!(step_times, Float64(days_since_start))
                push!(step_rates, 0.0)
                current_rate = 0.0
                println("DEBUG: Added step down to zero: day=$(days_since_start) for month $y-$m (no data)")
            else
                println("DEBUG: No data for month $y-$m, current rate already zero")
            end
        end
        
        # Move to next month
        current_month_date = next_m
    end
    
    
    
    # Ensure we end with zero if we had any non-zero rates
    if !isempty(step_rates) && step_rates[end] > 0
        days_total = (global_end_date - global_start_date).value + 1
        push!(step_times, Float64(days_total))
        push!(step_rates, 0.0)
    end
    
    # If we have no data, return empty arrays
    if isempty(step_times)
        #println("DEBUG: NO VALID STEP CHANGES FOUND. step_times and step_rates are empty.")
        return Float64[], Float64[]
    end
    
    # Print summary of steps
    #=
    println("DEBUG: Created $(length(step_times)) injection rate step changes")
    println("DEBUG: step_times: $(step_times)")
    println("DEBUG: step_rates: $(step_rates)")
    println("===== END DEBUG prepare_monthly_fsp_data =====\n")
    =#
    
    return step_times, step_rates
end

"""
    prepare_injection_tool_data(well_data::DataFrame, start_year::Int, end_year::Int, year_of_interest::Int, extrapolate::Bool)

Prepares injection data from injection tool format using monthly averages.
The injection tool data format has daily injection data with columns:
- "UIC Number" or "API Number" or "WellID" for well identification
- "Date of Injection" for the date
- "Volume Injected (BBLs)" for the daily injection volume
"""
function prepare_injection_tool_data(
    well_data::DataFrame, 
    start_year::Int, 
    inj_start_date::Date,
    end_year::Int, 
    inj_end_date::Date,
    year_of_interest::Int, 
    extrapolate::Bool=false,
    year_of_interest_date::Date=Date(year_of_interest-1, 12, 31)
)
    # If well starts after calculation end year, return empty arrays
    if inj_start_date > year_of_interest_date
        println("DEBUG: Well starts after calculation end year - start_year ($start_year) > end_year ($end_year)")
        return Float64[], Float64[]
    end
    
    # Check columns to determine data format
    if "Date of Injection" in names(well_data) && "Volume Injected (BBLs)" in names(well_data)
        # Standard injection tool data format with daily values
        return prepare_daily_injection_tool_data(well_data, start_year, inj_start_date, end_year, inj_end_date, extrapolate, year_of_interest_date)
    end
end

"""
    prepare_daily_injection_tool_data(well_data::DataFrame, start_year::Int, end_year::Int, extrapolate::Bool=false)

Processes standard injection tool data format with daily injection values.
The expected format has:
- "Date of Injection" column for dates
- "Volume Injected (BBLs)" column for daily injection volumes

Returns monthly step changes based on the average daily rates for each month.
"""
function prepare_daily_injection_tool_data(
    well_data::DataFrame, 
    start_year::Int, 
    inj_start_date::Date,
    end_year::Int, 
    inj_end_date::Date,
    extrapolate::Bool=false,
    year_of_interest_date::Date=Date(end_year-1, 12, 31)
)
    # Ensure we have the right columns
    date_col = "Date of Injection"
    volume_col = "Volume Injected (BBLs)"
    
    # Ensure dates are in Date format
    dates = Date[]
    
    # Check if dates are already Date objects
    if eltype(well_data[!, date_col]) <: Date
        dates = well_data[!, date_col]
    else
        # Need to parse from strings
        try
            # Try different date formats
            dates = Date.(well_data[!, date_col], dateformat"y-m-d")
        catch
            try
                dates = Date.(well_data[!, date_col], dateformat"m/d/y")
            catch
                try
                    dates = Date.(well_data[!, date_col], dateformat"m/d/yyyy")
                catch e
                    throw(ArgumentError("Could not parse dates in 'Date of Injection' column: $e"))
                end
            end
        end
    end
    
    if isempty(dates)
        @warn "No valid dates found in injection data"
        return Float64[], Float64[]
    end
    
    # Calculate earliest data date
    earliest_data_date = minimum(dates)
    
    # Use the actual earliest date from the data, but not earlier than start_year
    global_start_date = inj_start_date
    global_end_date = min(inj_end_date, year_of_interest_date)
    
    #println("DEBUG: Using actual start date from data: $global_start_date (earliest data: $earliest_data_date)")
    
    # Arrays for step changes
    step_times = Float64[]
    step_rates = Float64[]
    current_rate = 0.0
    
    # Helper function to get the first day of the next month
    function next_month(d::Date)
        y, m = year(d), month(d)
        if m < 12
            return Date(y, m+1, 1)
        else
            return Date(y+1, 1, 1)
        end
    end
    
    # Filter by date range
    # we filter for days that are on or after the start date and before the end date
    # findall() will return an array of indices that meet the condition
    filtered_indices = findall(d -> global_start_date <= d <= global_end_date, dates)
    if isempty(filtered_indices)
        @warn "No data found in the specified time range"
        return Float64[], Float64[]
    end
    
    filtered_data = well_data[filtered_indices, :]
    filtered_dates = dates[filtered_indices] # convert the filtered indices to dates
    
    # Group by month and calculate average daily rates
    monthly_totals = Dict{Tuple{Int, Int}, Tuple{Float64, Int}}()  # (year, month) => (total_volume, count)
    
    for i in 1:length(filtered_dates)
        d = filtered_dates[i]
        y, m = year(d), month(d)
        
        # Get volume for this day, default to 0 if missing or NaN
        daily_volume = 0.0
        try
            vol = filtered_data[i, volume_col]
            if !ismissing(vol) && !isnan(vol) && vol > 0
                daily_volume = vol
            end
        catch 
            # Keep default of 0
        end
        
        # Add to monthly totals
        if !haskey(monthly_totals, (y, m))
            monthly_totals[(y, m)] = (daily_volume, 1)  # (total_volume, count)
        else
            total, count = monthly_totals[(y, m)]
            monthly_totals[(y, m)] = (total + daily_volume, count + 1)
        end
    end
    
    # Convert monthly totals to average daily rates
    monthly_rates = Dict{Tuple{Int, Int}, Float64}()  # (year, month) => avg_daily_rate
    for key in keys(monthly_totals)
        total_volume, days_count = monthly_totals[key]
        monthly_rates[key] = total_volume / days_count
    end

    
    
    # Fill in missing months with 0 or extrapolated values
    if extrapolate && !isempty(monthly_rates)
        # Get the most recent known rate for extrapolation
        sorted_keys = sort(collect(keys(monthly_rates)))
        last_rate = monthly_rates[sorted_keys[end]]
        
        # Fill in all months in the range
        current_month_date = global_start_date
        while current_month_date < global_end_date
            y, m = year(current_month_date), month(current_month_date)
            if !haskey(monthly_rates, (y, m))
                monthly_rates[(y, m)] = last_rate
            end
            current_month_date = next_month(current_month_date)
        end
    end
    
    # Process month by month to create step changes
    current_month_date = global_start_date
    
    while current_month_date < global_end_date
        y, m = year(current_month_date), month(current_month_date)
        next_m = next_month(current_month_date)
        
        # Clamp if next_m > global_end_date
        if next_m >= global_end_date
            next_m = global_end_date
        end
        
        # Get average rate for this month
        month_daily_rate = get(monthly_rates, (y, m), 0.0)
        
        # If rate changed from current_rate, step up
        if month_daily_rate != current_rate
            push!(step_times, (current_month_date - global_start_date).value + 1)
            push!(step_rates, month_daily_rate)
            current_rate = month_daily_rate
        end
        
        # Check what's next
        if next_m == global_end_date
            # At final boundary, step down if still injecting
            if current_rate > 0
                push!(step_times, (global_end_date - global_start_date).value + 1)
                push!(step_rates, 0.0)
            end
        else
            # Check if next month's rate differs
            next_y, next_mo = year(next_m), month(next_m)
            next_month_rate = get(monthly_rates, (next_y, next_mo), 0.0)
            
            # Step down if rates change and current is non-zero
            if next_month_rate != current_rate && current_rate > 0
                push!(step_times, (next_m - global_start_date).value + 1)
                push!(step_rates, 0.0)
                current_rate = 0.0
            end
        end
        
        # Move to next month
        current_month_date = next_m
    end
    
    return step_times, step_rates
end





# creates a 2D grid of points over the field surface (in km)
# Used for the hydrology pressure field calculations
function create_spatial_grid_km(xmin, xmax, ymin, ymax, number_points)
    x_range_km = range(xmin, stop=xmax, length=number_points)
    y_range_km = range(ymin, stop=ymax, length=number_points)

    Xgrid_km, Ygrid_km = meshgrid(x_range_km, y_range_km)
    
    return Xgrid_km, Ygrid_km, x_range_km, y_range_km
end


# creates a lat/lon meshgrid
function create_spatial_grid_latlon(df::DataFrame, lat_column::String, lon_column::String, num_points::Int=50)
    # first we get the geographic bounds
    lat_min = minimum(df[!, lat_column])
    lat_max = maximum(df[!, lat_column])
    lon_min = minimum(df[!, lon_column])
    lon_max = maximum(df[!, lon_column])
    
    println("DEBUG: Input bounds: Lat [$lat_min, $lat_max], Lon [$lon_min, $lon_max]")
    
    # Validate input coordinates
    if isnan(lat_min) || isnan(lat_max) || isnan(lon_min) || isnan(lon_max)
        error("Invalid coordinates: NaN values detected in bounds")
    end
    
    # Add a small buffer to the bounds (in degrees)
    buffer_deg = 0.01  # Approximately 1 km at mid-latitudes
    lat_min -= buffer_deg
    lat_max += buffer_deg
    lon_min -= buffer_deg
    lon_max += buffer_deg
    
    println("DEBUG: Bounds with buffer (buffer is 0.01 degrees): Lat [$lat_min, $lat_max], Lon [$lon_min, $lon_max]")
    
    # Create 1D ranges for latitude and longitude
    lat_range = range(lat_min, stop=lat_max, length=num_points)
    lon_range = range(lon_min, stop=lon_max, length=num_points)
    
    # Create 2D meshgrids
    LAT_grid, LON_grid = meshgrid(lat_range, lon_range) 
    
    return LAT_grid, LON_grid, lat_range, lon_range
end

# Direct bounds version - takes lat/lon bounds directly
# This is the version that is used in the deterministic hydrology process
function create_spatial_grid_latlon(lat_min::Float64, lat_max::Float64, lon_min::Float64, lon_max::Float64, num_points::Int=50)
    println("DEBUG: Input bounds: Lat [$lat_min, $lat_max], Lon [$lon_min, $lon_max]")
    
    # Validate input coordinates
    if isnan(lat_min) || isnan(lat_max) || isnan(lon_min) || isnan(lon_max)
        error("Invalid coordinates: NaN values detected in bounds")
    end
    
    # Create 1D ranges for latitude and longitude
    # we added the buffer to the bounds before calling this function in the driver script
    lat_range = range(lat_min, stop=lat_max, length=num_points)
    lon_range = range(lon_min, stop=lon_max, length=num_points)
    
    # Create 2D meshgrids
    LAT_grid, LON_grid = meshgrid(lon_range, lat_range) 
    
    println("DEBUG: Created lat/lon grid with size: $(size(LAT_grid))")

    
    return LAT_grid, LON_grid, lat_range, lon_range
end


# function to create a uniform distribution for the probabilistic models
# currently, only the prob hydrology model uses this
function create_uniform_distribution(base_value::Union{Float64, Integer}, plus_minus::Union{Float64, Integer})
    # check for negative plus/minus values
    if plus_minus < 0
        throw(ArgumentError("Plus/minus value cannot be negative. Setting it to 0.0."))
        plus_minus = 0.0
    end

    # check that the plus_minus value is not greater than the base value
    if plus_minus > base_value
        throw(ArgumentError("Plus/minus value cannot be greater than the base value. Setting it to 0.0."))
        plus_minus = 0.0
    end

    min_value = base_value - plus_minus
    max_value = base_value + plus_minus
    return Uniform(min_value, max_value)
end

# This function will create a uniform distribution with proper bounds handling
# parameter_type can be "strike", "dip", "azimuth", or any other parameter
function create_bounded_uniform_distribution(base_value::Union{Float64, Integer}, uncertainty::Union{Float64, Integer}, parameter_type::String="default")
    # Check for negative uncertainty values
    if uncertainty < 0
        throw(ArgumentError("Uncertainty value for $parameter_type cannot be negative."))
    end
    
    # get the min and max values for the distribution before we call the uniform distribution
    # this way we prevent skewed distributions
    min_value = base_value - uncertainty
    max_value = base_value + uncertainty
    
    # Apply the parameter-specific bounds
    # the parameters strike, dip, and azimuth are all in degrees so we need to handle them differently
    # for strike and azimuth angles (0-360 degrees)
    # for dip angles (0-90 degrees)
    # for friction coefficient (0-1)
    # for porosity (0-1)
    # for permeability (must be positive)
    # TO DO: we also need to restrict the aphi value bounds (verify what the range should be)
    if parameter_type == "strike" || parameter_type == "azimuth" || parameter_type == "Strike" || parameter_type == "max_stress_azimuth"
        # For strike and azimuth angles (0-360 degrees)
        # No bounds needed as we'll sample within [min, max] and handle wrapping during sampling
        if max_value - min_value > 360.0
            # If the range exceeds 360 degrees, just use the full 0-360 range
            min_value = 0.0
            max_value = 360.0
        end
    elseif parameter_type == "dip" || parameter_type == "Dip"
        # For dip angles (0-90 degrees)
        min_value = max(min_value, 0.0)
        max_value = min(max_value, 90.0)
    elseif parameter_type == "friction_coefficient" || parameter_type == "FrictionCoefficient"
        # For friction coefficient (0-1)
        min_value = max(min_value, 0.0)
        max_value = min(max_value, 1.0)
    elseif parameter_type == "porosity" || parameter_type == "Porosity" # we get this from the user as a fraction 
        # For porosity (0-1)
        min_value = max(min_value, 0.0)
        max_value = min(max_value, 1.0)
    elseif parameter_type == "permeability" || parameter_type == "Permeability"
        # For permeability (must be positive)
        min_value = max(min_value, 0.0)
    else
        # For other parameters, apply general constraints
        # TO DO: check if we have any other special cases
        # Ensure min_value is not negative for parameters that should be positive
        if base_value > 0
            min_value = max(min_value, 0.0)
        end
    end
    
    # Create the uniform distribution with bounded values
    return Uniform(min_value, max_value)
end


# Converts a grid-based pressure dataset (with latitude, longitude, and pressure values) into 
# a new DataFrame where each grid cell is represented as a polygon (in WKT format) with an associated pressure value.
# Assumes input DataFrame `pressure_grid` has columns: Latitude, Longitude, Pressure_psi
# Requires the original 1D ranges used to create the grid to calculate cell boundaries
function reformat_pressure_grid_to_heatmap_data(pressure_grid::DataFrame, lat_range::AbstractVector, lon_range::AbstractVector)
    n_lat = length(lat_range)
    n_lon = length(lon_range)
    
    if nrow(pressure_grid) != n_lat * n_lon
        error("Pressure grid size does not match the product of latitude and longitude range lengths")
    end
    
    # Calculate half step sizes for latitude and longitude
    # Handle edge cases where range has only one point
    # If n_lat > 1, calculate the step size between consecutive latitude points (e.g., lat_range[2] - lat_range[1]) and divide by 2.
    # If n_lat == 1 (only one point), set half_lat_step to 0.0 (no width).
    half_lat_step = n_lat > 1 ? (lat_range[2] - lat_range[1]) / 2.0 : 0.0
    half_lon_step = n_lon > 1 ? (lon_range[2] - lon_range[1]) / 2.0 : 0.0

    # Initialize arrays for the new DataFrame
    shapes = String[]
    values = Float64[]

    # Iterate through the grid points (assuming row-major order from vec() in deterministic_hydrology_process.jl)
    # Note: The code of this process passes a LAT_grid containing longitude and a LON_grid containing latitude
    # We assume the input `pressure_grid` DataFrame respects this and has columns 'Latitude', 'Longitude'
    # where 'Latitude' corresponds to LON_grid (actual latitudes) and 'Longitude' to LAT_grid (actual longitudes)
    for row in eachrow(pressure_grid)
        lat = row.Latitude 
        lon = row.Longitude 
        pressure = row.Pressure_psi

        # corners of the polygon for this grid cell
        lon_min_cell = lon - half_lon_step
        lon_max_cell = lon + half_lon_step
        lat_min_cell = lat - half_lat_step
        lat_max_cell = lat + half_lat_step

        # polygon coordinates
        coordinates = [
            (lon_min_cell, lat_min_cell),
            (lon_max_cell, lat_min_cell),
            (lon_max_cell, lat_max_cell),
            (lon_min_cell, lat_max_cell),
            (lon_min_cell, lat_min_cell) 
        ]

        # Create a LibGEOS Polygon (check what the constructor expects)
        # Convert coordinates from Vector{Tuple{Float64, Float64}} to Vector{Vector{Float64}}
        shell_coords = [[lon, lat] for (lon, lat) in coordinates]
        
        polygon = LibGEOS.Polygon([shell_coords])


        # Convert the polygon to WKT format
        wkt_polygon = LibGEOS.writegeom(polygon)
        
        push!(shapes, wkt_polygon)
        push!(values, pressure)
    end

    heatmap_df = DataFrame(
        shape = shapes,
        value = values
    )

    # save the heatmap data to a csv file
    #CSV.write("output/heatmap_data.csv", heatmap_df)
    
    
    return heatmap_df
end

# New improved version of the function that correctly handles coordinate ordering
function reformat_pressure_grid_to_heatmap_data_v2(pressure_grid::DataFrame, lat_range::AbstractVector, lon_range::AbstractVector)
    n_lat = length(lat_range)
    n_lon = length(lon_range)
    
    if nrow(pressure_grid) != n_lat * n_lon
        error("Pressure grid size does not match the product of latitude and longitude range lengths")
    end
    
    # Calculate half step sizes for latitude and longitude
    half_lat_step = n_lat > 1 ? (lat_range[2] - lat_range[1]) / 2.0 : 0.0
    half_lon_step = n_lon > 1 ? (lon_range[2] - lon_range[1]) / 2.0 : 0.0

    # Initialize arrays for the new DataFrame
    shapes = String[]
    values = Float64[]

    # IMPORTANT: In create_spatial_grid_latlon, we use meshgrid(lon_range, lat_range)
    # This means in the resulting LAT_grid and LON_grid:
    # - The 'Longitude' column actually contains longitude values (from LAT_grid) 
    # - The 'Latitude' column actually contains latitude values (from LON_grid)
    # Both are correctly named in the DataFrame despite the confusing variable names in deterministic_hydrology_process.jl
    
    for row in eachrow(pressure_grid)
        # Get the actual lat/lon values - column names are correct
        lat = row.Latitude
        lon = row.Longitude
        pressure = row.Pressure_psi

        # Calculate corners of the polygon for this grid cell
        lat_min_cell = lat - half_lat_step
        lat_max_cell = lat + half_lat_step
        lon_min_cell = lon - half_lon_step
        lon_max_cell = lon + half_lon_step

        # Create polygon coordinates in correct (x,y) order for WKT
        # WKT format expects coordinates as (x,y) which corresponds to (longitude, latitude)
        coordinates = [
            (lon_min_cell, lat_min_cell),  # Bottom left
            (lon_max_cell, lat_min_cell),  # Bottom right
            (lon_max_cell, lat_max_cell),  # Top right
            (lon_min_cell, lat_max_cell),  # Top left
            (lon_min_cell, lat_min_cell)   # Close the polygon
        ]

        # Convert to the format LibGEOS expects: [[x1,y1], [x2,y2], ...]
        shell_coords = [[x, y] for (x, y) in coordinates]
        
        # Create the polygon and convert to WKT
        polygon = LibGEOS.Polygon([shell_coords])
        wkt_polygon = LibGEOS.writegeom(polygon)
        
        push!(shapes, wkt_polygon)
        push!(values, pressure)
    end

    return DataFrame(
        shape = shapes,
        value = values
    )
end




# function to get the date bounds from the inejction well data 
# supports all three FSP formats
function get_date_bounds(well_data::DataFrame)
    # check the format of the well data
    # if we have 'month' column, it's FSP Monthly
    if "Month" in names(well_data)
        inj_start_year = minimum(well_data[!, "Year"])
        inj_start_month = minimum(well_data[well_data[!, "Year"] .== inj_start_year, "Month"])
        inj_start_date = Date(inj_start_year, inj_start_month, 1)
        inj_end_year = maximum(well_data[!, "Year"])
        inj_end_month = maximum(well_data[well_data[!, "Year"] .== inj_end_year, "Month"])
        inj_end_date = Date(inj_end_year, inj_end_month, 1)
        inj_end_date = lastdayofmonth(inj_end_date)
    elseif "StartYear" in names(well_data)
        # FSP Annual format
        inj_start_year = first(well_data[!, "StartYear"])
        inj_end_year = first(well_data[!, "EndYear"])
        inj_start_date = Date(inj_start_year, 1, 1)
        inj_end_date = Date(inj_end_year-1, 12, 31)
    else
        # Injection tool data format
        dates = Date[]
        
        # Check if dates are already Date objects
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
                        throw(ArgumentError("Could not parse dates in 'Date of Injection' column: $e"))
                    end
                end
            end
        end
        
        if isempty(dates)
            throw(ArgumentError("No dates found in the well data"))
        end
        inj_start_date = minimum(dates)
        inj_end_date = maximum(dates)
    end

    return inj_start_date, inj_end_date
    
        
        
end

"""
    interpolate_cdf(x_values::Vector{Float64}, y_values::Vector{Float64}, x::Float64)

Interpolate a value on a CDF using Interpolations.jl for improved accuracy and performance.
Returns the probability corresponding to value x based on the provided CDF data points.

Parameters:
- x_values: Vector of x-coordinates (e.g., pressure values) on the CDF, must be sorted
- y_values: Vector of y-coordinates (probabilities) on the CDF, between 0 and 1
- x: The x value to interpolate at

Returns:
- The interpolated probability value
"""
function interpolate_cdf(x_values::Vector{Float64}, y_values::Vector{Float64}, x::Float64)
    # Handle edge cases
    if isempty(x_values) || isempty(y_values)
        @warn "Empty vectors provided to interpolate_cdf"
        return 0.0
    end
    
    # Make copies to avoid modifying the original data
    x_values_copy = copy(x_values)
    y_values_copy = copy(y_values)
    
    # Sort the values if they're not already sorted
    if !issorted(x_values_copy)
        p = sortperm(x_values_copy)
        x_values_copy = x_values_copy[p]
        y_values_copy = y_values_copy[p]
    end
    
    # Deduplicate knots explicitly to suppress warnings
    # This returns the indices of unique elements
    unique_indices = Interpolations.deduplicate_knots!(x_values_copy; move_knots=false)
    
    # Use only the unique x values and their corresponding y values
    if length(unique_indices) < length(y_values_copy)
        y_values_copy = y_values_copy[unique_indices]
    end
    
    # Create interpolation object with flat extrapolation behavior
    # This will return the endpoint values when x is outside the domain
    itp = LinearInterpolation(x_values_copy, y_values_copy, extrapolation_bc=Flat())
    
    # Return interpolated value
    return itp(x)
end


# Since we can have three possible formats, we need to reformat them all to the same format
function clean_up_well_data(well_data::DataFrame, injection_wells_format::String)
    # initialize a dataframe with a uniform format
    uniform_well_data = DataFrame(
        WellID = String[],
        Month = Int[],
        InjectionRate = Union{Float64, Int}[],
        Year = Int[],
        LastInjectionDate = Date[]
    )


    if injection_wells_format == "annual_fsp"
        # Here we have columns 'WellID', 'StartYear', 'EndYear', 'InjectionRate(bbl/day)'
        # We need to get the average monthly injection rates
        # get unique well ids and iterate over their rows
        well_ids = unique(well_data[!, "WellID"])
        for well_id in well_ids
            # get injection start and end dates
            inj_start_date = Date(minimum(well_data[well_data[!, "WellID"] .== well_id, "StartYear"]), 1, 1)
            inj_end_date = Date(maximum(well_data[well_data[!, "WellID"] .== well_id, "EndYear"]) - 1, 12, 31)

            # get the daily injection rate
            inj_rate = well_data[well_data[!, "WellID"] .== well_id, "InjectionRate(bbl/day)"]

            # find the monthly injection rates for all months (1-12)
            for month in 1:12
                # get the month's start and end dates
                month_start_date = Date(inj_start_date.year, month, 1)
                month_end_date = lastdayofmonth(Date(inj_start_date.year, month, 1))

                # get the average injection rate for the month
                avg_inj_rate = mean(inj_rate[inj_rate .>= month_start_date .&& inj_rate .<= month_end_date])

                # add the data to the uniform dataframe
                push!(uniform_well_data, (WellID = well_id, Month = month, InjectionRate = avg_inj_rate, Year = inj_start_date.year))
            end
            # add the last injection date
            push!(uniform_well_data, (WellID = well_id, LastInjectionDate = inj_end_date))
        end

    elseif injection_wells_format == "monthly_fsp"
        # Here we have columns 'WellID', 'Month', 'InjectionRate(bbl/month)', 'Year'
        # We already have 'WellID' 'Month' and 'Year' in the 'monthly_fsp' case
        # get unique well ids and iterate over their rows
        well_ids = unique(well_data[!, "WellID"])
        for well_id in well_ids
            # get injection start and end dates
            # find 
            
        end
        
    end
end



end # module

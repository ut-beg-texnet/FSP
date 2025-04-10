using CSV
using DataFrames
using Geodesy


function convert_latlon_to_easting_northing(df::DataFrame)
    required_columns = ["fault_id", "longitude", "latitude"]
    df_cols = lowercase.(names(df))

    missing_cols = setdiff(required_columns, df_cols)
    if !isempty(missing_cols)
        throw(ArgumentError("Missing required columns: $(missing_cols)"))
    end

    lat_column = findfirst(==(lowercase("latitude")), df_cols)
    lon_column = findfirst(==(lowercase("longitude")), df_cols)

    if any(lat -> !(-90 <= lat <= 90), df[:, lat_column])
        throw(ArgumentError("Latitude values must be in range [-90, 90] degrees"))
    end
    if any(lon -> !(-180 <= lon <= 180), df[:, lon_column])
        throw(ArgumentError("Longitude values must be in range [-180, 180] degrees"))
    end

    n_rows = nrow(df)
    eastings = Vector{Float64}(undef, n_rows)
    northings = Vector{Float64}(undef, n_rows)
    zones = Vector{String}(undef, n_rows)

    for i in 1:n_rows
        try
            # Retrieve latitude and longitude
            lat = df[i, lat_column]
            lon = df[i, lon_column]

            # Create an LLA object (altitude is set to 0.0)
            point_lla = LLA(lat, lon, 0.0)

            # Determine UTM zone and hemisphere
            zone = floor(Int, (lon + 180) / 6) + 1
            hemisphere = lat >= 0  # true for Northern Hemisphere

            # Define UTM transformation
            utm_transform = UTMfromLLA(zone, hemisphere, wgs84)

            # Convert to UTM coordinates
            point_utm = utm_transform(point_lla)

            # Extract easting and northing
            eastings[i] = point_utm.x
            northings[i] = point_utm.y

            # Extract zone
            zones[i] = "$(zone)$(hemisphere ? 'N' : 'S')"
        catch e
            throw(ArgumentError("Error converting latitude/longitude to easting/northing for row $i: $(e)"))
        end
    end

    # Add the new columns to the dataframe
    df_with_easting_northing = copy(df)
    df_with_easting_northing.easting = eastings
    df_with_easting_northing.northing = northings
    df_with_easting_northing.zone = zones

    return df_with_easting_northing
end




function main()
    # Create a dataframe with hardcoded data (fault_id, latitude, longitude)
    df = DataFrame(fault_id = [1, 2, 3],
                   latitude = [10, 25, 30.4],
                   longitude = [-20, 30, 40])

    # Convert the dataframe to include easting and northing
    df_with_easting_northing = convert_latlon_to_easting_northing(df)

    # Print the dataframe
    println(df_with_easting_northing)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end



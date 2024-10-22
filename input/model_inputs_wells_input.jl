module ModelInputsWellsInput

using CSV
using DataFrames
using JSON

# Define structures for Wells data
export Well, WellMonthly, load_wells_from_csv, convert_csv_to_json

# Well --> For manually entered Wells with constant injection rates
struct Well
    id::Int
    x::Float64
    y::Float64
    injection_rate::Float64  # bbl/day
    start_year::Int
    end_year::Int
end

# WellMonthly --> For wells with monthly injection data from the CSV
struct WellMonthly
    id::String
    x::Float64
    y::Float64
    year::Int
    month::Int
    volume::Float64  # bbl/month
end

# loads wells from a CSV file and returns appropriate structures ('Well' data structure for constant rates, 'WellMonthly' data structure for monthly data)
function load_wells_from_csv(file_path::String, header_lines::Int=1)
    try
        df = CSV.read(file_path, DataFrame; header=header_lines)

        if size(df, 2) == 6
            # Validate columns for monthly well data (6-column CSV)
            if names(df) != ["UniqueID/Name", "Easting (km)", "Northing (km)", "Year", "Month (1-12)", "InjectionVolume (bbl/month)"]
                throw(ArgumentError("CSV file must contain exactly 6 columns: UniqueID/Name, Easting (km), Northing (km), Year, Month (1-12), InjectionVolume (bbl/month)"))
            end

            wells = WellMonthly[]
            for row in eachrow(df)
                push!(wells, WellMonthly(
                    row[1],  # UniqueID/Name -> String
                    row[2],  # Easting (km) -> Float64
                    row[3],  # Northing (km) -> Float64
                    row[4],  # Year -> Int
                    row[5],  # Month (1-12) -> Int
                    row[6]   # InjectionVolume (bbl/month) -> Float64
                ))
            end
            return wells

        elseif size(df, 2) == 5
            # column validation for constant injection rates (5-column CSV)
            if names(df) != ["Xwell", "Ywell", "InjectionRate", "StartYear", "EndYear"]
                throw(ArgumentError("CSV file must contain exactly 5 columns: Xwell, Ywell, InjectionRate, StartYear, EndYear"))
            end

            wells = Well[]
            for (i, row) in enumerate(eachrow(df))
                push!(wells, Well(
                    i,         # Automatically assign an ID
                    row[1],    # Xwell -> Float64
                    row[2],    # Ywell -> Float64
                    row[3],    # InjectionRate -> Float64
                    row[4],    # StartYear -> Int
                    row[5]     # EndYear -> Int
                ))
            end
            return wells

        else
            throw(ArgumentError("CSV file must contain either 5 or 6 columns."))
        end
    catch e
        println("Error reading CSV file: ", e)
        rethrow(e)
    end
end

# converts the loaded wells data to JSON format and save it to a file
function convert_csv_to_json(csv_file_path::String, json_file_path::String)
    try
        wells = load_wells_from_csv(csv_file_path)
        json_data = JSON.json(wells)

        open(json_file_path, "w") do file
            write(file, json_data)
        end
        println("Data successfully saved to JSON at: $json_file_path")
        return json_file_path
    catch e
        println("Error during CSV to JSON conversion: ", e)
        rethrow(e)
    end
end

end  # End of module

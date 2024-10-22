module ModelInputsWellsOutput

using JSON


export reformat_json_data

function reformat_json_data(json_file_path::String)
    file_data = open(json_file_path, "r") do file
        JSON.parse(file)
    end

    # The reformatted data structure that D3.js will use
    reformatted_data = Dict{String, Vector{Dict{Symbol, Real}}}()

    for well_data in file_data
        # Well name/ID is the key
        well_id = string(well_data["id"])

        if haskey(well_data, "start_year") && haskey(well_data, "end_year") && haskey(well_data, "injection_rate")
            # Constant injection rate well (Well struct)
            start_year = well_data["start_year"]
            end_year = well_data["end_year"]
            injection_rate = well_data["injection_rate"]

            # Initialize empty vector for this well if it doesn't exist
            if !haskey(reformatted_data, well_id)
                reformatted_data[well_id] = []
            end

            # Add each year as an x-axis point with injection rate as y-axis
            for year in start_year:end_year
                push!(reformatted_data[well_id], Dict(:x => year, :y => injection_rate))
            end

        elseif haskey(well_data, "year") && haskey(well_data, "month") && haskey(well_data, "volume")
            # Monthly injection rate well (WellMonthly struct)
            year = well_data["year"]
            month = well_data["month"]
            volume = well_data["volume"]
            
            if !haskey(reformatted_data, well_id)
                reformatted_data[well_id] = []
            end

            # make sure 'year' is an int (no decimal places)
            push!(reformatted_data[well_id], Dict(:x => Int(year), :y => volume))
        else
            throw(ArgumentError("Unexpected well data format. Data must contain either (start_year, end_year, injection_rate) or (year, month, volume)."))
        end
    end

    return JSON.json(reformatted_data)
end

end  # End module

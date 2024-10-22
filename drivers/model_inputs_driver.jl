module DriverStep1

using JSON
using CSV

include("src/input/model_inputs_wells_input.jl")
include("src/output/model_inputs_wells_output.jl")

using .ModelInputsWellsInput
using .ModelInputsWellsOutput

# Driver logic for Step 1 line graphs
# helper script #1 reads the CSV and converts to JSON
# helper script #2 reformats the JSON for D3.js line graphs

# Step 1: Convert the CSV file to JSON using model_inputs_wells_input.jl
csv_file_path = "src/output/mock_wells.csv"
json_file_path = "src/output/well_data_line_graph.json"

# CSV --> JSON
json_file = convert_csv_to_json(csv_file_path, json_file_path)


# Step 2: Reformat the JSON data using model_inputs_wells_output.jl
reformatted_data = reformat_json_data(json_file)

# Step 3: return reformatted_data for use in the D3.js module
println("Reformatted Data for D3.js: ", reformatted_data)

# Step 4: Save the reformatted data to a new JSON file
output_json_file_path = "src/output/reformatted_well_data.json"
open(output_json_file_path, "w") do file
    write(file, reformatted_data)
end


end # End of driver module

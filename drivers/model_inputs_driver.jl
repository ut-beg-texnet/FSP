#module DriverStep1

using JSON
using CSV
using ArgParse
using DataFrames
using Dates
using Parameters
using StaticArrays
using Random
using LoggingExtras
using TypedTables
#using Serialization

#include("src/input/get_faults_input_csv.jl") # used to get fault data from csv
#include("src/input/get_well_input_csv.jl") # used to get well data from csv
#include("src/input/get_stress_input_csv.jl") # used to get stress data from csv
#include("src/input/get_hydrology_input_csv.jl") # used to get hydrology data from csv
#include("src/output/model_inputs_wells_output.jl")
#include("src/surfaceviz/surfaceviz.jl") # main data structure
#include("src/surfaceviz/setup_data.jl") # data initialization
#include("src/session/serialize_state.jl") # state serialization
include("src/core/fault_data_input.jl")
include("src/core/stress_data_input.jl")
include("src/core/hydrology_data_input.jl")

#using .GetFaultsInputCSV
#using .GetWellInputCSV
#using .GetStressDataCSV
#using .GetHydrologyInputCSV


#using .ModelInputsWellsOutput
#using .SurfaceViz
#using .SetupData



function parse_cli_args()
    s = ArgParseSettings(description="FSP Model Inputs Driver")

    @add_arg_table! s begin
        "--fault-data"
            help = "Path to fault data CSV"
            required = true
        "--stress-data"
            help = "Path to stress data CSV"
            required = true
        "--aphi"
            help = "A-phi value (required when using A-phi stress model)"
            arg_type = Float64
            required = false
        "--well-data"
            help = "Path to well data CSV (monthly or continuous format)"
            required = true
        "--friction-coefficient"
            help = "Friction coefficient value (applied to all faults)"
            arg_type = Float64
            required = true
        "--porosity"
            help = "Porosity (in %). Required for internal hydrology model"
            arg_type = Float64
            required = false
        "--permeability"
            help = "Permeability (in mD). Required for internal hydrology model"
            arg_type = Float64
            required = false
        "--aquifer-thickness"
            help = "Aquifer thickness (in feet). Required for internal hydrology model"
            arg_type = Float64
            required = false
        "--density"
            help = "Fluid density (in kg/m^3)."
            arg_type = Float64
            default = 1000.0
        "--dynamic_viscosity"
            help = "Dynamic viscosity (in Pa.s)."
            arg_type = Float64
            default = 0.0008
        "--fluid_compressibility"
            help = "Fluid compressibility (in 1/Pa)."
            arg_type = Float64
            default = 3.6e-10
        "--rock-compressibility"
            help = "Rock compressibility (in 1/Pa)."
            arg_type = Float64
            default = 1.08e-9
        "--output-json"
            help = "Path to output JSON file for this step's results"
            required = true
    end

    return parse_args(s)

end

# check for required arguments
function check_fault_data(args::Dict)
    #println("args[\"fault-data\"]: ", args["fault-data"])
    if haskey(args, "fault-data") && !isnothing(args["fault-data"])
        return read_fault_data(args["fault-data"], args["friction-coefficient"])
    else
        error("Must provide --faul-data argument")
    end
end

function check_stress_data(args::Dict)
    if haskey(args, "stress-data") && !isnothing(args["stress-data"]) && haskey(args, "aphi")
        return read_stress_data(args["stress-data"], args["aphi"])
    elseif haskey(args, "stress-data") && !isnothing(args["stress-data"] && !haskey(args, "aphi"))
        return read_stress_data(args["stress-data"])
    else 
        error("Must provide --stress-data argument")
    end
end

function check_well_data(args::Dict)
    if haskey(args, "well-data") && !isnothing(args["well-data"])
        return read_well_data(args["well-data"])
    else
        error("Must provide --well-data argument")
    end
end

# here we can expect two types of hydrology data
# 1. internal hydrology model (user provides porosity, permeability, aquifer thickness as cli arguments)
# 2. external hydrology model (user provides external hydrology model CSV file path as cli argument)
function check_hydrology_data(args)
    # check for external model
    if haskey(args, "hydrology-data") && !isnothing(args["hydrology-data"])
        if !isfile(args["hydrology-data"])
            error("Hydrology data file not found: $(args["hydrology-data"])")
        end
        return Dict(
            "hydrology_model_type" => "external",
            "file_source" => args["hydrology-data"],
            # read the CSV into a df and convert to a list of dictionaries
            "data" => CSV.read(args["hydrology-data"], DataFrame) |> DataFrame -> [Dict(pairs(row)) for row in eachrow(DataFrame)]
        )
    else
        # internal Model
        required_params = ["porosity", "permeability", "aquifer-thickness"]
        missing_params = filter(p -> !haskey(args, p) || isnothing(args[p]), required_params)

        if !isempty(missing_params)
            error("Missing required hydrology parameters: $(join(missing_params, ", "))")
        end

        return Dict(
            "hydrology_model_type" => "internal",
            "internal_hydrology_params" => Dict(
                "porosity" => args["porosity"],
                "permeability" => args["permeability"],
                "aquifer_thickness" => args["aquifer-thickness"]
                # might have to add compresibility, fluid density, dynamic viscosity etc here
            )
        )
    end

end

function get_session_metadata()
    return Dict{String,Any}(
        "timestamp" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        "julia_version" => string(VERSION),
        "os" => string(Sys.KERNEL),
        "machine" => string(Sys.MACHINE)
    )
end


function main()
    # parse CLI arguments
    args = parse_cli_args()

    # process input data
    fault_data = check_fault_data(args)
    stress_data = check_stress_data(args)
    well_data = check_well_data(args)
    hydrology_data = check_hydrology_data(args)

    # output directory for JSON file, create it if it doesn't exist
    output_dir = dirname(args["output-json"])
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    # output dict
    output = Dict{String,Any}(
        "metadata" => get_session_metadata(),
        "fault_data" => to_dict(fault_data),
        "stress_data" => to_dict(stress_data),
        "well_data" => hydrology_to_dict(well_data),
        "hydrology_data" => hydrology_data
    )

    # write output to JSON file
    open(args["output-json"], "w") do f
        JSON.print(f, output, 4)
    end
end



#end # step 1 end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end




#=
function main()
    # hardcode fault data CSV file path for testing
    #file_path = "/home/seiscomp/fsp_3/fsp_3/src/input/mock_fault_data_input.csv"
    file_path = "/home/seiscomp/fsp_3/fsp_3/src/input/Bluebonnet_Faults_for_FSP.csv"
    # load fault data from CSV
    faults = load_faults_from_csv(file_path)
    println("Faults loaded from CSV: ", faults)

    # convert fault data to JSON
    json_output_file = "step1_fault_data_output.json"
    fault_data_to_json(faults, json_output_file)

    open(json_output_file, "r") do file
        json_content = read(file, String)
        println("JSON content: ", json_content)
    end
    

    # Well data testing
    well_data_input_file_path = "/home/seiscomp/fsp_3/fsp_3/src/input/mock_well_data_input.csv"

    # load well data from CSV
    well_data = load_wells_from_csv(well_data_input_file_path)

    println("Wells loaded from CSV: ", well_data)

    # convert well data to JSON
    well_data_output_file = "step1_well_data_output.json"
    well_data_to_json(well_data, well_data_output_file)
    println("Well data JSON file created: ", well_data_output_file)

    

    # get stress data from csv testing 
    stress_input_csv_file = "/home/seiscomp/fsp_3/fsp_3/src/input/mock_stress_data_input.csv"
    stress_data = load_stress_from_csv(stress_input_csv_file)

    # convert stress data to JSON
    stress_data_json_output_file = "step1_stress_data_output.json"
    stress_data_to_json(stress_data, stress_data_json_output_file)
    println("Stress data: ", stress_data)

    

    # get hydrology data from csv testing
    hydrology_input_csv_file = "/home/seiscomp/fsp_3/fsp_3/src/input/mock_hydrology_data_input.csv"
    hydrology_data = load_hydrology_from_csv(hydrology_input_csv_file)

    # convert hydrology data to JSON
    hydrology_data_json_output_file = "step1_hydrology_data_output.json"
    hydrology_data_to_json(hydrology_data, hydrology_data_json_output_file)
    println("Hydrology data: ", hydrology_data)


end




end # step 1 end

DriverStep1.main()  # Call the main function to execute the driver step 1
=#

module DriverStep1

using JSON
using CSV
#using Serialization

include("src/input/get_faults_input_csv.jl") # used to get fault data from csv
include("src/input/get_well_input_csv.jl") # used to get well data from csv
include("src/input/get_stress_input_csv.jl") # used to get stress data from csv
include("src/input/get_hydrology_input_csv.jl") # used to get hydrology data from csv
#include("src/output/model_inputs_wells_output.jl")
#include("src/surfaceviz/surfaceviz.jl") # main data structure
#include("src/surfaceviz/setup_data.jl") # data initialization
#include("src/session/serialize_state.jl") # state serialization

using .GetFaultsInputCSV
using .GetWellInputCSV
using .GetStressDataCSV
using .GetHydrologyInputCSV
#using .ModelInputsWellsOutput
#using .SurfaceViz
#using .SetupData


#=
# Driver logic for Step 1 of FSP
mutable struct FSPState

    surfaceviz::SurfaceVizStruct

end


# create structure and set up default data
function initialize_state()::Union{FSPState, Nothing}
    try
        println("Initializing FSP State...")
        viz_state = SurfaceVizStruct()  # Create a new SurfaceViz structure
        setupdata(viz_state)  # Populate with default data
        println("FSP State initialized.")
        return FSPState(viz_state)
    catch e
        println("Error initializing FSP State: ", e)
        return nothing
    end
end

function run_step1(state::FSPState)
    try
        println("Running Step 1...")
        # Perform Step 1 operations on well data and fault data

        # Convert the CSV to JSON format for well and fault data
        print("Converting well data to JSON...")
        well_data_json = convert_csv_to_json("src/output/mock_wells.csv", "src/output/well_data.json")
        println("Well data converted to JSON.")
        print("Converting fault data to JSON...")
        fault_data_json = convert_csv_to_json("src/output/fault_data_test_input.csv", "src/output/fault_data.json")
        println("Fault data converted to JSON.")

        # get data from JSON file to update data structure
        println("Updating state with fault data...")
        num_faults, friction_coeff, strike_min, strike_max, dip_min, dip_max = get_fault_data_from_json("../output/fault_data_test_input.json")
        
        # Set fault data in SurfaceVizStruct
        set_fault_data(state.surfaceviz, num_faults, friction_coeff, strike_min, strike_max, dip_min, dip_max)

        println("Updated fault data in SurfaceVizStruct: ", state.surfaceviz.data[:fault])


        println("Reformatting well data for D3 visualization...")
        # Process and reformat data for D3 visualization
        reformatted_well_data = reformat_json_data(state.surfaceviz.data[:well_data])
        

        # Update state with reformatted data for visualization
        state.surfaceviz.plotdata[:well_data] = reformatted_well_data
        

        # Save updated plot data to JSON files for future use
        # Only well data gets reformatted (we are essentially renaming the axes for consistency within the D3 module)
        open("src/output/reformatted_well_data.json", "w") do file
            JSON.print(file, reformatted_well_data)
        end
        open("src/output/fault_data.json", "w") do file
            JSON.print(file, reformatted_fault_data)
        end

        # serialize the state so we can use it on step 2
        println("Serializing step 1 state...")
        serialize_state(state, "../output/fsp_step1_state.jls")
        println("Step 1 state serialized.")

        return state  # Returns an updated state
    catch e
        println("Error in Step 1: ", e)
        return nothing
    end
end
=#

function main()
    # hardcode fault data CSV file path for testing
    file_path = "/home/seiscomp/fsp_3/fsp_3/src/input/mock_fault_data_input.csv"
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

# OLD LOGIC (NOT BEING USED) --------------------------------------------------------------------------------------------


    #=
    # Initialize FSP state
    fsp_state = initialize_state()
    if fsp_state === nothing
        println("Failed to initialize state. Exiting Step 1.")
        return
    end
    

    new_fsp_state = run_step1(fsp_state)
    if new_fsp_state === nothing
        println("Step 1 encountered an error. Exiting.")
        return
    end
    
    println("Step 1 completed. Reformatted data saved for D3.js visualization")
    println("Reformatted well data: ", new_fsp_state.surfaceviz.plotdata[:well_data])
    println("Fault data: ", new_fsp_state.surfaceviz.plotdata[:fault_data])
    =#
end




end # step 1 end

DriverStep1.main()  # Call the main function to execute the driver step 1

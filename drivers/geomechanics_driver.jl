module GeomechanicsDriver

using JSON
using ArgParse
using DataFrames

include("src/core/DeterministicGeomechanicsCalculations.jl") 
using .DeterministicGeomechanicsCalculations




function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input-json"
            help = "Path to input JSON file from model_inputs step"
            arg_type = String
            required = true
        "--output-json"
            help = "Path to output JSON file"
            arg_type = String
            required = true
    end

    return parse_args(s)
end

# Function to load JSON data
function load_json_data(filepath::String)
    try
        return JSON.parsefile(filepath)
    catch e
        println("Error loading JSON data from ", filepath, ": ", e)
        return nothing
    end
end

# Function to save JSON data
function save_json_data(filepath::String, data::Dict{String,Any})
    open(filepath, "w") do f
        JSON.print(f, data, 4)
    end
end



# main driver function
function geomechanics_driver()


    println("\n=== Starting Deterministic Geomechanics Analysis ===")
    args = parse_commandline()

    # Read input JSON
    println("\nReading input data from: $(args["input-json"])")
    input_data = load_json_data(args["input-json"])
    if input_data === nothing
        println("Failed to load input data. Exiting.")
        return
    end
    
    println("Processing $(length(input_data["fault_data"])) faults")
    
    # Extract fault parameters
    fault_data = input_data["fault_data"]
    strikes = [fault["strike"] for fault in fault_data]
    dips = [fault["dip"] for fault in fault_data]
    muf = [fault["friction_coefficient"] for fault in fault_data]

    # extract stress data
    stress_data = input_data["stress_data"]
    SHdir = stress_data["max_horiz_direction"]
    dpth = stress_data["reference_depth"]


    # Calculate native pore pressure, ensuring it's wrapped as an array
    pp0 = [stress_data["pore_pressure_gradient"] * dpth]

    # Hardcoded Poisson's ratio and Biot coefficient
    nu = 0.5
    biot = 1.0

    # Parse CLI arguments for `aphi` and `min_hor_stress`
    #aphi_value, min_hor_stress = parse_commandline()
    

    # Prepare dp as an array of zeros for fault count, required for calculations
    dp = zeros(length(strikes))

    # get stress model
    stress_model = get(stress_data, "stress_model", nothing)

    # get min_hor_stress (if provided)
    min_hor_stress = get(stress_data, "min_horiz_gradient", nothing)

    # initialize empty stress vector with three elements
    sig = [get(stress_data, "vertical_gradient", 0.0), get(stress_data, "min_horiz_gradient", 0.0), get(stress_data, "max_horiz_gradient", 0.0)]

    # Determine which inputs to use based on CLI arguments and stress data availability
    if stress_model !== "complete"
        println("A-Phi value provided. Calculating SH and Sh from A-Phi...")
        # get aphi value (if provided)
        aphi_value = get(stress_data, "aphi_value", nothing)
        
        # Case 3: Both --aphi and --min_hor_stress are provided
        if get(stress_data, "min_horiz_gradient", nothing) !== nothing
            #println("Minimum horizontal stress provided: $min_horiz_gradient")
            #min_hor_stress = round(min_hor_stress * dpth, digits=2)
            # Calculate stress at reference depth
            sig[1] = round(stress_data["vertical_gradient"] * dpth, digits=2)
            sig[2] = round(stress_data["min_horiz_gradient"] * dpth, digits=2)

            geomechanics_calculations_inputs = (sig, pp0, strikes, dips, SHdir, dp, muf, biot, nu, aphi_value, min_hor_stress, stress_model)
        
        # Case 2: --aphi provided without --min_hor_stress (ignore min horizontal stress)
        else
            println("No minimum horizontal stress provided with A-Phi; calculating SH and Sh without min horizontal stress.")
            sig[1] = round(stress_data["vertical_gradient"] * dpth, digits=2)
            geomechanics_calculations_inputs = (
                sig, pp0, strikes, dips, SHdir, dp, muf, biot, nu, aphi_value, nothing, stress_model
            )
        end
    else
        # calculate all three stresses by multyplying the gradients by the depth
        sig[1] = round(stress_data["vertical_gradient"] * dpth, digits=2)
        sig[2] = round(stress_data["min_horiz_gradient"] * dpth, digits=2)
        sig[3] = round(stress_data["max_horiz_gradient"] * dpth, digits=2)
        geomechanics_calculations_inputs = (
            sig, pp0, strikes, dips, SHdir, dp, muf, biot, nu, nothing, min_hor_stress, stress_model
        )
    end

    # Run deterministic geomechanics calculations
    det_geomechanic_results = mohrs_3D(geomechanics_calculations_inputs)

    # Process and store results
    ppfail = det_geomechanic_results.outs["ppfail"]
    ppfail[ppfail .< 0] .= 0.0 # set negative values to zero
    println("ppfail results: ", ppfail)


    
    #save_json_data(args["output-json"], det_geomechanic_results.outs)

    # append the results to the JSON file from the previous step
    append_to_json(args["input-json"], det_geomechanic_results, args["output-json"])

    #=
    # Save results to JSON
    det_geomech_results_to_json(det_geomechanic_results, "src/output/deterministic_geomechanics_results.json")
    println("Geomechanics driver's results stored successfully!")

    # Serialize data for step 3
    serialize_data(fault_data, stress_data, "src/output/step2_user_input_data.json")
    =#
end


geomechanics_driver()

end # End of mmodule


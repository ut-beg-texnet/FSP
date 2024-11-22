module GeomechanicsDriver

using JSON

include("src/core/DeterministicGeomechanicsCalculations.jl") 
using .DeterministicGeomechanicsCalculations


function parse_cli_args()
    
    aphi_value = nothing
    min_hor_stress = nothing

    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--aphi"
            if i + 1 <= length(ARGS)
                i += 1
                aphi_value = try
                    parse(Float64, ARGS[i])
                catch e
                    error("Error: A-Phi value must be a valid Float64")
                end
            else
                error("Error: --aphi option requires a value")
            end
        elseif arg == "--min_hor_stress"
            if i + 1 <= length(ARGS)
                i += 1
                min_hor_stress = try
                    parse(Float64, ARGS[i])
                catch e
                    error("Error: min_hor_stress value must be a valid Float64")
                end
            else
                error("Error: --min_hor_stress option requires a value")
            end
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end

    return (aphi_value, min_hor_stress)
end


function load_json_data(file_path::String)
    try
        data = JSON.parsefile(file_path)
        println("Successfully loaded data from ", file_path)
        return data
    catch e
        println("Error loading data from ", file_path, ": ", e)
        return nothing
    end
end


# serialize user input data (to be used in next steps)
function serialize_data(fault_data, stress_data, file_path::String)
    user_data = Dict("fault_data" => fault_data, "stress_data" => stress_data)
    open(file_path, "w") do io
        JSON.print(io, user_data)
    end
    println("Step 2: User input data saved to ", file_path)
end




# main driver function
function geomechanics_driver()
    # Load fault data from step 1 output CSV
    fault_data = load_json_data("step1_fault_data_output.json")
    if fault_data === nothing
        println("Failed to load fault data. Exiting.")
        return
    end

    # Extract fault parameters
    thf = [fault["strike"] for fault in fault_data]
    dips = [fault["dip"] for fault in fault_data]
    muf = [fault["friction"] for fault in fault_data]

    # Load stress data from JSON file
    stress_data = load_json_data("step1_stress_data_output.json")
    if stress_data === nothing
        println("Failed to load stress data. Exiting.")
        return
    end

    # Extract parameters from stress data
    SHdir = stress_data["max_horizontal_stress_direction"]
    dpth = stress_data["reference_depth_ft"]

    # Calculate stress at reference depth
    sig = round.([
        stress_data["vertical_stress_gradient_psi_ft"],
        stress_data["min_horizontal_stress_gradient_psi_ft"],
        stress_data["max_horizontal_stress_gradient_psi_ft"]
    ] .* dpth, digits=2)

    # Calculate native pore pressure, ensuring it's wrapped as an array
    pp0 = [stress_data["initial_reservoir_pressure_gradient_psi_ft"] * dpth]

    # Hardcoded Poisson's ratio and Biot coefficient
    nu = 0.5
    biot = 1.0

    # Parse CLI arguments for `aphi` and `min_hor_stress`
    aphi_value, min_hor_stress = parse_cli_args()

    # Prepare dp as an array of zeros for fault count, required for calculations
    dp = zeros(length(thf))

    # Determine which inputs to use based on CLI arguments and stress data availability
    if aphi_value !== nothing
        println("A-Phi value provided. Calculating SH and Sh from A-Phi...")
        
        # Case 3: Both --aphi and --min_hor_stress are provided
        if min_hor_stress !== nothing
            println("Minimum horizontal stress provided via CLI: $min_hor_stress")
            min_hor_stress = round(min_hor_stress * dpth, digits=2)
            geomechanics_calculations_inputs = (
                sig, 0.00, pp0, thf, dips, SHdir, dp, muf, biot, nu, aphi_value, min_hor_stress
            )
        
        # Case 2: --aphi provided without --min_hor_stress (ignore min horizontal stress)
        else
            println("No minimum horizontal stress provided with A-Phi; calculating SH and Sh without min horizontal stress.")
            geomechanics_calculations_inputs = (
                sig, 0.00, pp0, thf, dips, SHdir, dp, muf, biot, nu, aphi_value
            )
        end
    else
        # Case 1: No CLI arguments; use stress data's min horizontal stress
        println("No A-Phi value provided. Calculating SH and Sh from stress data.")
        geomechanics_calculations_inputs = (
            sig, 0.00, pp0, thf, dips, SHdir, dp, muf, biot, nu
        )
    end

    # Run deterministic geomechanics calculations
    det_geomechanic_results = mohrs_3D(geomechanics_calculations_inputs)

    # Process and store results
    ppfail = det_geomechanic_results.outs["ppfail"]
    ppfail[ppfail .< 0] .= 0.0
    println("ppfail results: ", ppfail)

    # Save results to JSON
    det_geomech_results_to_json(det_geomechanic_results, "src/output/deterministic_geomechanics_results.json")
    println("Geomechanics driver's results stored successfully!")

    # Serialize data for step 3
    serialize_data(fault_data, stress_data, "src/output/step2_user_input_data.json")
end


geomechanics_driver()

end # End of mmodule

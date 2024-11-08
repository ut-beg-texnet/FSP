module GeomechanicsDriver

using JSON

include("src/core/DeterministicGeomechanicsCalculations.jl") 
using .DeterministicGeomechanicsCalculations


function parse_cli_args()
    # Initialize default values
    aphi_value = nothing
    min_hor_stress = nothing

    # Iterate over ARGS to parse command-line arguments
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


# main driver function
function geomechanics_driver()
    # load fault data from step 1 output CSV
    fault_data = load_json_data("step1_fault_data_output.json")
    if fault_data === nothing
        println("Failed to load fault data. Exiting.")
        return
    end

    #println("Fault data structure: ", typeof(fault_data)) # for debugging
    #println("Fault data content: ", fault_data)

    #println("Extracting fault parameters...")

    # get fault strikes 
    thf = [fault["strike"] for fault in fault_data]
    # get fault dips
    dips = [fault["dip"] for fault in fault_data]
    # get fault friction coefficients (the same one is used for all faults in the calculations)
    muf = [fault["friction"] for fault in fault_data]
    #println("Fault strikes (thf): ", thf)
    #println("Fault dips: ", dips)
    #println("Fault friction coefficients (muf): ", muf)

    stress_data = load_json_data("step1_stress_data_output.json")
    if stress_data === nothing
        println("Failed to load stress data. Exiting.")
        return
    end

    #println("Stress data structure: ", typeof(stress_data))
    #println("Stress data content: ", stress_data)
    println("stress data inputs:")
    println("Vertical stress gradient: ", stress_data["vertical_stress_gradient_psi_ft"])
    println("Minimum horizontal stress gradient: ", stress_data["min_horizontal_stress_gradient_psi_ft"])
    println("Maximum horizontal stress gradient: ", stress_data["max_horizontal_stress_gradient_psi_ft"])
    println("Maximum horizontal stress direction: ", stress_data["max_horizontal_stress_direction"])
    println("Reference depth: ", stress_data["reference_depth_ft"])
    println("Initial reservoir pressure gradient: ", stress_data["initial_reservoir_pressure_gradient_psi_ft"])

    
    # extract parameters from stress data
    # get maximum horizontal stress direction
    SHdir = stress_data["max_horizontal_stress_direction"]
    # get reference depth
    dpth = stress_data["reference_depth_ft"]

    # calculate stress at reference depth
    sig = [
        stress_data["vertical_stress_gradient_psi_ft"],
        stress_data["min_horizontal_stress_gradient_psi_ft"],
        stress_data["max_horizontal_stress_gradient_psi_ft"]
    ] .* dpth
    println("stress data after multiplying by depth:")
    println("Maximum horizontal stress direction (SHdir): ", SHdir)
    println("Reference depth (dpth): ", dpth)
    println("Vertical stress gradient (sig[1]): ", sig[1])
    println("Minimum horizontal stress gradient (sig[2]): ", sig[2])
    println("Maximum horizontal stress gradient (sig[3]): ", sig[3])
    

    #println("Scaled stress values (sig): ", sig)

    # Native pore pressure --> we get it by scaling the pore pressure gradient by depth
    pp0 = stress_data["initial_reservoir_pressure_gradient_psi_ft"] * dpth
    #println("Native pore pressure (pp0): ", pp0)

    
    nu = 0.5  #Poisson's ratio (will stay hardcoded for now)
    biot = 1.0 # Biot coefficient (hardcoded for testing and will be parsed through the web portal)


    # Parse CLI arguments
    aphi_value, min_hor_stress = parse_cli_args()

    if aphi_value !== nothing
        println("A-Phi value provided. Calculating SH and Sh from A-Phi...")
        if min_hor_stress !== nothing
            println("Minimum horizontal stress provided: $min_hor_stress")
            geomechanics_calculations_inputs = (
                sig, 0.00, pp0, thf, dips, SHdir, zeros(length(thf)), muf, biot, nu, aphi_value, min_hor_stress
            )
        else
            geomechanics_calculations_inputs = (
                sig, 0.00, pp0, thf, dips, SHdir, zeros(length(thf)), muf, biot, nu, aphi_value
            )
        end
    else
        println("No A-Phi value provided. Calculating SH and Sh from mu...")
        geomechanics_calculations_inputs = (
            sig, 0.00, pp0, thf, dips, SHdir, zeros(length(thf)), muf, biot, nu
        )
    end


    #println("running  mohrs_3D function...")
    det_geomechanic_results = mohrs_3D(geomechanics_calculations_inputs)


    ppfail = det_geomechanic_results.outs["ppfail"]
    ppfail[ppfail .< 0] .= 0.0
    println("ppfail results: ", ppfail)


    #println("Storing results...")
    det_geomech_results_to_json(det_geomechanic_results, "src/output/deterministic_geomechanics_results.json")
    println("Geomechanics driver's results stored successfully!")
end

geomechanics_driver()

end # End of mmodule

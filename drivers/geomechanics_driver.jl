module GeomechanicsDriver

using JSON


include("src/core/DeterministicGeomechanicsCalculations.jl") 
using .DeterministicGeomechanicsCalculations


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

    println("Fault data structure: ", typeof(fault_data)) # for debugging
    println("Fault data content: ", fault_data)

    # similarly, load stress data 
    stress_data = load_json_data("step1_stress_data_output.json")
    if stress_data === nothing
        println("Failed to load stress data. Exiting.")
        return
    end

    println("Stress data structure: ", typeof(stress_data))
    println("Stress data content: ", stress_data)

    
    println("Extracting fault parameters...")

    # get fault strikes 
    thf = [fault["strike"] for fault in fault_data]
    # get fault dips
    dips = [fault["dip"] for fault in fault_data]
    # get fault friction coefficients (the same one is used for all faults in the calculations)
    muf = [fault["friction"] for fault in fault_data]
    println("Fault strikes (thf): ", thf)
    println("Fault dips: ", dips)
    println("Fault friction coefficients (muf): ", muf)

    # extract parameters from stress data
    # get maximum horizontal stress direction
    SHdir = stress_data["max_horizontal_stress_direction"]
    # get reference depth
    dpth = stress_data["reference_depth_ft"]

    # multiply these stress values by depth
    sig = [
        stress_data["vertical_stress_gradient_psi_ft"],
        stress_data["min_horizontal_stress_gradient_psi_ft"],
        stress_data["max_horizontal_stress_gradient_psi_ft"]
    ] .* dpth
    println("Maximum horizontal stress direction (SHdir): ", SHdir)
    println("Calculation depth (dpth): ", dpth)

    println("Scaled stress values (sig): ", sig)

    # Native pore pressure --> we get it by scaling the pore pressure gradient by depth
    pp0 = stress_data["initial_reservoir_pressure_gradient_psi_ft"] * dpth
    println("Native pore pressure (pp0): ", pp0)

    # Hard-coded parameters (these will be input from the user parsed trough the web tools portal)
    nu = 0.25  #Poisson's ratio
    biot = 1.0 # Biot coefficient

    # the input array for mohrs_3D function
    inputCell = (sig, 0.00, pp0, thf, dips, SHdir, zeros(length(thf)), muf, biot, nu)

    println("Calling mohrs_3D function with input tuple ...")
    results = mohrs_3D(inputCell)

    println("Storing results...")
    store_results_to_json(results, "deterministic_geomechanics_results.json")
    println("Geomechanics driver's results stored successfully!")
end

geomechanics_driver()

end # End of mmodule

module GeomechanicsDriver

using JSON
using ArgParse


include("src/core/DeterministicGeomechanicsCalculations.jl") 
using .DeterministicGeomechanicsCalculations


function parse_cli_args()
    s = ArgParseSettings()

    # A-Phi optional argument
    @add_arg_table s begin
        "--aphi"
        help = "(Optional) A-Phi stress model calculates horizontal stresses based on stress criticality assumption."
        arg_type = Float64
        required = false 
    end

    # parse the arguments
    args = parse_args(s)


    # If the user provides --aphi, make sure it's a Float64
    if haskey(args, :aphi)
        try
            aphi_value = convert(Float64, args["aphi"]) # Ensures the input is Float64
            println("A-Phi value provided: ", aphi_value)
        catch e
            error("Error: A-Phi value must be a valid Float64")
        end
    else
        println("Skipping A-Phi value calculation...")
    end
    
    return args
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

    # similarly, load stress data 
    stress_data = load_json_data("step1_stress_data_output.json")
    if stress_data === nothing
        println("Failed to load stress data. Exiting.")
        return
    end

    #println("Stress data structure: ", typeof(stress_data))
    #println("Stress data content: ", stress_data)

    
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
    #println("Maximum horizontal stress direction (SHdir): ", SHdir)
    #println("Calculation depth (dpth): ", dpth)

    #println("Scaled stress values (sig): ", sig)

    # Native pore pressure --> we get it by scaling the pore pressure gradient by depth
    pp0 = stress_data["initial_reservoir_pressure_gradient_psi_ft"] * dpth
    #println("Native pore pressure (pp0): ", pp0)

    
    nu = 0.5  #Poisson's ratio (will stay hardcoded for now)
    biot = 1.0 # Biot coefficient (hardcoded for testing and will be parsed through the web portal)


    # flag to check if user provided A-Phi value
    # hardcoded for now
    aphi_check = false
    #aphi_value = 2.0 # for when we test it for true
    aphi_value = nothing
    
    # check for --aphi cli argument
    aphi = parse_cli_args()
    if haskey(aphi, :aphi)
        aphi_value = aphi[:aphi]
        aphi_check = true
    end


    if aphi_check == false
        println("No A-Phi value provided. Calculating SH and Sh from mu...")
        # (reminder) thf = fault strikes
        geomechanics_calculations_inputs = (sig, 0.00, pp0, thf, dips, SHdir, zeros(length(thf)), muf, biot, nu)
    else
        # in this case, the horizontal stress components in sig are being ignored
        # since we will recalculate them using the A-Phi model
        println("A-Phi value provided. Calculating SH and Sh from A-Phi...")
        geomechanics_calculations_inputs = (sig, 0.00, pp0, thf, dips, SHdir, zeros(length(thf)), muf, biot, nu, aphi_value)
    end

    # create grid for plotting

    println("Calling mohrs_3D function with input tuple ...")
    det_geomechanic_results = mohrs_3D(geomechanics_calculations_inputs)


    ppfail = det_geomechanic_results.outs["ppfail"]
    ppfail[ppfail .< 0] .= 0.0
    println("ppfail results: ", ppfail)


    println("Storing results...")
    det_geomech_results_to_json(det_geomechanic_results, "deterministic_geomechanics_results.json")
    println("Geomechanics driver's results stored successfully!")
end

geomechanics_driver()

end # End of mmodule

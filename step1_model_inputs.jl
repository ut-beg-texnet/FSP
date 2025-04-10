using JSON
using ArgParse


include("TexNetWebToolLauncherHelperJulia.jl")

using .TexNetWebToolLauncherHelperJulia


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--scratchPath"
            help="Path to the scratch directory"
            required=true
            arg_type=String
    end

    return parse_args(s)
end


function main()

    # read sratch path from the first argument the script was run with
    args = parse_commandline()

    scratchPath = args["scratchPath"]
    
    # instantiate the helper object
    helper = TexNetWebToolLaunchHelperJulia("scratchPath")



    # Read input JSON file
    input_file = "input/input_data.json"  # Specify the input JSON file path
    input_data = JSON.parsefile(input_file)

    # Debugging statements to check the values and types
    println("Parsed input data:")
    println(input_data)

    # Check stress state values
    println("Stress State Values:")
    # ... [existing debug prints]

    # Determine stress state based on input data
    stress_state = Dict(
        "reference_depth" => input_data["stress_state"]["reference_depth"],
        "vertical_stress" => input_data["stress_state"]["vertical_stress"],
        "min_horizontal_stress" => input_data["stress_state"]["min_horizontal_stress"] === nothing ? 0.0 : input_data["stress_state"]["min_horizontal_stress"],
        "max_horizontal_stress" => input_data["stress_state"]["max_horizontal_stress"] === nothing ? 0.0 : input_data["stress_state"]["max_horizontal_stress"],
        "pore_pressure" => input_data["stress_state"]["pore_pressure"],
        "max_stress_azimuth" => input_data["stress_state"]["max_stress_azimuth"],
        "aphi_value" => input_data["stress_state"]["aphi_value"] === nothing ? 0.0 : input_data["stress_state"]["aphi_value"],
        "model_type" => ""
    )

    
    if !haskey(input_data, "uncertainties")
        println("Error: No uncertainties provided in the input data.")
    end


    # Determine the model type based on the parameters provided
    if input_data["stress_state"]["aphi_value"] !== nothing
        if input_data["stress_state"]["min_horizontal_stress"] !== nothing
            stress_state["model_type"] = "aphi_min"
        else
            stress_state["model_type"] = "aphi_no_min"
        end
    else
        stress_state["model_type"] = "gradients"
    end

    # check for model_parameters and set default values any of its keys are missing
    if !haskey(input_data, "model_parameters")
        input_data["model_parameters"] = Dict(
            "fluid_compressibility" => 3.6e-10,
            "rock_compressibility" => 1.08e-10,
            "fluid_density" => 1000.0,
            "dynamic_viscosity" => 0.0008,
            "aphi_value" => 0.0
        )
    else
        if !haskey(input_data["model_parameters"], "fluid_compressibility")
            input_data["model_parameters"]["fluid_compressibility"] = 3.6e-10
        end
        if !haskey(input_data["model_parameters"], "rock_compressibility")
            input_data["model_parameters"]["rock_compressibility"] = 1.08e-10
        end
        if !haskey(input_data["model_parameters"], "fluid_density")
            input_data["model_parameters"]["fluid_density"] = 1000.0
        end
        if !haskey(input_data["model_parameters"], "dynamic_viscosity")
            input_data["model_parameters"]["dynamic_viscosity"] = 0.0008
        end
        if !haskey(input_data["model_parameters"], "aphi_value")
            input_data["model_parameters"]["aphi_value"] = 0.0
        end
    end

    # read the first occurence of friction_coefficient from faults and add that value to all faults
    if haskey(input_data, "faults")
        if length(input_data["faults"]) > 0
            friction_coefficient = input_data["faults"][1]["friction_coefficient"]
            for fault in input_data["faults"]
                fault["friction_coefficient"] = friction_coefficient
            end
        end
    end

    # the fault_id should be unique and will be an integer (we might want to make this accept strings as well)
    for (i, fault) in enumerate(input_data["faults"])
        fault["fault_id"] = i
    end

    
    hydrology_data = input_data["hydrology"]


    # output JSON structure
    output_data = Dict(
        "stress_state" => stress_state,
        "faults" => input_data["faults"],
        "injection_wells" => input_data["injection_wells"],
        "model_parameters" => input_data["model_parameters"],
        "uncertainties" => input_data["uncertainties"],
        "hydrology" => hydrology_data
    )
    #=
    # Process faults 
    for fault in input_data["faults"]
        fault_data = Dict(
            "fault_id" => fault["fault_id"],
            "strike" => fault["strike"],
            "dip" => fault["dip"],
            "friction_coefficient" => fault["friction_coefficient"],
            "length_km" => fault["length_km"],
            "easting" => fault["easting"],  
            "northing" => fault["northing"]  
        )
        push!(output_data["faults"], fault_data)
    end
    =#

    # Ensure output directory exists
    output_dir = "output"
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Write output JSON file
    output_file = "output/step1_output.json"  # Specify the output JSON file path
    try
        println("Writing output to $output_file...")
        open(output_file, "w") do f
            JSON.print(f, output_data, 4)  # Indent with 4 spaces
        end
        println("Output written to $output_file")
    catch e
        println("An error occurred while writing the output file: ", e)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
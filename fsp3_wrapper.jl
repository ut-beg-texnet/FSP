module FSP3Wrapper

# Include the core modules
include("core/bill_pfront.jl")
include("core/geomechanics_model.jl")
include("core/hydrology_calculations.jl")
include("core/utilities.jl")

include("graphs/julia_fsp_graphs.jl")

# Include the driver scripts
include("model_inputs_process.jl")
include("deterministic_geomechanics_process.jl")
include("probabilistic_geomechanics_process.jl")
include("deterministic_hydrology_process.jl")
include("probabilistic_hydrology_process.jl")
include("summary_process.jl")


function main(args)
    if length(args) != 1
        #println("Usage: fsp3_wrapper.jl <driver_file>")
        #println("Available driver files: model_inputs, deterministic_geomechanics, probabilistic_geomechanics, deterministic_hydrology, probabilistic_hydrology, summary")
        return
    end

    driver_file = args[1]
    scratch_path = args[2]

    driver_args = [scratch_path]

    if driver_file == "model_inputs"
        main(driver_args)
    elseif driver_file == "deterministic_geomechanics"
        main(driver_args)
    elseif driver_file == "probabilistic_geomechanics"
        main(driver_args)
    elseif driver_file == "deterministic_hydrology"
        main(driver_args)
    elseif driver_file == "probabilistic_hydrology"
        main(driver_args)
    elseif driver_file == "summary"
        main(driver_args)
    else
        #println("Invalid driver file: $driver_file")
        #println("Available driver files: model_inputs, deterministic_geomechanics, probabilistic_geomechanics, deterministic_hydrology, probabilistic_hydrology, summary")
    end

    main(ARGS)
    

end # module

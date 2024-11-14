module ProbabilisticGeomechanicsCalculations

using Random
using Base.Threads  

include("DeterministicGeomechanicsCalculations.jl")
using .DeterministicGeomechanicsCalculations: mohrs_3D

"""
# Inputs
- f: The function to run each Monte Carlo iteration on (calc_ppfail).
- in0: Nominal input values (center points of distribution).
- inSig: A vector of ranges for each element in in0 (uncertainty in input parameters)
- nruns: Number of iterations.

# Outputs
- out: Vector of output results from each Monte Carlo run.
- inj: Vector of input variations used in each run.
"""
function monte_carlo(f::Function, in0::Vector, inSig::Vector, nruns::Int)
    mc_out = Vector{Any}(undef, nruns)  # Stores outputs for each run
    # Stores the randomized inputs for each run
    inputs = Vector{Vector}(undef, nruns)  

    @threads for jj in 1:nruns
        # Generate perturbed inputs for this run
        perturbed_inputs = [in0[k] .+ inSig[k] .* (2 * rand() - 1) for k in 1:length(in0)]
        inputs[jj] = perturbed_inputs

        # Runs the function with those inputs
        mc_out[jj] = f(perturbed_inputs)
    end

    return mc_out, inputs
end

end # module

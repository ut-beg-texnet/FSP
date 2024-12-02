"""
FSP 3.0 Monte Carlo Simulation for Fault Stability Analysis
========================================================

This script performs Monte Carlo simulation to analyze fault stability under uncertainty.
It takes deterministic inputs and applies random variations based on specified uncertainties
to calculate a distribution of slip pressures for each fault.

Workflow:
1. Read input JSON containing:
   - Base stress state (gradients or A-phi model)
   - Fault geometries and properties
   - Uncertainties for each parameter

2. For each Monte Carlo simulation:
   a. Create modified stress state by applying random variations to:
      - Vertical stress gradient
      - Pore pressure gradient
      - Max horizontal stress azimuth
      - A-phi value (if using A-phi model)
      - Min horizontal stress (if using gradients or aphi_min)
      - Max horizontal stress (if using gradients)
   
   b. Create modified faults by applying random variations to:
      - Strike angles
      - Dip angles
      - Friction coefficients
   
   c. Calculate absolute stresses using appropriate model type
   
   d. Calculate slip pressure for each fault using modified parameters
   
3. Calculate statistics for each fault:
   - Mean slip pressure
   - Standard deviation
   - P10, P50, P90 percentiles
   - Min and max values

Usage:
    julia step3_monte_carlo.jl 
        --input-json path/to/input.json 
        --output-json path/to/output.json 
        --num-simulations 1000 
        --seed 42

Input JSON format:
{
    "stress_state": {...},
    "faults": [...],
    "uncertainties": {
        "vertical_stress_gradient": float,
        "initial_pore_pressure_gradient": float,
        "max_stress_azimuth": float,
        "aphi_value": float,  # only for aphi models
        "strike_angles": float,
        "dip_angles": float,
        "friction_coefficient": float
    }
}
"""

using JSON
using ArgParse
using Distributions
using Statistics
using Random
using LinearAlgebra

include("core/geomechanics_model.jl")
include("step2_deterministic_geomechanics.jl")

using .GeomechanicsModel

"""
Parse command line arguments
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input-json"
            help = "Input JSON file from previous step (step 2)"
            required = true
        "--output-json"
            help = "Output JSON file path"
            required = true
        "--num-simulations"
            help = "Number of Monte Carlo simulations to run"
            arg_type = Int
            default = 1000
        "--seed"
            help = "Random seed for reproducibility"
            arg_type = Int
            default = 42
    end

    return parse_args(s)
end


"""
Run Monte Carlo simulation
"""
function run_monte_carlo(input_data::Dict, n_sims::Int)
    # Extract base values and uncertainties
    stress_state = input_data["stress_state"]
    faults = input_data["faults"]
    uncertainties = get(input_data, "uncertainties", Dict())
    
    # Initialize results storage - store both pore pressure and slip pressure for each fault
    fault_results = [Dict(
        "pore_pressures" => Float64[],
        "slip_pressures" => Float64[]
    ) for _ in 1:length(faults)]
    
    # Run simulations
    for i in 1:n_sims
        # Create modified stress state for this simulation
        sim_stress = copy(stress_state)
        reference_depth = stress_state["reference_depth"]
        
        param_mapping = Dict()

        # check stress model type and map uncertainty parameters to stress state parameters
        if stress_state["model_type"] == "gradients"
            param_mapping = Dict(
                "vertical_stress_gradient" => "vertical_stress",
                "initial_pore_pressure_gradient" => "pore_pressure",
                "max_stress_azimuth" => "max_stress_azimuth",
                "max_horizontal_stress" => "max_horizontal_stress",
                "min_horizontal_stress" => "min_horizontal_stress"
            )
        elseif stress_state["model_type"] == "aphi_min"
            param_mapping = Dict(
                "vertical_stress_gradient" => "vertical_stress",
                "initial_pore_pressure_gradient" => "pore_pressure",
                "max_stress_azimuth" => "max_stress_azimuth",
                "aphi_value" => "aphi_value",
                "min_horizontal_stress" => "min_horizontal_stress"
            )
        elseif stress_state["model_type"] == "aphi_no_min"
            param_mapping = Dict(
                "vertical_stress_gradient" => "vertical_stress",
                "initial_pore_pressure_gradient" => "pore_pressure",
                "max_stress_azimuth" => "max_stress_azimuth",
                "aphi_value" => "aphi_value"
            )
        else
            error("Invalid stress model type: $(stress_state["model_type"])")
        end
        
        # Modify stress parameters if uncertainties exist
        for (uncertainty_param, stress_param) in param_mapping
            if haskey(uncertainties, uncertainty_param)
                base_value = stress_state[stress_param]
                if !isnothing(base_value)
                    uncertainty = uncertainties[uncertainty_param]
                    sim_stress[stress_param] = base_value + rand(Uniform(-uncertainty, uncertainty))
                end
            end
        end


        # Calculate absolute stresses using the appropriate model type
        stress_state_obj, initial_pressure = calculate_absolute_stresses(sim_stress, faults)

        # Update stress state with absolute values
        sim_stress["vertical_stress"] = stress_state_obj.principal_stresses[1]
        sim_stress["max_horizontal_stress"] = stress_state_obj.principal_stresses[3]
        sim_stress["min_horizontal_stress"] = stress_state_obj.principal_stresses[2]
        sim_stress["pore_pressure"] = initial_pressure

        
        # Create modified faults for this simulation
        sim_faults = []
        for fault in faults
            sim_fault = copy(fault)
            
            # Map uncertainty parameters to fault parameters
            fault_param_mapping = Dict(
                "strike_angles" => "strike",
                "dip_angles" => "dip",
                "friction_coefficient" => "friction_coefficient"
            )
            
            # Modify fault parameters if uncertainties exist
            for (uncertainty_param, fault_param) in fault_param_mapping
                if haskey(uncertainties, uncertainty_param)
                    uncertainty = uncertainties[uncertainty_param]
                    sim_fault[fault_param] = fault[fault_param] + rand(Uniform(-uncertainty, uncertainty))
                end
            end
            
            push!(sim_faults, sim_fault)
        end
        
        # Run geomechanical analysis with modified parameters
        stress_state_obj = StressState(
            [sim_stress["vertical_stress"], 
            sim_stress["min_horizontal_stress"], 
            sim_stress["max_horizontal_stress"]],
            sim_stress["max_stress_azimuth"]
        )
        
        # Calculate slip pressure for each fault
        for (fault_idx, fault) in enumerate(sim_faults)
            println("\nFault $fault_idx:")
            println("Strike: $(fault["strike"])")
            println("Dip: $(fault["dip"])")
            println("Friction: $(fault["friction_coefficient"])")
            println("Initial Pore Pressure: $(sim_stress["pore_pressure"])")
            
            sig_normal, tau_normal, s11, s22, s33, s12, n1, n2 = calculate_fault_effective_stresses(
                fault["strike"], 
                fault["dip"], 
                stress_state_obj, 
                sim_stress["pore_pressure"], 
                0.0  # dp is 0 for Monte Carlo
            )
            
            println("Normal Stress: $sig_normal")
            println("Shear Stress: $tau_normal")
            println("s11: $s11, s22: $s22, s33: $s33, s12: $s12")
            println("n1: $n1, n2: $n2")
            
            slip_pressure = calculate_slip_pressure(
                sig_normal,
                tau_normal,
                fault["friction_coefficient"],
                sim_stress["pore_pressure"],
                1.0,
                0.5,
                0.0,
                s11,
                s22,
                s33,
                s12,
                n1,
                n2
            )
            
            println("Slip Pressure: $slip_pressure")
            
            # Store both pore pressure and slip pressure
            push!(fault_results[fault_idx]["pore_pressures"], sim_stress["pore_pressure"])
            push!(fault_results[fault_idx]["slip_pressures"], slip_pressure)
        end
    end
    
    return fault_results
end

"""
Calculate statistics from Monte Carlo results
"""
function calculate_statistics(fault_results::Vector)
    stats = Dict()
    
    # For each fault
    for (i, fault_data) in enumerate(fault_results)
        stats["fault_$i"] = Dict(
            "slip_pressure" => Dict(
                "mean" => mean(fault_data["slip_pressures"]), # mean value
                "std" => std(fault_data["slip_pressures"]), # standard deviation
                "p50" => quantile(fault_data["slip_pressures"], 0.50), # median
                "min" => minimum(fault_data["slip_pressures"]), # min slip pressure
                "max" => maximum(fault_data["slip_pressures"]) # max slip pressure
            ),
            "iterations" => [
                Dict(
                    "pore_pressure" => p,
                    "slip_pressure" => s
                ) for (p, s) in zip(fault_data["pore_pressures"], fault_data["slip_pressures"])
            ]
        )
    end
    
    return stats
end

"""
Main function
"""
function main()
    # Parse command line arguments
    args = parse_commandline()
    
    # Set random seed for reproducibility
    Random.seed!(args["seed"])
    
    # Read input JSON from step 2
    input_data = JSON.parsefile(args["input-json"])
    
    # Get uncertainties from input JSON
    uncertainties = get(input_data, "uncertainties", Dict())
    
    println("\n=== Running Monte Carlo Simulation ===")
    println("Number of simulations: $(args["num-simulations"])")
    println("Using uncertainties:")
    for (param, value) in uncertainties
        println("  $param: ±$value")
    end
    
    # Run Monte Carlo simulation
    results = run_monte_carlo(input_data, args["num-simulations"])
    
    # Calculate statistics
    stats = calculate_statistics(results)
    
    # Add Monte Carlo results to input data
    input_data["monte_carlo_results"] = stats
    
    # Update metadata
    input_data["metadata"] = Dict(
        "step" => "monte_carlo",
        "input_file" => args["input-json"],
        "description" => "Monte Carlo simulation results",
        "num_simulations" => args["num-simulations"],
        "seed" => args["seed"]
    )
    
    # Write combined results to output JSON
    open(args["output-json"], "w") do f
        JSON.print(f, input_data, 4)  # indent with 4 spaces
    end
    
    println("\nResults written to $(args["output-json"])")
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

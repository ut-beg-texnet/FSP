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
using Distributed
using SharedArrays
using Plots
using ProgressMeter
gr()

include("core/geomechanics_model.jl")
include("step2_deterministic_geomechanics.jl")

using .GeomechanicsModel
using .GeomechanicsDriver

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


function plot_cdf(pps_to_slip::SharedArray{Float64, 2}, output_path::String)
    # Generate y-axis tick labels with '%' symbol
    y_ticks = 0:10:100
    y_labels = [string(x, "%") for x in y_ticks]

    # Initialize the plot
    p = plot(
        xlabel="Δ Pore Pressure to Slip (Psi)", 
        ylabel="Probability of Fault Slip (%)", 
        title="CDF of Pore Pressure to Slip for Each Fault",
        legend=:outerright,  # Move legend to the outer right
        linewidth=2,
        yticks=(y_ticks, y_labels),  # Add '%' symbol to y-axis tick labels
        legendfontsize=8  # Adjust legend font size if needed
    )

    # Loop through each fault and add its CDF to the plot
    for fault_idx in 1:size(pps_to_slip, 1)
        # Extract pore pressure data for the fault
        fault_data = pps_to_slip[fault_idx, :]

        # Sort the pore pressure values
        sorted_data = sort(fault_data)

        # Generate cumulative probabilities and convert to percentage
        cumulative_probabilities = range(0, stop=100, length=length(sorted_data))

        # Add the fault's CDF to the plot
        plot!(
            sorted_data, cumulative_probabilities,
            label="Fault $fault_idx"
        )
    end

    # Save the combined plot as a PNG
    combined_plot_file = joinpath(output_path, "prob_geo_cdf.png")
    savefig(p, combined_plot_file)
    println("Prob Geomechanics CDF plot saved at: $combined_plot_file")
    # display the figure on the screen
    display(p)
end


"""
Run Monte Carlo simulation
"""
function run_monte_carlo(input_data::Dict, n_sims::Int)
    # Extract base values and uncertainties
    stress_state = input_data["stress_state"]
    faults = input_data["faults"]
    uncertainties = get(input_data, "uncertainties", Dict())
    
    
    # Use SharedArray to store results (better for parallel processing)
    # rows are faults, columns are the pore pressure to slip for each simulation
    pps_to_slip = SharedArray{Float64, 2}(length(faults), n_sims)

    stress_param_mapping = Dict()

    # check stress model type and map uncertainty parameters to stress state parameters
    if stress_state["model_type"] == "gradients"
        stress_param_mapping = Dict(
            "vertical_stress_gradient" => "vertical_stress",
            "initial_pore_pressure_gradient" => "pore_pressure",
            "max_stress_azimuth" => "max_stress_azimuth",
            "max_horizontal_stress" => "max_horizontal_stress",
            "min_horizontal_stress" => "min_horizontal_stress"
        )
    elseif stress_state["model_type"] == "aphi_min"
        stress_param_mapping = Dict(
            "vertical_stress_gradient" => "vertical_stress",
            "initial_pore_pressure_gradient" => "pore_pressure",
            "max_stress_azimuth" => "max_stress_azimuth",
            "aphi_value" => "aphi_value",
            "min_horizontal_stress" => "min_horizontal_stress"
        )
    elseif stress_state["model_type"] == "aphi_no_min"
        stress_param_mapping = Dict(
            "vertical_stress_gradient" => "vertical_stress",
            "initial_pore_pressure_gradient" => "pore_pressure",
            "max_stress_azimuth" => "max_stress_azimuth",
            "aphi_value" => "aphi_value"
        )
    else
        error("Invalid stress model type: $(stress_state["model_type"])")
    end

    # Map uncertainty parameters to fault parameters
    fault_param_mapping = Dict(
        "strike_angles" => "strike",
        "dip_angles" => "dip",
        "friction_coefficient" => "friction_coefficient"
    )


    # Initialize progress meter
    progress = Progress(n_sims, 1, "Monte Carlo Simulations")
    

    # Run simulations
    Threads.@threads for i in 1:n_sims
        randomized_stress_values = Dict{String, Float64}()
        for (uncertainty_param, stress_param) in stress_param_mapping
            if haskey(uncertainties, uncertainty_param)
                base_value = stress_state[stress_param]
                if !isnothing(base_value)
                    uncertainty = uncertainties[uncertainty_param]
                    # Randomize stress values for this thread
                    randomized_stress_values[stress_param] = base_value + rand(Uniform(-uncertainty, uncertainty))
                end
            end
        end

        # Create modified stress state for this simulation
        sim_stress = deepcopy(stress_state)

        # Assign thread-local randomized values to stress parameters
        for (stress_param, randomized_value) in randomized_stress_values
            sim_stress[stress_param] = randomized_value
        end

        # Calculate absolute stresses using the appropriate model type
        stress_state_obj, initial_pressure = calculate_absolute_stresses(sim_stress, faults)

        # Update stress state with absolute values
        sim_stress["vertical_stress"] = stress_state_obj.principal_stresses[1]
        sim_stress["max_horizontal_stress"] = stress_state_obj.principal_stresses[3]
        sim_stress["min_horizontal_stress"] = stress_state_obj.principal_stresses[2]
        sim_stress["pore_pressure"] = initial_pressure

        
        # Create modified faults for this simulation
        sim_faults = deepcopy(faults)
        for fault in sim_faults
            for (uncertainty_param, fault_param) in fault_param_mapping
                if haskey(uncertainties, uncertainty_param)
                    uncertainty = uncertainties[uncertainty_param]
                    fault[fault_param] += rand(Uniform(-uncertainty, uncertainty))
                end
            end
        end
        
        # Run geomechanical analysis with modified parameters
        stress_state_obj = StressState(
            [sim_stress["vertical_stress"], 
            sim_stress["min_horizontal_stress"], 
            sim_stress["max_horizontal_stress"]],
            sim_stress["max_stress_azimuth"]
        )

        # Thread-local buffer to store results before writing to the shared array
        local_results = zeros(Float64, length(faults))
        
        # Calculate slip pressure for each fault
        for (fault_idx, fault) in enumerate(sim_faults)
            #println("\nFault $fault_idx:")
            #println("Strike: $(fault["strike"])")
            #println("Dip: $(fault["dip"])")
            #println("Friction: $(fault["friction_coefficient"])")
            #println("Initial Pore Pressure: $(sim_stress["pore_pressure"])")
            
            sig_normal, tau_normal, s11, s22, s33, s12, n1, n2 = calculate_fault_effective_stresses(
                fault["strike"], 
                fault["dip"], 
                stress_state_obj, 
                sim_stress["pore_pressure"], 
                0.0  # dp is 0 for Monte Carlo
            )
            
            #println("Normal Stress: $sig_normal")
            #println("Shear Stress: $tau_normal")
            #println("s11: $s11, s22: $s22, s33: $s33, s12: $s12")
            #println("n1: $n1, n2: $n2")
            
            pp_to_slip = calculate_slip_pressure(
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
            
            #println("Fault $fault_idx pore pressure to slip: $pp_to_slip")
            
            # Store results in shared arrays
            #pore_pressures[fault_idx, i] = sim_stress["pore_pressure"]
            local_results[fault_idx] = pp_to_slip
        end
        # Write thread-local results back to the shared array
        pps_to_slip[:, i] = local_results

        # Update progress meter
        next!(progress)

    end
    
    return pps_to_slip
end

"""
Calculate statistics from Monte Carlo results
"""
function calculate_statistics(pps_to_slip::SharedArray{Float64, 2})
    stats = Dict()

    # For each fault (row in SharedArray)
    for fault_idx in 1:size(pps_to_slip, 1)
        fault_data = pps_to_slip[fault_idx, :]  # Extract the data for this fault (row)

        stats["fault_$fault_idx"] = Dict(
            "slip_pressure" => Dict(
                "mean" => mean(fault_data),         # Mean value of slip pressures
                "std" => std(fault_data),           # Standard deviation
                "p50" => quantile(fault_data, 0.50), # Median (50th percentile)
                "min" => minimum(fault_data),       # Minimum slip pressure
                "max" => maximum(fault_data)        # Maximum slip pressure
            ),
            "iterations" => [                     # Iteration results
                Dict("slip_pressure" => sp) for sp in fault_data
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
    mc_pp_results = run_monte_carlo(input_data, args["num-simulations"])
    
    # Calculate statistics
    stats = calculate_statistics(mc_pp_results)
    
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

    output_plot_path = dirname(args["output-json"])
    println("\nGenerating CDF plots for each fault...")
    plot_cdf(mc_pp_results, output_plot_path)
    println("\nCDF plots saved in $(output_plot_path)")
    
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

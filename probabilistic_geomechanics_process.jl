"""
FSP 3.0 Monte Carlo Simulation for Fault Stability Analysis
========================================================

This script performs Monte Carlo simulation to analyze fault stability under uncertainty.
It takes deterministic inputs and applies random variations based on specified uncertainties
to calculate a distribution of slip pressures for each fault.
"""

# Set encoding to UTF-8 to avoid garbled text in console output
ENV["JULIA_UNICODE_INPUT"] = "true"
ENV["JULIA_UNICODE_OUTPUT"] = "true"

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
using DataFrames
using CSV
using Base.Threads
using PrettyTables
using StatsBase  # Add StatsBase for ecdf function
gr()

include("core/geomechanics_model.jl")
#include("step2_deterministic_geomechanics.jl")
include("deterministic_geomechanics_process.jl")
include("TexNetWebToolLauncherHelperJulia.jl")
include("graphs/julia_fsp_graphs.jl")



using .GeomechanicsModel
using .GeomechanicsDriver
using .TexNetWebToolLauncherHelperJulia
using .JuliaFSPGraphs

const ARGS_FILE_NAME = "args.json"
const RESULTS_FILE_NAME = "results.json"






"""
Runs Monte Carlo simulation
"""
function run_monte_carlo(stress_inputs::Dict, fault_inputs::DataFrame, uncertainties::Dict, n_sims::Int, stress_model_type::String, random_seed=nothing)
    n_faults = size(fault_inputs, 1)
    
    # Use SharedArray instead of regular Array for parallelization
    pps_to_slip = SharedArray{Float64}(n_faults, n_sims)
    
    # Initialize mapping dictionaries for parameter randomization
    stress_param_mapping = Dict()
    
    # We set the parameter mapping based on the stress model type
    if stress_model_type == "gradients" || stress_model_type == "all_gradients"
        stress_param_mapping = Dict(
            "vertical_stress_gradient_uncertainty" => "vertical_stress",
            "initial_pore_pressure_gradient_uncertainty" => "pore_pressure",
            "max_stress_azimuth_uncertainty" => "max_stress_azimuth",
            "max_horizontal_stress_uncertainty" => "max_horizontal_stress",
            "min_horizontal_stress_uncertainty" => "min_horizontal_stress"
        )
    elseif stress_model_type == "aphi_model" || stress_model_type == "aphi_min"
        stress_param_mapping = Dict(
            "vertical_stress_gradient_uncertainty" => "vertical_stress",
            "initial_pore_pressure_gradient_uncertainty" => "pore_pressure",
            "max_stress_azimuth_uncertainty" => "max_stress_azimuth",
            "aphi_value_uncertainty" => "aphi_value",
            "min_horizontal_stress_uncertainty" => "min_horizontal_stress"
        )
    elseif stress_model_type == "aphi_model_no_min" || stress_model_type == "aphi_no_min"
        stress_param_mapping = Dict(
            "vertical_stress_gradient_uncertainty" => "vertical_stress",
            "initial_pore_pressure_gradient_uncertainty" => "pore_pressure",
            "max_stress_azimuth_uncertainty" => "max_stress_azimuth",
            "aphi_value_uncertainty" => "aphi_value"
        )
    end
    
    # Map uncertainty parameters to fault parameters
    fault_param_mapping = Dict(
        "strike_angles_uncertainty" => "Strike",
        "dip_angles_uncertainty" => "Dip",
        "friction_coefficient_uncertainty" => "FrictionCoefficient"
    )
    
    # Pre-convert DataFrame to array of dictionaries
    fault_dicts = [Dict(name => row[name] for name in names(fault_inputs)) for row in eachrow(fault_inputs)]
    
    # Store actual fault IDs if present in the input dataframe
    actual_fault_ids = []
    if "FaultID" in names(fault_inputs)
        actual_fault_ids = string.(fault_inputs.FaultID)
    elseif "ID" in names(fault_inputs)
        actual_fault_ids = string.(fault_inputs.ID)
    else
        # If no ID column exists, use sequential numbers as strings
        actual_fault_ids = string.(1:n_faults)
    end
    
    # Set a base seed for deterministic results
    if isnothing(random_seed)
        # Use time-based seed for different results each run
        base_seed = time_ns()
    else
        # Use user-provided seed for reproducible results
        base_seed = random_seed  
    end
    
    # Create dictionaries to store randomized parameter values for each simulation
    # Structure: param_name => [value_sim_1, value_sim_2, ...]
    stress_param_values = Dict{String, Vector{Float64}}()
    for param in values(stress_param_mapping)
        stress_param_values[param] = Vector{Float64}(undef, n_sims)
    end
    
    # Structure: fault_idx => param_name => [value_sim_1, value_sim_2, ...]
    fault_param_values = Dict{Int, Dict{String, Vector{Float64}}}()
    for fault_idx in 1:n_faults
        fault_param_values[fault_idx] = Dict{String, Vector{Float64}}()
        for param in values(fault_param_mapping)
            fault_param_values[fault_idx][param] = Vector{Float64}(undef, n_sims)
        end
    end
    
    # Use threading to parallelize 
    # TO DO: read the threads documentaiton and configure this properly
    Threads.@threads for i in 1:n_sims
        # Create a thread-local random number generator to avoid race conditions
        local_rng = Random.MersenneTwister(base_seed + i)
        
        # Create a copy of the stress inputs for this simulation
        sim_stress = copy(stress_inputs)
        
        # Apply randomization to stress parameters
        for (uncertainty_param, stress_param) in stress_param_mapping
            if haskey(uncertainties, uncertainty_param)
                base_value = stress_inputs[stress_param]
                uncertainty_value = uncertainties[uncertainty_param]
                
                if !isnothing(base_value) && !isnothing(uncertainty_value)
                    # Get the uncertainty value
                    uncertainty = isa(uncertainty_value, Number) ? uncertainty_value : 
                                (isa(uncertainty_value, AbstractArray) && !isempty(uncertainty_value) ? uncertainty_value[1] : 0.0)
                    
                    
                    # Only apply randomization if uncertainty is greater than 0
                    if uncertainty > 0.0
                        sim_stress[stress_param] = base_value + rand(local_rng, Uniform(-uncertainty, uncertainty))
                    else
                        sim_stress[stress_param] = base_value
                    end
                    
                    # Store the randomized value
                    stress_param_values[stress_param][i] = sim_stress[stress_param]
                end
            end
        end
        
        # Prepare fault data for this simulation
        sim_faults = Vector{Dict{String, Any}}(undef, n_faults)
        
        for (idx, fault) in enumerate(fault_dicts)
            sim_fault = copy(fault)
            
            # Apply randomization to fault parameters
            for (uncertainty_param, fault_param) in fault_param_mapping
                if haskey(uncertainties, uncertainty_param)
                    uncertainty_value = uncertainties[uncertainty_param]
                    
                    if !isnothing(uncertainty_value)
                        # Get the uncertainty value
                        uncertainty = isa(uncertainty_value, Number) ? uncertainty_value : 
                                    (isa(uncertainty_value, AbstractArray) && !isempty(uncertainty_value) ? uncertainty_value[1] : 0.0)
                        
                        
                        # Only apply randomization if uncertainty is greater than 0
                        if uncertainty > 0.0
                            sim_fault[fault_param] += rand(local_rng, Uniform(-uncertainty, uncertainty))
                        end
                        
                        # Store the randomized value
                        fault_param_values[idx][fault_param][i] = sim_fault[fault_param]
                    end
                end
            end
            
            sim_faults[idx] = sim_fault
        end
        
        # Get friction coefficient from first fault 
        friction_coefficient = sim_faults[1]["FrictionCoefficient"]
        
        # Calculate absolute stresses using the first fault's friction coefficient
        stress_state_obj, initial_pressure = GeomechanicsModel.calculate_absolute_stresses(sim_stress, friction_coefficient, stress_model_type)
        
        # Update stress state with absolute values
        sim_stress["vertical_stress"] = stress_state_obj.principal_stresses[1]
        sim_stress["max_horizontal_stress"] = stress_state_obj.principal_stresses[3]
        sim_stress["min_horizontal_stress"] = stress_state_obj.principal_stresses[2]
        sim_stress["pore_pressure"] = initial_pressure
        
        # Create StressState object for slip calculations
        stress_state_obj = StressState(
            [sim_stress["vertical_stress"], 
            sim_stress["min_horizontal_stress"], 
            sim_stress["max_horizontal_stress"]],
            sim_stress["max_stress_azimuth"]
        )
        
        # Calculate slip pressure for each fault
        for (fault_idx, fault) in enumerate(sim_faults)
            sig_normal, tau_normal, s11, s22, s33, s12, n1, n2 = calculate_fault_effective_stresses(
                fault["Strike"], 
                fault["Dip"], 
                stress_state_obj, 
                sim_stress["pore_pressure"], 
                0.0  # dp is 0 for initial calculation
            )
            
            # Calculate critical pore pressure
            try
                pp_to_slip = ComputeCriticalPorePressureForFailure(
                    sig_normal,
                    tau_normal,
                    fault["FrictionCoefficient"],  # Use each fault's own friction coefficient (we apply the same friction coefficient to all faults so it doesn't matter)
                    sim_stress["pore_pressure"],
                    1.0,  # biot coefficient
                    0.5,  # Poisson's ratio
                    1.0   # dp
                )
                
                # Store result in the shared array
                pps_to_slip[fault_idx, i] = pp_to_slip
            catch e
                
                println("Warning: Error calculating slip pressure for fault $(fault_idx) in simulation $(i): $(e)")
                
            end
        end
    end
    
    # Create the result DataFrame
    # Build the DataFrame manually for clearer control
    result_rows = []
    
    for sim_id in 1:n_sims
        for fault_idx in 1:n_faults
            push!(result_rows, (
                SimulationID = sim_id,
                FaultID = actual_fault_ids[fault_idx],  # Use the actual fault ID instead of just the loop index
                SlipPressure = pps_to_slip[fault_idx, sim_id]
            ))
        end
    end
    
    result_df = DataFrame(result_rows)

    return result_df, stress_param_values, fault_param_values
end

"""
Calculate statistics from Monte Carlo results and return as a DataFrame
"""
function calculate_statistics(results::Union{DataFrame, SharedArray{Float64, 2}})
    # Create a DataFrame to store the statistics
    stats_df = DataFrame(
        FaultID = String[],  # Changed from Int to String
        Mean = Float64[],
        StdDev = Float64[],
        Median = Float64[],
        Min = Float64[],
        Max = Float64[]
    )
    
    if isa(results, DataFrame)
        # Process DataFrame results
        fault_ids = unique(results.FaultID)
        
        for fault_id in fault_ids
            # Filter data for this fault
            fault_data = results[results.FaultID .== fault_id, :SlipPressure]
            
            # Calculate statistics
            push!(stats_df, (
                fault_id,
                round(mean(fault_data), digits=3),
                round(std(fault_data), digits=3),
                round(quantile(fault_data, 0.50), digits=3),
                round(minimum(fault_data), digits=3),
                round(maximum(fault_data), digits=3)
            ))
        end
    else
        # Process SharedArray results
        for fault_idx in 1:size(results, 1)
            fault_data = results[fault_idx, :]
            
            # Calculate statistics
            push!(stats_df, (
                string(fault_idx),  # Convert to string
                mean(fault_data),
                std(fault_data),
                quantile(fault_data, 0.50),
                minimum(fault_data),
                maximum(fault_data)
            ))
        end
    end
    
    # Sort by FaultID for consistency
    sort!(stats_df, :FaultID)
    
    return stats_df
end

"""
Main function
"""
function main()
    
    
    scratchPath = ARGS[1]

    println("Scratch path: $scratchPath")

    helper = TexNetWebToolLaunchHelperJulia(scratchPath)

    # get the optional random_seed user input
    random_seed = get_parameter_value(helper, 3, "random_seed")

    println("Extracting stress state values from args.json...")
    stress_inputs = Dict(
        "reference_depth" => get_parameter_value(helper, 2, "reference_depth"),
        "vertical_stress" => get_parameter_value(helper, 2, "vertical_stress"),
        "min_horizontal_stress" => get_parameter_value(helper, 2, "min_horizontal_stress"),
        "max_horizontal_stress" => get_parameter_value(helper, 2, "max_horizontal_stress"),
        "pore_pressure" => get_parameter_value(helper, 2, "pore_pressure"),
        "max_stress_azimuth" => get_parameter_value(helper, 2, "max_stress_azimuth"),
        "aphi_value" => get_parameter_value(helper, 2, "aphi_value")
    )

    stress_model_type = get_parameter_value(helper, 2, "stress_model_type")

    #println("stress_inputs: $stress_inputs")

    
    println("stress_model_type: $stress_model_type")
    
    # create uncertainties dictionary
    uncertainties = Dict()

    if stress_model_type == "all_gradients" || stress_model_type == "gradients"
        uncertainties = Dict(
            "vertical_stress_gradient_uncertainty" => get_parameter_value(helper, 3, "vertical_stress_gradient_uncertainty"),
            "initial_pore_pressure_gradient_uncertainty" => get_parameter_value(helper, 3, "initial_pore_pressure_gradient_uncertainty"),
            "max_stress_azimuth_uncertainty" => get_parameter_value(helper, 3, "max_stress_azimuth_uncertainty"),
            "max_horizontal_stress_uncertainty" => get_parameter_value(helper, 3, "max_horizontal_stress_gradient_uncertainty"),
            "min_horizontal_stress_uncertainty" => get_parameter_value(helper, 3, "min_horizontal_stress_gradient_uncertainty")
        )
    elseif stress_model_type == "aphi_model" || stress_model_type == "aphi_min"
        uncertainties = Dict(
            "vertical_stress_gradient_uncertainty" => get_parameter_value(helper, 3, "vertical_stress_gradient_uncertainty"),
            "initial_pore_pressure_gradient_uncertainty" => get_parameter_value(helper, 3, "initial_pore_pressure_gradient_uncertainty"),
            "max_stress_azimuth_uncertainty" => get_parameter_value(helper, 3, "max_stress_azimuth_uncertainty"),
            "aphi_value_uncertainty" => get_parameter_value(helper, 3, "aphi_value_uncertainty"),
            "min_horizontal_stress_uncertainty" => get_parameter_value(helper, 3, "min_horizontal_stress_gradient_uncertainty")
        )
    elseif stress_model_type == "aphi_model_no_min" || stress_model_type == "aphi_min"
        uncertainties = Dict(
            "vertical_stress_gradient_uncertainty" => get_parameter_value(helper, 3, "vertical_stress_gradient_uncertainty"),
            "initial_pore_pressure_gradient_uncertainty" => get_parameter_value(helper, 3, "initial_pore_pressure_gradient_uncertainty"),
            "max_stress_azimuth_uncertainty" => get_parameter_value(helper, 3, "max_stress_azimuth_uncertainty"),
            "aphi_value_uncertainty" => get_parameter_value(helper, 3, "aphi_value_uncertainty"),
            "strike_angles_uncertainty" => get_parameter_value(helper, 3, "strike_angles_uncertainty"),
            "dip_angles_uncertainty" => get_parameter_value(helper, 3, "dip_angles_uncertainty"),
            "friction_coefficient_uncertainty" => get_parameter_value(helper, 3, "friction_coefficient_uncertainty")
        )
    end

    #println("uncertainties from the portal: $uncertainties")

    println("Extracting fault parameters from args.json...")
    fault_inputs_filepath = get_dataset_file_path(helper, 3, "faults")
    fault_inputs = CSV.read(fault_inputs_filepath, DataFrame)
    


    # get the number of MC simulations to run
    n_sims = get_parameter_value(helper, 3, "mc_iterations")

    
    println("Running Monte Carlo simulation with $(Threads.nthreads()) threads for $n_sims iterations")

    # create the shared array for the results (optimized for parallel processing)
    mc_pp_results, stress_param_values, fault_param_values = run_monte_carlo(stress_inputs, fault_inputs, uncertainties, n_sims, stress_model_type, random_seed)


    # save the monte carlo dataframe as a dataset to the portal
    #save_dataframe_as_parameter!(helper, 3, "prob_geomechanics_results", mc_pp_results)

    #println("prob_geomechanics_results:")
    #println(first(mc_pp_results, 10))


    # before we create the CDF graph data, we need to read the  colors from the previous step 
    # read the dataframe with the deterministic results
    deterministic_results_filepath = get_dataset_file_path(helper, 2, "det_geomechanics_results")
    deterministic_results_df = CSV.read(deterministic_results_filepath, DataFrame)
    
    
    
    # prepare the data for the d3.js CDF plot
    d3_cdf_data = prob_geomechanics_cdf(mc_pp_results, deterministic_results_df)

    # save the d3.js CDF data as a dataset to the portal
    save_dataframe_as_parameter!(helper, 3, "prob_geomechanics_cdf_graph_data", d3_cdf_data)

    # Generate input parameter distributions for histograms
    println("Generating input parameter distributions for histograms...")
    #=
    input_histograms_data = input_distribution_histograms_to_d3(
        mc_pp_results,
        stress_inputs,
        fault_inputs,
        uncertainties,
        stress_model_type,
        stress_param_values,
        fault_param_values
    )
        =#

    # Save the histogram data as a parameter
    #save_dataframe_as_parameter!(helper, 3, "prob_geomechanics_input_histograms", input_histograms_data)
    println("Input parameter distributions saved for visualization")

    # Calculate statistics
    stats = calculate_statistics(mc_pp_results)

    # save the statistics dataframe as a dataset to the portal
    save_dataframe_as_parameter!(helper, 3, "prob_geomechanics_stats", stats)
    
    # Generate sensitivity analysis (tornado chart) data for each fault
    println("Generating sensitivity analysis data for tornado charts...")
    
    # Add the fault parameters to the uncertainties dictionary if they're not already there
    if !haskey(uncertainties, "strike_angles_uncertainty")
        uncertainties["strike_angles_uncertainty"] = get_parameter_value(helper, 3, "strike_angles_uncertainty")
    end
    if !haskey(uncertainties, "dip_angles_uncertainty")
        uncertainties["dip_angles_uncertainty"] = get_parameter_value(helper, 3, "dip_angles_uncertainty")
    end
    if !haskey(uncertainties, "friction_coefficient_uncertainty")
        uncertainties["friction_coefficient_uncertainty"] = get_parameter_value(helper, 3, "friction_coefficient_uncertainty")
    end
    if !haskey(uncertainties, "min_horizontal_stress_uncertainty")
        uncertainties["min_horizontal_stress_uncertainty"] = get_parameter_value(helper, 3, "min_horizontal_stress_gradient_uncertainty")
    end
    if !haskey(uncertainties, "max_horizontal_stress_uncertainty")
        uncertainties["max_horizontal_stress_uncertainty"] = get_parameter_value(helper, 3, "max_horizontal_stress_gradient_uncertainty")
    end
    if !haskey(uncertainties, "initial_pore_pressure_gradient_uncertainty")
        uncertainties["initial_pore_pressure_gradient_uncertainty"] = get_parameter_value(helper, 3, "initial_pore_pressure_gradient_uncertainty")
    end
    if !haskey(uncertainties, "vertical_stress_gradient_uncertainty")
        uncertainties["vertical_stress_gradient_uncertainty"] = get_parameter_value(helper, 3, "vertical_stress_gradient_uncertainty")
    end
    

    # Calculate absolute stresses first to ensure we have the right values
    # Using the first fault's friction coefficient
    #friction_coefficient = fault_inputs[1, :FrictionCoefficient]
    #stress_state_obj, initial_pressure = calculate_absolute_stresses(stress_inputs, friction_coefficient, stress_model_type)
    
    

    # Generate uncertainty variability data
    uncertainty_variability_df = uncertainty_variability_inputs_to_d3(
        uncertainties, 
        stress_model_type, 
        stress_inputs, 
        fault_inputs
    )

    # print the uncertainty variability dataframe
    println("uncertainty_variability_df: $uncertainty_variability_df")
    
    # Save the dataframe as a parameter
    save_dataframe_as_parameter!(helper, 3, "uncertainty_variability_tornado_chart_data", uncertainty_variability_df)

    # Generate tornado chart data for each fault
    tornado_charts_data = Dict{String, DataFrame}()
    for (idx, row) in enumerate(eachrow(fault_inputs))
        fault_id = string(idx)
        println("Generating tornado chart data for Fault #$fault_id...")
        
        # Calculate absolute stresses first to ensure we have the right values for the tornado chart
        friction_coefficient = row.FrictionCoefficient
        stress_state_obj, initial_pressure = GeomechanicsModel.calculate_absolute_stresses(stress_inputs, friction_coefficient, stress_model_type)
        
        # MODIFIED: Use original gradient values for the tornado chart function
        # Instead of passing absolute stresses, we pass the original gradients
        local_stress_inputs = copy(stress_inputs)
        
        # We still need the reference values to calculate the baseline slip pressure
        # Store these as separate variables, not in the dictionary that gets passed to the function
        absolute_vertical_stress = stress_state_obj.principal_stresses[1]
        absolute_max_horizontal_stress = stress_state_obj.principal_stresses[3]
        absolute_min_horizontal_stress = stress_state_obj.principal_stresses[2]
        absolute_pore_pressure = initial_pressure

        # Create a separate stress state for baseline calculations using absolute values
        baseline_stress_state = StressState(
            [absolute_vertical_stress, 
             absolute_min_horizontal_stress, 
             absolute_max_horizontal_stress],
            local_stress_inputs["max_stress_azimuth"]
        )
        
        # Make sure all stress uncertainties are defined according to stress_model_type
        local_uncertainties = copy(uncertainties)
 
        println("local_uncertainties: ")
        for (key, value) in local_uncertainties
            println("  - $key: $value")
        end
        
        # Initialize the array for stress uncertainty parameters
        stress_uncertainty_params = String[]
        
        # Ensure we have all necessary uncertainty parameters based on stress_model_type
        if stress_model_type == "gradients" || stress_model_type == "all_gradients"
            
            stress_uncertainty_params = [
                "vertical_stress_gradient_uncertainty",
                "initial_pore_pressure_gradient_uncertainty", 
                "max_stress_azimuth_uncertainty",
                "max_horizontal_stress_uncertainty",
                "min_horizontal_stress_uncertainty"
            ]
        elseif stress_model_type == "aphi_model"
            
            stress_uncertainty_params = [
                "vertical_stress_gradient_uncertainty",
                "initial_pore_pressure_gradient_uncertainty", 
                "max_stress_azimuth_uncertainty",
                "aphi_value_uncertainty",
                "min_horizontal_stress_uncertainty"
            ]
        elseif stress_model_type == "aphi_model_no_min"
            
            stress_uncertainty_params = [
                "vertical_stress_gradient_uncertainty",
                "initial_pore_pressure_gradient_uncertainty", 
                "max_stress_azimuth_uncertainty",
                "aphi_value_uncertainty"
            ]
        end

        # the uncertainties for the faults are always the same
        fault_uncertainty_params = [
            "strike_angles_uncertainty",
            "dip_angles_uncertainty",
            "friction_coefficient_uncertainty"
        ]
        
        # MODIFIED: Add reference depth to the inputs so the function can convert gradients to absolute values
        # This ensures the tornado chart function has all needed information
        local_stress_inputs["absolute_vertical_stress"] = absolute_vertical_stress
        local_stress_inputs["absolute_min_horizontal_stress"] = absolute_min_horizontal_stress
        local_stress_inputs["absolute_max_horizontal_stress"] = absolute_max_horizontal_stress
        local_stress_inputs["absolute_pore_pressure"] = absolute_pore_pressure
        
        tornado_df = JuliaFSPGraphs.fault_sensitivity_tornado_chart_to_d3(
            local_stress_inputs,
            fault_inputs,
            local_uncertainties,
            idx,
            stress_model_type
        )
        
        # Store the tornado chart data for this fault
        tornado_charts_data[fault_id] = tornado_df
    end
    
    # Create and save the combined tornado chart data as a parameter
    println("Creating combined tornado chart data for all faults...")
    
    # Create vector to hold individual dataframes with fault_id added
    dfs_with_fault_id = []
    
    # Add fault_id to each DataFrame before combining them
    for (fault_id, df) in tornado_charts_data
        df_copy = copy(df)
        df_copy.fault_id = fill(fault_id, nrow(df))
        push!(dfs_with_fault_id, df_copy)
    end
    
    # Concatenate all DataFrames
    combined_tornado_df = vcat(dfs_with_fault_id...)
    
    
    
    # Add deterministic slip pressure to the combined dataframe
    # Create a mapping of fault IDs to deterministic slip pressures
    det_slip_pressures = Dict{String, Float64}()
    for row in eachrow(deterministic_results_df)
        # Extract fault ID - handle different column names
        fault_id = ""
        if "FaultID" in names(deterministic_results_df)
            fault_id = string(row.FaultID)
        elseif "ID" in names(deterministic_results_df)
            fault_id = string(row.ID)
        else
            # If no ID column, use index
            fault_id = string(findfirst(r -> r === row, eachrow(deterministic_results_df)))
        end
        
        # Extract slip pressure - assume column is named SlipPressure
        if "slip_pressure" in names(deterministic_results_df)
            det_slip_pressures[fault_id] = row.slip_pressure
        elseif "Slip_Pressure" in names(deterministic_results_df)
            det_slip_pressures[fault_id] = row.Slip_Pressure
        end

        
    end
    
    # Add the deterministic slip pressures to the combined dataframe
    # Don't remove fault_id column yet - use it to match with deterministic results
    combined_tornado_df.det_slip_pressure = map(row -> 
        haskey(det_slip_pressures, string(row.id)) ? 
        det_slip_pressures[string(row.id)] : NaN, 
        eachrow(combined_tornado_df))
    
    # Count how many non-NaN values were added
    non_nan_count = count(!isnan, combined_tornado_df.det_slip_pressure)
    #println("Added $(non_nan_count) non-NaN deterministic slip pressure values out of $(nrow(combined_tornado_df)) rows")
    
    # Look at some sample matches
    if nrow(combined_tornado_df) > 0
        for i in 1:min(10, nrow(combined_tornado_df))
            row = combined_tornado_df[i, :]
            println("Row $i: id=$(row.id), det_slip_pressure=$(row.det_slip_pressure)")
        end
    end
    
    # If no matches were found (or very few), try a different approach
    if non_nan_count < nrow(combined_tornado_df) / 8 # Less than 1/8th of rows have matches
        println("Few matches found. Trying alternative approach...")
        
        # Initialize column with NaN
        combined_tornado_df.det_slip_pressure .= NaN
        
        # Get unique fault IDs in combined dataframe
        unique_ids = unique(combined_tornado_df.id)
        
        # Directly set values for each fault ID
        for fault_id in unique_ids
            if haskey(det_slip_pressures, string(fault_id))
                combined_tornado_df[combined_tornado_df.id .== fault_id, :det_slip_pressure] .= det_slip_pressures[string(fault_id)]
                println("Set slip pressure $(det_slip_pressures[string(fault_id)]) for fault $fault_id")
            end
        end
        
        # Count how many non-NaN values after alternative approach
        non_nan_count_after = count(!isnan, combined_tornado_df.det_slip_pressure)
        println("After alternative approach: $(non_nan_count_after) non-NaN values")
    end
    
    # Now remove the redundant fault_id column
    select!(combined_tornado_df, Not(:fault_id))

    println("combined_tornado_df: $combined_tornado_df")
    
    # Debug: Print deterministic results dataframe
    println("deterministic_results_df: $deterministic_results_df")
    pretty_table(deterministic_results_df)
    
    # Save the combined DataFrame as the only tornado chart parameter
    save_dataframe_as_parameter!(helper, 3, "prob_geomechanics_fault_sensitivity_tornado_chart_data", combined_tornado_df)

    #error("Stop hereeeeeeeeeeeeee")
    
    # Also save as CSV file in the current directory for easy access
    CSV.write("tornado_chart_all_faults.csv", combined_tornado_df)
    println("Combined tornado chart data saved as parameter and CSV file")

    #output_plot_path = dirname(args["output-json"])
    #println("\nGenerating CDF plots for each fault...")
    #plot_cdf(mc_pp_results, output_plot_path)
    #println("\nCDF plots saved in $(output_plot_path)")

    write_final_args_file(helper, joinpath(helper.scratch_path, ARGS_FILE_NAME))

    # explicitly set this step's success state to true
    set_success_for_step_index!(helper, 3, true)

    write_results_file(helper)
    

    println("Probabilistic geomechanics process completed successfully.")
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

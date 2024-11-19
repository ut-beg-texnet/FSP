module ProbabilisticGeomechanicsDriver

using JSON
using CSV
using DataFrames    
using Statistics

#using UncertaintyQuantification
#using BenchmarkTools


include("src/core/ProbabilisticGeomechanicsCalculations.jl")
include("src/core/DeterministicGeomechanicsCalculations.jl")
include("src/graph_helpers/FaultSlipCDFPlot.jl")

using .ProbabilisticGeomechanicsCalculations: monte_carlo
using .DeterministicGeomechanicsCalculations: mohrs_3D, calc_ppfail
using .FaultSlipCDFPlot: plot_cdf


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


function load_uncertainties(file_path::String)
    try
        df = CSV.read(file_path, DataFrame)
        return df
    catch e
        println("Error loading uncertainties from ", file_path, ": ", e)
        return nothing
    end
end


# New FSP function----------------------------------------------------------
# function to create samples for Monte Carlo simulation using the user defined PDFs (assuming they are all uniform for now)
# Note that the input is just 1 row of the dataframes containing the use defined PDFs, where each row correspond to a different location
function create_mc_samples_from_pdfs(uncertainties_df::DataFrames.DataFrameRow{DataFrames.DataFrame, DataFrames.Index}, nSamples::Int64)
    # container that holds the user defined PDFs for sampling
    tmp = Union{UncertaintyQuantification.Parameter, UncertaintyQuantification.RandomVariable}[]
    for name in names(uncertainties_df)
        if typeof(uncertainties_df[name]) == Array{UncertaintyQuantification.Parameter,1}
            tmp = [tmp; uncertainties_df[name]]  
            
        elseif typeof(uncertainties_df[name]) == UncertaintyQuantification.Parameter
            tmp = [tmp; uncertainties_df[name]]  ## constant value assigned to this location
            
        elseif typeof(uncertainties_df[name]) == Int64 || typeof(uncertainties_df[name]) == Float64
            tmp = [tmp; UncertaintyQuantification.Parameter.(uncertainties_df[name], Symbol(name)) ]  ## constant value assigned to this location

        elseif typeof(uncertainties_df[name]) != String 
            tmp = [tmp; uncertainties_df[name]]  ## pdf assigned to this location            

        end
    end

    samples=UncertaintyQuantification.sample(tmp, nSamples)

    return samples

end



# checks if the uncertainties have a string as first row (column)
# makes sure all vlaues are Float64
function validate_uncertainties_data(df::DataFrame)
    # Check if the first row has strings (assuming headers in the first row)
    is_row_1_header = all(x -> occursin(r"^[A-Za-z]", string(x)), df[1, names(df)])
    start_row = is_row_1_header ? 2 : 1

    # Ensure all columns are Float64
    for col in names(df)
        if eltype(df[!, col]) <: AbstractString
            # Convert each element in the column to Float64 in-place
            df[!, col] .= parse.(Float64, df[!, col])
        end
    end
    return start_row
end


# validate uncertainties input data format
function check_uncertainties_format(uncertainties_df::DataFrame, deterministic_vals::Dict{String, Any}, stress_model_type::Symbol)
    # Define the expected columns
    columns_check = if stress_model_type == :aphi
        ["Vertical Stress Grad", "Initial PP Grad", "Strike Angles", "Dip Angles", 
         "Max Horiz. Stress Dir", "Friction Coeff Mu", "A Phi Parameter"]
    else
        ["Vertical Stress Grad", "Min Horiz. Grad", "Max Horiz. Grad", "Initial PP Grad", 
         "Strike Angles", "Dip Angles", "Max Horiz. Stress Dir", "Friction Coeff Mu"]
    end

    # Check for required columns
    for col in columns_check
        if !(col in names(uncertainties_df))
            println("Error: Column '$col' not found in uncertainties file.")
            return false
        end
        
        # Access the column values
        col_values = uncertainties_df[!, col]
        if any(col_values .< 0)
            println("Error: Column '$col' has negative values in uncertainties file.")
            return false
        end
        
        # Ensure deterministic value is scalar
        det_val = deterministic_vals[col]
        if isa(det_val, Vector)
            det_val = det_val[1]  # Use the first element if it's a vector
        end
        
        # Compare each uncertainty value with the deterministic value
        for val in col_values
            if val >= det_val
                println("Error: Uncertainty in '$col' ($val) exceeds or equals deterministic value ($det_val).")
                return false
            end
        end
    end
    return true
end





# calculates bounds based on uncertainties
function get_prob_bounds(deterministic_vals::Dict, uncertainties_df::DataFrame, stress_model_type::Symbol)

    depth = deterministic_vals[:Depth]

    bounds = Dict()

    # stress bounds
    sv_grad = deterministic_vals["Vertical Stress Grad"]
    sv_pm = uncertainties_df[1, "Vertical Stress Grad"]
    sv_bounds = [sv_grad - sv_pm, sv_grad, sv_grad + sv_pm] .* depth
    bounds["Vertical Stress Grad"] = sv_bounds
    
    # pressure bounds
    pp_grad = deterministic_vals["Initial PP Grad"]
    pp_pm = uncertainties_df[1, "Initial PP Grad"]
    pp_bounds = [pp_grad - pp_pm, pp_grad, pp_grad + pp_pm] .* depth
    bounds["Initial PP Grad"] = pp_bounds
    

    # model check
    if stress_model_type == :aphi
        aphi_value = deterministic_vals["A Phi Parameter"]
        aphi_pm = uncertainties_df[1, "A Phi Parameter"]
        aphi_bounds = [aphi_value - aphi_pm, aphi_value, aphi_value + aphi_pm]
    else #gradient stress model
        min_horiz_grad = deterministic_vals["Min Horiz. Grad"]
        max_horiz_grad = deterministic_vals["Max Horiz. Grad"]
        min_pm = uncertainties_df[1, "Min Horiz. Grad"]
        max_pm = uncertainties_df[1, "Max Horiz. Grad"]
        bounds["Min Horiz. Grad"] = [min_horiz_grad - min_pm, min_horiz_grad, min_horiz_grad + min_pm]
        bounds["Max Horiz. Grad"] = [max_horiz_grad - max_pm, max_horiz_grad, max_horiz_grad + max_pm]
    end

    return bounds
end


# WORK ON THIS
# GET FAULT AND STRESS DATA FROM STEP 2 JSON FILES (IMPLEMENT THE SAME FOR UNCERTAINTIES)
function probabilistic_geomechanics_driver(step2_data_file::String, uncertainties_data_file::String, det_geomechanics_results_file::String, nruns::Int)
    # load data from previous steps
    step_2_data = load_json_data(step2_data_file)
    if step_2_data === nothing
        println("Error loading data from step 2.")
        return nothing
    end
    fault_data = step_2_data["fault_data"]
    stress_data = step_2_data["stress_data"]


    det_geomechanics_data = load_json_data(det_geomechanics_results_file)
    if det_geomechanics_data === nothing
        println("Error loading deterministic geomechanics results.")
        return nothing
    end

    geomechanics_results = det_geomechanics_data["geomechanics_results"]

    # get 'aphi_use' and 'aphi_value' from deterministic geomechanics step
    aphi_use = geomechanics_results["aphi_use"]
    aphi_value = geomechanics_results["aphi_value"]
    min_hor_stress = get(geomechanics_results, "min_hor_stress", nothing)

    # determine stress model
    aphi_use = geomechanics_results["aphi_use"]
    aphi_value = geomechanics_results["aphi_value"]
    min_hor_stress = get(geomechanics_results, "min_hor_stress", nothing)
    model_type = if aphi_use
        min_hor_stress !== nothing ? :aphi_with_min_hor : :aphi_without_min_hor
    else
        :gradient
    end
    println("Stress model used: ", model_type)

    # create dictionary with deterministic values
    # used to check uncertainties' bounds
    deterministic_vals = Dict(
        "Vertical Stress Grad" => stress_data["vertical_stress_gradient_psi_ft"],
        "Min Horiz. Grad" => stress_data["min_horizontal_stress_gradient_psi_ft"],
        "Max Horiz. Grad" => stress_data["max_horizontal_stress_gradient_psi_ft"],
        "Initial PP Grad" => stress_data["initial_reservoir_pressure_gradient_psi_ft"],
        "Max Horiz. Stress Dir" => stress_data["max_horizontal_stress_direction"],
        "Strike Angles" => [fault["strike"] for fault in fault_data],
        "Dip Angles" => [fault["dip"] for fault in fault_data],
        "Friction Coeff Mu" => fault_data[1]["friction"], # assuming same friction coefficient for all faults
        "A Phi Parameter" => aphi_value,
        "Depth" => stress_data["reference_depth_ft"]
    )

    println("Deterministic values: ", deterministic_vals)

    # get uncertainties from CSV
    println("Loading uncertainties...")
    uncertainties_df = load_uncertainties(uncertainties_data_file)
    if uncertainties_df === nothing
        println("Error loading uncertainties data.")
        return nothing
    end
    println("Uncertainties loaded from file: ", uncertainties_df)

    start_row = validate_uncertainties_data(uncertainties_df)
    # check uncertainties format
    if !check_uncertainties_format(uncertainties_df, deterministic_vals, model_type)
        println("Error: Uncertainties data format is incorrect.")
        return nothing
    end
    println("Uncertainties values have been validated.")


    # prepare inputs for monte carlo simulation
    # in0 = Nominal input values (center points of distribution)
    # inSig = A vector of ranges or standard deviations for each element in in0
    if model_type == :gradient
        in0 = [
            deterministic_vals["Vertical Stress Grad"],
            deterministic_vals["Min Horiz. Grad"],
            deterministic_vals["Max Horiz. Grad"],
            deterministic_vals["Initial PP Grad"],
            deterministic_vals["Strike Angles"],
            deterministic_vals["Dip Angles"],
            deterministic_vals["Max Horiz. Stress Dir"],
            deterministic_vals["Friction Coeff Mu"]
        ]
        inSig = [
            uncertainties_df[start_row, "Vertical Stress Grad"],
            uncertainties_df[start_row, "Min Horiz. Grad"],
            uncertainties_df[start_row, "Max Horiz. Grad"],
            uncertainties_df[start_row, "Initial PP Grad"],
            uncertainties_df[start_row, "Strike Angles"],
            uncertainties_df[start_row, "Dip Angles"],
            uncertainties_df[start_row, "Max Horiz. Stress Dir"],
            uncertainties_df[start_row, "Friction Coeff Mu"]
        ]
    elseif model_type == :aphi_with_min_hor
        in0 = [
            deterministic_vals["Vertical Stress Grad"],
            deterministic_vals["Initial PP Grad"],
            deterministic_vals["Strike Angles"],
            deterministic_vals["Dip Angles"],
            deterministic_vals["Max Horiz. Stress Dir"],
            deterministic_vals["Friction Coeff Mu"],
            deterministic_vals["A Phi Parameter"]
        ]
        inSig = [
            uncertainties_df[start_row, "Vertical Stress Grad"],
            uncertainties_df[start_row, "Initial PP Grad"],
            uncertainties_df[start_row, "Strike Angles"],
            uncertainties_df[start_row, "Dip Angles"],
            uncertainties_df[start_row, "Max Horiz. Stress Dir"],
            uncertainties_df[start_row, "Friction Coeff Mu"],
            uncertainties_df[start_row, "A Phi Parameter"]
        ]
    elseif model_type == :aphi_without_min_hor
        in0 = [
            deterministic_vals["Vertical Stress Grad"],
            deterministic_vals["Initial PP Grad"],
            deterministic_vals["Strike Angles"],
            deterministic_vals["Dip Angles"],
            deterministic_vals["Max Horiz. Stress Dir"],
            deterministic_vals["Friction Coeff Mu"]
        ]
        inSig = [
            uncertainties_df[start_row, "Vertical Stress Grad"],
            uncertainties_df[start_row, "Initial PP Grad"],
            uncertainties_df[start_row, "Strike Angles"],
            uncertainties_df[start_row, "Dip Angles"],
            uncertainties_df[start_row, "Max Horiz. Stress Dir"],
            uncertainties_df[start_row, "Friction Coeff Mu"]
        ]
    end
    println("MONTE CARLO INPUTS-----------------------------")
    println("Monte Carlo input in0: ", in0)
    println("Monte Carlo input inSig (uncertainties): ", inSig)
    println("MONTE CARLO INPUTS-----------------------------")
    

    # call calc_ppfail within monte_carlo
    function mc_wrapper(inputs)
        #println("Running mohrs_wrapper. Inputs used: ", inputs)
        depth = deterministic_vals["Depth"]
        
        # Initialize variables
        Sig0 = Float64[]
        #p0 = Float64()
        biot = 1.0  # Biot coefficient (hardcoded)
        nu = 0.5    # Poisson's ratio (hardcoded)
        #mu = Float64()
        strikes = inputs[5]
        #dp = zeros(length(strikes))  # Adjust as necessary
        
        dips = inputs[6]
        #SHdir = Float64()
    
        if model_type == :gradient
            # Inputs Vector for :gradient:
            # [1] Vertical Stress Grad
            # [2] Min Horiz. Grad
            # [3] Max Horiz. Grad
            # [4] Initial PP Grad
            # [5] Strike Angles
            # [6] Dip Angles
            # [7] Max Horiz. Stress Dir
            # [8] Friction Coeff Mu
    
            # Map inputs
            Sig0 = [
                inputs[1] * depth,  # Vertical Stress Grad
                inputs[2] * depth,  # Min Horiz. Grad
                inputs[3] * depth   # Max Horiz. Grad
            ]
            p0 = inputs[4] * depth  # Initial PP Grad
    
            strikes = isa(inputs[5], Float64) ? [inputs[5]] : inputs[5]
            dips = isa(inputs[6], Float64) ? [inputs[6]] : inputs[6]
    
            # Access SHdir from inputs[7]
            if isa(inputs[7], Vector) && length(inputs[7]) == 1
                SHdir = Float64(inputs[7][1])
            elseif isa(inputs[7], Float64)
                SHdir = inputs[7]
            else
                error("Error: SHdir input is not a Float64 or a Vector with one element for :gradient model.")
            end
    
            mu = inputs[8]
            dp = zeros(length(strikes))  # Adjust dp as needed
    
        elseif model_type == :aphi_with_min_hor
            # Inputs Vector for :aphi_with_min_hor:
            # [1] Vertical Stress Grad
            # [2] Initial PP Grad
            # [3] Strike Angles
            # [4] Dip Angles
            # [5] Max Horiz. Stress Dir
            # [6] Friction Coeff Mu
            # [7] A Phi Parameter
    
            # Map inputs
            Sig0 = [inputs[1], deterministic_vals["Min Horiz. Grad"], deterministic_vals["Max Horiz. Grad"]]
            p0 = inputs[2] * depth  # Initial PP Grad
    
            strikes = isa(inputs[3], Float64) ? [inputs[3]] : inputs[3]
            dips = isa(inputs[4], Float64) ? [inputs[4]] : inputs[4]
    
            # Access SHdir from inputs[5]
            if isa(inputs[5], Vector) && length(inputs[5]) == 1
                SHdir = Float64(inputs[5][1])
            elseif isa(inputs[5], Float64)
                SHdir = inputs[5]
            else
                error("Error: SHdir input is not a Float64 or a Vector with one element for :aphi_with_min_hor model.")
            end
    
            mu = inputs[6]
            dp = zeros(length(strikes))  # Adjust dp as needed
    
            # aphi_param is present in inputs[7], but not used in calc_ppfail
            aphi_param = inputs[7]
    
        elseif model_type == :aphi_without_min_hor
            # Inputs Vector for :aphi_without_min_hor:
            # [1] Vertical Stress Grad
            # [2] Initial PP Grad
            # [3] Strike Angles
            # [4] Dip Angles
            # [5] Max Horiz. Stress Dir
            # [6] Friction Coeff Mu
    
            # Map inputs
            Sig0 = [inputs[1], deterministic_vals["Min Horiz. Grad"], deterministic_vals["Max Horiz. Grad"]]
            p0 = inputs[2] * depth  # Initial PP Grad
    
            strikes = isa(inputs[3], Float64) ? [inputs[3]] : inputs[3]
            dips = isa(inputs[4], Float64) ? [inputs[4]] : inputs[4]
    
            # Access SHdir from inputs[5]
            if isa(inputs[5], Vector) && length(inputs[5]) == 1
                SHdir = Float64(inputs[5][1])
            elseif isa(inputs[5], Float64)
                SHdir = inputs[5]
            else
                error("Error: SHdir input is not a Float64 or a Vector with one element for :aphi_without_min_hor model.")
            end
    
            mu = inputs[6]
            dp = zeros(length(strikes))  # Adjust dp as needed
    
        else
            error("Unsupported model_type: $model_type")
        end

        
        # Debugging Outputs
        println("Constructed calc_ppfail inputs:")
        println("  Sig0: ", Sig0)
        println("  az (SHdir): ", SHdir)
        println("  p0: ", p0)
        println("  biot: ", biot)
        println("  nu: ", nu)
        println("  mu: ", mu)
        println("  dp: ", dp)
        println("  strikes: ", strikes)
        println("  dips: ", dips)
        
    
        # Call calc_ppfail with Sig0, az, p0, biot, nu, mu, dp, strikes, dips
        failout, tau_fault, sig_fault = calc_ppfail(Sig0, SHdir, p0, biot, nu, mu, dp, strikes, dips)
    
        return failout
    end
    
    
    # run MC simulation
    mc_results, monte_carlo_inputs = monte_carlo(mc_wrapper, in0, inSig, nruns)

    println("Monte Carlo simulation completed. Number of iterations: $nruns")

    
    
    # Check if mc_results is empty
    if isempty(mc_results)
        println("Error: MC Results are empty.")
        return nothing
    end
    

    # Convert the MC results to a matrix (input to function that plots CDF)
    mc_results_matrix = hcat(mc_results...)


    # Additional debug: print size and a sample of mc_results_matrix
    println("mc_results_matrix size: ", size(mc_results_matrix))

    #=
    for fault_id in 1:size(mc_results_matrix, 1)
        fault_results = mc_results_matrix[fault_id, :]
        println("Fault $fault_id: min = $(minimum(fault_results)), max = $(maximum(fault_results)), mean = $(mean(fault_results))")
    end
    
    for fault_id in 1:size(mc_results_matrix, 1)
        fault_results = mc_results_matrix[fault_id, :]
        println("Fault $fault_id first 5 results: ", first(fault_results, 5))
    end
    =#

    # Plot the CDF
    plot_cdf(mc_results_matrix)


end



probabilistic_geomechanics_driver("src/output/step2_user_input_data.json", "src/input/uncertainties_gradient.csv", "src/output/deterministic_geomechanics_results.json", 1000)

end # module

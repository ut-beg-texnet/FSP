"""
Deterministic Geomechanics Analysis Driver Script

This script performs deterministic geomechanical analysis for fault stability assessment.
It takes stress state and fault data from step1 (Model Inputs) and calculates stability metrics
for each fault.

Workflow:
1. Read input data from step1 JSON file containing:
   - Stress state (vertical, min horizontal, max horizontal stresses)
   - Fault data (strike, dip, friction coefficient)
   - Reference depth and pore pressure

2. Calculate absolute stresses at reference depth:
   a. For gradient model: Use provided stress gradients directly
   b. For A-phi with min horizontal stress: Calculate max horizontal using modified model
   c. For A-phi without min horizontal: Calculate both horizontal stresses using standard model

3. For each fault:
   a. Transform principal stresses to fault coordinates
   b. Calculate normal and shear stresses on fault plane
   c. Calculate stability metrics:
      - Slip pressure: Pore pressure required for fault slip
      - Slip tendency: Ratio of shear to normal stress
      - Coulomb Failure Function (CFF)
      - Shear Capacity Utilization (SCU)

4. Output results to JSON file containing:
   - Fault stability metrics for each fault
   - Metadata about the analysis

Usage:
    julia step2_deterministic_geomechanics.jl --input-json path/to/input.json --output-json path/to/output.json

Dependencies:
    - JSON: For reading/writing JSON files
    - ArgParse: For command line argument parsing
    - LinearAlgebra: For matrix operations
    - GeomechanicsModel: Custom module for geomechanical calculations
"""

using JSON
using ArgParse
using LinearAlgebra

include("core/geomechanics_model.jl")
using .GeomechanicsModel

"""
Parse command line arguments
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input-json"
            help = "Path to input JSON file from model_inputs step"
            required = true
        "--output-json"
            help = "Path to output JSON file"
            required = true
    end

    return parse_args(s)
end

"""
Calculate n and Phi values from APhi
Following MATLAB implementation in getHorFromAPhi.m
"""
function calculate_n_phi(aphi::Float64)
    if aphi >= 0 && aphi < 1
        n = 0
    elseif aphi >= 1 && aphi < 2
        n = 1
    elseif aphi >= 2 && aphi <= 3
        n = 2
    else
        error("APhi value must be in range [0,3]. Got: $aphi")
    end
    
    phi = (aphi - (n + 0.5))/(-1)^n + 0.5
    return n, phi
end

"""
Calculate horizontal stresses using Modified A-Phi model (aphi_use == 12)
Following MATLAB implementation in getHorFromAPhi.m
"""
function calculate_modified_aphi_stresses(n::Int, phi::Float64, sv::Float64, sh::Float64, p0::Float64)
    # Calculate effective stresses
    sv_eff = sv - p0
    sh_eff = sh - p0
    
    # Calculate SH based on case
    sH = if n == 0
        # Case 0: Calculate SH directly
        phi * (sv_eff - sh_eff) + sh_eff + p0
    elseif n == 1
        # Case 1: Use provided sh and solve for SH
        (sv_eff - sh_eff + phi*sh_eff)/phi + p0
    elseif n == 2
        # Case 2: Use provided sh and solve for SH
        (sh_eff - sv_eff + phi*sv_eff)/phi + p0
    else
        error("Invalid n value for Modified A-Phi model: $n")
    end
    
    return sH, sh
end

"""
Calculate horizontal stresses using standard A-Phi model with friction
Following MATLAB implementation in getHorFromAPhi.m
"""
function calculate_standard_aphi_stresses(n::Int, phi::Float64, sv::Float64, p0::Float64, mu::Float64)
    if mu <= 0
        return sv, sv  # Both horizontal stresses equal vertical
    end
    
    k = (mu + sqrt(1 + mu^2))^2
    
    if n == 0
        # Case 0: Calculate Sh first, then SH
        sh = (sv - p0)/k + p0
        sH = phi * (sv - sh) + sh
    elseif n == 1
        # Case 1: Solve system of equations
        A = [1.0 -k; phi (1-phi)]
        b = [p0 - k*p0; sv]
        x = A \ b
        sH, sh = x[1], x[2]
    elseif n == 2
        # Case 2: Calculate SH first, then Sh
        sH = k * (sv - p0) + p0
        sh = phi * (sH - sv) + sv
    end
    
    
    return sH, sh
end

"""
Calculate absolute stresses at reference depth
Handles three stress models:
1. Gradients: All stresses provided directly
2. A-phi with min horizontal stress: Calculate max horizontal using modified model
3. A-phi without min horizontal stress: Calculate both horizontal stresses using standard model
"""
function calculate_absolute_stresses(stress_data::Dict, fault_data::Vector)
    # Extract common parameters
    reference_depth = stress_data["reference_depth"]
    vertical_gradient = stress_data["vertical_stress"]
    model_type = stress_data["model_type"]
    pore_pressure_gradient = stress_data["pore_pressure"]
    max_stress_azimuth = stress_data["max_stress_azimuth"]
    
    
    # Calculate absolute stresses at reference depth
    sV = round(vertical_gradient * reference_depth, digits=1)
    p0 = round(pore_pressure_gradient * reference_depth, digits=1)
    
    
    # Get friction coefficient from first fault (assume all faults have same friction)
    μ = fault_data[1]["friction_coefficient"]
    println("\nFriction coefficient from first fault: $(μ)")
    
    # Calculate horizontal stresses based on model type
    if model_type == "gradients"
        println("\nUsing gradients model (all stresses provided)")
        # Convert gradients to absolute stresses
        max_horizontal_gradient = stress_data["max_horizontal_stress"]
        min_horizontal_gradient = stress_data["min_horizontal_stress"]
        
        sH = round(max_horizontal_gradient * reference_depth, digits=1)
        sh = round(min_horizontal_gradient * reference_depth, digits=1)
        
    elseif model_type == "aphi_min" || model_type == "aphi_no_min"
        println("\nUsing A-phi model: $(model_type)")
        # Get A-phi value and calculate n and phi
        aphi = stress_data["aphi_value"]
        n, phi = calculate_n_phi(aphi)
        println("A-phi parameters:")
        println("  A-phi value = $(aphi)")
        println("  n = $(n)")
        println("  φ = $(phi)")
        
        if model_type == "aphi_min" && haskey(stress_data, "min_horizontal_stress") && !isnothing(stress_data["min_horizontal_stress"])
            
            # Convert min horizontal gradient to absolute stress
            sh = stress_data["min_horizontal_stress"] * reference_depth
            sH, _ = calculate_modified_aphi_stresses(n, phi, sV, sh, p0)
        else
            
            # Calculate both horizontal stresses using A-phi model
            sH, sh = calculate_standard_aphi_stresses(n, phi, sV, p0, μ)
        end
    else
        error("Invalid stress model type: $model_type")
    end
    
    
    # Create stress state object with principal stresses vector [Svert, shmin, sHmax]
    stress_state = StressState([sV, sh, sH], max_stress_azimuth)
    
    return stress_state, p0
end

"""
Process each fault and calculate geomechanical parameters
"""
function process_faults(fault_data::Vector, stress_state::StressState, initial_pressure::Float64)
    results = []
    
    for fault in fault_data
        # Calculate geomechanical parameters for each fault
        result = analyze_fault(
            fault["strike"],
            fault["dip"],
            fault["friction_coefficient"],
            stress_state,
            initial_pressure,
            0.0  # No pressure change for deterministic analysis
        )
        
        # Add fault metadata to result
        result["fault_id"] = fault["fault_id"]
        result["easting"] = fault["easting"]
        result["northing"] = fault["northing"]
        result["strike"] = fault["strike"]
        result["dip"] = fault["dip"]
        
        push!(results, result)
    end
    
    return results
end

function main()
    println("\n=== Starting Deterministic Geomechanics Analysis ===")
    
    # Parse command line arguments
    
    args = parse_commandline()
    
    # Read input JSON
    
    input_data = JSON.parsefile(args["input-json"])
    
    # Extract data
    
    fault_data = input_data["faults"]
    stress_data = input_data["stress_state"]
    
    # Calculate absolute stresses at reference depth
    
    stress_state, initial_pressure = calculate_absolute_stresses(stress_data, fault_data)
    
    # Process each fault
    
    results = process_faults(fault_data, stress_state, initial_pressure)
    
    # Prepare output by copying all data from step 1 and adding deterministic results
    println("\nPreparing output...")
    output = deepcopy(input_data)  # Deep copy to avoid modifying original
    
    # Add deterministic results
    output["metadata"] = Dict(
        "step" => "deterministic_geomechanics",
        "input_file" => args["input-json"],
        "description" => "Deterministic geomechanical analysis results"
    )
    output["det_geomechanics_results"] = results
    
    # Create output directory if it doesn't exist
    
    output_dir = dirname(args["output-json"])
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    # Write output JSON
    println("\nSaving results to: $(args["output-json"])")
    open(args["output-json"], "w") do f
        JSON.print(f, output, 4)  # indent with 4 spaces
    end
    
    println("\n=== Deterministic Geomechanics Analysis Complete ===\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

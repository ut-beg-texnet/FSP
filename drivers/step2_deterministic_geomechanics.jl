"""
Deterministic Geomechanics Analysis Driver Script

This script performs deterministic geomechanical analysis for fault stability assessment.
It takes stress state and fault data from step1 (Model Inputs) and calculates stability metrics
for each fault. Run with --threads argument to enable multi-threading and specify number of threads.

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
using Base.Threads

include("core/geomechanics_model.jl")
using .GeomechanicsModel

export process_faults, calculate_absolute_stresses

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
Calculate horizontal stresses using Modified A-Phi model
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
Calculate horizontal stresses using standard A-Phi model with friction (no APhi value provided)
MATLAB --> getHorFromAPhi.m
"""
function calculate_standard_aphi_stresses(n::Int, phi::Float64, sv::Float64, p0::Float64, mu::Float64)
    if mu <= 0
        return sv, sv  # Both horizontal stresses equal vertical
    end
    
    k = (mu + sqrt(1 + mu^2))^2
    
    if n == 0
        sh = (sv - p0)/k + p0
        sH = phi * (sv - sh) + sh

    elseif n == 1
        A = [1.0 -k; phi (1-phi)]
        b = [p0 - k*p0; sv]
        x = A \ b
        sH, sh = x[1], x[2]

    elseif n == 2
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
    # check if we have min horizontal stress, if not set to 0
    if !haskey(stress_data, "min_horizontal_stress") || stress_data["min_horizontal_stress"] === nothing
        stress_data["min_horizontal_stress"] = 0.0
    end

    if !haskey(stress_data, "max_horizontal_stress") || stress_data["max_horizontal_stress"] === nothing
        stress_data["max_horizontal_stress"] = 0.0
    end

    println("Fault data: $(fault_data)")

    #=
    println("Reference depth: $(reference_depth)")
    println("Vertical gradient: $(vertical_gradient)")
    println("Model type: $(model_type)")
    println("Pore pressure gradient: $(pore_pressure_gradient)")
    println("Max stress azimuth: $(max_stress_azimuth)")
    println("Min horizontal stress: $(stress_data["min_horizontal_stress"])")
    =#

    


    # check for aphi_value and set default value if missing
    if !haskey(stress_data, "aphi_value") || stress_data["aphi_value"] === nothing
        stress_data["aphi_value"] = 0.0
    end
    #println("aphi_value: $(stress_data["aphi_value"])")
    
    
    # Calculate absolute stresses at reference depth
    sV = round(vertical_gradient * reference_depth, digits=2)
    p0 = round(pore_pressure_gradient * reference_depth, digits=2)

    println("Vertical stress: $(sV)")
    
    
    # Get friction coefficient from first fault (assume all faults have same friction)
    μ = fault_data[1]["friction_coefficient"]
    #println("Friction coefficient from first fault: $(μ)")
    #println("\nFriction coefficient from first fault: $(μ)")
    
    # Calculate horizontal stresses based on model type
    if model_type == "gradients"
        println("\nUsing gradients model (all stresses provided)")
        # Convert gradients to absolute stresses
        max_horizontal_gradient = stress_data["max_horizontal_stress"]
        min_horizontal_gradient = stress_data["min_horizontal_stress"]

        println("Max horizontal gradient: $(max_horizontal_gradient)")
        println("Min horizontal gradient: $(min_horizontal_gradient)")
        
        sH = round(max_horizontal_gradient * reference_depth, digits=2)
        sh = round(min_horizontal_gradient * reference_depth, digits=2)
        
    elseif model_type == "aphi_min" || model_type == "aphi_no_min"
        println("\nUsing A-phi model: $(model_type)")
        # Get A-phi value and calculate n and phi
        aphi = stress_data["aphi_value"]
        n, phi = calculate_n_phi(aphi)
        #println("A-phi parameters:")
        println("  A-phi value = $(aphi)")
        println("  n = $(n)")
        println("  φ = $(phi)")
        

        
        if model_type == "aphi_min"
            println("Stress model type: A-phi with min horizontal stress")
            
            sh = stress_data["min_horizontal_stress"] * reference_depth
            sH, _ = calculate_modified_aphi_stresses(n, phi, sV, sh, p0)
        else
            println("Stress model type: A-phi without min horizontal stress")
            
            # Calculate both horizontal stresses using A-phi model
            sH, sh = calculate_standard_aphi_stresses(n, phi, sV, p0, μ)
        end
    else
        error("Invalid stress model type: $model_type")
    end
    
    
    # Create stress state object with principal stresses vector [Svert, shmin, sHmax]
    stress_state = StressState([sV, sh, sH], max_stress_azimuth)
    #println("\nStress state at reference depth:")
    # print stress state
    println("  Vertical stress (Sv) = $(sV) MPa") # VERIFY UNITS!!!!!!!!!!!!!!!!
    println("  Minimum horizontal stress (Sh) = $(sh) MPa")
    println("  Maximum horizontal stress (SH) = $(sH) MPa")
    println("  Max stress azimuth = $(max_stress_azimuth) degrees")

    
    return stress_state, p0
end


# This can replace calculate_absolute_stresses (not used yet)
function ComputeStressTensor_CS_Model(Sv::Real, SHmax::Real, Shmin::Real, Pp::Real, μ::Real, Aphi::Real)
    #=
       Compute stress tensor for the critically stressed model
    =#

    df=DataFrames.DataFrame([])
    
    if Aphi > 0 && Aphi <= 1 #normal stress regimen
        df=ComputeStressTensor_CS_Normal_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=0, Shmin=Shmin, SHmax=SHmax)
    elseif  Aphi > 1 && Aphi <= 2 #strike-slip stress regimen
        df=ComputeStressTensor_CS_Strike_Slip_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=1, Shmin=Shmin, SHmax=SHmax)
    elseif Aphi > 2 && Aphi <=3 #reverse stress regimen
        df=ComputeStressTensor_CS_Reverse_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=2, Shmin=Shmin, SHmax=SHmax)
    else
        throw(ErrorException("Error!!! Aphi parameter must be between 0 and 3"))
    end

    return df
end #ComputeStressTensor_CS_Model()

# Josimar code (not used yet)
function ComputeStressTensor_CS_Normal_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=0, Shmin=[], SHmax=[])
    #=
    Compute SHmax for the critically stressed faults for normal faults
    =#

    ## Positive values means compressive stress
    Sv = abs.(Sv)
    Pp = abs.(Pp)

    if Aphi < 0 || Aphi > 1
        throw(ErrorException("Error!!! The Aphi parameter for normal faults should be between 0 and 1"))
    end

    ## Compute the critically stressed fault factor
    B = FSP3D.ComputeFrictionalStressLimitFactor(μ)

    # Compute the ϕ parameters
    ϕ = FSP3D.ComputePhiParameter(Aphi, n)

    if Shmin == 0 && SHmax == 0
        Shmin = ((Sv .- Pp) ./ B) .+ Pp
        SHmax = ϕ .* (Sv .- Shmin) .+ Shmin
    elseif Shmin == 0 && SHmax != 0
        Shmin = (SHmax .- ϕ .* Sv) ./ (1 .- ϕ)
    elseif Shmin != 0 && SHmax == 0
        SHmax = ϕ .* (Sv .- Shmin) .+ Shmin        
    end

    if any(isempty([Sv, Pp, SHmax, Shmin, Aphi, μ]))
        throw(ErrorException("Error!!! Stress is empty"))
    end

    return DataFrames.DataFrame(:Sv => Sv, :Pp => Pp, :μ => μ, :Aphi => Aphi, :SHmax => SHmax, :Shmin => Shmin)
    
end #ComputeStressTensor_CS_Normal_Faults


# Josimar code (not used yet)
function ComputeStressTensor_CS_Strike_Slip_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=1, Shmin=[], SHmax = [])
    #=
      Compute SHmax for the critically stressed faults for strike-slip faults
    =#

    ## Positive values means compressive stress
    Sv = abs.(Sv)
    Pp = abs.(Pp)
    
    if Aphi < 1 || Aphi > 2
        throw(ErrorException("Error!!! The Aphi parameter for strike slip faults should be between 1 and 2"))
    end

    ## Compute the critically stressed fault factor
    B = FSP3D.ComputeFrictionalStressLimitFactor(μ)

    ## Computing the phi parameters
    ϕ = FSP3D.ComputePhiParameter(Aphi, n)

    if Shmin == 0 && SHmax == 0
        tmp1 = (Sv ./ ϕ) .+ (Pp ./ (ϕ .* B)) .- (Pp ./ ϕ) .- (Pp ./ B) .+ Pp
        tmp2 = 1 .+ (1 ./ (ϕ .* B)) .- (1 ./ B)
        SHmax = tmp1 ./ tmp2
        Shmin = ( (SHmax .- Pp) ./ B ) .+ Pp
    elseif Shmin == 0 && SHmax != 0
        Shmin = ( (SHmax .- Pp) ./ B ) .+ Pp
    elseif Shmin != 0 && SHmax == 0
        SHmax = B .* Shmin .- Pp .* (B .- 1)
    end

    if any(isempty([Sv, Pp, SHmax, Shmin, Aphi, μ]))
        throw(ErrorException("Error!!! Stress is empty"))
    end

    return DataFrames.DataFrame(:Sv => Sv, :Pp => Pp, :μ => μ, :Aphi => Aphi, :SHmax => SHmax, :Shmin => Shmin)

end #ComputeStressTensor_CS_Strike_Slip_Faults

# Josimar code (not used yet)
function ComputeStressTensor_CS_Reverse_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=2, Shmin = [], SHmax = [])
    #=
    Compute SHmax for the critically stressed faults for normal faults
    =#

    ## Positive values means compressive stress
    Sv = abs.(Sv)
    Pp = abs.(Pp)
    
    if Aphi < 2 || Aphi > 3
        throw(ErrorException("Error!!! The Aphi parameter for normal faults should be between 0 and 1"))
    end

    ## Compute the critically stressed fault factor
    B = FSP3D.ComputeFrictionalStressLimitFactor(μ)

    # Compute the ϕ parameters
    ϕ = FSP3D.ComputePhiParameter(Aphi, n)

    if Shmin == 0 && SHmax == 0
        SHmax = B .* (Sv .- Pp) .+ Pp
        Shmin = ϕ .* (SHmax .- Sv) .+ Sv
    elseif Shmin == 0 && SHmax != 0
        Shmin = ϕ .* (SHmax .- Sv) .+ Sv
    elseif Shmin != 0 && SHmax == 0
        SHmax = (1 ./ ϕ) * (Shmin .- Sv .* (1 - ϕ))
    end

    if any(isempty([Sv, Pp, SHmax, Shmin, Aphi, μ]))
        throw(ErrorException("Error!!! Stress is empty"))
    end

    return DataFrames.DataFrame(:Sv => Sv, :Pp => Pp, :μ => μ, :Aphi => Aphi, :SHmax => SHmax, :Shmin => Shmin)
    
end #ComputeStressTensor_CS_Reverse_Faults


"""
Process each fault and calculate geomechanical parameters
"""

function process_faults(fault_data::Vector, stress_state::StressState, initial_pressure::Float64)
    results = []
    
    #for fault in fault_data
    # analyze faults in parallel using threading
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

#=
# Monte Carlo version - only calculates slip pressure
function process_faults(faults::Vector, stress_state::StressState, initial_pressure::Float64, ::Val{:monte_carlo})
    results = []
    
    for fault in faults
        strike = fault["strike"]
        dip = fault["dip"]
        friction = fault["friction_coefficient"]
        
        # Calculate fault stresses
        fault_stresses = calculate_fault_effective_stresses(stress_state, strike, dip)
        
        # Calculate slip pressure
        slip_pressure = calculate_slip_pressure(fault_stresses, friction, initial_pressure)
        
        # Store only slip pressure for Monte Carlo
        push!(results, Dict(
            "slip_pressure" => slip_pressure
        ))
    end
    
    return results
end
=#

function main()
    println("\n=== Starting Deterministic Geomechanics Analysis ===")
    
    # Parse command line arguments
    
    args = parse_commandline()
    
    # Read input JSON
    
    input_data = JSON.parsefile(args["input-json"])
    
    # Extract data
    
    fault_data = input_data["faults"]
    stress_data = input_data["stress_state"]

    #println("fault_data: $(fault_data)")

    #println("Stress data: $(stress_data)")
    
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

# To do:
# Improve error handling
# Add efficient parallel processing that prevents data-races
# Add efficient parallel processing potentially making use of Threads.@spawn and :interactive

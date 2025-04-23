module GeomechanicsDriver

using JSON
using ArgParse
using LinearAlgebra
using Base.Threads
using DataFrames
using Geodesy
using CSV
using InlineStrings

include("graphs/julia_fsp_graphs.jl")
include("core/geomechanics_model.jl")
include("TexNetWebToolLauncherHelperJulia.jl")
include("core/utilities.jl")

using .GeomechanicsModel
using .JuliaFSPGraphs
using .TexNetWebToolLauncherHelperJulia
using .Utilities

const ARGS_FILE_NAME = "args.json"
const RESULTS_FILE_NAME = "results.json"

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
            required = false
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

    # check for aphi_value and set default value if missing
    if !haskey(stress_data, "aphi_value") || stress_data["aphi_value"] === nothing
        stress_data["aphi_value"] = 0.0
    end
    #println("aphi_value: $(stress_data["aphi_value"])")
    
    
    # Calculate absolute stresses at reference depth
    sV = round(vertical_gradient * reference_depth, digits=4)
    p0 = round(pore_pressure_gradient * reference_depth, digits=4)

    
    
    # Get friction coefficient from first fault (assume all faults have same friction)
    μ = fault_data[1]["friction_coefficient"]
    #println("Friction coefficient from first fault: $(μ)")
    #println("\nFriction coefficient from first fault: $(μ)")
    
    # Calculate horizontal stresses based on model type
    if model_type == "gradients"
        @assert haskey(stress_data, "max_horizontal_stress") "Max horizontal stress required for gradients model"
        @assert haskey(stress_data, "min_horizontal_stress") "Min horizontal stress required for gradients model"
        @assert haskey(stress_data, "max_stress_azimuth") "Max stress azimuth required for gradients model"

        #println("\nUsing gradients model (all stresses provided)")
        # Convert gradients to absolute stresses
        max_horizontal_gradient = stress_data["max_horizontal_stress"]
        min_horizontal_gradient = stress_data["min_horizontal_stress"]
        
        sH = round(max_horizontal_gradient * reference_depth, digits=2)
        sh = round(min_horizontal_gradient * reference_depth, digits=2)
        
    elseif model_type == "aphi_min" || model_type == "aphi_no_min"
        println("\nUsing A-phi model: $(model_type)")
        # Get A-phi value and calculate n and phi
        aphi = stress_data["aphi_value"]
        n, phi = calculate_n_phi(aphi)
        
        
        
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

    
    return stress_state, p0
end


#=
# PORTAL VERSION (accepts df for faults instead of vector)
function calculate_absolute_stresses(stress_data::Dict, friction_coefficient::Real, stress_model_type::String)
    # Extract common parameters
    reference_depth = stress_data["reference_depth"]
    vertical_gradient = stress_data["vertical_stress"]
    pore_pressure_gradient = stress_data["pore_pressure"]
    max_stress_azimuth = stress_data["max_stress_azimuth"]
    
    # Ensure optional values are set
    stress_data["min_horizontal_stress"] = get(stress_data, "min_horizontal_stress", 0.0)
    stress_data["max_horizontal_stress"] = get(stress_data, "max_horizontal_stress", 0.0)
    stress_data["aphi_value"] = get(stress_data, "aphi_value", 0.0)

    # Calculate absolute stresses at reference depth
    sV = round(vertical_gradient * reference_depth, digits=4)
    p0 = round(pore_pressure_gradient * reference_depth, digits=4)

    # Use the passed friction coefficient
    μ = friction_coefficient

    if stress_model_type == "all_gradients"
        #println("\nUsing gradients model (all stresses provided)")
        sH = round(stress_data["max_horizontal_stress"] * reference_depth, digits=2)
        sh = round(stress_data["min_horizontal_stress"] * reference_depth, digits=2)
    elseif stress_model_type == "aphi_model" || stress_model_type == "aphi_no_min"
        println("\nUsing A-phi model: $(stress_model_type)")
        aphi = stress_data["aphi_value"]
        n, phi = calculate_n_phi(aphi)
        
        if stress_model_type == "aphi_model"
            println("Stress model type: A-phi with min horizontal stress")
            sh = stress_data["min_horizontal_stress"] * reference_depth
            sH, _ = calculate_modified_aphi_stresses(n, phi, sV, sh, p0)
        else
            println("Stress model type: A-phi without min horizontal stress")
            sH, sh = calculate_standard_aphi_stresses(n, phi, sV, p0, μ)
        end
    else
        error("Invalid stress model type: $stress_model_type")
    end

    stress_state = StressState([sV, sh, sH], max_stress_azimuth)
    return stress_state, p0
end
=#


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

function process_faults(fault_data::Vector, stress_state::GeomechanicsModel.StressState, initial_pressure::Float64, friction_coefficient::Float64; tab::String = "det_geo", dp::Vector = Vector{Float64}())
    
    # if we run it as a Monte Carlo simulation, we need to store the results in a different way  
    if tab == "prob_geo"
        # for each fault, we have a vector of results (default: 1000 for each fault)
        num_faults = length(fault_data)
        num_iterations = 1000
        results = Vector{Vector{Dict{String, Any}}}(undef, num_faults)
        for i in 1:num_faults
            results[i] = Vector{Dict{String, Any}}(undef, num_iterations)
        end
    end

    
    if tab == "det_geo"
        println("\nProcessing faults for deterministic geomechanics analysis...")
        results = Vector{Dict{String, Any}}(undef, length(fault_data)) # pre-allocate results array
        
        @threads for i in eachindex(fault_data)
            fault = fault_data[i]
            results[i] = analyze_fault(
                Float64(fault["Strike"]),
                Float64(fault["Dip"]),
                friction_coefficient,
                stress_state,
                initial_pressure,
                0.0  # No pressure change for deterministic analysis
            )

            # Add fault metadata to the result
            for key in ["FaultID", "Latitude(WGS84)", "Longitude(WGS84)", "Strike", "Dip"]
                results[i][key] = fault[key]
            end
        end

    elseif tab == "det_hydro"
        println("\nProcessing faults for deterministic hydrology analysis...")
        results = Vector{Dict{String, Any}}(undef, length(fault_data))  # Preallocate results array
        
        
        @threads for i in eachindex(fault_data)
            fault = fault_data[i]
            
            # Use the fault index (i) to access dp
            dp_this_fault = dp[i]  # Ensure dp is a vector and aligned with fault_data

            # Calculate geomechanical parameters for each fault
            results[i] = analyze_fault_hydro(
                fault["strike"],
                fault["dip"],
                friction_coefficient,
                stress_state,
                initial_pressure,
                dp_this_fault  # Use the value from dp for this fault
            )

            #println("input of analyze_fault_hydro for fault $i: $(fault)")

            # Add fault metadata to the result
            for key in ["fault_id", "strike", "dip"]
                results[i][key] = fault[key]
            end
        end
    
    end

    return results
end


function get_stress_regime(Sv::Real, Shmin::Real, SHmax::Real)
    if abs(Sv) >= abs(SHmax) && abs(SHmax) >= abs(Shmin)
        return "Normal"
    elseif abs(SHmax) >= abs(Shmin) && abs(Shmin) >= abs(Sv)
        return "Reverse"
    elseif abs(SHmax) >= abs(Sv) && abs(Sv) >= abs(Shmin)
        return "Strike-Slip"
    else
        error("Error: Stress regime could not be determined")
    end
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

function get_stress_model_type(stress_inputs::Dict)
    if stress_inputs["aphi_value"] !== nothing
        if stress_inputs["min_horizontal_stress_gradient"] !== nothing
            return "aphi_min"
        else 
            return "aphi_no_min"
        end
    else
        return "gradients"
        
    end
end

function main()

    #=
    To test with portal, use these on FSP
    fault1: x=400, y=400, strike=10, dip=40, length=10
    fault2: x=398.9400, y=400.1700, strike=20, dip=60, length=8

    These correspond to the '2_faults_test_.csv' file
    =#
    
    println("\n=== Starting Deterministic Geomechanics Process ===")
    

    scratchPath = ARGS[1]

     
    helper = TexNetWebToolLaunchHelperJulia(scratchPath)
    

    
    println("Extracting stress state values from args.json...")
    stress_inputs = Dict(
        "reference_depth" => get_parameter_value(helper, 2, "reference_depth"),
        "vertical_stress" => get_parameter_value(helper, 2, "vertical_stress"),
        #"min_horizontal_stress" => get_parameter_value(helper, 2, "min_horizontal_stress") === nothing ? nothing : get_parameter_value(helper, 2, "min_horizontal_stress"),
        "min_horizontal_stress" => get_parameter_value(helper, 2, "min_horizontal_stress"),
        #"max_horizontal_stress" => get_parameter_value(helper, 2, "max_horizontal_stress") === nothing ? nothing : get_parameter_value(helper, 2, "max_horizontal_stress"),
        "max_horizontal_stress" => get_parameter_value(helper, 2, "max_horizontal_stress"),
        "pore_pressure" => get_parameter_value(helper, 2, "pore_pressure"),
        "max_stress_azimuth" => get_parameter_value(helper, 2, "max_stress_azimuth"),
        "aphi_value" => get_parameter_value(helper, 2, "aphi_value") === nothing ? nothing : get_parameter_value(helper, 2, "aphi_value"),
        "stress_field_mode" => get_parameter_value(helper, 2, "stress_field_mode"),
        "friction_coefficient" => get_parameter_value(helper, 2, "friction_coefficient")
    )

    println("friction_coefficient from the portal: $(stress_inputs["friction_coefficient"])")
    
    # REMOVE THIS
    #=
    if stress_inputs["max_horizontal_stress"] === nothing
        add_message_with_step_index!(helper, 2, "Max Horizontal Stress Gradient is not provided, using default value of 1.22", 2)
        stress_inputs["max_horizontal_stress"] = 1.22
    end
    =#

   

    # if both aphi_value and max_horizontal_stress are not nothing, then we need to throw an error
    if stress_inputs["aphi_value"] !== nothing && stress_inputs["max_horizontal_stress"] !== nothing
        add_message_with_step_index!(helper, 2, "Aphi value and Max Horizontal Stress Gradient cannot both be provided", 2)
        throw(ErrorException("Error: Aphi value and Max Horizontal Stress Gradient cannot both be provided"))
    end


    

    # make the 'stress_model_type' of the stress_inputs get the value from the 'get_stress_model_type' function
    #stress_model_type = get_stress_model_type(stress_inputs)


    
    
    #stress_model_type = get_stress_model_type(stress_inputs)

    #println("Using stress model type: $(stress_inputs["stress_field_mode"])")
    # this is the output parameter from the portal
    set_parameter_value!(helper, 2, "stress_model_type", stress_inputs["stress_field_mode"])


    # input validation
    
    @assert stress_inputs["reference_depth"] > 0 "Reference depth must be positive"
    @assert stress_inputs["vertical_stress"] > 0 "Vertical stress must be positive"
    if stress_inputs["min_horizontal_stress"] !== nothing
        @assert stress_inputs["min_horizontal_stress"] > 0 "Minimum horizontal stress must be positive"
    end
    if stress_inputs["max_horizontal_stress"] !== nothing
        @assert stress_inputs["max_horizontal_stress"] > 0 "Maximum horizontal stress must be positive"
    end
    @assert stress_inputs["pore_pressure"] > 0 "Pore pressure must be positive"
    @assert stress_inputs["max_stress_azimuth"] >= 0 && stress_inputs["max_stress_azimuth"] <= 360 "Max stress azimuth must be between 0 and 360"
    
    if stress_inputs["aphi_value"] !== nothing
        @assert stress_inputs["aphi_value"] >= 0 && stress_inputs["aphi_value"] <= 3 "Aphi value must be between 0 and 3"
    end



    # get friction coefficient from the first fault and apply it to all faults
    #println("Extracting fault data from the CSV at the scratch path...")
    faults_csv_filepath = get_dataset_file_path(helper, 2, "faults_model_inputs_output")
    if faults_csv_filepath !== nothing
        faults_inputs = CSV.read(faults_csv_filepath, DataFrame)
    else
        error("Error: No faults dataset provided.")
    end

    
    

    println("stress_field_mode: $(stress_inputs["stress_field_mode"])")

    # Calculate absolute stresses at reference depth
    stress_state, initial_pressure = GeomechanicsModel.calculate_absolute_stresses(stress_inputs, stress_inputs["friction_coefficient"], stress_inputs["stress_field_mode"])
    
    


    stress_regime = get_stress_regime(stress_state.principal_stresses[1], stress_state.principal_stresses[2], stress_state.principal_stresses[3])

   

    # convert faults df to a vector so we can use it in the process_faults function
    fault_data = collect(eachrow(faults_inputs))
    
    # Process each fault
    results = process_faults(fault_data, stress_state, initial_pressure, stress_inputs["friction_coefficient"])

    # Print information for Fault1
    for result in results
        if string(result["FaultID"]) == "1" || string(result["FaultID"]) == "Fault1"
            println("\nFault1 Results:")
            println("  Slip Pressure: $(round(result["slip_pressure"], digits=3)) psi")
            println("  Shear Capacity Utilization: $(round(result["shear_capacity_utilization"], digits=3))")
            println("  Coulomb Failure Function: $(round(result["coulomb_failure_function"], digits=3))")
        end
    end

    # create a CSV dataframe from the results
    # for each fault (row), we have the columns: fault_id, slip_pressure, coulomb_failure_function, shear_capacity_utilization, normal_stress, shear_stress
    
    # TO DO: remove this later, remove also from the portal
    process_output_df = DataFrame(
        :fault_id => [result["FaultID"] for result in results],
        :slip_pressure => [round(result["slip_pressure"], digits=3) for result in results],
        :coulomb_failure_function => [round(result["coulomb_failure_function"], digits=3) for result in results],
        :shear_capacity_utilization => [round(result["shear_capacity_utilization"], digits=3) for result in results],
        :normal_stress => [round(result["normal_stress"], digits=3) for result in results],
        :shear_stress => [round(result["shear_stress"], digits=3) for result in results]
    )

    # read the faults_model_inputs_output dataframe from step 1
    step1_faults_output_filepath = get_dataset_file_path(helper, 2, "faults_model_inputs_output")
    step1_faults_output = CSV.read(step1_faults_output_filepath, DataFrame)

    

    # Add the results columns to the existing dataframe from step 1
    step1_faults_output.slip_pressure = [round(result["slip_pressure"], digits=3) for result in results]
    step1_faults_output.coulomb_failure_function = [round(result["coulomb_failure_function"], digits=3) for result in results]
    step1_faults_output.shear_capacity_utilization = [round(result["shear_capacity_utilization"], digits=3) for result in results]
    step1_faults_output.normal_stress = [round(result["normal_stress"], digits=3) for result in results]
    step1_faults_output.shear_stress = [round(result["shear_stress"], digits=3) for result in results]

    # Save the updated dataframe for the next step
    step2_faults_output = step1_faults_output

    
    save_dataframe_as_parameter!(helper, 2, "det_geomechanics_results", step2_faults_output)

    #add_message_with_step_index!(helper, 2, "step2_faults_output: $(step2_faults_output)", 0)

    # extract tau and sigme effective from results (to use in the Mohr diagram)
    tau_effective_faults = [result["shear_stress"] for result in results]
    sigma_effective_faults = [result["normal_stress"] for result in results]
    slip_pressure_faults = [result["slip_pressure"] for result in results]
    
    # Extract fault IDs for the Mohr diagram
    fault_ids = [string(result["FaultID"]) for result in results]
    # the line above creates a vector of type Vector{InlineStrings.String7}
    # we need to convert it to a vector of type Vector{String}
    if typeof(fault_ids[1]) == String7
        fault_ids = String.(fault_ids)
    end

    # extract strikes from each fault and store in a vector
    strikes = Float64[fault["Strike"] for fault in fault_data]
    
    # Prepare output by copying all data from step 1 and adding deterministic results
    
    output = deepcopy(faults_inputs)  # Deep copy to avoid modifying original

    
   
    
    # get data for mohr diagram plot 
    arcsDF, slipDF, faultDF = mohr_diagram_data_to_d3_portal(stress_state.principal_stresses[2], stress_state.principal_stresses[3], stress_state.principal_stresses[1], tau_effective_faults, sigma_effective_faults, initial_pressure, 1.0, 0.5, 0.0, strikes, stress_inputs["friction_coefficient"], stress_regime, slip_pressure_faults, fault_ids)
    
    save_dataframe_as_parameter!(helper, 2, "arcsDF", arcsDF)
    save_dataframe_as_parameter!(helper, 2, "slipDF", slipDF)
    save_dataframe_as_parameter!(helper, 2, "faultDF", faultDF)

   

    #write_final_args_file(helper, joinpath(helper.scratch_path, ARGS_FILE_NAME))

    # explicitly set this step's success state to true
    set_success_for_step_index!(helper, 2, true)

    # write the results to the results.json file
    write_results_file(helper)
    
    #println("\n=== Deterministic Geomechanics Analysis Complete ===\n")

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end # module GeomechanicsDriver



# To do:
# Improve error handling
# Add efficient parallel processing that prevents data-races
# Add efficient parallel processing potentially making use of Threads.@spawn and :interactive

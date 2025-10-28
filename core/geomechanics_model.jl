"""
Module for deterministic geomechanical calculations

    Workflow:
    1. Read stress and fault data
    2. Calculate normal and shear stresses on fault plane
    3. Calculate pressure required for slip
    4. Calculate slip tendency
    
    More specifically:
    1. Read stress and fault data from step 1 (Model Inputs) JSON output
    2. Calculate absolute stresses at reference depth
    3. For each fault:
        - Transform stresses to fault coordinates
        - Calculate normal and shear stresses
        - Determine slip pressure and stability metrics
    4. Output results to JSON file (used as input to step 3)
"""
module GeomechanicsModel

using LinearAlgebra
using StaticArrays

export calculate_fault_effective_stresses, calculate_slip_tendency, StressState, analyze_fault, process_faults, calculate_slip_pressure, analyze_fault_hydro, ComputeCriticalPorePressureForFailure
export calculate_modified_aphi_stresses, calculate_standard_aphi_stresses, calculate_n_phi, calculate_absolute_stresses


# Constants
const DEG2RAD = π/180.0

"""
Structure to hold stress state
Matches MATLAB implementation order: [Svert, shmin, sHmax]
"""
struct StressState
    principal_stresses::Union{AbstractVector{Float64}, AbstractVector{Integer}}  # [Svert, shmin, sHmax]
    sH_azimuth::Union{Float64, Integer}    # Azimuth of maximum horizontal stress (degrees clockwise from North)
end


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






function calculate_n_phi(aphi::Union{Float64, Int64})
    println("APhi value provided: $aphi")
    
    aphi_float = Float64(aphi)  # Convert to Float64 for calculations
    println("APhi value converted to Float64: $aphi_float")
    if aphi_float >= 0 && aphi_float < 1
        n = 0
    elseif aphi_float >= 1 && aphi_float < 2
        n = 1
    elseif aphi_float >= 2 && aphi_float <= 3
        n = 2
    else
        error("APhi value must be in range [0,3]. Got: $aphi_float")
    end
    
    phi = (aphi_float - (n + 0.5))/(-1)^n + 0.5
    
    
    return n, phi
end



# PORTAL VERSION (accepts df for faults instead of vector)
function calculate_absolute_stresses(stress_data::Dict, friction_coefficient::Real, stress_model_type::String)
    # common stressparameters
    reference_depth = stress_data["reference_depth"]
    vertical_gradient = stress_data["vertical_stress"]
    pore_pressure_gradient = stress_data["pore_pressure"]
    max_stress_azimuth = stress_data["max_stress_azimuth"]

    
    
    # Ensure optional values are set
    stress_data["min_horizontal_stress"] = get(stress_data, "min_horizontal_stress", nothing)
    stress_data["max_horizontal_stress"] = get(stress_data, "max_horizontal_stress", nothing)
    stress_data["aphi_value"] = get(stress_data, "aphi_value", nothing)

    

    # Calculate absolute stresses at reference depth
    sV = round(vertical_gradient * reference_depth, digits=4)
    p0 = round(pore_pressure_gradient * reference_depth, digits=4)

    # Use the passed friction coefficient
    μ = friction_coefficient

    if stress_model_type == "all_gradients" || stress_model_type == "gradients"
        #println("\nUsing gradients model (all stresses provided), all_gradients")
        sH = round(stress_data["max_horizontal_stress"] * reference_depth, digits=2)
        sh = round(stress_data["min_horizontal_stress"] * reference_depth, digits=2)
    elseif stress_model_type == "aphi_model" || stress_model_type == "aphi_no_min" || stress_model_type == "aphi_min"
        #println("\nUsing A-phi model: $(stress_model_type)")
        aphi = stress_data["aphi_value"]
        n, phi = calculate_n_phi(aphi)
        
        if stress_data["min_horizontal_stress"] !== nothing
            #println("Stress model type: A-phi with min horizontal stress, aphi_min")
            sh = stress_data["min_horizontal_stress"] * reference_depth
            sH, _ = calculate_modified_aphi_stresses(n, phi, sV, sh, p0)
        else
            
            #println("Stress model type: A-phi without min horizontal stress, aphi_no_min")
            sH, sh = calculate_standard_aphi_stresses(n, phi, sV, p0, μ)
        end
    else
        error("Invalid stress model type: $stress_model_type")
    end

    stress_state = StressState([sV, sh, sH], max_stress_azimuth)
    return stress_state, p0
end





"""
Calculate the stress transformation matrix for given strike and dip
Implements fault normal vector calculation from MATLAB mohrs_3D.m
"""
function compute_unit_vectors(strike::Float64, dip::Float64)
    # Convert angles to radians
    strike_rad = deg2rad(strike)
    dip_rad = deg2rad(dip)
    
    # normal vector
    n = [-sin(strike_rad) * sin(dip_rad), cos(strike_rad) * sin(dip_rad), -cos(dip_rad)]
    # strike vector
    ns = [cos(strike_rad), sin(strike_rad), 0.0]
    # dip vector 
    nd = [-sin(strike_rad) * cos(dip_rad), cos(strike_rad) * cos(dip_rad), sin(dip_rad)]
    
    # Create transformation matrix Q
    return @SMatrix [n[1] n[2] n[3]; ns[1] ns[2] ns[3]; nd[1] nd[2] nd[3]]
end

"""
Calculate normal and shear stresses on fault plane
Direct implementation of MATLAB mohrs_3D.m calc_ppfail function
Returns sames results as AuxFunctionComputePcriticalSample from FSP3D.jl (sigmaN and tauN)
"""
function calculate_fault_effective_stresses(strike::Union{Float64, Integer}, dip::Union{Float64, Integer}, stress_state::StressState, p0::Union{Float64, Integer}, dp::Union{Float64, Integer})

    # extract azimuth from stress state
    az = stress_state.sH_azimuth

    # Convert angles from degrees to radians
    az_rad = deg2rad.(az)
    str_rad = deg2rad.(strike)
    dip_rad = deg2rad.(dip)

    # Cos and sin of azimuth
    cos_az = cos.(az_rad)
    sin_az = sin.(az_rad)

    # Cos and sin of strike and dip angles
    cs = cos.(str_rad)
    ss = sin.(str_rad)
    cd = cos.(dip_rad)
    sd = sin.(dip_rad)

    
    

    # Total stresses from input data
    Svert = stress_state.principal_stresses[1]
    shmin = stress_state.principal_stresses[2]
    sHmax = stress_state.principal_stresses[3]

    # Poisson's ratio effect factor (matches MATLAB exactly)
    biot = 1.0
    nu = 0.5

    # Factor based on nu for horizontal stress change
    # Biot is currently hardcoded to 1.0, and is probably good for fault stress analysis,
    # similar to Abaqus where effective stress is always total-pp even when Biot effect is included via porous bulk moduli option

    f = biot * nu / (1 - nu)


    # Effective stresses at current pressure
    """
    Note that factor f is used for dp because dp is pressure change due to
    injection upto present time - horizontal stresses will have Poisson's
    ratio effect (and Biot coeff, if other than 1; but see comments above)
    However, p0 is initial pressure before injection, so factor f is not needed
    """
    s11 = shmin .* cos_az.^2 + sHmax .* sin_az.^2 - p0 - f .* dp #VERIFIED
    s22 = shmin .* sin_az.^2 + sHmax .* cos_az.^2 - p0 - f .* dp #VERIFIED
    s33 = Svert - p0 - dp #VERIFIED
    s12 = (sHmax - shmin) .* cos_az .* sin_az #VERIFIED

    


    # Components of the unit normal vector to fault planes
    n1 = sd .* cs
    n2 = -sd .* ss
    n3 = cd  # not used below
    

    # Precompute repeated terms
    n1_sq = n1.^2
    n2_sq = n2.^2
    s11_s33 = s11 .- s33
    s22_s33 = s22 .- s33
    n1_n2 = n1 .* n2
    n1_cubed = n1.^3
    n2_cubed = n2.^3
    n1_sq_n2_sq = n1_sq .* n2_sq

    
    """
    Note that tau_fault is absolute magnitude, but sig_fault is signed
    value (which should be positive for most cases, but, with high pore
    pressure, resolved normal stress on fault can become negative any of 
    the principal stresses is negative)
    """
    # Calculate the expression inside sqrt first to avoid domain errors
    # sometimes this expression would yield really small negative numbers due to floating point precision and sqrt would throw domain errors
    sqrt_argument = n2_sq .* (s12.^2 .- (-1 .+ n2_sq) .* s22_s33.^2) .-
        n1_sq .* n1_sq .* s11_s33.^2 .+
        4 .* n1_cubed .* n2 .* s12 .* (-s11_s33) .+
        2 .* n1_n2 .* s12 .* (s11 .+ s22 .- 2 .* n2_sq .* s22 .+ 2 .* (-1 .+ n2_sq) .* s33) .+
        n1_sq .* (
            s11.^2 .+ (1 .- 4 .* n2_sq) .* s12.^2 .-
            2 .* s11 .* (n2_sq .* s22_s33 .+ s33) .+
            s33 .* (2 .* n2_sq .* s22_s33 .+ s33)
        )
    
    # Ensure sqrt argument is non-negative to avoid domain errors
    # If negative (usually due to numerical precision), clamp to zero
    sqrt_argument = max.(sqrt_argument, 0.0)
    
    # shear stress on fault plane
    tau_normal = sqrt.(sqrt_argument)

    # normal stress on fault plane (signed)
    sig_normal = 2 .* n1_n2 .* s12 .+ n1_sq .* s11_s33 .+ n2_sq .* s22_s33 .+ s33

    # if sig_normal or tau_normal are negative, set them to 0.0
    sig_normal = max.(sig_normal, 0.0)
    tau_normal = max.(tau_normal, 0.0)

    #println("sig_normal: ", sig_normal)
    #println("tau_normal: ", tau_normal)
    #println("s11: ", s11)
    #println("s22: ", s22)
    #println("s33: ", s33)
    #println("s12: ", s12)
    #println("n1: ", n1)
    #println("n2: ", n2)

    return sig_normal, tau_normal, s11, s22, s33, s12, n1, n2
end


"""
Calculate pore pressure required for fault slip using quadratic equation approach
    Inputs: sig_fault: Total normal stress
            tau_fault: Shear stress
            μ: Target friction coefficient
            p0: Initial pore pressure
            biot: Biot coefficient
            nu: Poisson's ratio
            dp: Current pressure perturbation
            s11, s22, s33, s12: Stress components
            n1, n2: on fault normal components (after being resolved) [x y coordinates on Mohr circle]
    Output: pore_pressure_to_slip: Pore pressure required for fault slip
"""

# Equivalent to Josimar's ComputeCriticalPorePressureForFailure
# Currently not using this
function calculate_slip_pressure(sig_fault::Float64, tau_fault::Float64, mu::Float64, p0::Float64, 
                               biot::Float64=1.0, nu::Float64=0.5, dp::Float64=0.0, s11::Float64=0.0, s22::Float64=0.0, s33::Float64=0.0, s12::Float64=0.0, n1::Float64=0.0, n2::Float64=0.0)
    

    # Calculate mobilized friction coefficient
    mobmu = sig_fault > 0.0 ? tau_fault/sig_fault : mu
    
    # Calculate quadratic coefficients for dp that brings fault to failure
    f = biot * nu/(1 - nu)  # Poisson effect factor


    # Coefficients of quadratic equation A*dp^2 + B*dp + C = 0
    # such that the mobilized friction coefficient on fault = mu
    # solving for chnage in pp that makes the fault point intersect with frictional slip line
    # ASK RALL ABOUT THIS
    C = -4 * (1 + mu^2) * n1^3 * n2 * s12 * (s11 - s33) -
        (1 + mu^2) * n1^4 * (s11 - s33)^2 -
        (1 + mu^2) * n2^4 * (s22 - s33)^2 -
        mu^2 * s33^2 +
        2 * n1 * n2 * s12 * (s11 + (1 - 2 * (1 + mu^2) * n2^2) * s22 + 2 * (1 + mu^2) * (-1 + n2^2) * s33) +
        n2^2 * (s12^2 + s22^2 - 2 * (1 + mu^2) * s22 * s33 + (1 + 2 * mu^2) * s33^2) +
        n1^2 * (s11^2 + (1 - 4 * (1 + mu^2) * n2^2) * s12^2 - 2 * (1 + mu^2) * s11 * (n2^2 * (s22 - s33) + s33) + s33 * (2 * (1 + mu^2) * n2^2 * (s22 - s33) + s33 + 2 * mu^2 * s33))

    B = 2 * (2 * (-1 + f) * (1 + mu^2) * n1^3 * n2 * s12 +
        2 * n1 * n2 * (-(1 + mu^2) * (-1 + n2^2) + f * (-1 + (1 + mu^2) * n2^2)) * s12 +
        (-1 + f) * (1 + mu^2) * n1^4 * (s11 - s33) +
        (-1 + f) * (1 + mu^2) * n2^4 * (s22 - s33) +
        mu^2 * s33 +
        n2^2 * ((1 - f + mu^2) * s22 + (-1 + f - 2 * mu^2 + f * mu^2) * s33) +
        n1^2 * ((-(1 + mu^2) * (-1 + n2^2) + f * (-1 + (1 + mu^2) * n2^2)) * s11 + (-1 + f) * (1 + mu^2) * n2^2 * (s22 - 2 * s33) + (-1 + f - 2 * mu^2 + f * mu^2) * s33))

    A = -mu^2 * (1 + (-1 + f) * n1^2 + (-1 + f) * n2^2)^2 -
        (-1 + f)^2 * (n1^4 + n2^2 * (-1 + n2^2) + n1^2 * (-1 + 2 * n2^2))

    
    Bsq_minus_4AC = B^2 - 4 * A * C

    # TO DO: check that this won't break anything else
    #=
    Here we had:
    ppfail1 = (-B - sqrt(Bsq_minus_4AC)) / (2 * A)
    ppfail2 = (-B + sqrt(Bsq_minus_4AC)) / (2 * A)

    However, Bsq_minus_4AC can be negative in some cases,
    which would cause the sqrt to throw an error. 

    Moved those two lines to the 'else' block below

    =#

    


    # For cases with no solution (Bsq_minus_4AC < 0), use horizontal distance
    if Bsq_minus_4AC < 0
        ppfail_horiz_dist = sig_fault - tau_fault / mu
        ppfail1 = -ppfail_horiz_dist
        ppfail2 = ppfail_horiz_dist
    else
        ppfail1 = (-B - sqrt(Bsq_minus_4AC)) / (2 * A)
        ppfail2 = (-B + sqrt(Bsq_minus_4AC)) / (2 * A)
    end

    # Select appropriate solution based on mobilized friction
    if mobmu < mu  # Fault below failure line
        if ppfail1 > 0 && ppfail2 > 0
            # If both roots positive, choose smaller one
            ppfail = min(ppfail1, ppfail2)
        elseif ppfail1 < 0 && ppfail2 < 0
            # Both roots negative - error condition
            error("Pressure to slip calculation error - no positive solution found")
        else
            # One positive, one negative - choose positive
            ppfail = ppfail1 > 0 ? ppfail1 : ppfail2
        end
    elseif mobmu > mu  # Fault above failure line
        if ppfail1 < 0 && ppfail2 < 0
            # Both negative - choose smaller magnitude
            ppfail = abs(ppfail1) < abs(ppfail2) ? ppfail1 : ppfail2
        elseif ppfail1 > 0 && ppfail2 > 0
            # Both positive (unusual case) - use horizontal distance
            ppfail = -(sig_fault - tau_fault / mu)
        else
            # One positive, one negative - choose negative
            ppfail = ppfail1 < 0 ? ppfail1 : ppfail2
        end
    else  # mobmu == mu
        ppfail = 0.0
    end

    
    return ppfail
end


# Josimar's simplified function (only run in geomechanics since hydrology needs dp as input)
# Original FSP equivalent: calculate_slip_pressure
# I have verified that it yields the same results as the original function
function ComputeCriticalPorePressureForFailure(
        sig_fault::Union{Float64, Integer},
        tau_fault::Union{Float64, Integer},
        mu::Union{Float64, Integer},
        p0::Union{Float64, Integer},
        biot::Union{Float64, Integer} = 1.0,
        nu::Union{Float64, Integer} = 0.5,
        dp::Union{Float64, Integer} = 1.0
    )

    # input validation
    if sig_fault < 0.0
        error("Error: Normal stress 'sig_fault' must be positive. Got: $sig_fault")
    end

    if mu <= 0.0
        error("Error: Friction coefficient 'mu' must be positive. Got: $mu")
    end

    # Biot coefficient is between 0.6 and 0.8 for most rocks
    if biot < 0.0 || biot > 1.0
        error("Error: Poisson's ratio 'nu' must be between 0 and 1. Got: $nu")
    end


    # Calculate mobilized friction coefficient
    mobmu = sig_fault > 0.0 ? abs(tau_fault/sig_fault) : mu

    if mobmu < mu # fault below failure line
        Pcritical = ((mu - mobmu) / mu) * abs(sig_fault)
    else
        # fault above failure line: already slipping 
        Pcritical = 0.0
    end

    return Pcritical
end
    






#Calculate slip tendency (tau_fault/sig_fault)
# How close is the fault to failure right now
function calculate_slip_tendency(sig_fault::Float64, tau_fault::Float64, p0::Float64)
    # sig_fault is already the effective stress
    # Avoid division by zero
    if sig_fault ≤ 0.0
        return 1.0  # Maximum slip tendency
    end
    
    return tau_fault/sig_fault
end


#Calculate shear capacity utilization (SCU)

function calculate_scu(sig_fault::Float64, tau_fault::Float64, μ::Float64)
    
    scu = tau_fault/(μ*sig_fault)
    
    # Avoid division by zero and cap at 1.0
    if !isfinite(scu) || scu > 1.0
        return 1.0
    end
    
    return scu
end


#Calculate Coulomb failure function (CFF)

function calculate_cff(sig_fault::Float64, tau_fault::Float64, μ::Float64)
    cff = tau_fault - μ*sig_fault
    return cff
end

"""
For a single fault, runs functions to get:
- pore pressure required for slip
- slip tendency
- Coulomb Failure Function
- Shear Capacity Utilization


"""
function analyze_fault(strike::Float64, dip::Float64, friction::Float64,
                      stress_state::StressState, p0::Float64, dp::Float64)
    
    # Calculate normal and shear stresses on fault plane (sig_normal, tau_normal)
    sig_fault, tau_fault, s11, s22, s33, s12, n1, n2 = calculate_fault_effective_stresses(strike, dip, stress_state, p0, dp) # AuxFunctionComputePcriticalSample --> FSP3D.jl

    # Calculate pp to bring fault to failure
    #slip_pressure = calculate_slip_pressure(sig_fault, tau_fault, friction, p0, 1.0, 0.5, dp, s11, s22, s33, s12, n1, n2) # compare to ComputeCriticalPorePressureForFailure --> FSP3D.jl

    # Josimar's simplified function
    slip_pressure = ComputeCriticalPorePressureForFailure(sig_fault, tau_fault, friction, p0, 1.0, 0.5)

    # Calculate slip tendency
    slip_tendency = calculate_slip_tendency(sig_fault, tau_fault, p0)
    
    # Calculate Coulomb Failure Function
    cff = calculate_cff(sig_fault, tau_fault, friction)
    
    # Calculate Shear Capacity Utilization using total stress
    scu = calculate_scu(sig_fault, tau_fault, friction)

    # REMOVE THIS
    #=
    if strike == 30.0 && dip == 75.0
        #println("Fault1 sig_fault: ", sig_fault)
        #println("Fault1 tau_fault: ", tau_fault)
        #println("Fault1 slip_pressure: ", slip_pressure)
        #println("Fault1 slip_tendency: ", slip_tendency)
        #println("Fault1 cff: ", cff)
        #println("Fault1 scu: ", scu)
    end
    =#

    # print the slip pressure, slip tendency, cff, and scu for each fault
    #println("slip pressure: ", slip_pressure)
    #println("slip tendency: ", slip_tendency)
    #println("cff: ", cff)
    #println("scu: ", scu)

    
    return Dict(
        "normal_stress" => sig_fault,
        "shear_stress" => tau_fault,
        "slip_pressure" => slip_pressure,
        "slip_tendency" => round(slip_tendency, digits=4),
        "coulomb_failure_function" => round(cff),
        "shear_capacity_utilization" => round(scu, digits=4)
    )
    

end


function analyze_fault_hydro(strike::Float64, dip::Float64, friction::Float64,
    stress_state::StressState, p0::Float64, dp::Float64)
    # Calculate normal and shear stresses on fault plane (sig_normal, tau_normal)
    sig_fault, tau_fault, s11, s22, s33, s12, n1, n2 = calculate_fault_effective_stresses(strike, dip, stress_state, p0, dp) # AuxFunctionComputePcriticalSample --> FSP3D.jl


    # Josimar's simplified function
    #slip_pressure = ComputeCriticalPorePressureForFailure(sig_fault, tau_fault, friction, p0, 1.0, 0.5)

    # Calculate pp to bring fault to failure (MATLAB approach)
    slip_pressure = calculate_slip_pressure(sig_fault, tau_fault, friction, p0, 1.0, 0.5, dp, s11, s22, s33, s12, n1, n2) # compare to ComputeCriticalPorePressureForFailure --> FSP3D.jl


    return Dict(
        "normal_stress" => sig_fault,
        "shear_stress" => tau_fault,
        "slip_pressure" => slip_pressure
)
end





#=

"""
Process fault and stress data to calculate geomechanical parameters (Monte Carlo version)
"""
function analyze_fault(strike::Float64, dip::Float64, friction::Float64,
                      stress_state::StressState, p0::Float64, dp::Float64, MC::Bool)
    # Calculate normal and shear stresses on fault plane (sig_normal, tau_normal)
    sig_fault, tau_fault, s11, s22, s33, s12, n1, n2 = calculate_fault_effective_stresses(strike, dip, stress_state, p0, dp)

    # Calculate pp to bring fault to failure
    slip_pressure = calculate_slip_pressure(sig_fault, tau_fault, friction, p0, 1.0, 0.5, dp, s11, s22, s33, s12, n1, n2)
       
    return Dict(
        "normal_stress" => round(sig_fault, digits=2),
        "shear_stress" => round(tau_fault, digits=2),
        "slip_pressure" => round(slip_pressure, digits=2)
    )
end




# Monte Carlo version - only calculates slip pressure
function process_faults(faults::Vector, stress_state::StressState, initial_pressure::Float64, ::Val{:monte_carlo})
    results = []
    
    for fault in faults
        strike = fault["strike"]
        dip = fault["dip"]
        friction = fault["friction_coefficient"]
        
        fault_stresses = calculate_fault_effective_stresses(stress_state, strike, dip)
        slip_pressure = calculate_slip_pressure(fault_stresses, friction, initial_pressure)
        
        push!(results, Dict(
            "slip_pressure" => slip_pressure
        ))
    end
    
    return results
end
=#

end # module

# Optimize tau_normal and sig_normal calculations


# NOTE: use strike to determine direction and start/end points for WKT format

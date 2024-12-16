"""
Deterministic Geomechanics Step

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
        - Determine slip pressure, scu, cff
    4. Output results to JSON file (used as input to step 3)
"""
module GeomechanicsModel

using LinearAlgebra
using StaticArrays

export calculate_fault_effective_stresses, calculate_slip_tendency, StressState, analyze_fault, process_faults, calculate_slip_pressure

# Constants
const DEG2RAD = π/180.0

"""
Structure to hold stress state
Matches MATLAB implementation order: [Svert, shmin, sHmax]
"""
struct StressState
    principal_stresses::Vector{Float64}  # [Svert, shmin, sHmax]
    sH_azimuth::Float64    # Azimuth of maximum horizontal stress (degrees clockwise from North)
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
function calculate_fault_effective_stresses(strike::Float64, dip::Float64, stress_state::StressState, p0::Float64, dp::Float64)

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
    f = biot * nu / (1 - nu)


    # Effective stresses at current pressure
    s11 = shmin .* cos_az.^2 + sHmax .* sin_az.^2 - p0 - f .* dp
    s22 = shmin .* sin_az.^2 + sHmax .* cos_az.^2 - p0 - f .* dp
    s33 = Svert - p0 - dp
    s12 = (sHmax - shmin) .* cos_az .* sin_az

    

    # Components of the unit normal vector to fault planes
    n1 = sd .* cs
    n2 = -sd .* ss
    n3 = cd  # n3 is not explicitly used below
    

    # Shear stress resolved on the fault
    tau_normal = sqrt.(
        n2.^2 .* (s12.^2 - (-1 .+ n2.^2) .* (s22 .- s33).^2) 
        - n1.^4 .* (s11 .- s33).^2 
        + 4 .* n1.^3 .* n2 .* s12 .* (-s11 .+ s33)
        + 2 .* n1 .* n2 .* s12 .* (s11 .+ s22 .- 2 .* n2.^2 .* s22 .+ 2 .* (-1 .+ n2.^2) .* s33)
        + n1.^2 .* (
            s11.^2 .+ (1 .- 4 .* n2.^2) .* s12.^2 
            - 2 .* s11 .* (n2.^2 .* (s22 .- s33) .+ s33) 
            + s33 .* (2 .* n2.^2 .* (s22 .- s33) .+ s33)
        )
    )

    # Normal stress resolved on the fault (signed)
    sig_normal = 2 .* n1 .* n2 .* s12 .+ n1.^2 .* (s11 .- s33) .+ n2.^2 .* (s22 .- s33) .+ s33


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
            n1, n2: Fault normal components
    Output: pore_pressure_to_slip: Pore pressure required for fault slip
"""
function calculate_slip_pressure(sig_fault::Float64, tau_fault::Float64, mu::Float64, p0::Float64, 
                               biot::Float64=1.0, nu::Float64=0.5, dp::Float64=0.0, s11::Float64=0.0, s22::Float64=0.0, s33::Float64=0.0, s12::Float64=0.0, n1::Float64=0.0, n2::Float64=0.0)
    

    # Calculate mobilized friction coefficient
    mobmu = sig_fault > 0.0 ? tau_fault/sig_fault : mu
    
    # Calculate quadratic coefficients for dp that brings fault to failure
    f = biot * nu/(1 - nu)  # Poisson effect factor
    
    # Coefficients of quadratic equation A*dp^2 + B*dp + C = 0
    # such that the mobilized friction coefficient on fault = mu
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

    ppfail1 = (-B - sqrt(Bsq_minus_4AC)) / (2 * A)
    ppfail2 = (-B + sqrt(Bsq_minus_4AC)) / (2 * A)


    # For cases with no solution (Bsq_minus_4AC < 0), use horizontal distance
    if Bsq_minus_4AC < 0
        ppfail_horiz_dist = sig_fault - tau_fault / mu
        ppfail1 = -ppfail_horiz_dist
        ppfail2 = ppfail_horiz_dist
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

    
    return round(ppfail, digits=2)
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
Process fault and stress data 

"""
function analyze_fault(strike::Float64, dip::Float64, friction::Float64,
                      stress_state::StressState, p0::Float64, dp::Float64)
    # Calculate normal and shear stresses on fault plane (sig_normal, tau_normal)
    sig_fault, tau_fault, s11, s22, s33, s12, n1, n2 = calculate_fault_effective_stresses(strike, dip, stress_state, p0, dp) # AuxFunctionComputePcriticalSample --> FSP3D.jl

    # Calculate pp to bring fault to failure
    slip_pressure = calculate_slip_pressure(sig_fault, tau_fault, friction, p0, 1.0, 0.5, dp, s11, s22, s33, s12, n1, n2) # compare to ComputeCriticalPorePressureForFailure --> FSP3D.jl
    
    # Calculate slip tendency
    slip_tendency = calculate_slip_tendency(sig_fault, tau_fault, p0)
    
    # Calculate Coulomb Failure Function
    cff = calculate_cff(sig_fault, tau_fault, friction)
    
    # Calculate Shear Capacity Utilization using total stress
    scu = calculate_scu(sig_fault, tau_fault, friction)
    
    return Dict(
        "normal_stress" => round(sig_fault, digits=2),
        "shear_stress" => round(tau_fault, digits=2),
        "slip_pressure" => round(slip_pressure, digits=2),
        "slip_tendency" => round(slip_tendency, digits=2),
        "coulomb_failure_function" => round(cff),
        "shear_capacity_utilization" => round(scu, digits=2)
    )
end

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



# Default process_faults for deterministic analysis
function process_faults(faults::Vector, stress_state::StressState, initial_pressure::Float64)
    results = []
    
    for fault in faults
        strike = fault["strike"]
        dip = fault["dip"]
        friction = fault["friction_coefficient"]
        
        fault_stresses = calculate_fault_effective_stresses(stress_state, strike, dip)
        slip_pressure = calculate_slip_pressure(fault_stresses, friction, initial_pressure)
        scu = calculate_scu(fault_stresses, friction)
        cff = calculate_cff(fault_stresses, friction)
        
        push!(results, Dict(
            "slip_pressure" => slip_pressure,
            "scu" => scu,
            "cff" => cff
        ))
    end
    
    return results
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

end # module

# Optimize tau_normal and sig_normal calculations

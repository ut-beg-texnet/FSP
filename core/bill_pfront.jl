module BillPFront

export pressureScenario_constant_rate, pressureScenario_Rall, pressureScenario_Rall_scalar

using LinearAlgebra
using SpecialFunctions
#using Profile
#using ProfileView



"""
  pressureScenario_constant_rate(bpds, days, r_meters, STRho)

Computes pressure change (in PSI) at distances `r_meters` from a well
that has been injecting at a **constant** rate. The Theis equation is used.

Inputs:
- bpds: Injection rates in bbl/day (all the same).
- days: The day-steps for injection. (We'll only use the max day.)
- r_meters: Distances from the well in meters.
- STRho: Tuple `(S, T, rho)` -> storativity, transmissivity, fluid density.

Output:
- A vector of pressure changes in PSI, one per distance in `r_meters`.
"""
function pressureScenario_constant_rate(bpds::Vector{Float64},
                                        days::Vector{Float64},
                                        r_meters::AbstractVector{Float64},
                                        STRho::Tuple{Float64, Float64, Float64})

    S, T, rho = STRho

    # 1) Find the final time (max day) in seconds
    maxDay = maximum(days)
    t_final_sec = maxDay * 24 * 3600

    # 2) If all rates are constant, we just pick the same injection rate
    #    or use the last entry. We'll assume they are all the same:
    Q_bpd = bpds[end]  # bbl/day
    # Convert bbl/day -> m^3/s
    Q_m3s = Q_bpd * 1.84013e-6

    if isempty(bpds)
        # no injection data, return zeros
        return zeros(length(r_meters))
    end

    if t_final_sec <= 0
        # no injection data, return zeros
        return zeros(length(r_meters))
    end

    # 3) Precompute: u(r) = r^2 * S / (4 T t_final)
    ppp = (r_meters .^ 2 .* S) ./ (4.0 * T * t_final_sec)

    # 4) Theis solution for head (Theis equation, single time t_final)
    #    head(r) = (Q/(4πT)) * expint(ppp)
    head = (Q_m3s / (4 * pi * T)) .* expint.(ppp)

    # 5) Convert hydraulic head to pressure in PSI:
    #    dp (Pa) = rho * g * head
    #    dp (psi) = dp (Pa) / 6894.76
    dp_pascals = head .* (rho * 9.81)
    dp_psi = dp_pascals ./ 6894.76

    # 6) Handle any infinities or NaN
    dp_psi[.!isfinite.(dp_psi)] .= 0.0

    return dp_psi
end





"""
Inputs:
- bpds: injection rates [bbl/day], length = N
- days: times [days from start], length = N
  (already truncated up to `year_of_interest` if you want 
   to stop at that year)
- r_meters: radial distances from the well [meters], can be 1D or 2D
- STRho = (S, T, rho): storativity [dimensionless], 
                     transmissivity [m²/s], 
                     fluid density [kg/m³]
- evaluation_days_from_start: Optional parameter to specify evaluation time in days from start
                              If not provided, defaults to maximum(days)

Output:
- A Vector or Matrix of pressure change in PSI 
  (same shape as r_meters if it was 2D).

Logic:
1) Convert final time = evaluation_days_from_start to seconds,
   or if not provided, use maximum(days) * 86400.
2) Step superposition:
   For i in 2..N:
     - dt = (t_final_sec - days[i-1]*86400)
     - If dt>0, compute dimensionless u = (r²*S)/(4*T*dt)
     - Theis well_func = expint.(u)
     - Add well_func * [ Q(i) - Q(i-1) ] to running sum
       where Q(i) = bpds[i]*1.84013e-6 [m³/s].
3) Convert final "head" to pressure in PSI.
4) Return as dp_psi, same shape as r_meters.
"""
function pressureScenario_Rall(
    bpds::Vector{Float64},
    days::Vector{Float64},
    r_meters::AbstractVector{Float64},
    STRho::Tuple{Float64, Float64, Float64},
    evaluation_days_from_start::Union{Float64, Nothing} = nothing
) :: Array{Float64}

    # 0) If there's no data, return zeros
    if isempty(bpds) || isempty(days)
        return zeros(size(r_meters))
    end

    
    S, T, rho = STRho
    g = 9.81  # gravitational acceleration

    # Possibly reshape r_meters to a column vector if needed
    original_size = size(r_meters)
    reshaped = false
    if ndims(r_meters) == 2 && original_size[2] != 1
        r_meters = reshape(r_meters, :, 1)
        reshaped = true
    end

    # 1) Convert final time: t_final_sec
    # Either use the provided evaluation time or the maximum injection time
    t_final_sec = isnothing(evaluation_days_from_start) ? 
                  maximum(days) * 86400.0 : 
                  evaluation_days_from_start * 86400.0

    # 2) Convert rates to m³/s
    Q_m3s = bpds .* 1.84013e-6


    # 3) Build ΔQ array, plus a time array for each step
    #    step_times[i] = days[i], step_dQ[i] = Q_m3s[i] - Q_m3s[i-1],
    #    with Q_m3s[-1] = 0 (initial jump).
    n = length(Q_m3s)
    step_times = Vector{Float64}(undef, n)
    step_dQ    = Vector{Float64}(undef, n)

    # First step: from Q=0 to Q_m3s[1] at time days[1]
    step_times[1] = days[1]
    step_dQ[1]    = Q_m3s[1]  # for example Q(1) - 0 for the first step

    for i in 2:n
        step_times[i] = days[i]
        step_dQ[i]    = Q_m3s[i] - Q_m3s[i-1]
    end

    # 4) Accumulate Theis contributions
    #    Preallocate buffers once to avoid per-step heap allocation inside the hot loop.
    nr = length(r_meters)
    tstep_sum = zeros(Float64, nr)
    u_buf     = Vector{Float64}(undef, nr)
    wf_buf    = Vector{Float64}(undef, nr)

    inv_4T = 1.0 / (4.0 * T)
    r2 = r_meters .^ 2  # computed once; safe because r_meters doesn't change inside the loop

    for i in 1:n
        dQ = step_dQ[i]
        if dQ == 0
            continue
        end

        t_i_sec = step_times[i] * 86400.0
        dt = t_final_sec - t_i_sec

        if dt <= 0
            continue
        end

        # Theis dimensionless parameter: u(r) = r² * S / (4 T dt)
        coeff = S * inv_4T / dt
        @. u_buf = r2 * coeff

        # Well function W(u) = expint(u) — computed in-place into wf_buf
        @. wf_buf = expint(u_buf)

        # Accumulate: tstep_sum += W(u) * dQ
        @. tstep_sum += wf_buf * dQ
    end

    # 5) Convert accumulated head to PSI in-place (reuse tstep_sum as dp_psi)
    inv_4piT        = 1.0 / (4π * T)
    rho_g_over_psi  = rho * g / 6894.76
    @. tstep_sum = tstep_sum * inv_4piT * rho_g_over_psi

    # Clean up NaN/Inf and negatives in-place
    @inbounds for k in eachindex(tstep_sum)
        v = tstep_sum[k]
        tstep_sum[k] = (isfinite(v) && v > 0.0) ? v : 0.0
    end

    dp_psi = tstep_sum

    # Reshape back if needed
    if reshaped
        dp_psi = reshape(dp_psi, original_size)
    end

    return dp_psi
end




"""
    pressureScenario_Rall_scalar(bpds, days, r_meters, STRho, evaluation_days_from_start) -> Float64

Scalar specialisation of `pressureScenario_Rall` for a single radial distance.
Avoids allocating a length-1 vector and running vector `expint` — used in the
fault-level hydrology path (MC and deterministic) where r_meters is a plain Float64.
"""
function pressureScenario_Rall_scalar(
    bpds::Vector{Float64},
    days::Vector{Float64},
    r_meters::Float64,
    STRho::Tuple{Float64, Float64, Float64},
    evaluation_days_from_start::Union{Float64, Nothing} = nothing
) :: Float64

    if isempty(bpds) || isempty(days)
        return 0.0
    end

    S, T, rho = STRho
    g = 9.81

    t_final_sec = isnothing(evaluation_days_from_start) ?
                  maximum(days) * 86400.0 :
                  evaluation_days_from_start * 86400.0

    n = length(bpds)
    inv_4T   = 1.0 / (4.0 * T)
    r2       = r_meters * r_meters
    tstep_sum = 0.0

    # dQ[1] = bpds[1] * 1.84013e-6, for i>1 dQ[i] = (bpds[i]-bpds[i-1]) * 1.84013e-6
    prev_Q = 0.0
    for i in 1:n
        cur_Q = bpds[i] * 1.84013e-6
        dQ    = cur_Q - prev_Q
        prev_Q = cur_Q

        if dQ == 0.0
            continue
        end

        t_i_sec = days[i] * 86400.0
        dt = t_final_sec - t_i_sec
        if dt <= 0.0
            continue
        end

        u = r2 * S * inv_4T / dt
        tstep_sum += expint(u) * dQ
    end

    inv_4piT       = 1.0 / (4.0 * π * T)
    rho_g_over_psi = rho * g / 6894.76
    dp_psi = tstep_sum * inv_4piT * rho_g_over_psi

    return (isfinite(dp_psi) && dp_psi > 0.0) ? dp_psi : 0.0
end


end # module




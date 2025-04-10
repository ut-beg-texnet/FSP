

function calculate_slip_pressure(sig_fault::Float64, tau_fault::Float64, mu::Float64, p0::Float64, 
    biot::Float64=1.0, nu::Float64=0.5; dp::Float64=0.0, s11::Float64=0.0, s22::Float64=0.0, s33::Float64=0.0, s12::Float64=0.0, n1::Float64=0.0, n2::Float64=0.0)


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


function ComputeCriticalPorePressureForFailure(sigma_n::Union{Array{Real}, Real}, tau::Union{Array{Real}, Real}; mu = 0.6)
    #=
        Compute the critical pore pressure necessary for failure
    =#

    mu1 = abs(tau ./ sigma_n)
    mu2 = mu
    Pcritical = ((mu2 .- mu1) ./ mu2) .* abs(sigma_n)

    return Pcritical
end #ComputeCriticalPorePressureForFailure




#=
(sig_fault::Float64, tau_fault::Float64, mu::Float64, p0::Float64, 
    biot::Float64=1.0, nu::Float64=0.5, dp::Float64=0.0, s11::Float64=0.0, s22::Float64=0.0, s33::Float64=0.0, s12::Float64=0.0, n1::Float64=0.0, n2::Float64=0.0)
=#

jos_function = ComputeCriticalPorePressureForFailure(8035.13, 1013.37; mu = 0.58)

rall_function = calculate_slip_pressure(-1.5, 0.6, 0.6, 1.0, 1.0, 0.5)

println(jos_function)
println(rall_function)

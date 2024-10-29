module MohrsCircle

using JSON

include("get_hor_from_APhi.jl")  # used to calculate horizontal stress
using .GetHorFromAPhi 

export mohrs_3D

# prepares JSON input for the Mohr diagram and calculates the required parameters
function mohrs_3D(indatacell, hDV::Dict{Symbol, Any})

    # Extracting inputs from indatacell
    Sig0 = indatacell[1]  # Initial total stresses 
    p0 = indatacell[3]    # Initial Reference pore pressure
    strike = indatacell[4] # Strikes of faults in degrees
    dip = indatacell[5]    # Dips of faults in degrees
    SHdir = indatacell[6]  # Maximum horizontal stress direction
    dp = indatacell[7]     # Pressure perturbation at each fault
    mu = indatacell[8]     # Friction coefficients for faults
    biot = indatacell[9]   # Biot coefficient
    nu = indatacell[10]    # Poisson's ratio

    # Handle A-phi model if 11th element exists
    # In this case, we call get_hor_from_APhi.jl to adjust initial stresses (Sig0)
    if length(indatacell) > 10
        APhi = indatacell[11]
        SH, Sh = getHorFromAPhi(APhi, mu[1], Sig0[2], Sig0[1], p0, 1)
        Sig0[2] = Sh
        Sig0[3] = SH
    end

    # Handling NaN in pressure perturbation (dp)
    dp = replace(dp, NaN => 0.0)  # Set NaN values in dp to 0

    # Call calc_ppfail function
    failout, tau_fault, sig_fault = calc_ppfail(Sig0, SHdir, p0, biot, nu, mu[1], dp, strike, dip)

    # Initialize outputs
    outs = Dict()
    outs["ppfail"] = failout
    outs["cff"] = tau_fault .- mu .* sig_fault
    outs["scu"] = tau_fault ./ (mu .* sig_fault)

    # Sorting total stresses in descending order
    Sig0sorted = sort(Sig0, rev=true)
    N = length(strike)
    ixSv = findall(x -> x == Sig0[1], Sig0sorted)[1]  # Find index of vertical stress

    # Compute poroelastic total stresses
    Sig = repeat(Sig0sorted', N, 1)  # Replicate Sig0sorted for each fault
    Ds = repeat(biot * (1 - 2 * nu) / (1 - nu) .* dp, 1, 3)
    Ds[:, ixSv] .= 0.0  # Set vertical stress row to 0
    Sig .= Sig .+ Ds

    # Compute Mohr Circle radii
    a = LinRange(0, π, 100)
    c = exp.(im * a)


    R1 = 0.5 .* (Sig[:, 1] .- Sig[:, 3])
    R2 = 0.5 .* (Sig[:, 2] .- Sig[:, 3])
    R3 = 0.5 .* (Sig[:, 1] .- Sig[:, 2])
    # compute centers of circles
    C1 = [R1[k] * c .+ (Sig[k, 1] + Sig[k, 3]) / 2 .- (p0 + dp[k]) for k in 1:N]
    C2 = [R2[k] * c .+ (Sig[k, 2] + Sig[k, 3]) / 2 .- (p0 + dp[k]) for k in 1:N]
    C3 = [R3[k] * c .+ (Sig[k, 1] + Sig[k, 2]) / 2 .- (p0 + dp[k]) for k in 1:N]

    # Prepare a dictionary for JSON export
    json_data = Dict(
        "C1" => [Dict("center" => real(mean(C1[k])), "radius" => real(R1[k])) for k in 1:N],
        "C2" => [Dict("center" => real(mean(C2[k])), "radius" => real(R2[k])) for k in 1:N],
        "C3" => [Dict("center" => real(mean(C3[k])), "radius" => real(R3[k])) for k in 1:N],
        "sig_fault" => sig_fault, # Effective normal stress
        "tau_fault" => tau_fault, # Effective shear stress
        "failout" => failout, 
        "mu" => mu[1],
        "outs" => outs
    )

    # Save to JSON file
    open("mohr_circle_inputs.json", "w") do file
        JSON.print(file, json_data)
    end

    return (failout, outs, C1, C2, C3, sig_fault, tau_fault)
end

# Function to calculate delta pore pressure to slip (failure) and the resolved shear and normal stresses
# friction coefficient mu is assumed to be constant for all faults
function calc_ppfail(Sig0, az, p0, biot, nu, mu, dp, str, dip)
    # Cos and sin of azimuth
    cos_az = cos(az * π / 180)
    sin_az = sin(az * π / 180)

    # Cos and sin of strike and dip angles (vectorized)
    cs = cos.(str * π / 180)  # Element-wise cos
    ss = sin.(str * π / 180)  # Element-wise sin
    cd = cos.(dip * π / 180)  # Element-wise cos
    sd = sin.(dip * π / 180)  # Element-wise sin

    # Total stresses from input data
    Svert = Sig0[1]
    shmin = Sig0[2]
    sHmax = Sig0[3]

    # Factor based on Poisson's ratio for horizontal stress change
    f = biot * nu / (1 - nu)

    # Effective stresses at current pressure (broadcast subtraction and addition)
    s11 = shmin * cos_az^2 .+ sHmax * sin_az^2 .- p0 .- f .* dp
    s22 = shmin * sin_az^2 .+ sHmax * cos_az^2 .- p0 .- f .* dp
    s33 = Svert .- p0 .- dp  # Corrected to use broadcast subtraction
    s12 = (sHmax - shmin) * cos_az * sin_az

    # Components of unit normal vector to fault planes
    n1 = sd .* cs
    n2 = -sd .* ss
    n3 = cd  # n3 is not used below

    # Shear and normal stresses resolved on fault 
    tau_fault = sqrt.(
        (n2.^2 .* (s12.^2 - (-1 .+ n2.^2) .* (s22 .- s33).^2))
        - n1.^4 .* (s11 .- s33).^2
        + 4 * n1.^3 .* n2 .* s12 .* (-s11 .+ s33)
        + 2 * n1 .* n2 .* s12 .* (s11 .+ s22 .- 2 * n2.^2 .* s22 .+ 2 * (-1 .+ n2.^2) .* s33)
        + n1.^2 .* (s11.^2 .+ (1 .- 4 * n2.^2) .* s12.^2 .- 2 * s11 .* (n2.^2 .* (s22 .- s33) .+ s33) .+ s33 .* (2 * n2.^2 .* (s22 .- s33) .+ s33))
    )

    sig_fault = 2 * n1 .* n2 .* s12 + n1.^2 .* (s11 .- s33) + n2.^2 .* (s22 .- s33) .+ s33  

    # Mobilized friction coefficient with present stress and pressure
    mobmu = tau_fault ./ sig_fault
    mobmu[isnan.(mobmu)] .= 99.99
    mobmu[mobmu .< 0] .= 99.99

    # Coefficients of quadratic equation A*dp^2 + B*dp + C = 0 for fault slip
    C = -4 * (1 + mu^2) * n1.^3 .* n2 .* s12 .* (s11 .- s33)
    C .+= - (1 + mu^2) * n1.^4 .* (s11 .- s33).^2
    C .+= - (1 + mu^2) * n2.^4 .* (s22 .- s33).^2
    C .+= - mu^2 * s33.^2
    C .+= 2 * n1 .* n2 .* s12 .* (s11 .+ (1 - 2 * (1 + mu^2) * n2.^2) .* s22 .+ 2 * (1 + mu^2) * (-1 + n2.^2) .* s33)

    C .+= n2.^2 .* (s12.^2 .+ s22.^2 .- 2 * (1 + mu^2) .* s22 .* s33 .+ (1 + 2 * mu^2) .* s33.^2)
    C .+= n1.^2 .* (s11.^2 .+ (1 - 4 * (1 + mu^2) * n2.^2) .* s12.^2 .- 2 * (1 + mu^2) .* s11 .* (n2.^2 .* (s22 .- s33) .+ s33) .+ s33 .* (2 * (1 + mu^2) * n2.^2 .* (s22 .- s33) .+ s33 + 2 * mu^2 .* s33))

    B = 2 * (2 * (-1 + f) * (1 + mu^2) * n1.^3 .* n2 .* s12)
    B .+= 2 * n1 .* n2 .* (-(1 + mu^2) * (-1 + n2.^2) + f * (-1 + (1 + mu^2) * n2.^2)) .* s12
    B .+= (-1 + f) * (1 + mu^2) * n1.^4 .* (s11 .- s33)
    B .+= (-1 + f) * (1 + mu^2) * n2.^4 .* (s22 .- s33)
    B .+= mu^2 * s33
    B .+= n2.^2 .* ((1 - f + mu^2) .* s22 .+ (-1 + f - 2 * mu^2 + f * mu^2) .* s33)
    B .+= n1.^2 .* ((-(1 + mu^2) * (-1 + n2.^2) + f * (-1 + (1 + mu^2) * n2.^2)) .* s11 .+ (-1 + f) * (1 + mu^2) * n2.^2 .* (s22 .- 2 * s33) .+ (-1 + f - 2 * mu^2 + f * mu^2) .* s33)

    A = -mu^2 .* (1 .+ (-1 + f) .* n1.^2 .+ (-1 + f) .* n2.^2).^2
    A .+= -(-1 + f)^2 .* (n1.^4 .+ n2.^2 .* (-1 + n2.^2) .+ n1.^2 .* (-1 + 2 * n2.^2))

    Bsq_minus_4AC = B.^2 .- 4 .* A .* C

    # Calculate dp for fault slip (two solutions to quadratic equation)
    ppfail1 = (-B .- sqrt.(Bsq_minus_4AC)) ./ (2 * A)
    ppfail2 = (-B .+ sqrt.(Bsq_minus_4AC)) ./ (2 * A)

    # Handling cases where there is no solution
    no_solution = Bsq_minus_4AC .< 0
    ppfail_horiz_dist = sig_fault .- tau_fault / mu
    ppfail1[no_solution] .= -ppfail_horiz_dist[no_solution]
    ppfail2[no_solution] .= ppfail_horiz_dist[no_solution]

    # Initialize ppfail
    ppfail = fill(NaN, length(dip))

    # Loop over faults
    for k in 1:length(dip)
        if mobmu[k] < mu
            if ppfail1[k] > 0 && ppfail2[k] > 0
                ppfail[k] = min(ppfail1[k], ppfail2[k])
            elseif ppfail1[k] < 0 && ppfail2[k] < 0
                error("Pressure to slip calculation error")
            else
                ppfail[k] = ppfail1[k] > 0 ? ppfail1[k] : ppfail2[k]
            end
        elseif mobmu[k] > mu
            if ppfail1[k] < 0 && ppfail2[k] < 0
                ppfail[k] = min(ppfail1[k], ppfail2[k])
            elseif ppfail1[k] > 0 && ppfail2[k] > 0
                ppfail[k] = sig_fault[k] - tau_fault[k] / mu
            else
                ppfail[k] = ppfail1[k] > 0 ? ppfail2[k] : ppfail1[k]
            end
        else
            ppfail[k] = 0.0
        end
    end

    return (ppfail, tau_fault, sig_fault)
end




# Testing the module with hardcoded mock data
if abspath(PROGRAM_FILE) == @__FILE__
    # Mock input data
    indatacell = [
        [6000.0, 4500.0, 7000.0],  # Initial total stresses (Sig0)
        nothing,                   # Placeholder for the removed parameter
        2500.0,                    # Initial Reference pore pressure (p0)
        [45.0, 30.0],              # Strikes of faults in degrees
        [20.0, 25.0],              # Dips of faults in degrees
        60.0,                      # Maximum horizontal stress direction (SHdir)
        [100.0, 200.0],            # Pressure perturbation at each fault (dp)
        [0.6],                     # Friction coefficients for faults (mu)
        0.8,                       # Biot coefficient (biot)
        0.25                       # Poisson's ratio (nu)
    ]

    # Mock dictionary for hDV
    hDV = Dict{Symbol, Any}(
        :stress => Dict{Symbol, Any}(
            :aphi => Dict{Symbol, Any}(:use => 1)
        )
    )

    # Run the test with the mock data
    println("Running test with mock data...")
    result = mohrs_3D(indatacell, hDV)
    println("Test result:")
    println(result)
end

end # module MohrsCircle

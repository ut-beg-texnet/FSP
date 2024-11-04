module DeterministicGeomechanicsCalculations

include("get_hor_from_APhi.jl")

using JSON
using LinearAlgebra

using .GetHorFromAPhi

export DeterministicGeomechanicsResults, mohrs_3D, det_geomech_results_to_json

struct DeterministicGeomechanicsResults
    failout::Vector{Float64} # dpp to failure (only output used in monte carlo)
    outs::Dict # Results of Mohr Circle calulations (includes the ppfail)
    C1::Matrix{ComplexF64}
    C2::Matrix{ComplexF64}
    C3::Matrix{ComplexF64}
    sig_fault::Vector{Float64}
    tau_fault::Vector{Float64}
end

# Mohr's 3D Calculation Function
function mohrs_3D(geomechanics_inputs::Tuple)
    # Unpack the input array and extract variables
    (Sig0, _, pp0, strikes, dips, SHdir, dp, mu, biot, nu, aphi_value...) = geomechanics_inputs
    
    # Debug: Print initial inputs
    println("Initial inputs - Sig0: ", Sig0, ", pp0: ", pp0, ", strikes: ", strikes, ", dips: ", dips)
    
    # If user enters A-Phi value, calculate SH and Sh from A-Phi
    # run getHorFromAPhi 
    if !isempty(aphi_value)
        APhi = aphi_value[1]
        println("Using A-Phi value: ", APhi)
        SH, Sh = getHorFromAPhi(APhi, mu[1], Sig0[2], Sig0[1], pp0, true)
        Sig0[2], Sig0[3] = Sh, SH
    end

    # Handle NaNs in dp by setting them to zero
    dp[isnan.(dp)] .= 0.0  

    failout, tau_fault, sig_fault = calc_ppfail(Sig0, Float64(SHdir), pp0, biot, nu, mu[1], dp, strikes, Float64.(dips))

    # Debug: Print failure output
    println("Failure Output - failout: ", failout, ", tau_fault: ", tau_fault, ", sig_fault: ", sig_fault)

    # If no A-Phi is provided, compute additional outputs
    if isempty(aphi_value)
        outs = Dict("ppfail" => failout, "cff" => tau_fault .- mu .* sig_fault, "scu" => tau_fault ./ (mu .* sig_fault))
    else
        outs = Dict("ppfail" => failout)
    end

    # Calculate Mohr Circle parameters for plotting
    N = length(strikes)
    sorted_Sig0 = sort(Sig0, rev=true)
    ixSv = findfirst(x -> x == Sig0[1], sorted_Sig0)
    
    # Pass all necessary parameters, including `biot` and `nu`
    sig_fault_matrices, C1, C2, C3 = calculate_mohr_circle_parameters(sorted_Sig0, pp0, dp, ixSv, N, biot, nu)

    # Return structured results
    return DeterministicGeomechanicsResults(failout, outs, C1, C2, C3, sig_fault, tau_fault)
end

# Helper function for Mohr Circle parameters (used for plotting)
function calculate_mohr_circle_parameters(Sig0sorted, p0, dp, ixSv, N, biot, nu)
    a = range(0, stop=pi, length=100)
    c = exp.(1im .* a)

    # Debug: Print initial parameters for Mohr circle
    println("Mohr Circle - Sig0sorted: ", Sig0sorted, ", ixSv: ", ixSv, ", biot: ", biot, ", nu: ", nu)

    # Initialize `Sig` with each row matching the shape of `dp`
    Sig = hcat([fill(Sig0sorted[i], N) for i in 1:3]...)
    
    # Initialize `Ds` as a matrix matching the dimensions of `Sig`
    # We calculate `Ds` based on biot and nu, but only set the specified column (ixSv) to zero
    Ds = [biot * (1 - 2 * nu) / (1 - nu) * dp[i] for i in 1:N, _ in 1:3]
    Ds[:, ixSv] .= 0  # Set column `ixSv` to zero in all rows

    # Broadcasting addition for `Sig += Ds`
    Sig = Sig .+ Ds

    R1, R2, R3 = 0.5 .* (Sig[:,1] .- Sig[:,3]), 0.5 .* (Sig[:,2] .- Sig[:,3]), 0.5 .* (Sig[:,1] .- Sig[:,2])
    C1 = R1 .* c' .+ (Sig[:,1] .+ Sig[:,3]) / 2 .- (p0 .+ dp)
    C2 = R2 .* c' .+ (Sig[:,2] .+ Sig[:,3]) / 2 .- (p0 .+ dp)
    C3 = R3 .* c' .+ (Sig[:,1] .+ Sig[:,2]) / 2 .- (p0 .+ dp)

    # Debug: Print Mohr Circle matrices
    #println("Mohr Circle - C1: ", C1, ", C2: ", C2, ", C3: ", C3)

    return Sig, C1, C2, C3
end

# Function to calculate pore pressure to slip (failure) and shear and normal stresses on faults
# Friciton coefficient mu is the same for all faults and is a scalar
function calc_ppfail(Sig0::Vector{Float64}, az::Float64, p0::Float64, biot::Float64, nu::Float64, mu::Float64, dp::Vector{Float64}, strikes::Vector{Float64}, dips::Vector{Float64})
    # cos and sin of azimuth
    cos_az, sin_az = cosd(az), sind(az)
    # cos and sin of strike and dip angles 
    cs, ss, cd, sd = cosd.(strikes), sind.(strikes), cosd.(dips), sind.(dips)

    # total stresses from input data
    Svert, shmin, sHmax = Sig0[1], Sig0[2], Sig0[3]

    # factor based on nu for horizontal stress change
    # Biot is currently hardcoded to 1
    f = biot * nu / (1 - nu)


    # effective stresses at current pressure, compresive positive notation
    # dp is also an array, which makes s11, s22, s33, s12 also arrays
    # Note that factor f is used for dp because dp is pressure change due to
    # injection upto present time - horizontal stresses will have Poisson's
    # ratio effect (and Biot coeff, if other than 1; but see comments above)
    # However, p0 is initial pressure before injection, so factor f is not needed
    s11 = shmin * cos_az^2 .+ sHmax * sin_az^2 .- p0 .- f .* dp
    s22 = shmin * sin_az^2 .+ sHmax * cos_az^2 .- p0 .- f .* dp
    s33 = Svert .- p0 .- dp
    s12 = (sHmax - shmin) .* cos_az .* sin_az  
    
    println("Effective stresses at current pressure: s11: ", s11, ", s22: ", s22, ", s33: ", s33, ", s12: ", s12)
    

    # components of unit normal vector to fault planes
    n1 = sd .* cs
    n2 = -sd .* ss
    n3 = cd # n3 is not used below, since by def n3 = sqrt(1 - n1^2 - n2^2)

    # Shear and normal stresses resolved on fault
    # Note that tau_fault is absolute magnitude, but sig_fault is signed
    # value (which should be positive for most cases, but, with high pore
    # pressure, resolved normal stress on fault can become negative if any of 
    # the principal stresses is negative
    tau_fault = sqrt.((n2.^2 .* (s12.^2 .- (-1 .+ n2.^2) .* (s22 .- s33).^2) .-
              n1.^4 .* (s11 .- s33).^2  .+
              4 .* n1.^3 .* n2 .* s12 .* (-s11 .+ s33) .+ 
              2 .* n1 .* n2 .* s12 .* (s11 .+ s22 .- 2 .* n2.^2 .* s22 .+ 2 .* (-1 .+ n2.^2) .* s33) .+
              n1.^2 .* (s11.^2 .+ (1 .- 4 .* n2.^2) .* s12.^2 .- 2 .* s11 .* (n2.^2 .* (s22 .- s33) .+ s33) .+ s33 .* (2 .* n2.^2 .* (s22 .- s33) .+ s33))))

    sig_fault = 2 .* n1 .* n2 .* s12 .+ n1.^2 .* (s11 .- s33) .+ n2.^2 .* (s22 .- s33) .+ s33  # signed sig_fault  # signed sig_fault
    println("Resolved stresses on fault: tau_fault (shear stress): ", tau_fault, ", sig_fault (normal stress): ", sig_fault)


    # mobilized friction coefficient with present stress and pressure used below
    # mob_mu determines if pressure to slip should be positive or negative
    mob_mu = tau_fault ./ sig_fault
    mob_mu[isnan.(mob_mu)] .= 99.99 # set a high number if sig_fault is zero and mob_mu is NaN
    mob_mu[mob_mu .< 0] .= 99.99 # set a high number if sig_fault is zero and mob_mu is negative


    # Coefficients of quadratic equation A*dp^2 + B*dp + C = 0
    # such that the mobilized friction coefficient on fault = mu
    # i.e. dp that will bring the fault stresses to failure point
    C = -4 .* (1 .+ mu^2) .* n1.^3 .* n2 .* s12 .* (s11 .- s33) .- 
    (1 .+ mu^2) .* n1.^4 .* (s11 .- s33).^2 .-
    (1 .+ mu^2) .* n2.^4 .* (s22 .- s33).^2 .-
    mu^2 .* s33.^2 .+
    2 .* n1 .* n2 .* s12 .* (s11 .+ (1 .- 2 .* (1 .+ mu^2) .* n2.^2) .* s22 .+ 2 .* (1 .+ mu^2) .* (-1 .+ n2.^2) .* s33) .+
    n2.^2 .* (s12.^2 .+ s22.^2 .- 2 .* (1 .+ mu^2) .* s22 .* s33 .+ (1 .+ 2 .* mu^2) .* s33.^2) .+
    n1.^2 .* (s11.^2 .+ (1 .- 4 .* (1 .+ mu^2) .* n2.^2) .* s12.^2 .- 2 .* (1 .+ mu^2) .* s11 .* (n2.^2 .* (s22 .- s33) .+ s33) .+ s33 .* (2 .* (1 .+ mu^2) .* n2.^2 .* (s22 .- s33) .+ s33 .+ 2 .* mu^2 .* s33))


    B = 2 .* (2 .* (-1 .+ f) .* (1 .+ mu^2) .* n1.^3 .* n2 .* s12 .+
    2 .* n1 .* n2 .* (.-(1 .+ mu^2) .* (-1 .+ n2.^2) .+ f .* (-1 .+ (1 .+ mu^2) .* n2.^2)) .* s12 .+
    (-1 .+ f) .* (1 .+ mu^2) .* n1.^4 .* (s11 .- s33) .+
    (-1 .+ f) .* (1 .+ mu^2) .* n2.^4 .* (s22 .- s33) .+
    mu^2 .* s33 .+
    n2.^2 .* ((1 .- f .+ mu^2) .* s22 .+ (-1 .+ f .- 2 .* mu^2 .+ f .* mu^2) .* s33) .+
    n1.^2 .* (.-(1 .+ mu^2) .* (-1 .+ n2.^2) .+ f .* (-1 .+ (1 .+ mu^2) .* n2.^2)) .* s11 .+
    (-1 .+ f) .* (1 .+ mu^2) .* n2.^2 .* (s22 .- 2 .* s33) .+
    (-1 .+ f .- 2 .* mu^2 .+ f .* mu^2) .* s33)


    A = -mu^2 .* (1 .+ (-1 .+ f) .* n1.^2 .+ (-1 .+ f) .* n2.^2).^2 .- (-1 .+ f).^2 .* (n1.^4 .+ n2.^2 .* (-1 .+ n2.^2) .+ n1.^2 .* (-1 .+ 2 .* n2.^2))


    Bsq_minus_4AC = B .* B .- 4 .* A .* C

    ppfail1 = (-B .- sqrt.(Complex.(Bsq_minus_4AC))) ./ (2 .* A)
    ppfail2 = (-B .+ sqrt.(Complex.(Bsq_minus_4AC))) ./ (2 .* A)




    # For some cases, with low Poisson's ratio, there may be no solution
    # to above equation, indicated by Bsq_minus_4AC < 0, since its 
    # sqrt will be imaginary
    #
    # for such cases set ppfail to the horizontal distance from failure line 
    # (i.e. value for nu=0.5)
    # Assign this to ppfail1 and ppfail2 with opposite signs; mobmu
    # relative to mu below will determine correct sign (since there is no
    # solution this is just a placeholder anyway)
    #
    # Cannot set all these to same value such as zero because probability
    # distribution calculation (outside this function) has some issues
    
    no_solution = Bsq_minus_4AC .< 0
    ppfail_horiz_dist = sig_fault .- tau_fault ./ mu
    ppfail1[no_solution] .= -ppfail_horiz_dist[no_solution]
    ppfail2[no_solution] .= ppfail_horiz_dist[no_solution]

    # initialize ppfail to NaN
    ppfail = fill(NaN, length(dips))

    # loop over the faults 

    for k in 1:length(dips)
        println("Processing fault index ", k)
        println("  mobilized friction coefficient (mob_mu): ", mob_mu[k])
        println("  specified fault friction coefficient (mu): ", mu)
        println("  computed ppfail1 (real part): ", real(ppfail1[k]), ", ppfail2 (real part): ", real(ppfail2[k]))
    
        if mob_mu[k] < mu # if fault is below failure line
            # if both roots are positive, choose the smaller one
            if real(ppfail1[k]) > 0 && real(ppfail2[k]) > 0
                ppfail[k] = min(real(ppfail1[k]), real(ppfail2[k]))
                println("Choosing the smaller positive root as ppfail: ", ppfail[k])
            # if both roots are negative, there must be an error
            # we need at least one positive root (at some pressure, all Mohr circles will move enough to the left)
            elseif real(ppfail1[k]) < 0 && real(ppfail2[k]) < 0
                error("Pressure to slip calculation error. Both roots are negative.")
            else
                # if one root is negative, choose the positive one
                ppfail[k] = real(ppfail1[k]) > 0 ? real(ppfail1[k]) : real(ppfail2[k])
                println("Choosing the positive root among mixed signs as ppfail: ", ppfail[k])

            end
        elseif mob_mu[k] > mu # if fault is above failure line
            if ppfail1[k] < 0 && ppfail2[k] < 0 
                # both roots are negative, choose the one with the smaller 'magnitude'
                ppfail[k] = abs(real(ppfail1[k])) < abs(real(ppfail2[k])) ? real(ppfail2[k]) : real(ppfail1[k])
                println("Both roots are negative; choosing one with smaller magnitude as ppfail: ", ppfail[k])
            elseif ppfail1[k] > 0 && ppfail2[k] > 0
                # both roots are positive (even though fault is initially above failure line)
                # happens in a few cases with low nu due to nonlinear stress path
                # in this case, return horizontal distance
                ppfail[k] = sig_fault[k] - tau_fault[k] / mu
                println("Both roots positive but above failure line; setting ppfail to horizontal distance: ", ppfail[k])
            else
                # one should be negative and one positive, choose the negative one
                ppfail[k] = real(ppfail1[k]) > 0 ? real(ppfail2[k]) : real(ppfail1[k])
                println("Choosing the negative root among mixed signs as ppfail: ", ppfail[k])
            end
        
        else # mob_mu == mu (very unlikely due to finite numberical precision)
            ppfail[k] = 0
            println("Setting ppfail to 0 due to finite precision issue with mob_mu == mu.")
        end
    end
    return ppfail, tau_fault, sig_fault
end

# Function to export results to JSON
function det_geomech_results_to_json(results::DeterministicGeomechanicsResults, output_file::String)
    results_data = Dict(
        "outs" => results.outs,
        "C1" => results.C1,
        "C2" => results.C2,
        "C3" => results.C3,
        "sig_fault" => results.sig_fault,
        "tau_fault" => results.tau_fault
    )
    json_data = JSON.json(results_data)
    open(output_file, "w") do io
        write(io, json_data)
    end
    println("Deterministic Geomechanics Results exported to JSON: ", output_file)
end

end # End of module

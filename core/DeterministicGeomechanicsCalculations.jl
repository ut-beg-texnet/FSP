module DeterministicGeomechanicsCalculations

include("get_hor_from_APhi.jl")

using JSON
using LinearAlgebra
using OrderedCollections # to preserve order of dictionary keys in JSON output
using LoopVectorization
using Base.Threads
using BenchmarkTools

using .GetHorFromAPhi

export DeterministicGeomechanicsResults, mohrs_3D, det_geomech_results_to_json

struct DeterministicGeomechanicsResults
    failout::Vector{Float64} # deterministic pore pressure to failure (only output used in monte carlo)
    outs::Dict # Results of Mohr Circle calulations (includes the ppfail)
    C1::Tuple{Float64, Float64}  
    C2::Tuple{Float64, Float64}  
    C3::Tuple{Float64, Float64}  
    R1::Vector{Float64}
    R2::Vector{Float64}
    R3::Vector{Float64}
    center_C1::Float64
    center_C2::Float64
    center_C3::Float64
    sig_fault::Vector{Float64}
    tau_fault::Vector{Float64}
    frictional_slip_line::Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}} # contains start and end points of the frictional slip line
end

# Det Geomechanics main calculation function
function mohrs_3D(geomechanics_inputs::Tuple)
    (Sig0, _, pp0, strikes, dips, SHdir, dp, mu, biot, nu, rest...) = geomechanics_inputs
    # Sig0: Initial stress state (vertical, min horizontal, max horizontal)
    # pp0: Initial pore pressure
    # strikes: Fault strike angles
    # dips: Fault dip angles
    # SHdir: max horizontal stress direction
    # dp: pore pressure change
    # mu: friction coefficient
    # biot: Biot coefficient
    # nu: Poisson's ratio
    # aphi_value: A-Phi value (optional for APhi model use)

    aphi_value = nothing
    min_hor_stress = nothing

    # Process the rest of the arguments
    if length(rest) >= 1
        aphi_value = rest[1]
        if length(rest) >= 2
            min_hor_stress = rest[2]
        end
    end
    
    # keep for debuggin: Print initial inputs
    println("Initial inputs - Sig0: ", Sig0, ", pp0: ", pp0, ", strikes: ", strikes, ", dips: ", dips)
    
    # If user enters A-Phi value, calculate SH and Sh from A-Phi
    # this is an optional cli argument (--aphi) that expects a Float64
    # in this case, run getHorFromAPhi function to get SH and Sh
    if aphi_value !== nothing
        aphi_use = true
    else
        aphi_use = false
    end

    if aphi_value !== nothing
        println("Using A-Phi value: $aphi_value")
        if min_hor_stress !== nothing
            println("Using provided minimum horizontal stress: $min_hor_stress")
            # Use min_hor_stress in calculations
            SH, Sh = getHorFromAPhi(aphi_value, mu[1], min_hor_stress, Sig0[1], pp0, true, min_hor_stress)
        else
            # Calculate SH and Sh without min_hor_stress
            SH, Sh = getHorFromAPhi(aphi_value, mu[1], Sig0[2], Sig0[1], pp0, true)
        end
        Sig0[2], Sig0[3] = Sh, SH
    end

    # Handle NaNs in dp by setting them to zero
    dp[isnan.(dp)] .= 0.0  

    failout, tau_fault, sig_fault = calc_ppfail(Sig0, Float64(SHdir), pp0, biot, nu, mu[1], dp, strikes, Float64.(dips))

    # Debug print
    println("Failure Output - failout: ", failout, ", tau_fault: ", tau_fault, ", sig_fault: ", sig_fault)

    # If no A-Phi is provided, compute additional outputs
    
    if !aphi_use
        outs = Dict("ppfail" => failout, "cff" => tau_fault .- mu .* sig_fault, "scu" => tau_fault ./ (mu .* sig_fault))
    else
        outs = Dict("ppfail" => failout, "cff" => [], "scu" => [])
    end

    # Calculate Mohr Circle parameters for plotting
    N = length(strikes)
    Sig0sorted = sort(Sig0, rev=true)
    ixSv = findfirst(x -> x == Sig0[1], Sig0sorted) # index of vertical stress, use only one index if there's multiple
    
    # Pass all necessary parameters, including biot and nu
    sig_fault_matrices, C1, C2, C3 , R1, R2, R3 , center_C1, center_C2, center_C3, frictional_slip_line = calculate_mohr_circle_parameters(Sig0sorted, pp0, dp, ixSv, N, biot, nu, mu[1])

    # Return structured results
    return DeterministicGeomechanicsResults(failout, outs, C1, C2, C3, R1, R2, R3, center_C1, center_C2, center_C3, sig_fault, tau_fault, frictional_slip_line)
end

# Helper function for Mohr Circle parameters (used for plotting)
function calculate_mohr_circle_parameters(Sig0sorted, p0, dp, ixSv, N, biot, nu, mu)
    
    #= ORIGINAL CODE ---------------------------------------------------------------
    a = range(0, stop=pi, length=100) # range of 100 points from 0 to pi
    c = exp.(1im .* a) # complex exponential for each point in a
    # real part of complex number --> x-axis (normal stress)
    # imaginary part of complex number --> y-axis (shear stress)
    =# 
    #---------------------------------------------------------------------------


    # OPTIMIZED CODE (WITHOUT USING COMPLEX NUMBERS)----------
    a = range(0, stop=pi, length=100)
    cos_a = cos.(a)
    sin_a = sin.(a)
    


    # Debug: Print initial parameters for Mohr circle
    #println("Mohr Circle - Sig0sorted: ", Sig0sorted, ", ixSv: ", ixSv, ", biot: ", biot, ", nu: ", nu)

    # Initialize Sig with each row matching the shape of dp
    Sig = hcat([fill(Sig0sorted[i], N) for i in 1:3]...) # --> N x 3 matrix that repeats each element of Sig0sorted
    
    # Ds --> matrix matching the dimensions of Sig
    # calculate Ds based on biot and nu, but only set the specified column (ixSv) to zero
    #= ORIGINAL CODE ----------------------------------------------------------------------
    Ds = [biot * (1 - 2 * nu) / (1 - nu) * dp[i] for i in 1:N, _ in 1:3]
    Ds[:, ixSv] .= 0  # Set column ixSv to zero in all rows
    =#
    # -------------------------------------------------------------------------------------

    # OPTIMIZED CODE USING LOOP VECTORIZATION ----------------------
    Ds = Matrix{Float64}(undef, N, 3)
    @tturbo for i in 1:N, j in 1:3
        Ds[i, j] = biot * (1 - 2 * nu) / (1 - nu) * dp[i]
    end
    Ds[:, ixSv] .= 0  # Set column ixSv to zero in all rows
    # ---------------------------------------------------------------
    
    
    Sig = Sig .+ Ds

    #=
    # FRICTIONAL SLIP LINE PARAMS (make sure they are correct)
    # Generate the x values for the frictional slip line using the max stress value
    max_sig = maximum(Sig0sorted)  # Max value for normal stress in the frictional slip line
    σ_values = range(0, stop=max_sig, length=100)  # Generate 100 points from 0 to max(σ)

    # Calculate y-values (τ) for the frictional slip line
    τ_values = [μ * σ for σ in σ_values]  # τ = μ * σ
    =#

    # ALTERNATIVE METHOD FOR FRICTIONAL SLIP LINE (FASTER, WITHOUT CALCULATING INDIVIDUAL POINTS)
    # Calculate maximum normal stress (σ_max) for the frictional slip line
    σ_max = maximum(Sig0sorted)

    # Calculate the second point on the frictional slip line
    τ_max = mu * σ_max  # This is the maximum τ corresponding to σ_max

    # The frictional slip line can be represented by the two points (0, 0) and (σ_max, τ_max)
    frictional_slip_line = ((0.0, 0.0), (σ_max, τ_max))

    # radii of circles
    #= ORIGINAL CODE ---------------------------------------------------------------
    R1, R2, R3 = 0.5 .* (Sig[:,1] .- Sig[:,3]), 0.5 .* (Sig[:,2] .- Sig[:,3]), 0.5 .* (Sig[:,1] .- Sig[:,2])
    #---------------------------------------------------------------------------
    =#

    # OPTIMIZED CODE -------------------------------
    R1 = @turbo 0.5 .* (Sig[:,1] .- Sig[:,3])
    R2 = @turbo 0.5 .* (Sig[:,2] .- Sig[:,3])
    R3 = @turbo 0.5 .* (Sig[:,1] .- Sig[:,2])
    # ---------------------------------------------


    center_C1 = (Sig0sorted[1] + Sig0sorted[3]) / 2  # largest circle
    center_C2 = (Sig0sorted[2] + Sig0sorted[3]) / 2  # middle 
    center_C3 = (Sig0sorted[1] + Sig0sorted[2]) / 2  # small
    
    #= ORIGINAL CODE-------------------------------------------------------------------------
    # calculate points on each Mohr circle (do we really need that?) --> play around with 'a' range (line 85)
    #(replaced loop with broadcasting)
    C1 = R1 .* c' .+ (Sig[:,1] .+ Sig[:,3]) / 2 .- (p0 .+ dp)

    C2 = R2 .* c' .+ (Sig[:,2] .+ Sig[:,3]) / 2 .- (p0 .+ dp)
    C3 = R3 .* c' .+ (Sig[:,1] .+ Sig[:,2]) / 2 .- (p0 .+ dp)
    =#
    # ---------------------------------------------------------------------------------------

    # OPTIMIZED CODE (WITHOUT USING COMPLEX NUMBERS)----------
    C1_x = @turbo R1 .* cos_a' .+ (Sig[:,1] .+ Sig[:,3]) / 2 .- (p0 .+ dp)
    C2_x = @turbo R2 .* cos_a' .+ (Sig[:,2] .+ Sig[:,3]) / 2 .- (p0 .+ dp)
    C3_x = @turbo R3 .* cos_a' .+ (Sig[:,1] .+ Sig[:,2]) / 2 .- (p0 .+ dp)
    
    #=
    C1_min = minimum(real(C1))
    C1_max = maximum(real(C1))
    C2_min = minimum(real(C2))
    C2_max = maximum(real(C2))
    C3_min = minimum(real(C3))
    C3_max = maximum(real(C3))
    =#

    C1_min, C1_max = minimum(C1_x), maximum(C1_x)
    C2_min, C2_max = minimum(C2_x), maximum(C2_x)
    C3_min, C3_max = minimum(C3_x), maximum(C3_x)
    
    # Debug: Mohr Circle matrices
    #println("Mohr Circle - C1: ", C1, ", C2: ", C2, ", C3: ", C3)

    return Sig, (C1_min, C1_max), (C2_min, C2_max), (C3_min, C3_max), R1, R2, R3, center_C1, center_C2, center_C3, frictional_slip_line
end



# Function to calculate pore pressure to slip (failure) and shear and normal stresses on faults
# Friciton coefficient mu is the same for all faults and is a scalar
function calc_ppfail(Sig0::Vector{Float64}, az::Float64, p0::Float64, biot::Float64, nu::Float64, mu::Float64, dp::Vector{Float64}, strikes::Vector{Float64}, dips::Vector{Float64})
    
    println("CALC_PPFAIL_INPUTS: Sig0: ", Sig0, ", az: ", az, ", p0: ", p0, ", biot: ", biot, ", nu: ", nu, ", mu: ", mu, ", dp: ", dp, ", strikes: ", strikes, ", dips: ", dips)
    
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
    s33 = Svert .- p0 .- dp # this is σV in the Mohr 
    s12 = (sHmax - shmin) .* cos_az .* sin_az  
    
    #println("Effective stresses at current pressure: s11: ", s11, ", s22: ", s22, ", s33: ", s33, ", s12: ", s12)
    

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
    #println("Effective stresses projected onto fault: tau_fault (shear stress): ", tau_fault, ", sig_fault (normal stress): ", sig_fault)


    # mobilized friction coefficient with present stress and pressure used below
    # mob_mu determines if pressure to slip should be positive or negative
    mob_mu = tau_fault ./ sig_fault
    mob_mu[isnan.(mob_mu)] .= 99.99 # set a high number if sig_fault is zero and mob_mu is NaN
    mob_mu[mob_mu .< 0] .= 99.99 # set a high number if sig_fault is zero and mob_mu is negative


    # Coefficients of quadratic equation A*dp^2 + B*dp + C = 0
    # such that the mobilized friction coefficient on fault = mu
    # here we find the dp that will bring the fault stresses to failure point
    # PRE COMPUTE REPEATED TERMS (TO DO)
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
    # to above equation, indicated by Bsq_minus_4AC < 0, since its sqrt will be imaginary
    # For such cases, we set ppfail to the horizontal distance from failure line 
    
    # Assign this to ppfail1 and ppfail2 with opposite signs; mobmu relative to mu below will determine correct sign (since there is no
    # solution this is just a placeholder anyway)

    # Cannot set all these to same value such as zero because probability
    # distribution calculation (outside this function) has some issues
    
    no_solution = Bsq_minus_4AC .< 0
    ppfail_horiz_dist = sig_fault .- tau_fault ./ mu
    ppfail1[no_solution] .= -ppfail_horiz_dist[no_solution]
    ppfail2[no_solution] .= ppfail_horiz_dist[no_solution]

    # initialize ppfail (NaN)
    ppfail = fill(NaN, length(dips))

    # iterate through the faults (length of dips indicates how many faults we have)
    @threads for k in 1:length(dips)
        #println("Processing fault number ", k)
        #println("  mobilized friction coefficient (mob_mu): ", mob_mu[k])
        #println("  specified fault friction coefficient (mu): ", mu)
        #println("  computed ppfail1 (real part): ", real(ppfail1[k]), ", ppfail2 (real part): ", real(ppfail2[k]))
    
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

# Function to export calculation results to JSON

function det_geomech_results_to_json(results::DeterministicGeomechanicsResults, output_file::String)
    # Combine all outputs into one structured dictionary
    results_data = OrderedDict(
        "geomechanics_results" => OrderedDict(
            "ppfail" => results.outs["ppfail"],
            "cff" => results.outs["cff"],
            "scu" => results.outs["scu"],
            "tau_fault" => results.tau_fault,
            "sig_fault" => results.sig_fault
        ),
        "mohr_circle_data" => OrderedDict(
            "frictional_slip_line" => Dict(
                "start" => results.frictional_slip_line[1],
                "end" => results.frictional_slip_line[2]
            ),
            "R1" => results.R1,
            "R2" => results.R2,
            "R3" => results.R3,
            "center_C1" => results.center_C1,
            "center_C2" => results.center_C2,
            "center_C3" => results.center_C3,
            "C1" => Dict("min" => results.C1[1], "max" => results.C1[2]),
            "C2" => Dict("min" => results.C2[1], "max" => results.C2[2]),
            "C3" => Dict("min" => results.C3[1], "max" => results.C3[2])
        )
    )
    
    # Convert dictionary to JSON format and write to file
    json_data = JSON.json(results_data)
    open(output_file, "w") do io
        write(io, json_data)
    end
    println("Geomechanical results exported to JSON: ", output_file)
end



end # End of module

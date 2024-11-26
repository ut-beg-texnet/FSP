module DeterministicGeomechanicsCalculations

include("get_hor_from_APhi.jl")

using JSON
using LinearAlgebra
using OrderedCollections # to preserve order of dictionary keys in JSON output
using LoopVectorization
using Base.Threads
using BenchmarkTools

using .GetHorFromAPhi

export DeterministicGeomechanicsResults, mohrs_3D, det_geomech_results_to_json, calc_ppfail, calculate_mohr_diagram_params, append_to_json

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
    # Unpack input tuple
    (Sig0, pp0, strikes, dips, SHdir, dp, mu, biot, nu, aphi_value, min_hor_stress, stress_model) = geomechanics_inputs
    # Sig0: Initial stress state (vertical, min horizontal, max horizontal)
    # pp0: Initial pore pressure
    # strikes: Fault strike angles
    # dips: Fault dip angles
    # SHdir: max horizontal stress direction
    # dp: pore pressure change
    # mu: friction coefficient
    # biot: Biot coefficient
    # nu: Poisson's ratio
    # rest: additional inputs for optional A-Phi model parameters
    println("aphi value: ", aphi_value)

    # DEBUGGING FOR TYPE MISMATCH WHEN CALLED THROUGH PROB GEOMECHANICS STEP
    println("mohrs_3D inputs:")
    println("  Sig0 (Initial stress state): ", Sig0, " (Type: ", typeof(Sig0), ")")
    println("  pp0 (Initial pore pressure): ", pp0, " (Type: ", typeof(pp0), ")")
    println("  strikes (Fault strike angles): ", strikes, " (Type: ", typeof(strikes), ")")
    println("  dips (Fault dip angles): ", dips, " (Type: ", typeof(dips), ")")
    println("  SHdir (Max horizontal stress direction): ", SHdir, " (Type: ", typeof(SHdir), ")")
    println("  dp (Pore pressure change): ", dp, " (Type: ", typeof(dp), ")")
    println("  mu (Friction coefficient): ", mu, " (Type: ", typeof(mu), ")")
    println("  biot (Biot coefficient): ", biot, " (Type: ", typeof(biot), ")")
    println("  nu (Poisson's ratio): ", nu, " (Type: ", typeof(nu), ")")
    #println("  rest (Additional A-Phi parameters): ", rest, " (Type: ", typeof(rest), ")")

    # Initialize optional variables for A-Phi model
    #aphi_value = nothing
    #min_hor_stress = nothing
    #stress_model = nothing  # will be assigned later based on inputs

    # Ensure `dp` and `pp0` are wrapped in arrays for consistent indexing if they are scalars
    dp = isa(dp, Float64) ? [dp] : dp
    pp0 = isa(pp0, Float64) ? [pp0] : pp0

    #println("  dp: ", dp, " (Type: ", typeof(dp), ")")
    #println("  pp0: ", pp0, " (Type: ", typeof(pp0), ")")

    
    # Process the rest of the arguments
    # If provided, `aphi_value` and `min_hor_stress` are extracted from `rest`
    #=
    if length(rest) >= 1
        println("inside lenght!!!")
        aphi_value = rest[1]
        if length(rest) >= 2
            min_hor_stress = rest[2]
        end
    end
    =#

    #=
    # Ensure `aphi_value` and `min_hor_stress` are scalars, not tuples, if they are set
    if aphi_value !== nothing && isa(aphi_value, Tuple)
        aphi_value = aphi_value[1]
    end
    if min_hor_stress !== nothing && isa(min_hor_stress, Tuple)
        min_hor_stress = min_hor_stress[1]
    end

    println("Optional parameters:")
    println("  aphi_value: ", aphi_value, " (Type: ", typeof(aphi_value), ")")
    println("  min_hor_stress: ", min_hor_stress, " (Type: ", typeof(min_hor_stress), ")")
    

    # Determine the stress model based on the presence of A-Phi parameters
    if aphi_value !== nothing
        aphi_use = true
    else
        aphi_use = false
        stress_model = "gradient"
    end
    =#
    # If A-Phi model is used, compute `SH` and `Sh` based on `aphi_value` and optional `min_hor_stress`
    println("Determining get_hor_from_APhi case...")
    if aphi_value !== nothing
        if min_hor_stress !== nothing
            # Use min_hor_stress in calculations
            SH, Sh = getHorFromAPhi(aphi_value, mu[1], Sig0[1], pp0[1], min_hor_stress)
            Sig0[2], Sig0[3] = Sh, SH
        else
            SH, Sh = getHorFromAPhi(aphi_value, mu[1], Sig0[1], pp0[1])
            Sig0[2], Sig0[3] = Sh, SH
        end
    end

    # Handle NaNs in dp by setting them to zero
    dp[isnan.(dp)] .= 0.0  

    # Calculate pore pressure to slip, and shear and normal stresses on faults using `calc_ppfail`
    # `pp0[1]` is used for scalar expectations in `calc_ppfail`
    failout, tau_fault, sig_fault = calc_ppfail(Sig0, Float64(SHdir), pp0[1], biot, nu, mu[1], dp, strikes, Float64.(dips))

    # Debug print for failure output results
    println("Failure Output - failout: ", failout, ", tau_fault: ", tau_fault, ", sig_fault: ", sig_fault)

    # Compile results in dictionary `outs`
    outs = Dict(
        "stress_model" => stress_model,
        "ppfail" => failout,
        "cff" => tau_fault .- mu .* sig_fault,
        "scu" => tau_fault ./ (mu .* sig_fault)
    )

    # Calculate Mohr Circle parameters for plotting
    N = length(strikes)
    Sig0sorted = sort(Sig0, rev=true)
    ixSv = findfirst(x -> x == Sig0[1], Sig0sorted)  # index of vertical stress

    # Call `calculate_mohr_circle_parameters` to get Mohr Circle parameters
    sig_fault_matrices, C1, C2, C3, R1, R2, R3, center_C1, center_C2, center_C3, frictional_slip_line = calculate_mohr_circle_parameters(Sig0sorted, pp0[1], dp, ixSv, N, biot, nu, mu[1])

    # Return results encapsulated in `DeterministicGeomechanicsResults`
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


# OLD CALC_PPFAIL -----------------------------------------------------------------
# Function to calculate pore pressure to slip (failure) and shear and normal stresses on faults
# Friciton coefficient mu is the same for all faults and is a scalar
function calc_ppfail(Sig0::Vector{Float64}, az::Float64, p0::Float64, biot::Float64, nu::Float64, mu::Float64, dp::Vector{Float64}, strikes::Vector{Float64}, dips::Vector{Float64})
    
    #println("CALC_PPFAIL_INPUTS: Sig0: ", Sig0, ", az: ", az, ", p0: ", p0, ", biot: ", biot, ", nu: ", nu, ", mu: ", mu, ", dp: ", dp, ", strikes: ", strikes, ", dips: ", dips)
    #=
    # cos and sin of azimuth
    cos_az, sin_az = cosd(az), sind(az)
    # cos and sin of strike and dip angles 
    cs, ss, cd, sd = cosd.(strikes), sind.(strikes), cosd.(dips), sind.(dips)
    =#

    @assert !isempty(strikes) "Calc_ppfail: Strike angles cannot be empty."
    @assert !isempty(dips) "Calc_ppfail: Dip angles cannot be empty."

    
    # Cosine and sine of azimuth angle (az)
    cos_az = cos(az * π / 180)
    sin_az = sin(az * π / 180)

    cs = cos.(strikes .* π ./ 180)  # cosine of strike angles
    ss = sin.(strikes .* π ./ 180)  # sine of strike angles
    cd = cos.(dips .* π ./ 180)  # cosine of dip angles
    sd = sin.(dips .* π ./ 180)  # sine of dip angles
    
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
                #println("Choosing the smaller positive root as ppfail: ", ppfail[k])
            # if both roots are negative, there must be an error
            # we need at least one positive root (at some pressure, all Mohr circles will move enough to the left)
            elseif real(ppfail1[k]) < 0 && real(ppfail2[k]) < 0
                
                println("Pressure to slip calculation error. Both roots are negative.")
                #error("Pressure to slip calculation error. Both roots are negative.")
            else
                # if one root is negative, choose the positive one
                ppfail[k] = real(ppfail1[k]) > 0 ? real(ppfail1[k]) : real(ppfail2[k])
                #println("Choosing the positive root among mixed signs as ppfail: ", ppfail[k])

            end
        elseif mob_mu[k] > mu # if fault is above failure line
            if real(ppfail1[k]) < 0 && real(ppfail2[k]) < 0 
                # both roots are negative, choose the one with the smaller 'magnitude'
                ppfail[k] = abs(real(ppfail1[k])) < abs(real(ppfail2[k])) ? real(ppfail2[k]) : real(ppfail1[k])
                #println("Both roots are negative; choosing one with smaller magnitude as ppfail: ", ppfail[k])
            elseif real(ppfail1[k]) > 0 && real(ppfail2[k]) > 0
                # both roots are positive (even though fault is initially above failure line)
                # happens in a few cases with low nu due to nonlinear stress path
                # in this case, return horizontal distance
                ppfail[k] = sig_fault[k] - tau_fault[k] / mu
                #println("Both roots positive but above failure line; setting ppfail to horizontal distance: ", ppfail[k])
            else
                # one should be negative and one positive, choose the negative one
                ppfail[k] = real(ppfail1[k]) > 0 ? real(ppfail2[k]) : real(ppfail1[k])
                #println("Choosing the negative root among mixed signs as ppfail: ", ppfail[k])
            end
        
        else # mob_mu == mu (very unlikely due to finite numberical precision)
            ppfail[k] = 0
            #println("Setting ppfail to 0 due to finite precision issue with mob_mu == mu.")
        end
    end
    ppfail = abs.(ppfail)  # ensure ppfail is positive
    return ppfail, tau_fault, sig_fault
end


# function to append the results of this step with another JSON file passed as argument
function append_to_json(prev_step_file::String, results::DeterministicGeomechanicsResults, output_file::String)
    # Read the existing JSON file
    json_data = JSON.parsefile(prev_step_file)

    # check for 'step2_output' key
    if !haskey(json_data, "step2_output")
        json_data["step2_output"] = Dict()
    end

    # add geomechanics results under the 'step2_output' key
    json_data["step2_output"]["geomechanics_results"] = Dict(
        "ppfail" => results.outs["ppfail"],
        "cff" => results.outs["cff"],
        "scu" => results.outs["scu"]
    )

    # same with Mohr Circle data
    json_data["step2_output"]["mohr_circle_data"] = Dict(
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
    
    
    # Write the updated JSON data to the output file
    open(output_file, "w") do io
        write(io, JSON.json(json_data))
    end
    println("Results appended to JSON file: ", output_file)
end


# Function to export calculation results to JSON
function det_geomech_results_to_json(results::DeterministicGeomechanicsResults, output_file::String)
    # Combine all outputs into one structured dictionary
    results_data = OrderedDict(
        "geomechanics_results" => OrderedDict(
            "ppfail" => round.(results.outs["ppfail"], digits=2),
            "cff" => round.(results.outs["cff"], digits=2),
            "scu" => round.(results.outs["scu"], digits=2),
            "tau_fault" => results.tau_fault,
            "sig_fault" => results.sig_fault
        ),
        "mohr_circle_data" => OrderedDict(
            "frictional_slip_line" => Dict(
                "start" => (round(results.frictional_slip_line[1][1]), round(results.frictional_slip_line[1][2])),
                "end" => (round(results.frictional_slip_line[2][1]), round(results.frictional_slip_line[2][2]))
            ),
            "R1" => results.R1,
            "R2" => results.R2,
            "R3" => results.R3,
            "center_C1" => round(results.center_C1),
            "center_C2" => round(results.center_C2),
            "center_C3" => round(results.center_C3),
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





























#=
function calc_ppfail(principal_stresses, SHdir, p0, biot, nu, friction, dp, strike, dip)
    # Convert angles to radians
    SHdir_rad = deg2rad(SHdir)
    
    # Cos and sin of azimuth
    cos_az = cos(SHdir_rad)
    sin_az = sin(SHdir_rad)
    
    # Extract principal stresses in correct order
    sHmax = principal_stresses[3]  # Maximum horizontal stress
    shmin = principal_stresses[2]  # Minimum horizontal stress
    Svert = principal_stresses[1]  # Vertical stress
    
    # Poisson's ratio effect factor (matches MATLAB exactly)
    f = convert(Float64, biot * nu / (1 - nu))
    
    # Effective stresses at current pressure (using total stresses)
    s11 = shmin*cos_az^2 + sHmax*sin_az^2 - p0 - f*dp
    s22 = shmin*sin_az^2 + sHmax*cos_az^2 - p0 - f*dp
    s33 = Svert - p0 - dp
    s12 = (sHmax - shmin)*cos_az*sin_az
    
    # Convert fault orientation to radians
    str_rad = deg2rad(strike)
    dip_rad = deg2rad(dip)
    
    # Components of unit normal vector to fault plane
    n1 = sin(dip_rad) * cos(str_rad)
    n2 = -sin(dip_rad) * sin(str_rad)
    n3 = cos(dip_rad)
    
    # Calculate shear and normal stresses on fault (matches MATLAB)
    tau_fault = sqrt(
        n2^2 * (s12^2 - (-1 + n2^2)*(s22-s33)^2) - 
        n1^4 * (s11-s33)^2 + 
        4*n1^3 * n2 * s12*(-s11+s33) + 
        2*n1 * n2 * s12*(s11 + s22 - 2*n2^2 * s22 + 2*(-1+n2^2) * s33) + 
        n1^2 * (s11^2 + (1-4*n2^2)*s12^2 - 2*s11*(n2^2 * (s22-s33) + s33) + 
                s33*(2*n2^2 * (s22-s33) + s33))
    )
    
    sig_fault = 2*n1*n2*s12 + n1^2*(s11-s33) + n2^2*(s22-s33) + s33
    
    # Mobilized friction coefficient (matches MATLAB)
    mobmu = tau_fault/sig_fault
    if isnan(mobmu) || mobmu < 0
        mobmu = 99.99  # Same as MATLAB for NaN or negative cases
    end
    
    # Coefficients of quadratic equation A*dp^2 + B*dp + C = 0
    # Ensure all terms are Float64 for consistent precision with MATLAB
    C = convert(Float64, -4*(1+friction^2)*n1^3 * n2*s12*(s11-s33) -
        (1+friction^2)*n1^4 * (s11-s33)^2 -
        (1+friction^2)*n2^4 * (s22-s33)^2 -
        friction^2 * s33^2 +
        2*n1 * n2*s12 * (s11 + (1-2*(1+friction^2)*n2^2)*s22 + 2*(1+friction^2)*(-1+n2^2)*s33) +
        n2^2 * (s12^2 + s22^2 - 2*(1+friction^2)*s22*s33 + (1+2*friction^2)*s33^2) +
        n1^2 * (s11^2 + (1-4*(1+friction^2)*n2^2)*s12^2 - 2*(1+friction^2)*s11*(n2^2 *(s22-s33) +s33) + 
                s33*(2*(1+friction^2)*n2^2 *(s22-s33) +s33 +2*friction^2 *s33)))
    
    B = convert(Float64, 2*( 2*(-1+f)*(1+friction^2)*n1^3 * n2*s12 +
            2*n1*n2 * (-(1+friction^2)*(-1+n2^2) + f*(-1+(1+friction^2)*n2^2))*s12 +
            (-1+f)*(1+friction^2)*n1^4 *(s11-s33) +
            (-1+f)*(1+friction^2)*n2^4 *(s22-s33) +
            friction^2 *s33 +
            n2^2 *((1-f+friction^2)*s22 + (-1 +f -2*friction^2 +f*friction^2)*s33) +
            n1^2 * ((-(1+friction^2)*(-1+n2^2) + f*(-1+(1+friction^2)*n2^2))*s11 + 
                    (-1 +f)*(1+friction^2)*n2^2 *(s22-2*s33) +(-1+f-2*friction^2 +f*friction^2)*s33) ))
    
    A = convert(Float64, -friction^2 *(1+(-1+f)*n1^2 +(-1+f)*n2^2)^2 -
        (-1+f)^2 *(n1^4 + n2^2 *(-1+n2^2) + n1^2 *(-1+2*n2^2)))
    
    # Calculate discriminant with high precision
    Bsq_minus_4AC = B^2 - 4*A*C
    
    # Handle no-solution case (matches MATLAB)
    if Bsq_minus_4AC < 0
        ppfail = sig_fault - tau_fault/friction
        return ppfail, tau_fault, sig_fault
    end
    
    # Calculate both solutions with high precision
    ppfail1 = (-B - sqrt(Bsq_minus_4AC))/(2*A)
    ppfail2 = (-B + sqrt(Bsq_minus_4AC))/(2*A)
    
    # Choose appropriate solution based on mobilized friction (matches MATLAB)
    if mobmu < friction
        # Below failure line
        if ppfail1 > 0 && ppfail2 > 0
            ppfail = min(ppfail1, ppfail2)
        elseif ppfail1 < 0 && ppfail2 < 0
            @warn "Pressure to slip calculation error: Both roots negative"
            ppfail = sig_fault - tau_fault/friction
        else
            ppfail = max(ppfail1, ppfail2)  # Choose positive root
        end
    else
        # Above failure line
        if ppfail1 < 0 && ppfail2 < 0
            ppfail = abs(ppfail1) < abs(ppfail2) ? ppfail1 : ppfail2
        elseif ppfail1 > 0 && ppfail2 > 0
            ppfail = sig_fault - tau_fault/friction
        else
            ppfail = min(ppfail1, ppfail2)  # Choose negative root
        end
    end
    
    return ppfail, tau_fault, sig_fault
end
=#

#=
function calculate_mohr_diagram_params(vertical_stress, min_horizontal_stress, max_horizontal_stress, p0)
    # Calculate effective stresses
    σv_eff = vertical_stress - p0
    σh_eff = min_horizontal_stress - p0
    σH_eff = max_horizontal_stress - p0
    
    # For APhi calculations, we need to be careful about the order
    # σh is always the minimum stress
    # σH is always the maximum stress
    # σv is in between for normal faulting (n=0)
    # σv is maximum for reverse faulting (n=2)
    # σv varies for strike-slip (n=1)
    
    # Calculate the three circles:
    # 1. Large circle: from σh to σH
    # 2. Middle circle: from σh to σv
    # 3. Small circle: from σv to σH
    
    # Circle 1 (largest): σh to σH
    R1 = abs(σH_eff - σh_eff) / 2
    C1 = (σH_eff + σh_eff) / 2
    
    # Circle 2 (middle): σh to σv
    R2 = abs(σv_eff - σh_eff) / 2
    C2 = (σv_eff + σh_eff) / 2
    
    # Circle 3 (smallest): σv to σH
    R3 = abs(σH_eff - σv_eff) / 2
    C3 = (σH_eff + σv_eff) / 2
    
    # Create dictionary for JSON output
    mohr_diagram_params = Dict(
        "effective_stresses" => Dict(
            "σh" => σh_eff,
            "σv" => σv_eff,
            "σH" => σH_eff
        ),
        "circles" => [
            Dict("center" => C1, "radius" => R1, "label" => "σH-σh"),  # Largest circle
            Dict("center" => C2, "radius" => R2, "label" => "σv-σh"),  # Middle circle
            Dict("center" => C3, "radius" => R3, "label" => "σH-σv")   # Smallest circle
        ],
        "max_stress" => maximum([σh_eff, σv_eff, σH_eff]),
        "max_shear" => maximum([R1, R2, R3])
    )
    
    # Print debug information
    println("\n=== Mohr Circle Parameters ===")
    println("Effective Stresses:")
    println("- σh = $(round(σh_eff, digits=2)) psi")
    println("- σv = $(round(σv_eff, digits=2)) psi")
    println("- σH = $(round(σH_eff, digits=2)) psi")
    println("\nCircles:")
    println("1. σH-σh: center = $(round(C1, digits=2)), radius = $(round(R1, digits=2))")
    println("2. σv-σh: center = $(round(C2, digits=2)), radius = $(round(R2, digits=2))")
    println("3. σH-σv: center = $(round(C3, digits=2)), radius = $(round(R3, digits=2))")
    
    return mohr_diagram_params
end

=#

#=

function calc_mohr_circles(stress_data, fault_data)
    # Extract principal stresses in correct order
    σ1 = stress_data["vertical_stress"]  # Svert
    σ2 = stress_data["min_horizontal_stress"]  # shmin 
    σ3 = stress_data["max_horizontal_stress"]  # sHmax
    
    # Calculate circle parameters
    circles = []
    
    # Calculate parameters for each circle
    for i in 1:2
        center = (σ1 + σ2) / 2
        radius = abs(σ1 - σ2) / 2
        
        push!(circles, Dict(
            "center" => center,
            "radius" => radius,
            "sigma1" => σ1,
            "sigma2" => σ2
        ))
        
        # Update stresses for second circle
        if i == 1
            σ1 = σ2
            σ2 = σ3
        end
    end
    
    return circles
end

=#





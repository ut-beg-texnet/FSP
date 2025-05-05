module JuliaFSPGraphs

export plot_pressure_distance_graph, plot_pressure_grid_heatmap, plot_mohr_diagram_geo, plot_mohr_diagram_hydro, plot_injection_rate_line_chart, plot_fault_surface_map,
plot_cdf_det_hydro, plot_prob_hydro_combined_cdf, fault_surface_map_data_to_d3, mohr_diagram_data_to_d3_portal, injection_rate_data_to_d3, prob_geomechanics_cdf, fault_sensitivity_tornado_chart_to_d3,
uncertainty_variability_inputs_to_d3, prob_hydrology_cdf, input_distribution_histograms_to_d3, mohr_diagram_hydro_data_to_d3_portal, hydro_input_distribution_histograms_to_d3

#using Plots
#using JSON
#using Colors
using LinearAlgebra
using SpecialFunctions
using Statistics
using Dates
using FileIO
using DataFrames
using CSV
using Printf
using Random
using Distributions
using StatsBase  # Added for ecdf function
using PrettyTables
include("../TexNetWebToolLauncherHelperJulia.jl")
include("../core/geomechanics_model.jl")

using .TexNetWebToolLauncherHelperJulia
using .GeomechanicsModel



# Hydrology: Plot pressure vs distance
# For the well(s) of interest, plot the pressure vs distance graph
# Using pressureScenario_constant_rate function for pressure calculation
# This graph does not superimpose them
function plot_pressure_distance_graph(
    radial_data::Vector{Tuple{String,Vector{Float64},Vector{Float64}}}
)
    if isempty(radial_data)
        println("No radial data to plot. Skipping pressure-distance graph.")
        return
    end

    # get maximum pressure for all curves
    max_pressure = maximum([maximum(pressure_psi) for (_, _, pressure_psi) in radial_data])

    p = plot(
        xlabel="Distance [km]",
        ylabel="Pressure [psi]",
        title="Pressure vs. Distance",
        ylim=(0, max_pressure + 100),
        legend=:topleft
    )

    for (wid, r_km, p_psi) in radial_data
        plot!(p, r_km, p_psi, label="Well $wid")
    end

    savefig(p, "pressure_distance_graph.png")
    println("Created pressure_distance_graph.png")
end







function plot_mohr_diagram_hydro(
    sh::Float64, 
    sH::Float64, 
    sV::Float64, 
    tau_effective::Union{Float64, AbstractVector{Float64}}, 
    sigma_effective::Union{Float64, AbstractVector{Float64}}, 
    p0::Float64, 
    biot::Float64, 
    nu::Float64, 
    dp::Union{AbstractVector{Float64}, AbstractVector{Integer}}, 
    strike::Vector{Float64},
    mu::Union{Float64, Integer}
    )

    N = length(strike)

    # Ensure dp is a vector (for the hydrology step) or broadcast scalar dp across all faults
    dp = dp isa Float64 ? fill(dp, N) : dp

    # create the stress vector
    Sig0sorted = sort([sh, sH, sV], rev=true)
    
    # find index of the vertical stress sV
    ixSv = findall(x -> x == sV, Sig0sorted)
    ixSv = ixSv[1]

    # Initialize storage for stress changes
    Ds = zeros(N, 3)  # Each row corresponds to a fault

    # Compute stress changes due to pore pressure for each fault
    for i in 1:N
        Ds[i, :] .= biot * (1 - 2 * nu) / (1 - nu) * dp[i]  # Apply dp[i] per fault
        Ds[i, ixSv] = 0.0  # Vertical stress remains unchanged
    end

    # Initialize Sig (each row corresponds to a fault)
    Sig = repeat(Sig0sorted', N, 1)  # Repeat sorted stresses N times as rows

    # Add stress changes to the initial stress state
    Sig .= Sig .+ Ds 

    # create angular variable for plotting the semicircles
    a = range(0, stop=pi, length=100) # use 100 points (can increase for quality)
    c = exp.(1im .* a) # complex exponential

    # Initialize Mohr circle centers
    C1 = zeros(Complex{Float64}, N, length(a))
    C2 = zeros(Complex{Float64}, N, length(a))
    C3 = zeros(Complex{Float64}, N, length(a))

    # radii of circles
    R1 = 0.5 .* (Sig[:, 1] .- Sig[:, 3])
    R2 = 0.5 .* (Sig[:, 2] .- Sig[:, 3])
    R3 = 0.5 .* (Sig[:, 1] .- Sig[:, 2])

    # Compute Mohr circles for each fault (apply dp[i] individually)
    for k in 1:N
        C1[k, :] = R1[k] .* c .+ (Sig[k, 1] + Sig[k, 3]) / 2 .- (p0 + dp[k])
        C2[k, :] = R2[k] .* c .+ (Sig[k, 2] + Sig[k, 3]) / 2 .- (p0 + dp[k])
        C3[k, :] = R3[k] .* c .+ (Sig[k, 1] + Sig[k, 2]) / 2 .- (p0 + dp[k])
    end

    # plot the Mohr circles
    p = plot(legend=:topright, xlabel="σ (Effective Normal Stress)", ylabel="τ (Shear Stress)", title="Mohr Diagram Hydrology")
    for k in 1:N
        plot!(real(C1[k, :]), imag(C1[k, :]), label="Circle 1 - Fault $k")
        plot!(real(C2[k, :]), imag(C2[k, :]), label="Circle 2 - Fault $k")
        plot!(real(C3[k, :]), imag(C3[k, :]), label="Circle 3 - Fault $k")
    end

    # Add fault locations
    scatter!(p, [tau_effective], [sigma_effective], label="Fault locations", color=:red)

    # frictional slip line
    sigma_range = range(0, stop=maximum(Sig[:, 1]), length=100)
    tau_line = mu .* sigma_range
    plot!(p, sigma_range, tau_line, label="Frictional Slip Line", color=:black, lw=2)

    savefig(p, "mohr_diagram_det_hydro.png")
    println("Mohr diagram saved as mohr_diagram_det_hydro.png")

end





# function to prepare the graph data for D3.js to plot a Mohr diagram
function mohr_diagram_data_to_d3(
    sh::Float64, 
    sH::Float64, sV::Float64, 
    tau_effective::Union{Float64, AbstractVector{Float64}}, 
    sigma_effective::Union{Float64, AbstractVector{Float64}}, 
    p0::Float64, 
    biot::Float64, 
    nu::Float64, 
    dp::Float64, 
    strike::Vector{Float64},
    mu::Union{Float64, Integer}
)


    N = length(strike)

    # create the stress tensor
    Sig0sorted = sort([sh, sH, sV], rev=true)

    # find index of the vertical stress sV
    ixSv = findall(x -> x ==sV, Sig0sorted)
    ixSv = ixSv[1]

    # change in principal stresses due to pore pressure
    Ds = kron(biot * (1-2*nu) / (1-nu)*dp, [1.0, 1.0, 1.0])

    # set vertical stress change to 0
    Ds[ixSv] = 0.0

    # initialize Sig
    Sig = repeat(Sig0sorted', N, 1) # Repeat the sorted stresses N times as rows

    # Expand Ds to match the dimensions of Sig
    Ds_expanded = repeat(Ds', N, 1)  # Repeat Ds as rows to match Sig

    # Add stress changes to the initial stress state
    Sig .= Sig .+ Ds_expanded

    # create angular variable for plotting the semicircles
    a = range(0, stop=pi, length=100) # use 100 points (can increase for quality)
    c = exp.(1im .* a) # complex exponential

    # C1, C2, C3 are the centers of the three semicircles
    C1 = zeros(Complex{Float64}, N, length(a))
    C2 = zeros(Complex{Float64}, N, length(a))
    C3 = zeros(Complex{Float64}, N, length(a))

    # radii of circles
    R1 = 0.5 .* (Sig[:, 1] .- Sig[:, 3])
    R2 = 0.5 .* (Sig[:, 2] .- Sig[:, 3])
    R3 = 0.5 .* (Sig[:, 1] .- Sig[:, 2])

    #preallocate arrays for rows
    rows = DataFrame(series=String[], x=Float64[], y=Float64[], stress_regime=String[])

    # ask why we need the loop for
    for k in 1:N
        # Semicircle 1
        circle1 = R1[k] .* c .+ (Sig[k, 1] + Sig[k, 3]) / 2 - (p0 + dp)
        for pt in circle1
            push!(rows, ("circle1", real(pt), imag(pt), stress_regime))
        end
        # Semicircle 2
        circle2 = R2[k] .* c .+ (Sig[k, 2] + Sig[k, 3]) / 2 - (p0 + dp)
        for pt in circle2
            push!(rows, ("circle2", real(pt), imag(pt), stress_regime))
        end
        # Semicircle 3
        circle3 = R3[k] .* c .+ (Sig[k, 1] + Sig[k, 2]) / 2 - (p0 + dp)
        for pt in circle3
            push!(rows, ("circle3", real(pt), imag(pt), stress_regime))
        end
    end

    faults = Dict("x" => tau_effective, "y" => sigma_effective)
    for (xval, yval) in zip(faults["x"], faults["y"])
        push!(rows, ("faults", xval, yval))
    end

    # Frictional slip line: sample points along the line
    sigma_range = range(0, stop=maximum(Sig[:, 1]), length=100)
    tau_line = mu .* sigma_range
    for (xval, yval) in zip(sigma_range, tau_line)
        push!(rows, ("frictional_slip_line", xval, yval))
    end

    CSV.write("graphs/mohr_diagram_data_d3.csv", rows)


end


function mohr_diagram_data_to_d3_portal(
    sh::Float64, 
    sH::Float64, 
    sV::Float64, 
    tau_effective::Union{Float64, AbstractVector{Float64}}, 
    sigma_effective::Union{Float64, AbstractVector{Float64}}, 
    p0::Float64, 
    biot::Float64, 
    nu::Float64, 
    dp::Float64, 
    strike::Vector{Float64},
    mu::Union{Float64, Integer},
    stress_regime::String,
    slip_pressure::Union{Float64, AbstractVector{Float64}},
    fault_ids::Union{Vector{String}, Vector{Int}, Vector{Any}}=["fault"]
)

    # check if tau_effective or sigma_effective is a vector 
    # if they are, make sure we set all negative values to 0
    if isa(tau_effective, AbstractVector)
        tau_effective = max.(tau_effective, 0.0)
    end
    if isa(sigma_effective, AbstractVector)
        sigma_effective = max.(sigma_effective, 0.0)
    end

    # chekc if tau_effective or shear_stress is a float, and if it's negative, set it to 0
    if isa(tau_effective, Float64)
        if tau_effective < 0.0
            tau_effective = 0.0
        end
    end

    if isa(sigma_effective, Float64)
        if sigma_effective < 0.0
            sigma_effective = 0.0
        end
    end
    


    N = length(strike)

    # Assign principal stresses based on the regime.
    # sh = min horizontal, sH = max horizontal, sV = vertical.
    if stress_regime == "Normal"
        S1 = sV    
        S2 = sH    
        S3 = sh   
        arc1_label = "σV"      # S1 = σV
        arc2_label = "σH"      # S2 = σH
        arc3_label = "σh"      # S3 = σh 
        vertical_index = 1  
    elseif stress_regime == "Reverse"
        S1 = sH    
        S2 = sh    
        S3 = sV    
        arc1_label = "σH"      # S1 = σH
        arc2_label = "σh"      # S2 = σh
        arc3_label = "σV"      # S3 = σV
        vertical_index = 3  
    elseif stress_regime == "Strike-Slip"
        S1 = sH    
        S2 = sV    
        S3 = sh    
        arc1_label = "σH"      # S1 = σH
        arc2_label = "σV"      # S2 = σV
        arc3_label = "σh"      # S3 = σh
        vertical_index = 2  
    else
        error("Invalid stress regime: $stress_regime. Must be one of 'Normal', 'Reverse', or 'Strike-Slip'.")
    end

    stress_array = [S1, S2, S3]
    
    # Pore pressure correction factor:
    Δ = biot * (1 - 2*nu) / (1 - nu) * dp
    Ds = Δ * ones(3)
    # Do not modify σV
    Ds[vertical_index] = 0.0
    stress_array .= stress_array .+ Ds

    S1, S2, S3 = stress_array[1], stress_array[2], stress_array[3]

    # Define arc parameters.
    # Arc 1: between S1 and S3
    centerX1 = (S1 + S3)/2 - (p0 + dp)
    radius1  = max(0, 0.5 * (S1 - S3)) # prevent negative radius
    
    # Arc 2: between S2 and S3
    centerX2 = (S2 + S3)/2 - (p0 + dp)
    radius2  = max(0, 0.5 * (S2 - S3))
    
    # Arc 3: between S1 and S2
    centerX3 = (S1 + S2)/2 - (p0 + dp)
    radius3  = max(0, 0.5 * (S1 - S2))
    
    # Calculate effective stress values for label positioning
    S1_eff = S1 - (p0 + dp)
    S2_eff = S2 - (p0 + dp)
    S3_eff = S3 - (p0 + dp)

    # --- Create Arc DataFrame ---
    arcsDF = DataFrame(
        label = String[],
        centerX = Float64[], 
        centerY = Float64[], 
        radius = Float64[],
        stress_regime = String[],
        start_angle = Float64[],
        end_angle = Float64[],
        labelPosX = Float64[], # we add 300 to the effective stress to avoid overlap with the arc lines
        labelPosY = Float64[]
    )
    push!(arcsDF, (arc1_label, centerX1, 0.0, radius1, stress_regime, 0.0, π, S1_eff+300, 0.0))
    push!(arcsDF, (arc2_label, centerX2, 0.0, radius2, stress_regime, 0.0, π, S2_eff+300, 0.0))
    push!(arcsDF, (arc3_label, centerX3, 0.0, radius3, stress_regime, 0.0, π, S3_eff+300, 0.0))
    
    # --- Compute Frictional Slip Line (common for all faults) ---
    sigma_range = collect(range(0, stop=S1, length=100))
    tau_line = mu .* sigma_range
    slipDF = DataFrame(label=String[], x=Float64[], y=Float64[])
    # start point for the slip line
    push!(slipDF, ("frictional_slip_line", 0.0, 0.0))
    # end point for the slip line
    max_sigma = S1
    max_tau = mu * max_sigma
    push!(slipDF, ("frictional_slip_line", max_sigma, max_tau))

    # --- Compute Fault Markers ---
    faultDF = DataFrame(fault_id=String[], x=Float64[], y=Float64[], pore_pressure_slip=Float64[])
    if isa(tau_effective, AbstractVector) && isa(sigma_effective, AbstractVector) && isa(slip_pressure, AbstractVector)
        # Make sure fault_ids has the right length, default to "fault#" if not provided
        if length(fault_ids) == 1 && length(tau_effective) > 1
            fault_ids = ["fault$i" for i in 1:length(tau_effective)]
        end
        
        for i in eachindex(tau_effective)
            te = tau_effective[i]
            se = sigma_effective[i]
            sp = slip_pressure[i]
            # Use actual fault ID instead of generic "fault"
            fault_id = string(fault_ids[min(i, length(fault_ids))])
            push!(faultDF, (fault_id, se, te, sp))
        end
    else
        # Single fault case
        push!(faultDF, (string(fault_ids[1]), sigma_effective, tau_effective, slip_pressure))
    end
    
    
    
    return (arcsDF, slipDF, faultDF)
end

# helper function for Mohr diagram (maps color to fault proximity)
function get_fault_color_mohr(tau, sigma, mu, min_psi_tolerance=0.0, max_psi_tolerance=4000.0)

    # find distance from fault to the frictional slip line
    distance = abs(tau - mu * sigma)

    # clamp to non-negative values (if faults is below the frictional slip line)
    distance = max(distance, 0.0)

    # tolerance range
    tolerance_range = max_psi_tolerance - min_psi_tolerance

    # normalize the distance to the tolerance range
    normalized_distance = (distance - min_psi_tolerance) / tolerance_range


    # normalize the distance to a value between 0 and 1
    normalized_distance = clamp(normalized_distance, 0.0, 1.0)

    # interpolate between green and red
    # when normalized_distance = 0, color = red (255,0,0)
    # when normalized_distance = 1, color = green (0,255,0)
    red = round(Int, 255 * normalized_distance)
    green = round(Int, 255 * (1.0 - normalized_distance))
    blue = 0
    
    # return the hex color
    return @sprintf("#%02X%02X%02X", red, green, blue)

end


# Plots a line graph for injection well
# x axis: time (years)
# y axis: injection rate (bbl/day)

function plot_injection_rate_line_chart(
    well_data::Dict{String,Any},
    extrapolate::Bool=false
)
    if isempty(well_data)
        println("No well data to plot. Skipping injection rate line chart.")
        return
    end

    injection_wells = well_data["injection_wells"]

    # Initialize the plot
    p = plot(
        xlabel="Time [years]",
        ylabel="Injection Rate [bbl/day]",
        title="Injection Rate vs. Time",
        legend=:topleft
    )

    # D3.js module data for line graph
    d3_line_graph_data = Dict("data" => Dict(), "config" => Dict(
        "xlabel" => "Time [years]",
        "ylabel" => "Injection Rate [bbl/day]",
        "title" => "Injection Rate vs. Time",
        "colors" => Dict()
    ))

    # Iterate through wells
    for (well_id, well_info) in injection_wells
        injection_rate_info = well_info["injection_rate"]
        rate_mode = injection_rate_info["rate_mode"]
        start_year = injection_rate_info["start_year"]
        end_year = injection_rate_info["end_year"]

        if rate_mode == "constant"
            injection_rate = injection_rate_info["rate_bbl_day"]
            years = collect(start_year:end_year)
            rates = [injection_rate for _ in years]

            # Plot constant rates
            plot!(p, years, rates, label="Well $well_id")

            # add the data to the graph JSON
            d3_line_graph_data["data"]["Well $well_id"] = [
                Dict("x" => year, "y" => rate) for year in years
            ]
            d3_line_graph_data["config"]["colors"]["Well $well_id"] = "auto"
            

        elseif rate_mode == "monthly"
            monthly_bbl = injection_rate_info["monthly_bbl"]
            time_steps = []
            rates = []

            for year in start_year:end_year
                if haskey(monthly_bbl, string(year))
                    for month in 1:12
                        month_start = year + (month - 1) / 12
                        month_end = year + month / 12

                        push!(time_steps, month_start)
                        push!(rates, monthly_bbl[string(year)][string(month)] / day(Dates.lastdayofmonth(Date(year, month, 1))))

                        # Add step at the end of the month
                        push!(time_steps, month_end)
                        push!(rates, monthly_bbl[string(year)][string(month)] / day(Dates.lastdayofmonth(Date(year, month, 1))))
                    end
                elseif extrapolate && !isempty(rates)
                    for month in 1:12
                        month_start = year + (month - 1) / 12
                        month_end = year + month / 12

                        push!(time_steps, month_start)
                        push!(rates, rates[end])

                        # Add step at the end of the month
                        push!(time_steps, month_end)
                        push!(rates, rates[end])
                    end
                end
            end

            # Plot variable rates with step changes
            plot!(p, time_steps, rates, label="Well $well_id", lw=2, linestyle=:solid)

            # add the data to the graph JSON
            d3_line_graph_data["data"]["Well $well_id"] = [
                Dict("x" => time_steps[i], "y" => rates[i]) for i in 1:length(time_steps)
            ]
            d3_line_graph_data["config"]["colors"]["Well $well_id"] = "auto"
        end
    end

    # prepare the JSON for D3 module


    # Save the plot
    println("Saving graph as injection_rate_line_chart.png")
    savefig(p, "injection_rate_line_chart.png")

    # output JSON file for D3.js module
    open("graphs/injection_rate_line_chart.json", "w") do file
        write(file, JSON.json(d3_line_graph_data))
    end
    println("Saved D3.js data as injection_rate_line_chart.json")
end



# TO DO: Make this read the CSV directly
function injection_rate_line_chart_data_to_d3(well_data::Vector{Dict{String,Any}},
    extrapolate::Bool=false
    )

    if isempty(well_data)
        println("No well data to plot. Skipping injection rate line chart.")
        return
    end

    injection_wells = well_data["injection_wells"]

    d3_line_graph_data = Dict("data" => Dict(), "config" => Dict(
        "xlabel" => "Time [years]",
        "ylabel" => "Injection Rate [bbl/day]",
        "title" => "Injection Rate vs. Time",
        "colors" => Dict()
    ))

    # Iterate through wells
    for (well_id, well_info) in injection_wells
        
        start_year = injection_rate["start_year"]
        end_year = injection_rate["end_year"]

        injection_rate = well_info["rate_bbl_day"]

    end




end


# STEP 1: Model Inputs --> Fault Surface Map (2D)
function plot_fault_surface_map(
    input_json::Dict{String,Any}, 
    graph_output_dir::String = "graphs"
)

    if !isdir(graph_output_dir)
        mkdir(graph_output_dir)
    end


    # Extract the fault data
    fault_data = input_json["faults"]
    # Extract the well data
    well_data = input_json["injection_wells"]

    # initialize plot
    p = plot(
        xlabel="X Easting [km]",
        ylabel="Y Northing [km]",
        title="Fault Surface Map",
        legend=:topright,
        aspect_ratio=:equal,
        grid=true
    )


    #prepare data for D3 module
    d3_fault_surface_map_data = Dict(
        "faults" => [],
        "wells" => []
    )

    # add faults to the graph
    for fault in fault_data
        x = fault["easting"]
        y = fault["northing"]
        length = fault["length_km"]
        strike = fault["strike"]
        fault_id = fault["fault_id"]

        # calculate the end points of the fault
        half_length = length / 2
        angle_radians = deg2rad(strike)

        dx = half_length * cos(angle_radians)
        dy = half_length * sin(angle_radians)

        x_start = x - dx
        x_end = x + dx
        y_start = y - dy
        y_end = y + dy

        plot!([x_start, x_end], [y_start, y_end], label="Fault $(fault["fault_id"])", lw=2, color=:black)

        # add the fault to the D3 module data
        push!(d3_fault_surface_map_data["faults"], Dict(
            "id" => fault_id,
            "start" => Dict("x" => x_start, "y" => y_start),
            "end" => Dict("x" => x_end, "y" => y_end)
        ))
    
    end

    # add wells to the graph
    for (well_id, well_info) in well_data
        x = well_info["location"]["easting_km"]
        y = well_info["location"]["northing_km"]

        scatter!([x], [y], label="Well $well_id", color=:red, marker=:circle, ms=5)
        annotate!(x+1, y+1, text("Well $well_id", :blue, 10))

        # add the well to the D3 module data
        push!(d3_fault_surface_map_data["wells"], Dict(
            "id" => well_id,
            "location" => Dict("x" => x, "y" => y),
            "color" => "red"
        ))
    end

    # save plot
    output_path = joinpath(graph_output_dir, "fault_surface_map.png")
    println("Saving fault surface map to $output_path")
    savefig(p, output_path)

    # save D3 module data
    d3_output_path = joinpath(graph_output_dir, "fault_surface_map.json")
    open(d3_output_path, "w") do file
        write(file, JSON.json(d3_fault_surface_map_data))
    end
    println("Saved D3.js data as $d3_output_path")

end




# prepares the data for the fault surface map in model Inputs
# read DF for faults and wells 
# used to parse the data from the portal
# Helper function for degrees-to-radians conversion
deg2rad(deg) = deg * (π / 180)

function fault_surface_map_data_to_d3(faults_df::DataFrame, wells_df::DataFrame)
    try
        if faults_df === nothing || wells_df === nothing
            throw(ArgumentError("Input dataframes cannot be 'nothing'"))
        end

        # Helper function to normalize column names
        normalize_colname(col) = strip(lowercase(string(col)))
        
        # Get normalized column names
        fault_cols = normalize_colname.(names(faults_df))
        well_cols = normalize_colname.(names(wells_df))

        # Column validation for faults (case insensitive)
        required_fault_cols = ["faultid", "lengthkm", "strike", "latitude(wgs84)", "longitude(wgs84)", "frictioncoefficient", "dip"]
        missing_fault_cols = filter(col -> !in(normalize_colname(col), fault_cols), required_fault_cols)
        if !isempty(missing_fault_cols)
            throw(ArgumentError("Faults dataframe is missing columns: $(join(missing_fault_cols, ", "))"))
        end

        # Column validation for wells (case insensitive)
        required_well_cols = ["apinumber", "latitude(wgs84)", "longitude(wgs84)"]
        missing_well_cols = filter(col -> !in(normalize_colname(col), well_cols), required_well_cols)
        if !isempty(missing_well_cols)
            throw(ArgumentError("Wells dataframe is missing columns: $(join(missing_well_cols, ", "))"))
        end

        # Create output tables
        fault_table = DataFrame(
            "FaultID" => String[],
            "Latitude(WGS84)" => Float64[],
            "Longitude(WGS84)" => Float64[],
            "Color" => String[]
        )

        # Use distinguishable_colors for fault coloring
        num_faults = nrow(faults_df)
        palette = distinguishable_colors(num_faults)

        # Process faults
        for (i, row) in enumerate(eachrow(faults_df))
            try
                fault_id = string(row[findfirst(==(normalize_colname("FaultID")), fault_cols)])
                lat = row[findfirst(==(normalize_colname("Latitude(WGS84)")), fault_cols)]
                lon = row[findfirst(==(normalize_colname("Longitude(WGS84)")), fault_cols)]
                length_km = row[findfirst(==(normalize_colname("LengthKm")), fault_cols)]
                strike_deg = row[findfirst(==(normalize_colname("Strike")), fault_cols)]

                # Calculate start/end points based on length and strike
                # For lat/lon calculations, we need to account for Earth's curvature
                # Using approximation: 1 degree of latitude ≈ 111 km, 1 degree of longitude ≈ 111*cos(lat) km
                half_length = length_km / 2
                angle_radians = deg2rad(90 - strike_deg)  # Convert from azimuth to math angle
                
                # Convert distance to degrees (approximate)
                lat_km_per_degree = 111.0
                lon_km_per_degree = 111.0 * cos(deg2rad(lat))
                
                dlat = (half_length * cos(angle_radians)) / lat_km_per_degree
                dlon = (half_length * sin(angle_radians)) / lon_km_per_degree

                # Calculate endpoints
                lat_start = lat - dlat
                lon_start = lon - dlon
                lat_end = lat + dlat
                lon_end = lon + dlon

                fault_color = "#" * hex(palette[i])

                # Push start point
                push!(fault_table, (
                    fault_id,
                    lat_start,
                    lon_start,
                    fault_color
                ))

                # Push end point
                push!(fault_table, (
                    fault_id,
                    lat_end,
                    lon_end,
                    fault_color
                ))

            catch e
                @warn "Error processing fault $i: $(sprint(showerror, e))"
                continue
            end
        end

        # Create wells table
        well_table = DataFrame(
            "APINumber" => String[],
            "Latitude(WGS84)" => Float64[],
            "Longitude(WGS84)" => Float64[],
            "Color" => String[]
        )

        well_color = "#" * hex(RGB(1,0.5,0))

        # Process wells
        for (i, row) in enumerate(eachrow(wells_df))
            try
                well_id = string(row[findfirst(==(normalize_colname("APINumber")), well_cols)])
                lat = Float64(row[findfirst(==(normalize_colname("Latitude(WGS84)")), well_cols)])
                lon = Float64(row[findfirst(==(normalize_colname("Longitude(WGS84)")), well_cols)])

                push!(well_table, (
                    well_id,
                    lat,
                    lon,
                    well_color
                ))

            catch e
                @warn "Error processing well $i: $(sprint(showerror, e))"
                continue
            end
        end

        return fault_table, well_table

    catch e
        throw(ErrorException("Failed to process fault/well data: $(sprint(showerror, e))"))
    end
end


function year_to_date(year::Int)
    return Date(year, 1, 1)
end


function date_to_js_timestamp(date::Date)::BigInt
    return Dates.datetime2unix(DateTime(date)) * 1000 # seconds to milliseconds
end





function injection_rate_data_to_d3(well_df::DataFrame, injection_data_type::String)

    # possible types: 'annual_fsp', 'monthly_fsp', 'injection_tool_data'

    # Helper function to normalize column names
    normalize_colname(col) = strip(lowercase(string(col)))
    
    # get the columns names of the df
    colnames = lowercase.(strip.(names(well_df)))

    #println("DEBUG: Input well dataframe columns: $(colnames)")
    try
        
        if injection_data_type == "annual_fsp"

            # well data for each row
            wells_reformatted = DataFrame(
                "WellID" => String[],
                "InjectionRate(bbl/month)" => Float64[],
                "Year" => Int[],
                "Month" => Int[],
                "Date" => String[],
                "Timestamp" => Float64[]
            )

            

            for (i, row) in enumerate(eachrow(well_df))
                println("Iterating over rows...")
                
                try
                    well_id = string(row["WellID"])
                    injection_rate_daily = row["InjectionRate(bbl/day)"]
                    start_year = row["StartYear"]
                    end_year = row["EndYear"]
                    
                    # For each year in the range
                    for year in start_year:end_year
                        # For each month in the year
                        for month in 1:12
                            # Define the date for the start of the month
                            start_date = Date(year, month, 1)

                            # Skip months before start date if it's the first year
                            # (Adjust logic slightly to check against start_date)
                            # Assuming StartYear/EndYear define the range inclusively at the year level
                            # If more precise start/end dates are needed, this logic might need refinement
                            if year == start_year && start_year != end_year # Add check if start/end year are different
                                # Simplification: Include all months if start/end is the same year or range spans multiple years
                            elseif year < start_year || year > end_year
                                continue
                            end

                            # Calculate days in the month for conversion
                            days_in_month = Dates.daysinmonth(start_date)
                            
                            # Convert daily rate to monthly rate
                            injection_rate_monthly = injection_rate_daily * days_in_month
                            
                            # Get end date of the month
                            end_date = Dates.lastdayofmonth(start_date)

                            # Create date string (month/day/year) - use start date for representation
                            date_string = Dates.format(start_date, "m/d/Y")

                            # convert start and end dates to unix timestamps
                            start_timestamp = date_to_js_timestamp(start_date)
                            end_timestamp = date_to_js_timestamp(end_date)
                            
                            # Add start-of-month data point
                            push!(wells_reformatted, (
                                String(well_id),
                                injection_rate_monthly,
                                year,
                                month,
                                date_string, # Representational date string
                                start_timestamp
                            ))
                            # Add end-of-month data point with the same rate
                             push!(wells_reformatted, (
                                String(well_id),
                                injection_rate_monthly,
                                year,
                                month,
                                date_string, # Representational date string
                                end_timestamp
                            ))
                        end
                    end
                    
                    println("Processed well: $well_id, converted daily rate $injection_rate_daily to monthly for $start_year-$end_year")
                    
                catch e
                    error("Error processing well $i: $(sprint(showerror, e))")
                end
            end

            
            
            return wells_reformatted

        # monthly injection rates (FSP format)
        elseif injection_data_type == "monthly_fsp"

            println("we have fsp monthly injection data in injection_rate_data_to_d3")

            # well data for each row
            wells_reformatted = DataFrame(
                "WellID" => String[],
                "InjectionRate(bbl/month)" => Float64[],
                "Month" => Int[],
                "Year" => Int[],
                "Date" => String[],
                "Timestamp" => Float64[]
            )

            # Get unique Well IDs numbers
            unique_well_ids = unique(string.(well_df.WellID))
            num_unique_wells = length(unique_well_ids)
            
            
            
            for (i, row) in enumerate(eachrow(well_df))
                try
                    well_id = string(row["WellID"])
                    injection_rate = row["InjectionRate(bbl/month)"]
                    month = row["Month"]
                    year = row["Year"]
                    start_date = Date(year, month, 1)
                    end_date = Dates.lastdayofmonth(start_date)
                    
                    # Format date as "month/day/year" - use start date for representation
                    date_formatted = Dates.format(start_date, "m/d/Y")

                    # convert start and end dates to unix timestamps
                    start_timestamp = date_to_js_timestamp(start_date)
                    end_timestamp = date_to_js_timestamp(end_date)

                    # Add start-of-month data point
                    push!(wells_reformatted, (
                        String(well_id),
                        injection_rate,
                        month,
                        year,
                        date_formatted, # Representational date string
                        start_timestamp
                    ))
                    # Add end-of-month data point with the same rate
                    push!(wells_reformatted, (
                        String(well_id),
                        injection_rate,
                        month,
                        year,
                        date_formatted, # Representational date string
                        end_timestamp
                    ))
                catch e
                    error("Error processing well $i: $(sprint(showerror, e))")
                end
            end

            
            
            

            
            return wells_reformatted

        elseif injection_data_type == "injection_tool_data"

            # First, ensure "Date of Injection" is a Date object
            # TO DO: try parsing it as Date fallback option
            if eltype(well_df[!, "Date of Injection"]) <: AbstractString
                well_df[!, "Date of Injection"] = Date.(well_df[!, "Date of Injection"])
            end

            # Get unique API numbers
            unique_api_numbers = unique(string.(well_df[!, "API Number"]))

            # print the unique api numbers
            println("Unique API numbers: $unique_api_numbers")
            
            
            # well data for each row
            wells_reformatted = DataFrame(
                "WellID" => String[],
                "InjectionRate(bbl/month)" => Float64[],
                "Year" => Int[],
                "Month" => Int[],
                "Date" => String[],
                "Timestamp" => Float64[]
            )

            # Process each well separately to calculate monthly averages
            for (i, well_id) in enumerate(unique_api_numbers)
                println("Processing well $well_id...")
                
                try
                    # Filter data for this well
                    well_data = well_df[well_df[!, "API Number"] .== well_id, :]
                    
                    # Create temporary DataFrame with date and injection rate
                    temp_df = DataFrame(
                        "Date" => well_data[!, "Date of Injection"],
                        "InjectionRate" => well_data[!, "Volume Injected (BBLs)"]
                    )
                    
                    # Add year and month columns
                    temp_df[!, :Year] = year.(temp_df.Date)
                    temp_df[!, :Month] = month.(temp_df.Date)
                    
                    # Group by year and month, calculate sum (total monthly volume)
                    monthly_totals = combine(groupby(temp_df, [:Year, :Month]), 
                        :InjectionRate => sum => :MonthlyInjectionRate)
                    
                    
                    
                    # Add each monthly data point to our result DataFrame
                    for row in eachrow(monthly_totals)
                        year_val = row.Year
                        month_val = row.Month
                        monthly_rate = row.MonthlyInjectionRate
                        start_date_obj = Date(year_val, month_val, 1)
                        end_date_obj = Dates.lastdayofmonth(start_date_obj)
                        
                        date_formatted = Dates.format(start_date_obj, "m/d/Y") # Use start date for representation
                        
                        # convert start and end dates to unix timestamps
                        start_timestamp = date_to_js_timestamp(start_date_obj)
                        end_timestamp = date_to_js_timestamp(end_date_obj)
                        
                        # Add start-of-month data point
                        push!(wells_reformatted, (
                            String(well_id),
                            monthly_rate,
                            year_val,
                            month_val,
                            date_formatted, # Representational date string
                            start_timestamp
                        ))
                        # Add end-of-month data point with the same rate
                        push!(wells_reformatted, (
                            String(well_id),
                            monthly_rate,
                            year_val,
                            month_val,
                            date_formatted, # Representational date string
                            end_timestamp
                        ))
                    end
                    
                    println("Processed well $well_id with $(nrow(monthly_totals)) monthly data points")
                    
                catch e
                    error("Error processing well $well_id: $(sprint(showerror, e))")
                end
            end

            

            return wells_reformatted
        end

    catch e
        @warn "Error processing injection data: $(sprint(showerror, e))"
        return nothing
    end
end

        







#=
function plot_cdf_det_hydro(prob_geo_values, det_hydro_values)
    plot()
    for (fault_id, fault_data) in prob_geo_values
        pressures = sort(fault_data["pressures"])
        probabilities = range(0, stop=1, length=length(pressures))

        # CDF of prob geomechanics
        plot!(pressures, probabilities, label="Fault $fault_id", linestyle=:dash, lw=2)

        # plot the vertical slicing lines (deterministic hydrology pressure values)
        if haskey(det_hydro_values, fault_id)
            pp_value = det_hydro_values[fault_id]
            plot!([pp_value, pp_value], [0, 1], label="", color=:blue, lw=2)
        end

        xlabel!("Δ Pore Pressure on fault [psi]")
        ylabel!("Probability")
        title!("CDF of Pore Pressure on Fault")
        savefig("graphs/cdf_det_hydro.png")
        println("Prob hydrology CDF saved as cdf_det_hydro.png")
    end
    
end
=#

function plot_cdf_det_hydro(prob_geo_values, det_hydro_values)
    # Generate distinguishable colors based on number of faults
    n_faults = length(prob_geo_values)
    colors = distinguishable_colors(n_faults)
    
    plot()
    
    for (fault_index, (fault_id, fault_data)) in enumerate(prob_geo_values)
        pressures = sort(fault_data["pressures"])
        probabilities = range(0, stop=1, length=length(pressures))

        # Get the distinguishable color for this fault
        current_color = colors[fault_index]

        # CDF of prob geomechanics
        plot!(pressures, probabilities, label="Fault $fault_id", linestyle=:dash, lw=2, color=current_color)

        # Plot the vertical slicing lines with matching color
        if fault_index <= length(det_hydro_values)
            pp_value = det_hydro_values[fault_index]
            plot!([pp_value, pp_value], [0, 1], label="", color=current_color, lw=2)
        end

        xlabel!("Δ Pore Pressure on fault [psi]")
        ylabel!("Probability")
        title!("CDF of Pore Pressure on Fault")
    end
    
    savefig("graphs/cdf_det_hydro.png")
    println("Prob hydrology CDF saved as cdf_det_hydro.png")
end




function plot_prob_hydro_combined_cdf(prob_geo_values, prob_hydro_results)
    plot()
    
    for (fault_id, fault_data) in prob_geo_values
        pressures = sort(fault_data["pressures"])
        probabilities = range(0, stop=1, length=length(pressures))
        plot!(pressures, probabilities, label="Geomechanics - $(fault_id)", lw=2, linestyle=:dash)
    end

    # Handle probabilistic hydrology results
    for (fault_idx, hydro_data) in enumerate(prob_hydro_results)
        if isa(hydro_data, Vector)
            pressures = sort(hydro_data)
            probabilities = range(0, stop=1, length=length(pressures))
            plot!(pressures, probabilities, label="Hydrology - Fault $(fault_idx+1)", lw=2, linestyle=:solid)
        else
            println("Warning: Skipping fault $(fault_idx+1) as hydrology data is not a vector.")
        end
    end

    xlabel!("Δ Pore Pressure on fault [psi]")
    ylabel!("Probability")
    title!("Combined CDF Curves for Faults")
    savefig("graphs/prob_hydro_combined_cdf.png")
end



########################## PORTAL GRAPHS #########################################

function prob_geomechanics_cdf(mc_results_df::DataFrame, det_geomechanics_results_df::DataFrame)
    points_df = DataFrame(
        "slip_pressure" => Float64[],
        "probability" => Float64[],
        "ID" => String[],
        "det_slip_pressure" => Float64[]
    )

    # check if the FaultID column of the det_geomechanics_results_df is a string, if not, convert it to a string
    if !(eltype(det_geomechanics_results_df[!, "FaultID"]) <: AbstractString)
        det_geomechanics_results_df[!, "FaultID"] = string.(det_geomechanics_results_df[!, "FaultID"])
    end

    # get the unique fault IDs
    fault_ids = unique(mc_results_df.FaultID)

    

    for fault_id in fault_ids
        
        # Get the 'slip pressure' column data for this fault
        fault_data = mc_results_df[mc_results_df.FaultID .== fault_id, :SlipPressure]


        # also get the deterministic slip pressure for this fault
        det_slip_pressure_df = det_geomechanics_results_df[det_geomechanics_results_df.FaultID .== fault_id, :slip_pressure]
       

        
        if isempty(det_slip_pressure_df)
            # Skip if no deterministic data for this fault
            continue
        end
        
        # Extract the single value (assuming one result per fault in deterministic data)
        det_slip_pressure = det_slip_pressure_df[1]

        # Skip if no data for this fault
        number_of_points = length(fault_data)
        if number_of_points == 0
            continue 
        end

        # Convert to vector of Float64 if needed (to handle Any type arrays)
        fault_data_float = convert(Vector{Float64}, fault_data)
        
        # Use StatsBase.ecdf to create the empirical CDF function
        ecdf_func = ecdf(fault_data_float)
        
        # Get sorted pressure values for creating the CDF points
        sorted_pressure_values = sort(fault_data_float)
        
        # Evaluate the ECDF function at each sorted pressure value
        # Note: StatsBase's ecdf returns probabilities in [0,1]
        for pressure in sorted_pressure_values
            probability = ecdf_func(pressure)
            push!(points_df, (pressure, probability, fault_id, det_slip_pressure))
        end
        
        #= Original implementation (commented out)
        # Sort the values
        sorted_pressure_values = sort(fault_data)

        # Cumulative probabilities (0 to 100%)
        number_of_points = length(sorted_pressure_values)
        if number_of_points == 0
            continue # Skip if no data for this fault
        end
        probability_values = range(0, stop=100, length=number_of_points)
        
        # Add points for this fault to the DataFrame
        for (pressure, probability) in zip(sorted_pressure_values, probability_values)
            push!(points_df, (pressure, probability, fault_id, det_slip_pressure))
        end
        =#
    end

    

    return points_df
end


"""
Creates a DataFrame representing the probability of exceedance curve
for probabilistic hydrology pore pressure results.
Takes raw Monte‑Carlo hydrology results (a DataFrame with one row per iteration per fault, 
and turns it into an exceedance‑probability curve for each fault.
"""
function prob_hydrology_cdf(prob_hydro_results_df::DataFrame)

    points_df = DataFrame(
        "slip_pressure" => Float64[],
        "probability" => Float64[],
        "ID" => String[]
    )

    # Get the unique fault IDs
    fault_ids = unique(prob_hydro_results_df.ID)
    
    for fault_id in fault_ids
        # Get the 'Pressure' column data for this fault from the probabilistic hydrology results
        fault_data = prob_hydro_results_df[prob_hydro_results_df.ID .== fault_id, :Pressure]

        # Skip if no data for this fault
        number_of_points = length(fault_data)
        if number_of_points == 0
            continue 
        end

        # Convert to vector of Float64 if needed 
        fault_data_float = convert(Vector{Float64}, fault_data)
        
        # create the empirical CDF function
        ecdf_func = ecdf(fault_data_float)
        
        # sort the values so we can create the exceedance curve
        sorted_pressure_values = sort(fault_data_float)
        
        
        # Note: StatsBase's ecdf returns probabilities in [0,1]
        for pressure in sorted_pressure_values
            # Convert CDF to exceedance probability (1-CDF)
            exceedance_probability = (1.0 - ecdf_func(pressure)) 
            push!(points_df, (pressure, exceedance_probability, fault_id))
        end
    end

    

    return points_df
end




function uncertainty_variability_inputs_to_d3(
    uncertainties::Dict,
    stress_model_type::String,
    stress_inputs::Dict,
    fault_inputs::DataFrame
)
    #println("\n====== DEBUG: Starting uncertainty_variability_inputs_to_d3 ======")
    #println("Stress Model Type: $stress_model_type")
    
    # Create mapping for display names
    stress_param_mapping = Dict()
    
    # Set the parameter mapping based on the stress model type
    if stress_model_type == "gradients" || stress_model_type == "all_gradients"
        stress_param_mapping = Dict(
            "vertical_stress_gradient_uncertainty" => "Vert Stress Grad",
            "initial_pore_pressure_gradient_uncertainty" => "Pore Press Grad",
            "max_stress_azimuth_uncertainty" => "SHmax Azimuth",
            "max_horizontal_stress_uncertainty" => "SHmax Gradient",
            "min_horizontal_stress_uncertainty" => "SHmin Gradient"
        )
    elseif stress_model_type == "aphi_model" && stress_inputs["min_horizontal_stress"] !== nothing
        stress_param_mapping = Dict(
            "vertical_stress_gradient_uncertainty" => "Vert Stress Grad",
            "initial_pore_pressure_gradient_uncertainty" => "Pore Press Grad",
            "max_stress_azimuth_uncertainty" => "SHmax Azimuth",
            "aphi_value_uncertainty" => "APhi Value",
            "min_horizontal_stress_uncertainty" => "SHmin Gradient"
        )
    elseif stress_model_type == "aphi_model" && stress_inputs["min_horizontal_stress"] === nothing
        stress_param_mapping = Dict(
            "vertical_stress_gradient_uncertainty" => "Vert Stress Grad",
            "initial_pore_pressure_gradient_uncertainty" => "Pore Press Grad",
            "max_stress_azimuth_uncertainty" => "SHmax Azimuth",
            "aphi_value_uncertainty" => "APhi Value"
        )
    end
    
    # Map uncertainty parameters to fault parameters
    fault_param_mapping = Dict(
        "strike_angles_uncertainty" => "Strike of fault",
        "dip_angles_uncertainty" => "Dip of fault",
        "friction_coefficient_uncertainty" => "Friction Coeff"
    )
    
    # Combine all parameters
    parameter_mapping = merge(stress_param_mapping, fault_param_mapping)
    
    #println("Combined Parameter Mapping Keys: $(sort(collect(keys(parameter_mapping))))")
    
    # Create DataFrame to store results with consistent column order
    uncertainty_df = DataFrame(
        label = String[],
        min = Float64[],
        max = Float64[],
        id = Int[]
    )
    
    # Process each parameter in the mapping
    id_counter = 1
    
    for (uncertainty_param, display_name) in parameter_mapping
        #println("Processing parameter: $uncertainty_param => $display_name")
        
        # Skip if uncertainty is not defined or is 0
        if !haskey(uncertainties, uncertainty_param) || isnothing(uncertainties[uncertainty_param]) || uncertainties[uncertainty_param] == 0
            println("  SKIPPING: Uncertainty not defined or is 0")
            continue
        end
        
        uncertainty_value = uncertainties[uncertainty_param]
        
        # Get the base parameter value for min/max calculation
        base_value = 0.0
        
        if haskey(stress_param_mapping, uncertainty_param)
            # Stress parameter - get corresponding base value
            if uncertainty_param == "vertical_stress_gradient_uncertainty"
                base_value = get(stress_inputs, "vertical_stress", 0.0)
            elseif uncertainty_param == "initial_pore_pressure_gradient_uncertainty"
                base_value = get(stress_inputs, "pore_pressure", 0.0)
            elseif uncertainty_param == "max_stress_azimuth_uncertainty"
                base_value = get(stress_inputs, "max_stress_azimuth", 0.0)
            elseif uncertainty_param == "max_horizontal_stress_uncertainty"
                base_value = get(stress_inputs, "max_horizontal_stress", 0.0)
            elseif uncertainty_param == "min_horizontal_stress_uncertainty"
                base_value = get(stress_inputs, "min_horizontal_stress", 0.0)
            elseif uncertainty_param == "aphi_value_uncertainty"
                base_value = get(stress_inputs, "aphi_value", 0.0)
            end
        else
            # Fault parameter - use average from fault_inputs
            if uncertainty_param == "strike_angles_uncertainty"
                if haskey(fault_inputs[1, :], "Strike")
                    base_value = mean(fault_inputs[:, "Strike"])
                else
                    println("  WARNING: 'Strike' column not found in fault_inputs")
                    base_value = 0.0
                end
            elseif uncertainty_param == "dip_angles_uncertainty"
                if haskey(fault_inputs[1, :], "Dip")
                    base_value = mean(fault_inputs[:, "Dip"])
                else
                    println("  WARNING: 'Dip' column not found in fault_inputs")
                    base_value = 0.0
                end
            elseif uncertainty_param == "friction_coefficient_uncertainty"
                if haskey(fault_inputs[1, :], "FrictionCoefficient")
                    base_value = mean(fault_inputs[:, "FrictionCoefficient"])
                else
                    println("  WARNING: 'FrictionCoefficient' column not found in fault_inputs")
                    base_value = 0.0
                end
            end
        end
        
        #println("  Base value: $base_value, Uncertainty: $uncertainty_value")
        
        # Skip if base value is 0 to avoid division by zero
        if base_value == 0.0
            println("  SKIPPING: Base value is zero, cannot calculate percentage")
            continue
        end
        
        # Calculate percentage deviation (not absolute min/max)
        percent_deviation = (uncertainty_value / base_value) * 100
        
        # Round to 2 decimal places
        percent_deviation = round(percent_deviation, digits=2)
        
        #println("  Calculated percentage deviation: $percent_deviation%")
        
        # Add to DataFrame with negative and positive percentage deviations
        push!(uncertainty_df, (
            label = display_name,
            min = -percent_deviation,  # Negative percentage deviation
            max = percent_deviation,   # Positive percentage deviation
            id = id_counter
        ))
        
        id_counter += 1
    end
    
    #println("Created uncertainty variability data with $(nrow(uncertainty_df)) parameters")
    #println("====== DEBUG: Ending uncertainty_variability_inputs_to_d3 ======\n")
    
    return uncertainty_df
end

# Helper function to map parameter names from stress model
function stress_model_mapping(param_name::String, stress_model_type::String)
    if param_name == "vertical_stress_gradient"
        return "vertical_stress"
    elseif param_name == "initial_pore_pressure_gradient"
        return "pore_pressure"
    elseif param_name == "max_stress_azimuth"
        return "max_stress_azimuth"
    elseif param_name == "aphi_value"
        return "aphi_value"
    elseif param_name == "max_horizontal_stress"
        return "max_horizontal_stress"
    elseif param_name == "min_horizontal_stress"
        return "min_horizontal_stress"
    else
        return param_name
    end
end

function fault_sensitivity_tornado_chart_to_d3(
    base_stress_inputs::Dict, 
    base_fault_inputs::DataFrame, 
    uncertainties::Dict, 
    fault_idx::Int,
    stress_model_type::String,
    friction_coefficient::Float64
)
    # Setup the parameters to test
    stress_param_mapping = Dict()

    
    
    # Set the parameter mapping based on the stress model type
    if stress_model_type == "gradients" || stress_model_type == "all_gradients"
        stress_param_mapping = Dict(
            "vertical_stress_gradient_uncertainty" => "Vert Stress Grad",
            "initial_pore_pressure_gradient_uncertainty" => "Pore Press Grad",
            "max_stress_azimuth_uncertainty" => "SHmax Azimuth",
            "max_horizontal_stress_uncertainty" => "SHmax Gradient",
            "min_horizontal_stress_uncertainty" => "SHmin Gradient"
        )
    elseif stress_model_type == "aphi_model" && base_stress_inputs["min_horizontal_stress"] !== nothing
        stress_param_mapping = Dict(
            "vertical_stress_gradient_uncertainty" => "Vert Stress Grad",
            "initial_pore_pressure_gradient_uncertainty" => "Pore Press Grad",
            "max_stress_azimuth_uncertainty" => "SHmax Azimuth",
            "aphi_value_uncertainty" => "APhi Value",
            "min_horizontal_stress_uncertainty" => "SHmin Gradient"
        )
    elseif stress_model_type == "aphi_model" && base_stress_inputs["min_horizontal_stress"] === nothing
        stress_param_mapping = Dict(
            "vertical_stress_gradient_uncertainty" => "Vert Stress Grad",
            "initial_pore_pressure_gradient_uncertainty" => "Pore Press Grad",
            "max_stress_azimuth_uncertainty" => "SHmax Azimuth",
            "aphi_value_uncertainty" => "APhi Value"
        )
    end
    
    #println("Stress Parameter Mapping Keys: $(sort(collect(keys(stress_param_mapping))))")
    
    # Map uncertainty parameters to fault parameters
    fault_param_mapping = Dict(
        "strike_angles_uncertainty" => "Strike of fault",
        "dip_angles_uncertainty" => "Dip of fault",
        "friction_coefficient_uncertainty" => "Friction Coeff"
    )
    
    # Combine all parameters
    parameter_mapping = merge(stress_param_mapping, fault_param_mapping)
    
    #println("Combined Parameter Mapping Keys: $(sort(collect(keys(parameter_mapping))))")
    
    # Convert the specific fault to a dictionary for easier manipulation
    fault_dict = Dict(name => base_fault_inputs[fault_idx, name] for name in names(base_fault_inputs))
    #friction_coefficient = fault_dict["FrictionCoefficient"]
    
    # Get fault ID
    fault_id = string(fault_dict["FaultID"])
    
    # MODIFIED: Check if absolute stress values are provided
    # Use them for baseline calculations if available
    if haskey(base_stress_inputs, "absolute_vertical_stress") && 
       haskey(base_stress_inputs, "absolute_min_horizontal_stress") && 
       haskey(base_stress_inputs, "absolute_max_horizontal_stress") && 
       haskey(base_stress_inputs, "absolute_pore_pressure")
       
        # Use the absolute stress values provided
        baseline_stress_state = StressState(
            [base_stress_inputs["absolute_vertical_stress"], 
             base_stress_inputs["absolute_min_horizontal_stress"], 
             base_stress_inputs["absolute_max_horizontal_stress"]],
            base_stress_inputs["max_stress_azimuth"]
        )
        
        # Initial pore pressure from inputs (absolute value)
        initial_pressure = base_stress_inputs["absolute_pore_pressure"]
        
        
    else
        # Fallback to calculating absolute stresses from gradients
        # This is the original behavior
        println("Absolute stresses not provided, calculating from gradients")
        
        # Calculate absolute stresses from gradients for baseline
        reference_depth = base_stress_inputs["reference_depth"]
        friction = fault_dict["FrictionCoefficient"]
        stress_state, pore = calculate_absolute_stresses(base_stress_inputs, friction, stress_model_type)
        
        baseline_stress_state = stress_state
        initial_pressure = pore
    end
    
    sig_normal, tau_normal, s11, s22, s33, s12, n1, n2 = calculate_fault_effective_stresses(
        fault_dict["Strike"], 
        fault_dict["Dip"], 
        baseline_stress_state, 
        initial_pressure, 
        0.0  
    )
    
    # Calculate baseline critical pore pressure (we use that to compate the results from calculations
    # with lower and upper bounds)
    baseline_pp_to_slip = ComputeCriticalPorePressureForFailure(
        sig_normal,
        tau_normal,
        friction_coefficient,
        initial_pressure,
        1.0,  # biot coefficient
        0.5,  # Poisson's ratio
        1.0   # dp
    )
    
    #println("Baseline PP to Slip: $baseline_pp_to_slip")
    
    # Create DataFrame to store results
    fault_sensitivity_df = DataFrame(
        parameter = String[],
        display_name = String[],
        percent_deviation = Float64[],
        baseline = Float64[],
        lower_bound = Float64[],
        upper_bound = Float64[],
        delta_pressure = Float64[]
    )
    
    # For each parameter, calculate slip pressure at lower and upper bounds
    for (uncertainty_param, display_name) in parameter_mapping
        #println("\nProcessing parameter: $uncertainty_param => $display_name")
        
        if !haskey(uncertainties, uncertainty_param) || isnothing(uncertainties[uncertainty_param]) || uncertainties[uncertainty_param] == 0
            #println("  SKIPPING: Uncertainty not defined or is 0")
            continue  # Skip if uncertainty is not defined or is 0
        end
        
        uncertainty_value = uncertainties[uncertainty_param]
        #println("  Uncertainty value: $uncertainty_value")
        
        # Process separately for stress and fault parameters
        if haskey(stress_param_mapping, uncertainty_param)
            #println("  This is a STRESS parameter")
            # Modify the affected parameters directly without using calculate_absolute_stresses
            try
                #println("  Calculating with negative change (-$uncertainty_value)")
                lower_pp_to_slip = calculate_with_direct_parameter(
                    base_stress_inputs, 
                    fault_dict, 
                    uncertainty_param, 
                        -uncertainty_value,
                        stress_model_type
                )
                
                    #println("  Calculating with positive change (+$uncertainty_value)")
                upper_pp_to_slip = calculate_with_direct_parameter(
                    base_stress_inputs, 
                    fault_dict, 
                    uncertainty_param, 
                        uncertainty_value,
                        stress_model_type
                    )
                
                #println("  Results: Lower PP: $lower_pp_to_slip, Upper PP: $upper_pp_to_slip")
            catch e
                #println("  ERROR calculating slip pressure: $e")
                continue
            end
        else
            #println("  This is a FAULT parameter")
            # Fault parameter
            fault_param = fault_param_mapping[uncertainty_param]
            if fault_param == "Strike of fault"
                param_key = "Strike"
            elseif fault_param == "Dip of fault"
                param_key = "Dip"
            else
                param_key = "FrictionCoefficient"
            end
            
            # Calculate with modified fault parameters
            lower_fault_dict = copy(fault_dict)
            lower_fault_dict[param_key] = fault_dict[param_key] - uncertainty_value
            
            upper_fault_dict = copy(fault_dict)
            upper_fault_dict[param_key] = fault_dict[param_key] + uncertainty_value
            
            #println("  Original value: $(fault_dict[param_key]), Lower: $(lower_fault_dict[param_key]), Upper: $(upper_fault_dict[param_key])")
            
            # Calculate slip pressure with lower bound
            sig_normal_lower, tau_normal_lower, _, _, _, _, _, _ = calculate_fault_effective_stresses(
                lower_fault_dict["Strike"], 
                lower_fault_dict["Dip"], 
                baseline_stress_state, 
                initial_pressure, 
                0.0
            )
            
            lower_pp_to_slip = ComputeCriticalPorePressureForFailure(
                sig_normal_lower,
                tau_normal_lower,
                lower_fault_dict["FrictionCoefficient"],
                initial_pressure,
                1.0,
                0.5,
                1.0
            )
            
            # Calculate slip pressure with upper bound
            sig_normal_upper, tau_normal_upper, _, _, _, _, _, _ = calculate_fault_effective_stresses(
                upper_fault_dict["Strike"], 
                upper_fault_dict["Dip"], 
                baseline_stress_state, 
                initial_pressure, 
                0.0
            )
            
            upper_pp_to_slip = ComputeCriticalPorePressureForFailure(
                sig_normal_upper,
                tau_normal_upper,
                upper_fault_dict["FrictionCoefficient"],
                initial_pressure,
                1.0,
                0.5,
                1.0
            )
            
            #println("  Results: Lower PP: $lower_pp_to_slip, Upper PP: $upper_pp_to_slip")
        end
        
        # Calculate percent deviation
        lower_deviation = ((lower_pp_to_slip - baseline_pp_to_slip) / baseline_pp_to_slip) * 100
        
        upper_deviation = ((upper_pp_to_slip - baseline_pp_to_slip) / baseline_pp_to_slip) * 100
        
        
        # find the max absolute deviation (we use that for sorting, since we want a tornado-like graph)
        max_deviation = max(abs(lower_deviation), abs(upper_deviation))
        
        # Calculate absolute pressure change
        delta_pressure = max(abs(lower_pp_to_slip - baseline_pp_to_slip), abs(upper_pp_to_slip - baseline_pp_to_slip))
        
        #println("  Deviations: Lower: $lower_deviation%, Upper: $upper_deviation%, Max: $max_deviation%")
        #println("  Delta pressure: $delta_pressure")
        
        # Add to results
        push!(fault_sensitivity_df, (
            uncertainty_param,
            display_name,
            max_deviation,
            baseline_pp_to_slip,
            lower_pp_to_slip,
            upper_pp_to_slip,
            delta_pressure
        ))
    end
    
    # Sort by absolute pressure change (delta_pressure) in descending order
    sort!(fault_sensitivity_df, :delta_pressure, rev=true)
    
    #println("\nFinal fault_sensitivity_df ($(nrow(fault_sensitivity_df)) rows):")
    #for row in eachrow(fault_sensitivity_df)
        #println("  $(row.display_name): delta_pressure=$(row.delta_pressure)")
    #end
    
    # Create the final DataFrame format suitable for D3.js visualization
    tornado_df = DataFrame(
        label = String[],
        min = Float64[],
        max = Float64[],
        id = String[]
    )
    
    for row in eachrow(fault_sensitivity_df)
        # Calculate absolute pressure changes
        pressure_min = row.lower_bound - row.baseline
        pressure_max = row.upper_bound - row.baseline
        
        # Ensure the min value is always less than max value
        if pressure_min > pressure_max
            pressure_min, pressure_max = pressure_max, pressure_min
        end
        
        push!(tornado_df, (
            row.display_name,  # label
            round(pressure_min, digits=2),      # min (pressure change at lower bound)
            round(pressure_max, digits=2),      # max (pressure change at upper bound)
            fault_id           # id (fault ID)
        ))
    end
    
    #println("\nFinal tornado_df ($(nrow(tornado_df)) rows):")
    #for row in eachrow(tornado_df)
        #println("  $(row.label): min=$(row.min), max=$(row.max), id=$(row.id)")
    #end
    
    #println("Created tornado chart data for fault ID: $fault_id with $(nrow(tornado_df)) parameters")
    #println("====== DEBUG: Ending fault_sensitivity_tornado_chart_to_d3 ======\n")
    
    return tornado_df
end

"""
Modified helper function to calculate slip pressure with modified parameters directly
instead of using the calculate_absolute_stresses function
"""
function calculate_with_direct_parameter(
    base_stress_inputs::Dict, 
    fault_dict::Dict, 
    uncertainty_param::String, 
    param_change::Float64,
    stress_model_type::String = "gradients"
)
    #println("\n---- DEBUG: calculate_with_direct_parameter ----")
    #println("Parameter: $uncertainty_param with change: $param_change")
    #println("Stress model type: $stress_model_type")
    
    # Create a copy of the inputs
    modified_stress = copy(base_stress_inputs)
    # stress_model_type is now passed as a parameter, not read from the dictionary
    
    # Get reference depth for calculations
    reference_depth = modified_stress["reference_depth"]
    #println("Reference depth: $reference_depth")
    
    # Apply changes to the gradient values (not absolute values)
    # We're working with gradient inputs from probabilistic_geomechanics_process.jl
    if uncertainty_param == "vertical_stress_gradient_uncertainty"
        modified_stress["vertical_stress"] += param_change
        #println("Changed vertical_stress gradient to: $(modified_stress["vertical_stress"])")
    elseif uncertainty_param == "initial_pore_pressure_gradient_uncertainty"
        modified_stress["pore_pressure"] += param_change
        #println("Changed pore_pressure gradient to: $(modified_stress["pore_pressure"])")
    elseif uncertainty_param == "max_stress_azimuth_uncertainty"
        modified_stress["max_stress_azimuth"] += param_change
        #println("Changed max_stress_azimuth to: $(modified_stress["max_stress_azimuth"])")
    elseif uncertainty_param == "max_horizontal_stress_uncertainty"
        modified_stress["max_horizontal_stress"] += param_change
        #println("Changed max_horizontal_stress gradient to: $(modified_stress["max_horizontal_stress"])")
    elseif uncertainty_param == "min_horizontal_stress_uncertainty"
        modified_stress["min_horizontal_stress"] += param_change
        #println("Changed min_horizontal_stress gradient to: $(modified_stress["min_horizontal_stress"])")
    elseif uncertainty_param == "aphi_value_uncertainty"
        # If we're varying aphi, we need to check if it's available in the inputs
        if haskey(modified_stress, "aphi_value")
            modified_stress["aphi_value"] += param_change
            #println("Changed aphi_value to: $(modified_stress["aphi_value"])")
        else
            println("WARNING: aphi_value not found in inputs, no change applied")
        end
    else
        println("WARNING: Unknown parameter: $uncertainty_param")
    end

    # Get the friction coefficient from the fault
    friction_coefficient = fault_dict["FrictionCoefficient"]
    #println("Friction coefficient: $friction_coefficient")
    
    # Calculate absolute stress values based on the changed parameters
    # This ensures we properly convert from gradients to absolute values
    #println("Calculating absolute stresses with calculate_absolute_stresses")
    stress_state_obj, initial_pressure = calculate_absolute_stresses(modified_stress, friction_coefficient, stress_model_type)
    
    # Extract the calculated absolute stresses
    vertical_stress = stress_state_obj.principal_stresses[1]
    min_horizontal_stress = stress_state_obj.principal_stresses[2]
    max_horizontal_stress = stress_state_obj.principal_stresses[3]
    pore_pressure = initial_pressure
    
    #println("Final absolute stresses calculated from gradients:")
    #println("  vertical_stress: $vertical_stress")
    #println("  min_horizontal_stress: $min_horizontal_stress")
    #println("  max_horizontal_stress: $max_horizontal_stress")
    #println("  pore_pressure: $pore_pressure")
    
    # Create StressState object with calculated values
    stress_state_obj = StressState(
        [vertical_stress, 
         min_horizontal_stress, 
         max_horizontal_stress],
        modified_stress["max_stress_azimuth"]
    )
    
    # Calculate fault stresses
    sig_normal, tau_normal, _, _, _, _, _, _ = calculate_fault_effective_stresses(
        fault_dict["Strike"], 
        fault_dict["Dip"], 
        stress_state_obj, 
        pore_pressure, 
        0.0  # dp is 0 for initial calculation
    )
    
    # Calculate critical pore pressure
    pp_to_slip = ComputeCriticalPorePressureForFailure(
        sig_normal,
        tau_normal,
        friction_coefficient,
        pore_pressure,
        1.0,  # biot coefficient
        0.5,  # Poisson's ratio
        1.0   # dp
    )
    
    #println("Calculated pp_to_slip: $pp_to_slip")
    return pp_to_slip
end




"""
Generates the data for the shifted Mohr diagram
For the arcs, it shifts the circles to the left by the mean dp value
For the slip line, it doesn't need to be shifted
For the faults, it adds the dp value to the existing pore pressure
If those dataframes are not provided, it will calculate from scratch using the existing code
"""

function mohr_diagram_hydro_data_to_d3_portal(
    sh::Float64, 
    sH::Float64, 
    sV::Float64, 
    tau_effective::Union{Float64, AbstractVector{Float64}}, 
    sigma_effective::Union{Float64, AbstractVector{Float64}}, 
    p0::Float64, 
    biot::Float64, 
    nu::Float64, 
    dp::Union{AbstractVector{Float64}, AbstractVector{Integer}}, 
    strike::Vector{Float64},
    mu::Union{Float64, Integer},
    stress_regime::String,
    slip_pressure::Union{Float64, AbstractVector{Float64}},
    fault_ids::Union{Vector{String}, Vector{Int}, Vector{Any}}=["fault"],
    geo_arcs_df::DataFrame=DataFrame(),
    geo_faults_df::DataFrame=DataFrame(),
    geo_slip_df::DataFrame=DataFrame(),
    fault_inputs::DataFrame=DataFrame()
)
    N = length(strike)

    # Ensure dp is a vector and calculate the mean for Mohr circles
    dp = dp isa AbstractVector ? dp : fill(dp, N)
    dp_mean = mean(dp)
    
    # Check if we have valid geomechanics dataframes to use
    has_geo_data = !isempty(geo_arcs_df) && !isempty(geo_faults_df) && !isempty(geo_slip_df)
    
    if has_geo_data
        # --- Modify Arcs DataFrame with pressure shift ---
        arcsDF = copy(geo_arcs_df)
        
        # Shift the circles to the left by dp_mean using string indexing
        arcsDF[!, "centerX"] .-= dp_mean
        arcsDF[!, "labelPosX"] .-= dp_mean
        
        # --- Keep Slip Line DataFrame unchanged ---
        slipDF = copy(geo_slip_df)
        
        # --- Create new Faults DataFrame with same structure ---
        faultsDF = DataFrame(
            [name => similar(geo_faults_df[!, name], 0) for name in names(geo_faults_df)]
        )
        
        # For each fault, use the pre-calculated effective stresses
        for i in 1:min(size(fault_inputs, 1), length(dp))
            # Get fault properties
            fault_id = fault_inputs[i, "FaultID"]
            
            # Use the precalculated effective stresses passed to the function
            sig_fault = sigma_effective[i]
            tau_fault = tau_effective[i]
            fault_slip_pressure = slip_pressure[i]
            
            println("Using calculated stresses for fault: $(fault_id)")
            println("  sig_fault: $(sig_fault), tau_fault: $(tau_fault)")
            
            # Create a new row for this fault with updated coordinates
            new_row = Dict(name => geo_faults_df[i, name] for name in names(geo_faults_df))
            
            # Update the x and y coordinates with the stresses
            new_row["x"] = sig_fault
            new_row["y"] = tau_fault
            new_row["pore_pressure_slip"] = fault_slip_pressure

            # Add the row to the output dataframe
            push!(faultsDF, new_row)
        end
        
        # Return dataframes with consistent structure
        #=
        println("MODIFIED DATAFRAMES FOR SHIFTED MOHR DIAGRAM:")
        println("arcsDF: ")
        pretty_table(arcsDF)
        println("faultsDF: ")
        pretty_table(faultsDF)
        =#
        
        return (arcsDF, slipDF, faultsDF)
    else
        # Calculate from scratch (existing logic)
        # ...
        error("No geomechanics dataframes provided for the hydrology Mohr diagram, cannot calculate from scratch")
        # Return all three DataFrames
        println("Created hydrology Mohr diagram data from scratch with mean dp = $(dp_mean) psi")
        return (arcsDF, slipDF, faultDF)
    end
end

"""
Returns a DataFrame with data for Histograms for input parameter distributions.
param_values: Dictionary mapping fault index to dictionary of parameter name => vector of samples
ex. fault1 => Dict("strike" => [10, 20, 30...], "dip" => [30, 40, 50...])
nbins: number of histogram bins (we will keep 25)
"""
# TO DO: also pass the stress parameters
# strike, dip, and slip pressure are the only ones that are different for each fault
function input_distribution_histograms_to_d3(
    fault_param_values::Dict{String, Dict{String, Vector{Float64}}},
    stress_values::Dict{String, Vector{Float64}},
    stress_field_model::String;
    nbins::Int=25
)

   

    histogram_d3_data = DataFrame(id=String[], subgraph=String[], count=Int[], bar_index=String[])
    # label mapping for all parameters
    # Added some additional labels to catch all possible mappings
    label_map = Dict(
        # Lowercase keys 
        "strike" => "Strike",
        "dip" => "Dip",
        "slip_pressure" => "Pore Pressure to Slip (PSI)",
        "vertical_stress_gradient_uncertainty" => "Sv (PSI/ft)",
        "initial_pore_pressure_gradient_uncertainty" => "Initial PP Grad (PSI/ft)",
        "max_stress_azimuth_uncertainty" => "SHmax Azimuth (deg)",
        "max_horizontal_stress_uncertainty" => "SHmax Grad (PSI/ft)",
        "min_horizontal_stress_uncertainty" => "Shmin Grad (PSI/ft)",
        "aphi_value_uncertainty" => "A-Phi Param",
        
        # Capitalized keys 
        "Strike" => "Strike",
        "Dip" => "Dip",
        "FrictionCoefficient" => "Friction Coefficient",
        "SlipPressure" => "Pore Pressure to Slip (PSI)",
        
        # Actual stress parameter names from run_monte_carlo
        "vertical_stress" => "Sv (PSI/ft)",
        "pore_pressure" => "Initial PP Grad (PSI/ft)",
        "max_stress_azimuth" => "SHmax Azimuth (deg)",
        "max_horizontal_stress" => "SHmax Grad (PSI/ft)",
        "min_horizontal_stress" => "Shmin Grad (PSI/ft)",
        "aphi_value" => "A-Phi Parameter",
        
        # Hydrology parameters
        "aquifer_thickness" => "Aquifer Thickness (ft)",
        "porosity" => "Porosity",
        "permeability" => "Permeability (mD)",
        "fluid_density" => "Fluid Density (kg/m³)",
        "dynamic_viscosity" => "Dynamic Viscosity (Pa·s)",
        "fluid_compressibility" => "Fluid Compressibility (Pa⁻¹)",
        "rock_compressibility" => "Rock Compressibility (Pa⁻¹)"
    )

    # For each fault, histogram its inputs and then the selected stress uncertainties
    for (fid, pdict) in fault_param_values
        # fault-specific parameters
        for (pname, vals) in pdict
            h = fit(Histogram, vals, nbins=nbins)
            edges = h.edges[1]; centers = [(edges[i]+edges[i+1])/2 for i in 1:length(h.weights)]
            for i in eachindex(h.weights)
                #println("pname: $pname, centers[i]: $(centers[i])")
                push!(histogram_d3_data, (
                    id=fid,
                    subgraph=label_map[pname],
                    count=Int(h.weights[i]),
                    bar_index=string(centers[i])
                ))
            end
        end
        
        # this shouldn't be empty but we should check
        if !isempty(stress_values)
            
            for (pname, vals) in stress_values
                if haskey(label_map, pname)
                    h = fit(Histogram, vals, nbins=nbins)
                    edges = h.edges[1]; centers = [(edges[i]+edges[i+1])/2 for i in 1:length(h.weights)]
                    for i in eachindex(h.weights)
                        push!(histogram_d3_data, (
                            id=fid,
                            subgraph=label_map[pname],
                            count=Int(h.weights[i]),
                            bar_index=string(centers[i])
                        ))
                    end
                else
                    error("Warning: No label mapping found for stress parameter '$pname'")
                end
            end
        else
            error("No stress values available - skipping stress uncertainty histograms")
        end
    end

    return histogram_d3_data
end


function hydro_input_distribution_histograms_to_d3(
    hydro_param_values::Dict{String, Vector{Float64}};
    nbins::Int=25
)

    histogram_d3_data = DataFrame(subgraph=String[], count=Int[], bar_index=String[])

    # label mapping for all hydrology parameters
    label_map = Dict(
        "aquifer_thickness" => "Injection Formation Thickness (ft)",
        "porosity" => "Porosity (%)",
        "permeability" => "Permeability (mD)",
        "fluid_density" => "Fluid Density (kg/m³)",
        "dynamic_viscosity" => "Dynamic Viscosity (Pa.s)",
        "fluid_compressibility" => "Fluid Compressibility (Pa⁻¹)",
        "rock_compressibility" => "Rock Compressibility (Pa⁻¹)"
    )

    # in the hydrology case, these histograms are the same for all faults
    # so we don't need to iterate over faults
    if !isempty(hydro_param_values)
        for (pname, vals) in hydro_param_values
            if haskey(label_map, pname)
                h = fit(Histogram, vals, nbins=nbins)
                edges = h.edges[1]; centers = [(edges[i]+edges[i+1])/2 for i in 1:length(h.weights)]
                for i in eachindex(h.weights)
                    push!(histogram_d3_data, (
                        subgraph=label_map[pname],
                        count=Int(h.weights[i]),
                        bar_index=string(centers[i])
                    ))
                end
            else
                error("Warning: No label mapping found for hydrology parameter '$pname'")
            end
        end
    else
        error("No randomly sampled hydrology parameters were provided for the histograms")
    end

    return histogram_d3_data
end


end # module
module DetHydrologyDriver

using JSON
using ArgParse
using LinearAlgebra
using MeshGrid
using PrettyTables
using Dates
using Plots
using Profile
using ProfileView




#using Dates # might use later for dynamic dates


include("core/Hydrology_model.jl")
include("core/hydrology_calculations.jl")
include("core/well_model.jl")
include("graphs/julia_fsp_graphs.jl")
include("core/bill_pfront.jl")
include("step2_deterministic_geomechanics.jl")
include("core/geomechanics_model.jl")


using .GeomechanicsDriver
using .WellModel
using .HydrologyModel
using .HydroCalculations
using .JuliaFSPGraphs
using .BillPFront
using .GeomechanicsModel
using .GeomechanicsModel: StressState


export prepare_well_data_for_pressureScenario, prepare_constant_well_data_for_pressureScenario, prepare_monthly_well_data_for_pressureScenario


# Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input-json"
            help = "Path to input JSON file from step3"
            required = true
        "--output-graph"
            help = "Path to output JSON file for D3.js graph data"
            required = false
        "--output-json"
            help = "Path to output JSON file"
            required = false
    end
    return parse_args(s)
end




# Converts a year (integer format) to a MATLAB datenum (float64 format) 1st of January
function year_to_datenum(year::Int)::Float64
    return Dates.value(Date(year, 1, 1))
end


#=
# creates the days and bpds arrays for a constant rate well (used in the pressureScenario_constant_rate function) 
function prepare_constant_well_data_for_pressureScenario(inj_per_day, start_year, end_year)
    # Validate input parameters
    @assert start_year < end_year "Start year must be less than end year"
    @assert inj_per_day > 0 "Injection rate must be greater than 0"

    # Define start and end dates
    start_date = Date(start_year, 1, 1)
    end_date = Date(end_year, 1, 1)

    # Generate the range of daily timestamps
    date_range = start_date:(Day(1)):(end_date - Day(1)) # Exclude the end_date

    # Preallocate arrays for days and rates
    num_days = length(date_range)
    days = Vector{Float64}(undef, num_days)
    #rates = fill(Float64(well["injection_rate"]["rate_bbl_day"]), num_days)
    rates = fill(inj_per_day, num_days)

    # Populate the days array (days since start_date)
    for (i, date) in enumerate(date_range)
        days[i] = (date - start_date).value + 1  # Convert to days since start_date
    end

    return days, rates
end
=#

function well_id_to_string(well_id::Union{String, Integer})
    return string(well_id)
end



# crates the days and bpds arrays for a well (used in the pressureScenario function)
# check for constant or monthly injection rates
function prepare_well_data_for_pressureScenario(well, start_year, end_year, extrapolate=false)
    # checks for constant or monthly injection rates (defaults to constant if the field is missing)
    rate_mode = get(well["injection_rate"], "rate_mode", "constant") 

    if rate_mode == "constant"
        return prepare_constant_well_data_for_pressureScenario(well, start_year, end_year)
    elseif rate_mode == "monthly"
        if extrapolate
            return prepare_monthly_well_data_for_pressureScenario(well, start_year, end_year, true)
        else
            return prepare_monthly_well_data_for_pressureScenario(well, start_year, end_year)
        end
    else
        error("Invalid rate mode: $rate_mode. Must be 'constant' or 'monthly'")
    end
end



"""
    prepare_monthly_well_data_for_pressureScenario(well, start_year, end_year)

For a "monthly_bbl" field, each key = year, subkey=month => total bbl.
We produce step changes:
  - Step up at the start of each month to (total_monthly / days_in_that_month)
  - Step down at the end if the next month differs or we reach end_year
"""
function prepare_monthly_well_data_for_pressureScenario(
    well::Dict, 
    start_year::Int, 
    end_year::Int,
    extrapolate::Bool=false
    )
    monthly_data = get(well["injection_rate"], "monthly_bbl", Dict()) 
    @assert start_year <= end_year "Start year must be <= end year"

    # If start_year == end_year => no valid interval
    if start_year == end_year
        return (Float64[], Float64[])
    end

    # Global date boundaries
    global_start_date = Date(start_year, 1, 1)
    global_end_date   = Date(end_year,   1, 1)

    # Arrays for step changes
    step_times = Float64[]
    step_rates = Float64[]
    current_rate = 0.0

    # function to retrieve total monthly volume or 0 if missing
    function monthly_volume(y, m)
        y_str, m_str = string(y), string(m)
        if haskey(monthly_data, y_str) && haskey(monthly_data[y_str], m_str)
            return monthly_data[y_str][m_str]  # total bbl that month
        elseif extrapolate
            # extrapolate from last known month
            if !isempty(monthly_data)
                last_year = maximum(parse(Int, k) for k in keys(monthly_data))
                last_month = maximum(parse(Int, k) for k in keys(monthly_data[string(last_year)]))
                return monthly_data[string(last_year)][string(last_month)]

            else 
                # no data to extrapolate from
                return 0.0
            end
        else
            return 0.0
        end
    end

    # Helper to get the first day of the next month
    function next_month(d::Date)
        y, m = year(d), month(d)
        if m < 12
            return Date(y, m+1, 1)
        else
            return Date(y+1, 1, 1)
        end
    end

    current_month_date = global_start_date

    while current_month_date < global_end_date
        y, m = year(current_month_date), month(current_month_date)

        # next month's 1st
        next_m = next_month(current_month_date)

        # clamp if next_m > global_end_date
        if next_m >= global_end_date
            next_m = global_end_date
        end

        # how many days in [current_month_date, next_m)
        ndays = (next_m - current_month_date).value

        # total barrels for that month from JSON
        vol_bbl = monthly_volume(y, m)
        # daily rate (if ndays>0)
        month_daily_rate = ndays > 0 ? (vol_bbl / ndays) : 0.0

        # If rate changed from current_rate, step up
        if month_daily_rate != current_rate
            push!(step_times, (current_month_date - global_start_date).value + 1)
            push!(step_rates, month_daily_rate)
            current_rate = month_daily_rate
        end

        # Check what's next. If next month's rate is different, or we've reached the end, we step down
        if next_m == global_end_date
            # final boundary => step down if still injecting
            if current_rate > 0
                push!(step_times, (global_end_date - global_start_date).value + 1)
                push!(step_rates, 0.0)
                current_rate = 0.0
            end
        else
            # see if next month rate differs
            next_y, next_mo = year(next_m), month(next_m)
            next_vol_bbl = monthly_volume(next_y, next_mo)
            # days in that next month segment
            next_ndays = (next_month(next_m) - next_m).value
            next_rate = next_ndays > 0 ? (next_vol_bbl / next_ndays) : 0.0

            if next_rate != current_rate && current_rate > 0
                # step down at next_m
                push!(step_times, (next_m - global_start_date).value + 1)
                push!(step_rates, 0.0)
                current_rate = 0.0
            end
        end

        # move to next month
        current_month_date = next_m
    end

    return step_times, step_rates
end



"""
    prepare_constant_well_data_for_pressureScenario(well, start_year, end_year)

Creates a minimal 2-step array:
 - Step up from 0 -> rate_bbl_day at day=1
 - Step down to 0 at day=(#days_in_interval + 1)

Thus, you don't create a full daily array of repeated rates.
"""
function prepare_constant_well_data_for_pressureScenario(
    well::Dict, 
    start_year::Int, 
    end_year::Int
    )
    @assert start_year <= end_year "Start year must be <= end year"

    # 1) If there's no actual injection time (start_year == end_year),
    #    return empty. Or handle partial day logic if you prefer.
    if start_year == end_year
        return (Float64[], Float64[])
    end

    # Extract the constant rate
    rate_bbl_day = well["injection_rate"]["rate_bbl_day"]
    @assert rate_bbl_day > 0 "Injection rate must be > 0"

    # 2) Define the start and end dates
    start_date = Date(start_year, 1, 1)
    end_date   = Date(end_year,   1, 1)

    # number_of_days = difference in days
    # If e.g. 2015..2016 => there's exactly 365 or 366 days (leap year).
    ndays = (end_date - start_date).value

    if ndays <= 0
        # no valid interval
        return (Float64[], Float64[])
    end

    # 3) Build a 2-step array:
    #    1) At day=1: rate = rate_bbl_day
    #    2) At day=ndays+1: rate=0 => we shut off
    # Note: This matches the superposition approach: injection from day=1..ndays

    days = Float64[]
    rates = Float64[]

    # step up at day=1
    push!(days, 1.0)
    push!(rates, rate_bbl_day)

    # step down at day=ndays+1
    push!(days, ndays + 1.0)
    push!(rates, 0.0)

    return days, rates
end


# creates a 2D grid of points over the field surface (in km)
function create_spatial_grid_km(xmin, xmax, ymin, ymax, number_points)
    x_range_km = range(xmin, stop=xmax, length=number_points)
    y_range_km = range(ymin, stop=ymax, length=number_points)

    Xgrid_km, Ygrid_km = meshgrid(x_range_km, y_range_km)
    
    return Xgrid_km, Ygrid_km, x_range_km, y_range_km
end

function prompt_for_well_highlight(injection_wells::Dict{String,Any})
    all_ids = collect(keys(injection_wells))
    println("Available wells: ", join(all_ids, ", "))
    println("Type 'all' for all wells, 'none' for none, or a single well ID.")
    print("Which wells do you want to highlight? > ")
    user_in = readline()
    user_in_list = strip.(split(lowercase(user_in), ","))
    
    return user_in_list
end


function main()
    println("\n=== Starting Hydrology Analysis ===")

    args = parse_commandline()
    input_data = JSON.parsefile(args["input-json"])

    hydrology_data = input_data["hydrology"]
    model_parameters = input_data["model_parameters"]
    injection_wells = input_data["injection_wells"]



    # plot fault surface map (testing)
    plot_fault_surface_map(input_data)

    




    # 2) User input for year of interest and extrapolate option
    println("Enter the year of interest (YYYY): ")
    year_of_interest = parse(Int, readline())

    println("Do you want to extrapolate monthly rates if any are missing? (y/n): ")
    extrapolate_rates = lowercase(strip(readline())) == "y"

    #=
    println("Plotting injection rate line chart...")
    if extrapolate_rates
        println("Extrapolating missing monthly rates...")
        plot_injection_rate_line_chart(input_data, true)
    else
        plot_injection_rate_line_chart(input_data)
    end
    =#


    # 3) Calculate storativity and transmissivity
    h_feet      = hydrology_data["aquifer_thickness"]
    phi         = hydrology_data["porosity"]
    kap_md      = hydrology_data["permeability_md"]
    rho         = model_parameters["fluid_density"]
    mu          = model_parameters["dynamic_viscosity"]
    beta        = model_parameters["fluid_compressibility"]
    alphav      = model_parameters["rock_compressibility"]
    g           = 9.81


    S, T = calcST(
        h_feet, phi, kap_md, rho,
        mu, 9.81, beta, alphav
    )
    STRho = (S, T, rho) # tuple for pressureScenario function
    println("Storativity: $S")
    println("Transmissivity: $T")

    

    # 4) Create the 2D grid for the pressure field (in km)
    local xmin, xmax = 2, 18
    local ymin, ymax = 2, 18
    local number_points = 50
    Xgrid_km, Ygrid_km, x_range_km, y_range_km = create_spatial_grid_km(
        xmin, xmax, ymin, ymax, number_points
    )
    @assert size(Xgrid_km) == (50, 50)

    local total_pressure_2d = zeros(size(Xgrid_km))

    # store well information for plotting
    local wells_for_plot = Vector{Dict{String, Any}}()


    # 5) Superposition pressure of all wells at the year of interest
    for well_id in keys(injection_wells)
        well_data = injection_wells[well_id]
        inj_start_year = well_data["injection_rate"]["start_year"]
        inj_end_year = well_data["injection_rate"]["end_year"]
        well_x_km = well_data["location"]["easting_km"]
        well_y_km = well_data["location"]["northing_km"]

        # check if well is injecting at year of interest
        if year_of_interest < inj_start_year
            println("Well $well_id is not injecting at year $year_of_interest")
            continue
        end
        local actual_end = min(inj_end_year, year_of_interest)
        if actual_end < inj_start_year
            println("Well $well_id is not injecting at year $year_of_interest")
            continue
        end

        # prepare days and injection rate arrays for pressure front calculation
        
        days, bpds = prepare_well_data_for_pressureScenario(
            well_data, inj_start_year, actual_end
        )
        if isempty(days)
            println("No valid injection data for well $well_id")
            continue
        end

        if well_data["injection_rate"]["rate_mode"] == "constant"
            println("Well $well_id has constant injection rates")
            # Calculate pressure field for the well
            """
            pfield_this_well = pfieldcalc_constant_rates(
                Xgrid_km, Ygrid_km, STRho, days, bpds,
                well_x_km, well_y_km, Date(year_of_interest, 1, 1)
            )
            """
            pfield_this_well = pfieldcalc_all_rates(
                Xgrid_km, Ygrid_km, STRho, days, bpds,
                well_x_km, well_y_km, Date(year_of_interest, 1, 1)
            )
            total_pressure_2d .+= pfield_this_well
        else
            println("Well $well_id has monthly injection rates")
            # Calculate pressure field for the well
            pfield_this_well = pfieldcalc_all_rates(
                Xgrid_km, Ygrid_km, STRho, days, bpds,
                well_x_km, well_y_km, Date(year_of_interest, 1, 1)
            )
            total_pressure_2d .+= pfield_this_well
        
        end

        

        # Store well data for plotting in D3.js
        push!(wells_for_plot, Dict(
            "well_id" => well_id,
            "x_easting_km" => well_x_km,
            "y_northing_km" => well_y_km
        ))

    end


    # plot the heatmap
    plot_pressure_grid_heatmap(collect(x_range_km), collect(y_range_km), total_pressure_2d, wells_for_plot)

    # 6) Compute radial curves 
    local radial_info = Vector{Tuple{String, Vector{Float64}, Vector{Float64}}}()
    local r_km = range(0.5, stop=20.0, length=50)
    local r_m = r_km .* 1000


    # plot the radial curves for all wells
    for well_id in keys(injection_wells)
        well_data = injection_wells[well_id]
        inj_start_year = well_data["injection_rate"]["start_year"]
        inj_end_year = well_data["injection_rate"]["end_year"]

        if year_of_interest < inj_start_year
            println("Well $well_id is not injecting at year $year_of_interest")
            continue
        end
        local actual_end = min(inj_end_year, year_of_interest)
        if actual_end <= inj_start_year
            println("Well $well_id is not injecting at year $year_of_interest")
            continue
        end

        days, bpds = prepare_well_data_for_pressureScenario(
            well_data, inj_start_year, actual_end
        )
        if isempty(days)
            println("No valid injection data for well $well_id")
            continue
        end

        local pressure_psi = pressureScenario_Rall(bpds, days, r_m, STRho)
        push!(radial_info, (well_id, collect(r_km), collect(pressure_psi)))
    end

    plot_pressure_distance_graph(radial_info)


    # Prepare graph data JSON
    local graph_data = Dict(
        "heatmap" => Dict(
            "x_range_km" => collect(x_range_km),
            "y_range_km" => collect(y_range_km),
            "pressure_grid" => total_pressure_2d,
            "wells_data" => wells_for_plot
        ),
        "radial" => Dict("lines" => radial_info)
    )

    # Save graph data JSON
    if haskey(args, "output-graph") && !isnothing(args["output-graph"])
        JSON.open(args["output-graph"], "w") do io
            JSON.print(io, graph_data, 2)
        end
        println("Saved graph data to ", args["output-graph"])
    else
        println("No graph_data output directory specified, saving to current directory...")
        # save to current directory
        JSON.open("graph_data.json", "w") do io
            JSON.print(io, graph_data, 2)
        end
    end

    # Run the hydrology modeling on faults
    # prepare fault data
    strikes = Float64[fault["strike"] for fault in input_data["faults"]]
    dips = Float64[fault["dip"] for fault in input_data["faults"]]
    SHdir = Float64(input_data["stress_state"]["max_stress_azimuth"])
    # friction coefficient (this will be a scalar since we are using the same one for all faults)
    mu = Float64(input_data["faults"][1]["friction_coefficient"])
    number_of_faults = length(input_data["faults"])

    # column vector to store pore pressure on each fault
    ppOnFault = zeros(number_of_faults)

    # loop over faults
    println("\nCalculating pore pressure on faults...")
    for i in 1:number_of_faults
        # loop over wells
        for well_id in keys(input_data["injection_wells"])
            # well data structure
            well_data = injection_wells[well_id]
            # get injection dates from well
            inj_start_year = well_data["injection_rate"]["start_year"]
            inj_end_year = well_data["injection_rate"]["end_year"]
            # fault coordinates
            x_fault = input_data["faults"][i]["easting"]
            y_fault = input_data["faults"][i]["northing"]

            # well coordinates
            well_x_km = well_data["location"]["easting_km"]
            well_y_km = well_data["location"]["northing_km"]

            if year_of_interest < inj_start_year
                println("Well $well_id is not injecting at year $year_of_interest")
                continue
            end
            # actual_end might be affected by the year of interest the user selected
            local actual_end = min(inj_end_year, year_of_interest)
            if actual_end <= inj_start_year
                println("Well $well_id is not injecting at year $year_of_interest")
                continue
            end

            # create the timestamps for this well
            days, bpds = prepare_well_data_for_pressureScenario(
                well_data, inj_start_year, actual_end
            )

            if isempty(days)
                println("No valid injection data for well $well_id")
                continue
            end


            
            ppOnFault[i] += pfieldcalc_all_rates(x_fault, y_fault, STRho, days, bpds, well_x_km, well_y_km, Date(year_of_interest, 1, 1))


        end
    end

    fault_data = input_data["faults"]
    stress_data = input_data["stress_state"]
    reference_depth = input_data["stress_state"]["reference_depth"]
    pp_gradient = input_data["stress_state"]["pore_pressure"]
    stress_state, pp0 = calculate_absolute_stresses(stress_data, fault_data)
    nu = 0.5 # poisson's ratio
    biot = 1.0 # biot's coefficient
    mu = fault_data[1]["friction_coefficient"]

    strikes = [fault["strike"] for fault in fault_data]

    hydro_results = process_faults(fault_data, stress_state, pp0; tab="det_hydro", dp=ppOnFault)


    # extract tau_effective and sigma_n from the results
    tau_effective_faults = [result["normal_stress"] for result in hydro_results]
    sigma_n_faults = [result["shear_stress"] for result in hydro_results]

    # plot the shifted Mohr diagram for testing
    plot_mohr_diagram_hydro(stress_state.principal_stresses[2], stress_state.principal_stresses[3], stress_state.principal_stresses[1], tau_effective_faults, sigma_n_faults, pp0, 1.0, 0.5, ppOnFault, strikes, mu)


    # save ppOnFault values to a JSON file
    open("output/hydro_ppOnFault.json", "w") do file
        JSON.print(file, ppOnFault)
    end




    println("Pore presssure on faults results: ", hydro_results)

    println("\n=== Hydrology Analysis Complete ===")

end
   

   

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


end # module
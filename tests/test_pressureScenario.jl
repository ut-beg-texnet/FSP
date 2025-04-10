using LinearAlgebra
using Plots
using CSV
using DataFrames
using Plots
using MeshGrid
using Dates
using PrettyTables

include("../core/bill_pfront.jl")
include("../core/hydrology_calculations.jl")
include("../graphs/julia_fsp_graphs.jl")

using .BillPFront
using .HydroCalculations
using .JuliaFSPGraphs

function prepare_constant_well_data_for_pressureScenario(inj_per_day, start_year, end_year)
    # Validate input parameters
    @assert start_year < end_year "Start year must be less than end year"

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


function plot_fsp_pressure_grid_heatmap(filepath::String)
    # load CSV to a DataFrame
    df = CSV.read(filepath, DataFrame)

    # Remove leading/trailing whitespaces from column names
    rename!(df, Symbol.(strip.(String.(names(df)))))
    # extract unique X and Y coordinates
    x = unique(df.X_Easting_km)
    y = unique(df.Y_Northing_km)

    # check if product of unique X and Y coordinates equals the data length
    if length(x) * length(y) != size(df, 1)
        error("Data cannot be reshaped due to irregularities in the grid structure")
    end

    # Reshape additional pressure data
    pressure_grid = reshape(df.additionalPressure_PSI, length(y), length(x))
    # create heatmap
    try
        heatmap(x, y, pressure_grid, xlabel="X Easting (km)", ylabel="Y Northing (km)", title="Pressure Field FSP", color=:viridis)
        savefig("pressure_field_fsp.png")
        println("Heatmap for FSP CSV created successfully.")
    catch e
        println("Error creating heatmap for FSP CSV: $e")
    end
end


# creates a 2D grid of points over the field surface (in meters)
function create_spatial_grid_km(xmin, xmax, ymin, ymax, number_points)
    x_range_km = range(xmin, stop=xmax, length=number_points)
    y_range_km = range(ymin, stop=ymax, length=number_points)

    Xgrid_km, Ygrid_km = meshgrid(x_range_km, y_range_km)
    
    return Xgrid_km, Ygrid_km, x_range_km, y_range_km
end


function main()

    # Hardcoded variables
    density = 1000.0          # rho
    dynamic_viscosity = 0.0008  # mu
    fluid_compressibility = 3.6e-10 # beta
    rock_compressibility = 1.08e-09 # alpha
    aquifer_thickness = 100.0 # h
    porosity = 0.1            # phi
    permeability = 200.0      # kap

    # Well data
    inj_start_year = 2015
    inj_end_year = 2016
    inj_per_day = 27000.0
    well_x_km = 10.0
    well_y_km = 10.0


    # ---------------------- HEATMAP ---------------------------------------

    # 1) Create the grid
    xmin, xmax = 2, 18
    ymin, ymax = 2, 18
    number_points = 50
    Xgrid_km, Ygrid_km, x_range_km, y_range_km = create_spatial_grid_km(
        xmin, xmax, ymin, ymax, number_points
    )

    # 2) Calculate S and T
    S, T = calcST(
        aquifer_thickness, porosity, permeability, density,
        dynamic_viscosity, 9.81, fluid_compressibility, rock_compressibility
    )
    STRho = (S, T, density)

    # 3) Get the day/rate arrays
    days, bpds = prepare_constant_well_data_for_pressureScenario(
        inj_per_day, inj_start_year, inj_end_year
    )

    # 4) Compute the 2D pressure field
    pressure_field_2D = pfieldcalc(
        Xgrid_km, Ygrid_km, STRho, days, bpds,
        well_x_km, well_y_km, Date(2024, 1, 1)
    )

    # 5) Plot the 2D heatmap
    plot_pressure_grid_heatmap(x_range_km, y_range_km, pressure_field_2D)

    # ---------------------- PRESSURE-DISTANCE GRAPH ------------------------

    # radial distances for pressure-distance graph
    r_km = range(0.5, stop=20, length=50)

    # radial pressure calculation
    radial_pressures = pressureScenario_constant_rate(bpds, days, r_km * 1000, STRho) # in meters
    plot_pressure_distance_graph(r_km, radial_pressures)


    # Plot the pressure field from FSP CSV
    #plot_fsp_pressure_grid_heatmap("hydro_grid_2016_1_well_2_years.csv")
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

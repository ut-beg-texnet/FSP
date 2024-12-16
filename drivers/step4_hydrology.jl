using JSON
using ArgParse
using LinearAlgebra
using MeshGrid
#using Dates # might use later for dynamic dates


include("core/Hydrology_model.jl")
include("core/hydrology_calculations.jl")
using .HydrologyModel
using .HydroCalculations


# Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input-json"
            help = "Path to input JSON file from step3"
            required = true
        "--output-json"
            help = "Path to output JSON file"
            required = false
    end
    return parse_args(s)
end

# Create datenum_barrels_per_day for each well (follows MATLAB format)
function generate_datenum_barrels(injection_wells, year_range)
    datenum_barrels_per_day = []
    for well_id in keys(injection_wells)
        well = injection_wells[well_id]
        start_year = well["injection_rate"]["start_year"]
        end_year = well["injection_rate"]["end_year"]
        rate_per_year = well["injection_rate"]["rate_bbl_day"]

        # Populate rates across year_range
        for year in year_range
            rate = (year >= start_year && year <= end_year) ? rate_per_year : 0.0
            push!(datenum_barrels_per_day, [year, rate])
        end
    end
    return hcat(datenum_barrels_per_day...)  # Combine into a matrix
end





function main()
    println("\n=== Starting Hydrology Analysis ===")

    args = parse_commandline()

    # Load input JSON file
    input_data = JSON.parsefile(args["input-json"])

    # Extract hydrology data from input JSON
    hydrology_data = input_data["hydrology"]
    porosity = hydrology_data["porosity"]
    aquifer_thickness = hydrology_data["aquifer_thickness"]
    fluid_density = hydrology_data["fluid_density"]
    fluid_compressibility = hydrology_data["fluid_compressibility"]
    dynamic_viscosity = hydrology_data["dynamic_viscosity"]
    rock_compressibility = hydrology_data["rock_compressibility"]
    permeability = hydrology_data["permeability_md"] # in mD (will be converted to m^2 in calcST.jl)

    # get the number of wells
    injection_wells = input_data["injection_wells"]
    
    #println("Type of injection wells: ", typeof(injection_wells))
    number_wells = length(keys(injection_wells))

    # number of points for the grid
    number_points = 50



    # need to implement check for internal calculation (based on whether the user entered an external model)

    # hardcode xmin, xmax, ymin, ymax for now
    xmin = 0
    xmax = 100
    ymin = 0
    ymax = 100

    # create a 2D grid of points over the field surface
    x = range(xmin, stop=xmax, length=number_points)
    y = range(ymin, stop=ymax, length=number_points)
    X, Y = meshgrid(x, y)

    # define radial distances for pressure curves (in km)
    # a linearly spaced vector from 0.5 to 20 with 'number_points' points
    r = range(0.5, 20, length=number_points)


    # this is read from the slider in matlab (year)
    # here we dynamically extract the minimum and maximum year from the input JSON file
    #year_range = range(1900, 3000, step=1)
    #println(typeof(ts))
    min_start_year = Inf
    max_end_year = -Inf
    for well in values(injection_wells)
        start_year = well["injection_rate"]["start_year"]
        end_year = well["injection_rate"]["end_year"]
        if start_year < min_start_year
            min_start_year = start_year
        end
        if end_year > max_end_year
            max_end_year = end_year
        end
    end
    year_range = range(min_start_year, max_end_year, step=1)
    
    #println("Year range: ", year_range)

    # bounds to the year range
    if length(year_range) > 100
        println("Year range is too large. Setting to 100 years.")
        year_range = range(min_start_year, stop=min_start_year + 100, step=1)
    end
    if length(year_range) == 0
        error("FSP expects at least one start and end year for hydrology.")
    end

    # read the years and injection rates from JSON file in this format: datenum_barrels_per_day = [2018 1000.0; 2019 1500.0; 2020 2000.0]
    datenum_barrels_per_day = generate_datenum_barrels(injection_wells, year_range)
    # make this a matrix of Float64
    datenum_barrels_per_day = convert(Matrix{Float64}, datenum_barrels_per_day)
    #println("datenum_barrels_per_day: ", datenum_barrels_per_day)
    #println(typeof(datenum_barrels_per_day))

    # pdfieldcalc calculation
    pField = HydroCalculations.pdfieldcalc(injection_wells, X, Y, number_wells, aquifer_thickness, porosity, permeability, fluid_density, dynamic_viscosity, fluid_compressibility, rock_compressibility, year_range, datenum_barrels_per_day)
    println("pField min: ", minimum(pField))
    println("pField max: ", maximum(pField))


end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end




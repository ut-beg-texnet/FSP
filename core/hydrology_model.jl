"""
Module for handling hydrology-related calculations and validations
"""
module HydrologyModel

using DataFrames
using ArgParse

export validate_hydrology_inputs, create_hydrology_data, HydrologyParams, parse_hydrology_args



"""
Structure to hold hydrology parameters
"""
struct HydrologyParams
    porosity::Union{Float64, Int64}              # in %
    permeability::Union{Float64, Int64}      # [mD]
    aquifer_thickness::Union{Float64, Int64}   # [ft]
    density::Union{Float64, Int64}           # [kg/m³]
    dynamic_viscosity::Union{Float64, Int64}  # [Pa.s]
    fluid_compressibility::Union{Float64, Int64}  # [1/Pa]
    rock_compressibility::Union{Float64, Int64}   # [1/Pa]
end

"""
Validate hydrology parameters
"""
function validate_hydrology_inputs(params::Dict)
    # Check porosity is between 0 and 1
    if params["porosity"] <= 0 || params["porosity"] >= 1
        error("Porosity must be between 0 and 1")
    end
    
    # Check permeability is positive
    if params["permeability"] <= 0
        error("Permeability must be positive")
    end
    
    # Check aquifer thickness is positive
    if params["aquifer_thickness"] <= 0
        error("Aquifer thickness must be positive")
    end
    
    # Check fluid properties are positive
    if params["fluid_density"] <= 0
        error("Fluid density must be positive")
    end
    
    if params["dynamic_viscosity"] <= 0
        error("Dynamic viscosity must be positive")
    end
    
    if params["fluid_compressibility"] <= 0
        error("Fluid compressibility must be positive")
    end
    
    if params["rock_compressibility"] <= 0
        error("Rock compressibility must be positive")
    end
end

"""
Create hydrology data from CLI arguments
"""
function create_hydrology_data(args::Dict)
    # Validate hydrology inputs
    validate_hydrology_inputs(args)
    
    # Convert permeability from mD to m²
    permeability = md_to_m2(args["permeability"])
    
    # Convert aquifer thickness from ft to m
    aquifer_thickness = args["aquifer_thickness"] * 0.3048
    
    # Create DataFrame with single row of hydrology data
    df = DataFrame(
        porosity = [args["porosity"]],
        permeability_md = [args["permeability"]],
        aquifer_thickness = [aquifer_thickness],
        fluid_density = [args["fluid_density"]],
        dynamic_viscosity = [args["dynamic_viscosity"]],
        fluid_compressibility = [args["fluid_compressibility"]],
        rock_compressibility = [args["rock_compressibility"]]
    )
    
    return df
end

"""
Parse command line arguments
"""
function parse_hydrology_args(args)
    # Define CLI arguments
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--porosity"
            arg_type = Float64
            required = true
        "--permeability"
            arg_type = Float64
            required = true
        "--aquifer-thickness"
            arg_type = Float64
            required = true
        "--fluid-density"
            arg_type = Float64
            required = true
        "--dynamic-viscosity"
            arg_type = Float64
            required = true
        "--fluid-compressibility"
            arg_type = Float64
            required = true
        "--rock-compressibility"
            arg_type = Float64
            required = true
    end
    
    # Parse CLI arguments
    parsed_args = parse_args(args, s)
    
    return parsed_args
end

"""
Main function
"""
function main()
    args = parse_hydrology_args(ARGS)
    hydro_data = create_hydrology_data(args)
    params = HydrologyParams(
        porosity = hydro_data.porosity[1],
        permeability = md_to_m2(hydro_data.permeability_md[1]),
        aquifer_thickness = hydro_data.aquifer_thickness[1],
        density = hydro_data.fluid_density[1],
        dynamic_viscosity = hydro_data.dynamic_viscosity[1],
        fluid_compressibility = hydro_data.fluid_compressibility[1],
        rock_compressibility = hydro_data.rock_compressibility[1]
    )
    return params
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end # module

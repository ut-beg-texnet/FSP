"""
Model Inputs Step

This script processes input data for the FSP analysis:
1. Reads fault data from CSV
2. Creates stress state from CLI arguments
3. Creates hydrology parameters from CLI arguments
4. Processes injection well data from CSV
5. Gets uncertainty values from CLI
6. Outputs all data to JSON for subsequent steps
"""

using ArgParse
using JSON
using CSV
using DataFrames

# Include the modules
include("core/stress_model.jl")
include("core/fault_model.jl")
include("core/hydrology_model.jl")

# Import the modules
using .StressModel
using .FaultModel
using .HydrologyModel

"""
Parse command line arguments for the Model Inputs step
"""
function parse_commandline()
    s = ArgParseSettings(description="FSP 3.0 Model Inputs Step")
    
    # Input/output files
    @add_arg_table! s begin
        "--fault-data"
            help = "Path to fault data CSV file"
            required = true
        "--injection-well-data"
            help = "Path to injection well data CSV file"
            required = true
    end
    
    # Stress state parameters
    @add_arg_table! s begin
        "--vertical-stress"
            help = "Vertical stress gradient [psi/ft]"
            arg_type = Float64
            required = true
        "--max-stress-azimuth"
            help = "Maximum horizontal stress azimuth [degrees]"
            arg_type = Float64
            required = true
        "--pore-pressure"
            help = "Pore pressure gradient [psi/ft]"
            arg_type = Float64
            required = true
        "--reference-depth"
            help = "Reference depth [ft]"
            arg_type = Float64
            required = true
    end

    # Hydrology parameters
    @add_arg_table! s begin
        "--porosity"
            help = "Porosity [fraction]"
            arg_type = Float64
            required = true
        "--permeability"
            help = "Permeability [mD]"
            arg_type = Float64
            required = true
        "--aquifer-thickness"
            help = "Aquifer thickness [ft]"
            arg_type = Float64
            required = true
        "--fluid-density"
            help = "Fluid density [kg/m³]"
            arg_type = Float64
            default = 1000.0
        "--dynamic-viscosity"
            help = "Dynamic viscosity [Pa·s]"
            arg_type = Float64
            default = 0.0008
        "--fluid-compressibility"
            help = "Fluid compressibility [1/Pa]"
            arg_type = Float64
            default = 3.6e-10
        "--rock-compressibility"
            help = "Rock compressibility [1/Pa]"
            arg_type = Float64
            default = 1.08e-09
    end

    # Model-specific parameters
    @add_arg_table! s begin
        "--model-type"
            help = "Stress model type (gradients, aphi_min, aphi_no_min)"
            required = true
        "--max-horizontal-stress"
            help = "Maximum horizontal stress gradient [psi/ft]"
            arg_type = Float64
        "--min-horizontal-stress"
            help = "Minimum horizontal stress gradient [psi/ft]"
            arg_type = Float64
        "--aphi-value"
            help = "A-phi value for stress model"
            arg_type = Float64
    end

    # Uncertainty parameters
    @add_arg_table! s begin
        "--unc-vertical-stress"
            help = "Uncertainty in vertical stress gradient [psi/ft]"
            arg_type = Float64
            required = false
        "--unc-pore-pressure"
            help = "Uncertainty in pore pressure gradient [psi/ft]"
            arg_type = Float64
            required = false
        "--unc-strike"
            help = "Uncertainty in strike angles [degrees]"
            arg_type = Float64
            required = false
        "--unc-dip"
            help = "Uncertainty in dip angles [degrees]"
            arg_type = Float64
            required = false
        "--unc-max-stress-azimuth"
            help = "Uncertainty in maximum stress azimuth [degrees]"
            arg_type = Float64
            required = false
        "--unc-friction"
            help = "Uncertainty in friction coefficient"
            arg_type = Float64
            required = false
        "--unc-aphi"
            help = "Uncertainty in A-phi value"
            arg_type = Float64
            required = false
    end

    return parse_args(s)
end

"""
Read input data from CSV files and CLI arguments
"""
function read_input_data(args)
    fault_data = FaultModel.read_fault_data(args["fault-data"])
    injection_data = CSV.read(args["injection-well-data"], DataFrame)
    
    # Create stress state from CLI arguments
    stress_data = create_stress_state(args)
    
    # Create hydrology data from CLI arguments
    hydro_data = Dict(
        "porosity" => args["porosity"],
        "permeability_md" => args["permeability"],
        "aquifer_thickness" => args["aquifer-thickness"],
        "fluid_density" => args["fluid-density"],
        "dynamic_viscosity" => args["dynamic-viscosity"],
        "fluid_compressibility" => args["fluid-compressibility"],
        "rock_compressibility" => args["rock-compressibility"]
    )
    
    # Get uncertainties from CLI arguments
    uncertainties = get_uncertainties(args)
    
    return fault_data, stress_data, hydro_data, injection_data, uncertainties
end

"""
Get uncertainties from command line arguments
"""
function get_uncertainties(args::Dict)
    uncertainties = Dict{String, Float64}()
    
    # Map CLI argument names to output names
    uncertainty_mapping = Dict(
        "unc-vertical-stress" => "vertical_stress_gradient",
        "unc-pore-pressure" => "initial_pore_pressure_gradient",
        "unc-strike" => "strike_angles",
        "unc-dip" => "dip_angles",
        "unc-max-stress-azimuth" => "max_stress_azimuth",
        "unc-friction" => "friction_coefficient",
        "unc-aphi" => "aphi_value"
    )
    
    # Add non-null uncertainties to output
    for (arg_name, output_name) in uncertainty_mapping
        if !isnothing(args[arg_name])
            uncertainties[output_name] = args[arg_name]
        end
    end
    
    return uncertainties
end

"""
Convert DataFrame to Dictionary for JSON serialization
"""
function df_to_dict(df::DataFrame)
    # Convert each row to a dictionary with column names as keys
    rows = [Dict(name => row[name] for name in names(df)) for row in eachrow(df)]
    return rows
end

"""
Process injection well data into a more convenient format
"""
function process_injection_data(df::DataFrame)
    # Group by well_id and create a dictionary of time series data
    wells = Dict{Int, Dict{String, Any}}()
    
    # Check if we have constant rate data (start_year, end_year) or monthly data (year, month)
    is_constant_rate = "start_year" in names(df)
    
    for well_id in unique(df.well_id)
        well_data = filter(row -> row.well_id == well_id, df)
        first_row = first(eachrow(well_data))
        
        wells[well_id] = Dict(
            "location" => Dict(
                "easting_km" => first_row.easting_km,
                "northing_km" => first_row.northing_km
            )
        )
        
        if is_constant_rate
            # For constant rate data, store rate and time period
            wells[well_id]["injection_rate"] = Dict(
                "rate_bbl_day" => first_row.injection_rate_bbl_day,
                "start_year" => first_row.start_year,
                "end_year" => first_row.end_year
            )
        else
            # For monthly data, store injection history
            wells[well_id]["injection_history"] = [
                Dict(
                    "year" => row.year,
                    "month" => row.month,
                    "volume_bbl" => row.injection_volume_bbl
                ) for row in eachrow(well_data)
            ]
        end
    end
    
    return wells
end

"""
Create stress state from CLI arguments
"""
function create_stress_state(args)
    # Validate model type and required parameters
    model_type = args["model-type"]
    if !(model_type in ["gradients", "aphi_min", "aphi_no_min"])
        error("Invalid model type. Must be one of: gradients, aphi_min, aphi_no_min")
    end
    
    # Common parameters for all models
    stress_state = Dict(
        "model_type" => model_type,
        "vertical_stress" => Float64(args["vertical-stress"]),
        "max_stress_azimuth" => Float64(args["max-stress-azimuth"]),
        "pore_pressure" => Float64(args["pore-pressure"]),
        "reference_depth" => Float64(args["reference-depth"])
    )
    
    # Add model-specific parameters
    if model_type == "gradients"
        if isnothing(args["max-horizontal-stress"]) || isnothing(args["min-horizontal-stress"])
            error("gradients model requires both max-horizontal-stress and min-horizontal-stress")
        end
        if !isnothing(args["aphi-value"])
            error("gradients model should not include aphi-value")
        end
        stress_state["max_horizontal_stress"] = Float64(args["max-horizontal-stress"])
        stress_state["min_horizontal_stress"] = Float64(args["min-horizontal-stress"])
        
    elseif model_type == "aphi_min"
        if isnothing(args["min-horizontal-stress"]) || isnothing(args["aphi-value"])
            error("aphi_min model requires both min-horizontal-stress and aphi-value")
        end
        if !isnothing(args["max-horizontal-stress"])
            error("aphi_min model should not include max-horizontal-stress")
        end
        stress_state["min_horizontal_stress"] = Float64(args["min-horizontal-stress"])
        stress_state["aphi_value"] = Float64(args["aphi-value"])
        stress_state["max_horizontal_stress"] = nothing
        
    else  # aphi_no_min
        if isnothing(args["aphi-value"])
            error("aphi_no_min model requires aphi-value")
        end
        if !isnothing(args["max-horizontal-stress"]) || !isnothing(args["min-horizontal-stress"])
            error("aphi_no_min model should not include max-horizontal-stress or min-horizontal-stress")
        end
        stress_state["aphi_value"] = Float64(args["aphi-value"])
        stress_state["max_horizontal_stress"] = nothing
        stress_state["min_horizontal_stress"] = nothing
    end
    
    # Validate stress state
    StressModel.validate_stress_state(stress_state)
    
    return stress_state
end

"""
Main function for the Model Inputs step
"""
function main()
    # Parse command line arguments
    args = parse_commandline()
    
    # Read input data files
    fault_data, stress_state, hydro_data, injection_data, uncertainties = read_input_data(args)
    
    # Create output directory if it doesn't exist
    output_dir = "output"
    mkpath(output_dir)
    
    # Prepare output data
    output_data = Dict(
        "model_parameters" => Dict(
            "fluid_density" => args["fluid-density"],
            "dynamic_viscosity" => args["dynamic-viscosity"],
            "fluid_compressibility" => args["fluid-compressibility"],
            "rock_compressibility" => args["rock-compressibility"]
        ),
        "faults" => df_to_dict(fault_data),
        "stress_state" => stress_state,
        "hydrology" => hydro_data,
        "injection_wells" => process_injection_data(injection_data)
    )
    
    # Add aphi value to model parameters if using an aphi model
    if stress_state["model_type"] in ["aphi_min", "aphi_no_min"]
        output_data["model_parameters"]["aphi_value"] = args["aphi-value"]
    end
    
    # Add uncertainties if any were provided
    if !isempty(uncertainties)
        output_data["uncertainties"] = uncertainties
    end
    
    # Write output to JSON file
    output_file = joinpath(output_dir, "model_inputs.json")
    open(output_file, "w") do f
        JSON.print(f, output_data, 4)  # 4 spaces for indentation
    end
    
    println("Model Inputs step completed. Output written to $output_file")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

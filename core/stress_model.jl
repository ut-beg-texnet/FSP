"""
Module for handling stress-related calculations and validations
"""
module StressModel

using DataFrames
using CSV

export validate_stress_inputs, read_stress_data, detect_stress_model, calculate_aphi_stresses, validate_stress_state

"""
Detect stress model type from input data
"""
function detect_stress_model(stress_data::DataFrame)
    # Check if we have all required columns for gradients model
    gradient_cols = ["vertical_stress", "max_horizontal_stress", "min_horizontal_stress", 
                    "max_stress_azimuth", "pore_pressure", "reference_depth"]
    
    has_all_cols = all(col -> col in names(stress_data), gradient_cols)
    if !has_all_cols
        error("Missing required columns for stress model")
    end
    
    # Check if max_horizontal_stress is empty
    if ismissing(stress_data.max_horizontal_stress[1])
        # If min_horizontal_stress is also empty, it's aphi_no_min model
        if ismissing(stress_data.min_horizontal_stress[1])
            return "aphi_no_min"
        else
            return "aphi_min"  # Only max_horizontal_stress is empty
        end
    else
        return "gradients"  # All stresses are provided
    end
end

"""
Read and validate stress data from CSV file
"""
function read_stress_data(file_path::String)
    # Define column types, allowing missing values
    col_types = Dict(
        :vertical_stress => Union{Float64, Missing},
        :max_horizontal_stress => Union{Float64, Missing},
        :min_horizontal_stress => Union{Float64, Missing},
        :max_stress_azimuth => Union{Float64, Missing},
        :pore_pressure => Union{Float64, Missing},
        :reference_depth => Union{Float64, Missing}
    )
    
    # Read CSV file with header (we know our example files have headers)
    stress_data = CSV.read(file_path, DataFrame, types=col_types)
    
    # Validate the data
    validate_stress_inputs(stress_data)
    
    return stress_data
end

"""
Validate stress input data
"""
function validate_stress_inputs(stress_data::DataFrame)
    # Required columns for stress data
    required_cols = ["vertical_stress", "max_horizontal_stress", "min_horizontal_stress", 
                    "max_stress_azimuth", "pore_pressure", "reference_depth"]
    
    # Check if all required columns exist
    for col in required_cols
        if !(col in names(stress_data))
            error("Missing required column: $col")
        end
    end
    
    # Validate data types and ranges
    for row in eachrow(stress_data)
        # Check if vertical stress is positive
        if !ismissing(row.vertical_stress) && row.vertical_stress <= 0
            error("vertical_stress must be positive")
        end
        
        # Check if max stress azimuth is between 0 and 360
        if !ismissing(row.max_stress_azimuth) && (row.max_stress_azimuth < 0 || row.max_stress_azimuth > 360)
            error("max_stress_azimuth must be between 0 and 360 degrees")
        end
        
        # Check if pore pressure is non-negative
        if !ismissing(row.pore_pressure) && row.pore_pressure < 0
            error("pore_pressure must be non-negative")
        end
        
        # Check if reference depth is positive
        if !ismissing(row.reference_depth) && row.reference_depth <= 0
            error("reference_depth must be positive")
        end
        
    end
end



"""
Validate stress state parameters
"""
function validate_stress_state(stress_state::Dict)
    # Check if vertical stress is positive
    if stress_state["vertical_stress"] <= 0
        error("vertical_stress must be positive")
    end
    
    # Check if max stress azimuth is between 0 and 360
    if stress_state["max_stress_azimuth"] < 0 || stress_state["max_stress_azimuth"] > 360
        error("max_stress_azimuth must be between 0 and 360 degrees")
    end
    
    # Check if pore pressure is non-negative
    if stress_state["pore_pressure"] < 0
        error("pore_pressure must be non-negative")
    end
    
    # Check if reference depth is positive
    if stress_state["reference_depth"] <= 0
        error("reference_depth must be positive")
    end
    
    
    # Validate based on stress model type
    model_type = stress_state["model_type"]::String
    
    if model_type == "gradients"
        # All stresses must be provided and valid
        if isnothing(stress_state["max_horizontal_stress"]) || stress_state["max_horizontal_stress"] <= 0
            error("max_horizontal_stress must be positive for gradients model")
        end
        if isnothing(stress_state["min_horizontal_stress"]) || stress_state["min_horizontal_stress"] <= 0
            error("min_horizontal_stress must be positive for gradients model")
        end
        
        # Check stress magnitude ordering
        if stress_state["max_horizontal_stress"] < stress_state["min_horizontal_stress"]
            error("max_horizontal_stress must be greater than min_horizontal_stress")
        end
        
        # Should not have aphi value
        if haskey(stress_state, "aphi_value")
            error("gradients model should not include aphi_value")
        end
        
    elseif model_type == "aphi_min"
        # Only min horizontal stress should be provided
        if !isnothing(stress_state["max_horizontal_stress"])
            error("max_horizontal_stress should not be provided for aphi_min model")
        end
        if isnothing(stress_state["min_horizontal_stress"]) || stress_state["min_horizontal_stress"] <= 0
            error("min_horizontal_stress must be positive for aphi_min model")
        end
        
        # Must have aphi value
        if !haskey(stress_state, "aphi_value") || isnothing(stress_state["aphi_value"])
            error("aphi_min model requires aphi_value")
        end
        if stress_state["aphi_value"] < 0 || stress_state["aphi_value"] > 3
            error("aphi_value must be between 0 and 3")
        end
        
    elseif model_type == "aphi_no_min"
        # Neither horizontal stress should be provided
        if !isnothing(stress_state["max_horizontal_stress"])
            error("max_horizontal_stress should not be provided for aphi_no_min model")
        end
        if !isnothing(stress_state["min_horizontal_stress"])
            error("min_horizontal_stress should not be provided for aphi_no_min model")
        end
        
        # Must have aphi value
        if !haskey(stress_state, "aphi_value") || isnothing(stress_state["aphi_value"])
            error("aphi_no_min model requires aphi_value")
        end
        if stress_state["aphi_value"] < 0 || stress_state["aphi_value"] > 3
            error("aphi_value must be between 0 and 3")
        end
    else
        error("Invalid stress model type: $model_type")
    end
end

end # module

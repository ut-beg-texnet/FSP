"""
Module for handling fault-related calculations and validations
"""
module FaultModel

using DataFrames
using CSV

export validate_fault_inputs, read_fault_data

"""
Read and validate fault data from CSV file
"""
function read_fault_data(file_path::String)
    # Define column types
    col_types = Dict(
        :fault_id => Int64,
        :easting => Float64,
        :northing => Float64,
        :strike => Float64,
        :dip => Float64,
        :length_km => Float64,
        :friction_coefficient => Float64
    )
    
    # Read CSV file with header (we know our example files have headers)
    fault_data = CSV.read(file_path, DataFrame, types=col_types)
    
    # Validate the data
    validate_fault_inputs(fault_data)
    
    return fault_data
end

"""
Validate fault input data
"""
function validate_fault_inputs(fault_data::DataFrame)
    # Required columns for fault data
    required_cols = ["fault_id", "easting", "northing", "strike", "dip", "length_km", "friction_coefficient"]
    
    # Check if all required columns exist
    for col in required_cols
        if !(col in names(fault_data))
            error("Missing required column: $col")
        end
    end
    
    # Validate data types and ranges
    for row in eachrow(fault_data)
        # Check if fault_id is a positive number
        if row.fault_id <= 0
            error("fault_id must be a positive number")
        end
        
        # Check if coordinates are numeric
        if !isa(row.easting, Number) || !isa(row.northing, Number)
            error("easting and northing must be numeric")
        end
        
        # Check if strike is between 0 and 360
        if !isa(row.strike, Number) || row.strike < 0 || row.strike > 360
            error("strike must be between 0 and 360 degrees")
        end
        
        # Check if dip is between 0 and 90
        if !isa(row.dip, Number) || row.dip < 0 || row.dip > 90
            error("dip must be between 0 and 90 degrees")
        end
        
        # Check if length is positive
        if !isa(row.length_km, Number) || row.length_km <= 0
            error("length_km must be positive")
        end
        
        # Check if friction coefficient is between 0 and 1
        if !isa(row.friction_coefficient, Number) || row.friction_coefficient < 0 || row.friction_coefficient > 1
            error("friction_coefficient must be between 0 and 1")
        end
    end
end

end # module

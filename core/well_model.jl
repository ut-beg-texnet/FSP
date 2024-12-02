"""
Module for handling injection well data and calculations
"""
module WellModel

using DataFrames
using CSV
using Dates

export validate_well_inputs, read_well_data, detect_well_data_format

"""
Detect the format of well data (constant rate vs monthly volumes)
"""
function detect_well_data_format(well_data::DataFrame)
    required_base_cols = ["well_id", "easting_km", "northing_km"]
    
    # Check if base columns exist
    missing_base_cols = setdiff(required_base_cols, names(well_data))
    if !isempty(missing_base_cols)
        error("Missing required base columns: $(join(missing_base_cols, ", "))")
    end
    
    # Check for monthly volume format columns
    has_monthly = all(col -> col in names(well_data), 
                     ["year", "month", "injection_volume_bbl"])
    
    # Check for constant rate format columns
    has_constant = all(col -> col in names(well_data), 
                      ["injection_rate_bbl_day", "start_year", "end_year"])
    
    if has_monthly && has_constant
        error("Data contains both monthly volume and constant rate columns. " *
              "Please provide only one format.")
    elseif has_monthly
        return "monthly"
    elseif has_constant
        return "constant"
    else
        error("Unable to determine well data format. " *
              "For constant rate: provide injection_rate_bbl_day, start_year, and end_year. " *
              "For monthly volumes: provide year, month, and injection_volume_bbl.")
    end
end

"""
Detect if a CSV file has a header row
"""
function detect_csv_header(file_path::String)
    # Read first few rows of the file
    rows = CSV.Rows(file_path, limit=2)
    first_row = first(rows)
    
    # Check if first row contains any non-numeric values (likely headers)
    header_candidates = [field for field in first_row]
    
    # If any field can't be parsed as a float, it's likely a header
    has_header = any(field -> try
                                parse(Float64, field)
                                false
                            catch
                                true
                            end,
                    header_candidates)
    
    return has_header
end

"""
Validate well model inputs from CSV data
"""
function validate_well_inputs(well_data::DataFrame)
    # Detect format type
    format_type = detect_well_data_format(well_data)
    
    # Validate well IDs are unique
    if length(unique(well_data.well_id)) != length(well_data.well_id)
        error("Well IDs must be unique")
    end
    
    # Validate coordinates
    if any(well_data.easting_km .< 0)
        error("Easting coordinates must be non-negative")
    end
    
    if any(well_data.northing_km .< 0)
        error("Northing coordinates must be non-negative")
    end
    
    # Format-specific validations
    if format_type == "constant"
        # Validate injection rate
        if any(well_data.injection_rate_bbl_day .< 0)
            error("Injection rate cannot be negative")
        end
        
        # Validate years
        if any(well_data.start_year .>= well_data.end_year)
            error("Start year must be before end year")
        end
        
        if any(well_data.start_year .< 1900) || any(well_data.end_year .> 2100)
            error("Years must be between 1900 and 2100")
        end
        
    elseif format_type == "monthly"
        # Validate monthly volumes
        if any(well_data.injection_volume_bbl .< 0)
            error("Monthly injection volume cannot be negative")
        end
        
        # Validate year and month
        if any(well_data.year .< 1900) || any(well_data.year .> 2100)
            error("Years must be between 1900 and 2100")
        end
        
        if any(well_data.month .< 1) || any(well_data.month .> 12)
            error("Months must be between 1 and 12")
        end
        
        # Check chronological order
        dates = Date.(well_data.year, well_data.month, 1)
        if !issorted(dates)
            error("Data must be in chronological order")
        end
    end
    
    return format_type
end

"""
Read well data from CSV file with header detection
"""
function read_well_data(file_path::String)
    has_header = detect_csv_header(file_path)
    
    # Read CSV file with appropriate header setting
    well_data = if has_header
        CSV.read(file_path, DataFrame)
    else
        # If no header, provide column names for all possible columns
        CSV.read(file_path, DataFrame, 
                header=["well_id", "easting_km", "northing_km",
                       # Monthly format columns
                       "year", "month", "injection_volume_bbl",
                       # Constant rate format columns
                       "injection_rate_bbl_day", "start_year", "end_year"])
    end
    
    # Validate the data and get format type
    format_type = validate_well_inputs(well_data)
    
    return well_data, format_type
end

end # module

#using JSON
using CSV
using DataFrames
#using Shapefile
using Geodesy 
#using PrettyTables
using InlineStrings


include("TexNetWebToolLauncherHelperJulia.jl")
include("graphs/julia_fsp_graphs.jl")
include("core/utilities.jl")

using .TexNetWebToolLauncherHelperJulia
using .JuliaFSPGraphs
using .Utilities

const ARGS_FILE_NAME = "args.json"
const RESULTS_FILE_NAME = "results.json"


# TO DO: read this from the 'stress_field_mode' parameter
# it will be either 'all_gradients' or 'aphi_model'
function get_stress_model_type(stress_inputs::Dict)
    if stress_inputs["aphi_value"] !== nothing
        if stress_inputs["min_horizontal_stress_gradient"] !== nothing
            return "aphi_min"
        else 
            return "aphi_no_min"
        end
    else
        return "gradients"
        
    end
end




function get_injection_dataset_path(helper::TexNetWebToolLaunchHelperJulia, step_index::Int)
    #println("DEBUG: get_injection_dataset_path called with step_index = $step_index")
    for param_name in ["injection_wells_annual", "injection_wells_monthly", "injection_tool_data"]
        #println("DEBUG: Trying to get file path for param_name = $param_name")
        filepath = get_dataset_file_path(helper, step_index, param_name)
        #println("DEBUG: filepath = $filepath, type = $(typeof(filepath))")
        if filepath !== nothing
            if param_name == "injection_wells_annual"
                injection_data_type = "annual_fsp"
                #println("DEBUG: Returning filepath = $filepath, type = $(typeof(filepath)), injection_data_type = $injection_data_type")
                return filepath, injection_data_type
            elseif param_name == "injection_wells_monthly"
                injection_data_type = "monthly_fsp"
                #println("DEBUG: Returning filepath = $filepath, type = $(typeof(filepath)), injection_data_type = $injection_data_type")
                return filepath, injection_data_type
            elseif param_name == "injection_tool_data"
                injection_data_type = "injection_tool_data"
                #println("DEBUG: Returning filepath = $filepath, type = $(typeof(filepath)), injection_data_type = $injection_data_type")
                return filepath, injection_data_type
            end
        end
    end
    
    error("No injection dataset found.")
    return nothing, nothing
end


function get_fault_dataset_path(helper::TexNetWebToolLaunchHelperJulia, step_index::Int)
    for param_name in ["faults", "FaultDataShapefile"]
        filepath = get_dataset_file_path(helper, step_index, param_name)
        if filepath !== nothing
            if param_name == "faults"
                fault_type = "fsp_native"
                return filepath, fault_type
            elseif param_name == "FaultDataShapefile"
                fault_type = "shapefile"
                return filepath, fault_type
            end
        end
    end
    error("No fault dataset found.")
    return nothing, nothing

end


function main()

    #println("\n=== Starting Deterministic Geomechanics Process ===")

    

    scratchPath = ARGS[1]

    helper = TexNetWebToolLaunchHelperJulia(scratchPath)    

    #=
    # for testing: define the path to the faults shapefile here
    faults_shapefile_path = "C:/Users/bakirtzisn/Desktop/FSP_dev_test/FSP_3/tests/shapefile_example.csv"

    # test the shapefile_to_fsp_csv function
    faults_df = shapefile_to_fsp_csv(faults_shapefile_path)
    # save to a CSV file
    CSV.write("test_shapefile_converted.csv", faults_df)
    error("stop Here")
    =#

    faults_csv_filepath, fault_csv_type = get_fault_dataset_path(helper, 1)
    if faults_csv_filepath !== nothing
        if fault_csv_type == "fsp_native"
            faults_df = CSV.read(faults_csv_filepath, DataFrame, types=Dict("FaultID" => String))
            
        elseif fault_csv_type == "shapefile"
            faults_df = shapefile_to_fsp_csv(faults_csv_filepath)
        end
    else
        error("No fault dataset provided.")
    end

    println("fault CSV type: $fault_csv_type")

    

    


    # get file path for faults dataset using the helper
    #faults_csv_filepath = get_dataset_file_path(helper, 1, "faults")
    if faults_csv_filepath !== nothing
        #faults_df = CSV.read(faults_csv_filepath, DataFrame, types=Dict("FaultID" => String))


        #get all the unique faults from the FaultID column
        unique_faults_num = unique(faults_df[!, "FaultID"])
        if length(unique_faults_num) > 10000
            add_message_with_step_index!(helper, 1, "Number of faults provided is greater than 10000. Please provide a smaller number of faults.", 2)

            error("Number of faults provided is greater than 1000. Please provide a smaller number of faults.")
            
        end

        #=
        # Ensure the column exists and is of the correct type
        # we create it if it doesn't exist
        if !("FrictionCoefficient" in names(faults_df))
            faults_df[!, "FrictionCoefficient"] = Vector{Union{Missing, Float64}}(missing, nrow(faults_df))
        end
        =#
        
        #=
        # TO DO: make this a single input parameter from the portal
        # get friction coefficient from the first fault
        mu = faults_df[1, "FrictionCoefficient"]
        # assign this friction coefficient to all faults
        faults_df[!, "FrictionCoefficient"] .= mu
        =#

        
        

        # convert lat/lon to wkt format and add to the dataframe
        latlon_to_wkt(faults_df)

        # append the new columns with no data (will be filled in later processes)
        faults_df[!, :slip_pressure] = Vector{Union{Missing, Float64}}(missing, nrow(faults_df))
        faults_df[!, :coulomb_failure_function] = Vector{Union{Missing, Float64}}(missing, nrow(faults_df))
        faults_df[!, :summary_end_year] = Vector{Union{Missing, Int64}}(missing, nrow(faults_df))
        faults_df[!, :summary_fsp] = Vector{Union{Missing, Float64}}(missing, nrow(faults_df))
        faults_df[!, :summary_pressure] = Vector{Union{Missing, Float64}}(missing, nrow(faults_df))
        faults_df[!, :prob_hydro_fsp] = Vector{Union{Missing, Float64}}(missing, nrow(faults_df))



        save_dataframe_as_parameter!(helper, 1, "faults_model_inputs_output", faults_df)

    else
        error("No faults dataset provided.")
    end

    
    
    

    

    # read injection wells CSV into a DataFrame
    #println("Extracting injection well data from the CSV at the scratch path...")
    injection_wells_csv_filepath, injection_data_type = get_injection_dataset_path(helper, 1)
    if injection_wells_csv_filepath !== nothing
        

        # if we have injection tool data format, read the 'API Number' column as a string
        # TO DO: we might have to do this for the rest of the formats
        # parsing it as a string allows us to keep any leading zeros in the API numbers. If we parse it as an integer, the leading zeros are stripped.
        if injection_data_type == "injection_tool_data"
            injection_wells_df = CSV.read(injection_wells_csv_filepath, DataFrame; types=Dict("API Number" => String))
        else
            injection_wells_df = CSV.read(injection_wells_csv_filepath, DataFrame, types=Dict("WellID" => String), pool=false)
            #=
            #println("column types:")
            for col in names(injection_wells_df)
                #println("$col: $(eltype(injection_wells_df[!, col]))")
            end
            =#
        end

        
        
        
        # Prepare standardized dataframe for d3 visualization of injection rate over time
        if injection_data_type == "annual_fsp"

            # get unique values from the 'WellID' column
            unique_well_ids = unique(injection_wells_df[!, "WellID"])
            if length(unique_well_ids) > 200
                add_message_with_step_index!(helper, 1, "Number of wells provided is greater than 200. Please provide a smaller number of wells.", 2)

                error("Number of wells provided is greater than 500. Please provide a smaller number of wells.")
            end


            #injection_rate_data = injection_rate_data_to_d3(injection_wells_df, injection_data_type)
            injection_rate_data = injection_rate_data_to_d3_bbl_day(injection_wells_df, injection_data_type)
            save_dataframe_as_parameter!(helper, 1, "injection_rate_d3_data", injection_rate_data)
            #println("DEBUG: Saved annual injection rate data for visualization")
            save_dataframe_as_parameter!(helper, 1, "injection_wells_annual_output", injection_wells_df)
        elseif injection_data_type == "monthly_fsp"

            # get unique values from the 'WellID' column
            unique_well_ids = unique(injection_wells_df[!, "WellID"])
            if length(unique_well_ids) > 200
                add_message_with_step_index!(helper, 1, "Number of wells provided is greater than 200. Please provide a smaller number of wells.", 2)

                error("Number of wells provided is greater than 200. Please provide a smaller number of wells.")
            end

            
            #injection_rate_data = injection_rate_data_to_d3(injection_wells_df, injection_data_type)
            injection_rate_data = injection_rate_data_to_d3_bbl_day(injection_wells_df, injection_data_type)
            save_dataframe_as_parameter!(helper, 1, "injection_rate_d3_data", injection_rate_data)
            #println("DEBUG: Saved monthly injection rate data for visualization")

            # print the type of each column in the injection_rate_data dataframe
            #println("Type of each column in the injection_rate_data dataframe:")
            


            save_dataframe_as_parameter!(helper, 1, "injection_wells_monthly_output", injection_wells_df)
        elseif injection_data_type == "injection_tool_data"


            # get unique values from the 'WellID' column
            unique_well_ids = unique(injection_wells_df[!, "API Number"])
            if length(unique_well_ids) > 200
                add_message_with_step_index!(helper, 1, "Number of wells provided is greater than 200. Please provide a smaller number of wells.", 2)

                error("Number of wells provided is greater than 200. Please provide a smaller number of wells.")
            end
            #injection_rate_data = injection_rate_data_to_d3(injection_wells_df, injection_data_type)
            injection_rate_data = injection_rate_data_to_d3_bbl_day(injection_wells_df, injection_data_type)
            save_dataframe_as_parameter!(helper, 1, "injection_rate_d3_data", injection_rate_data)

            

            # we also need to filter this so we can create a map layer
            # keep only the first row for each unique API Number
            injection_wells_df_filtered = unique(injection_wells_df, "API Number")
            # only keep the columns we need
            injection_wells_df_filtered = injection_wells_df_filtered[!, ["API Number", "Surface Latitude", "Surface Longitude"]]
            # rename "API Number" to "UWI"
            rename!(injection_wells_df_filtered, "API Number" => :UWI)
            # rename "Surface Latitude" to "Latitude(WGS84)"
            rename!(injection_wells_df_filtered, "Surface Latitude" => :"Latitude(WGS84)")
            # rename "Surface Longitude" to "Longitude(WGS84)"
            rename!(injection_wells_df_filtered, "Surface Longitude" => :"Longitude(WGS84)")
            
            save_dataframe_as_parameter!(helper, 1, "injection_tool_data_filtered_map_layer", injection_wells_df_filtered)

            save_dataframe_as_parameter!(helper, 1, "injection_tool_data_output", injection_wells_df)



            #println("DEBUG: Saved injection tool data for visualization")
        end
        
    else
        error("No injection wells dataset provided.")
    end

    

    #println("Generating D3.js visualization data...")
    #faults_csv, wells_csv = fault_surface_map_data_to_d3(faults_df, injection_wells_df)
    # CONTINUE FROM HERE FIX THIS (DEOSN'T READ CORRECT FILE)
    #   save_dataframe_as_parameter!(helper, 1, "fault_surface_map_faults", faults_csv)
    #   save_dataframe_as_parameter!(helper, 1, "fault_surface_map_injection_wells", wells_csv)

    # Write the modified args_data back to file.
    write_final_args_file(helper, joinpath(helper.scratch_path, ARGS_FILE_NAME))


    
    

    # explicitly set this step's success state to true
    set_success_for_step_index!(helper, 1, true)

    write_results_file(helper)

    #println("Model Inputs process: Completed successfully.")




end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
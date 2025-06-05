module TexNetWebToolLauncherHelperJulia


export get_parameter_state, loadOrigArgsData, get_parameter_value, get_dataset_file_path,
    get_dataset_contents_as_dataframe, set_parameter_value!, save_dataframe_as_parameter!,
    write_final_args_file, write_results_file, TexNetWebToolLaunchHelperJulia,
    set_success_for_step_index!, add_message_with_step_index!


using JSON3
using DataFrames
using CSV 
using FilePathsBase


const PARAMETER_TYPE_INT = 0
const PARAMETER_TYPE_FLOAT = 1
const PARAMETER_TYPE_BOOLEAN = 2
const PARAMETER_TYPE_STRING = 3
const PARAMETER_TYPE_USER_DATASET = 4
const PARAMETER_TYPE_DATE = 5
const PARAMETER_TYPE_USER_DATASET_ROW = 6
const ARGS_FILE_NAME = "args.json"
const RESULTS_FILE_NAME = "results.json"





mutable struct TexNetWebToolLaunchHelperJulia
    scratch_path::String
    orig_args_data::Dict{String, Any}
    args_data::Dict{String, Any}

    # Constructor (equivalent to Python's __init__)
    function TexNetWebToolLaunchHelperJulia(scratch_path::String)

        if !isdir(scratch_path)
            throw(ArgumentError("Scratch path directory does not exist: $scratch_path"))
        end

        args_file_path = joinpath(scratch_path, "args.json")
        if !isfile(args_file_path)
            throw(ArgumentError("args.json file does not exist in the scratch path: $args_file_path"))
        end

        #println("args.json file read: $args_file_path")
        #println("contents of args.json file: $(read(args_file_path, String))")

        #=
        json_contents = open(args_file_path) do file
            read(file, String)
        end

        orig_args_data = JSON3.read(json_contents, Dict)
        new(scratch_path, orig_args_data, orig_args_data)
        =#

        try
            json_contents = read(args_file_path, String)
            orig_args_data = JSON3.read(json_contents, Dict)

            args_data = deepcopy(orig_args_data)
            #println("args_data: $args_data")

            return new(scratch_path, orig_args_data, args_data)
        catch e
            throw(ArgumentError("Failed to initialize helper: $(e.msg)"))
        end
    end

end




# Returns the parameter state object (as a dictionary) with the given name from the step at the given index
function get_parameter_state(helper::TexNetWebToolLaunchHelperJulia, step_index::Int, param_name::String)
    if haskey(helper.args_data, "SessionState") &&
       haskey(helper.args_data["SessionState"], "StepState") &&
       step_index <= length(helper.args_data["SessionState"]["StepState"])
       
        step = helper.args_data["SessionState"]["StepState"][step_index]
        for param in step["InputParameterStates"]
            if param["ProcessStepParameterName"] == param_name
                return param
            end
        end

        for param in step["OutputParameterStates"]
            if param["ProcessStepParameterName"] == param_name
                return param
            end
        end
    end
    return nothing
end



function loadOrigArgsData(helper::TexNetWebToolLaunchHelperJulia)
    args_file_path = joinpath(helper.scratch_path, ARGS_FILE_NAME)
    json_contents = ""
    open(args_file_path) do file
        json_contents = read(file, String)
    end
    helper.orig_args_data = JSON3.read(json_contents, Dict)
end


# gets the value of the parameter with the given name from the step at the given index
function get_parameter_value(helper::TexNetWebToolLaunchHelperJulia, step_index::Int, param_name::String)
    
    param = get_parameter_state(helper, step_index, param_name)
    if param !== nothing && haskey(param, "Value")
        value = param["Value"]

        # if the value is an empty string or explicitly nothing, return nothing
        if isa(value, AbstractString) && isempty(strip(value)) || value === nothing
            #println("parameter $param_name not found in step $step_index")
            return nothing
        end

        data_type = param["DataType"]
        if data_type == PARAMETER_TYPE_INT  # 0 for integer
            return isa(value, AbstractString) ? parse(Int, value) : value
        elseif data_type == PARAMETER_TYPE_FLOAT  # 1 for float
            return isa(value, AbstractString) ? parse(Float64, value) : value
        elseif data_type == PARAMETER_TYPE_BOOLEAN  # 2 for boolean
            return string(value) == "true"
        else
            return value
        end
        
    end
    #println("parameter $param_name not found in step $step_index")
    return nothing
end


function get_dataset_file_path(helper::TexNetWebToolLaunchHelperJulia, step_index::Int, param_name::String)
    param_value = get_parameter_value(helper, step_index, param_name)
    
    
    if param_value !== nothing
        key_str = string(param_value)
        dataset_paths = helper.args_data["DatasetPaths"]
        if haskey(dataset_paths, key_str)
            path = dataset_paths[key_str]
            return isa(path, AbstractString) ? path : string(path)
        elseif haskey(dataset_paths, param_value)
            path = dataset_paths[param_value]
            return isa(path, AbstractString) ? path : string(path)
        elseif tryparse(Int, key_str) !== nothing && haskey(dataset_paths, tryparse(Int, key_str))
            path = dataset_paths[tryparse(Int, key_str)]
            return isa(path, AbstractString) ? path : string(path)
        else
            #println("DEBUG: DatasetPaths does not contain key '$key_str'. Available keys: $dataset_paths. Param_value: $param_value")
            error("DatasetPaths does not contain key '$key_str'. Available keys: $(keys(dataset_paths))")
        end
    end
    
    return nothing
end




function get_dataset_contents_as_dataframe(helper::TexNetWebToolLaunchHelperJulia, step_index::Int, param_name::String)
    file_path = get_dataset_file_path(helper, step_index, param_name)
    if file_path !== nothing
        return CSV.read(file_path, DataFrame)
    end
    return nothing
end


# sets the value of the parameter with the given name from the step at the given index
function set_parameter_value!(helper::TexNetWebToolLaunchHelperJulia, step_index::Int, param_name::String, value)
    param = get_parameter_state(helper, step_index, param_name)
    #println("get_parameter_state returned: $param")
    if param !== nothing
        if param["DataType"] == PARAMETER_TYPE_USER_DATASET_ROW
            param["Value"] = value
        else
            param["Value"] = string(value)
        end
    else
        error("Parameter $param_name not found in step $step_index")
    end
end


function save_dataframe_as_parameter!(helper::TexNetWebToolLaunchHelperJulia, step_index::Int, param_name::String, df::DataFrame)
    path = ""
    tries = 0

    while path == "" || isfile(path)
        filename = tries == 0 ? "$param_name.csv" : "$param_name-$tries.csv"
        path = joinpath(helper.scratch_path, filename)
        tries += 1
    end

    CSV.write(path, df)
    set_parameter_value!(helper, step_index, param_name, path)
end


# writes the modified arguments back to file
function write_final_args_file(helper::TexNetWebToolLaunchHelperJulia, output_path::String)
    try
        json_contents = JSON3.write(helper.args_data)
        open(output_path, "w") do file
            write(file, json_contents)
        end
    catch e
        throw(ArgumentError("Failed to write final args file: $(e.msg)"))
    end
end

# saves results to the results.json file
function write_results_file(helper::TexNetWebToolLaunchHelperJulia)
    output_path = joinpath(helper.scratch_path, "results.json")
    json_contents = JSON3.write(helper.args_data)
    open(output_path, "w") do file
        write(file, json_contents)
    end
end


# getter function
function Base.getproperty(obj::TexNetWebToolLaunchHelperJulia, sym::Symbol)
    if sym === :scratchPath
        return getfield(obj, :scratch_path)
    elseif sym === :origArgsData
        return getfield(obj, :orig_args_data)
    elseif sym === :argsData
        return getfield(obj, :args_data)
    else
        return getfield(obj, sym)
    end
end

# setter function
function Base.setproperty!(obj::TexNetWebToolLaunchHelperJulia, sym::Symbol, val)
    if sym === :argsData
        setfield!(obj, :args_data, val)
    else
        error("Cannot modify read-only property: $sym")
    end
end

# Define which fields can be accessed
function Base.propertynames(obj::TexNetWebToolLaunchHelperJulia)
    return (:scratchPath, :origArgsData, :argsData)
end

function set_success_for_step_index!(helper::TexNetWebToolLaunchHelperJulia, step_index::Int, new_success_state::Bool)
    if haskey(helper.args_data, "SessionState") &&
       haskey(helper.args_data["SessionState"], "StepState") &&
       step_index <= length(helper.args_data["SessionState"]["StepState"])
        
        step = helper.args_data["SessionState"]["StepState"][step_index]
        step["Success"] = new_success_state

        #println("DEBUG: Step $step_index success state set to $new_success_state")
    else
        error("Step index $step_index not found")
    end
end




# adds a message to the step at the given index
# message level can be 0 (info), 1 (warning), 2 (error)
function add_message_with_step_index!(helper::TexNetWebToolLaunchHelperJulia, step_index::Int, message_content::String, message_level::Int)
    if haskey(helper.args_data, "SessionState") &&
       haskey(helper.args_data["SessionState"], "StepState") &&
       step_index <= length(helper.args_data["SessionState"]["StepState"])
        
        step = helper.args_data["SessionState"]["StepState"][step_index]
        
        message = Dict{String, Any}(
            "MessageContent" => message_content,
            "MessageLevel" => message_level
        )
        
        if !haskey(step, "Messages")
            step["Messages"] = []
        end
        
        push!(step["Messages"], message)
    else
        error("Step index $step_index not found")
    end
end





end # module


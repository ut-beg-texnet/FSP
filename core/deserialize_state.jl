module DeserializeState

using Serialization

export deserialize_state


# Input: filepath (string) to the serialized state file
# Output: deserialized state object, or nothing if an error occurs
function deserialize_state(filepath::String)
    try
        state = nothing
        open(filepath, "r") do file
            state = deserialize(file)
        end
        println("State deserialized from $filepath")
        return state
    catch e
        println("Error deserializing state: ", e)
        return nothing
    end
end



end
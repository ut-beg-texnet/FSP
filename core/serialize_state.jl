module SerializeState

using Serialization

export serialize_state

# Input: state object to serialize and filepath (string) to save the serialized state
# Output: true if serialization is successful, false otherwise
function serialize_state(state, filepath::String)
    try
        open(filepath, "w") do file
            serialize(file, state)
        end
        println("State serialized to $filepath")
        return true
    catch e
        println("Error serializing state: ", e)
        return false
    end
end



end # module SerializeState
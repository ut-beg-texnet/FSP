# WellDataModule.jl
module WellModel

export WellData  # Make WellData available outside the module

# Define the WellData struct
struct WellData
    start_year::Float64
    end_year::Float64
    injection_rate::Float64
end

end # module

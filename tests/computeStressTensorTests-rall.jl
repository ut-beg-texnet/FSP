dir=@__FILE__
#mainDir="/" .* joinpath(split(dir,"/")[1:end-1])*"/"

additionalPath = "C:/Users/bakirtzisn/Desktop/FSP_dev_test/FSP_3/"
push!(LOAD_PATH, additionalPath)




include("../step2_deterministic_geomechanics.jl")
using .GeomechanicsModel
using LinearAlgebra
import JSON

# test inputs for calculate_absolute_stresses

# create test stress_data dictionary
stress_data_aphi_min = Dict(
    "vertical_stress" => 60.0,
    "max_horizontal_stress" => nothing,
    "min_horizontal_stress" => 20.0,
    "max_stress_azimuth" => 45.0,
    "pore_pressure" => 10.0,
    "reference_depth" => 11000,
    "aphi_value" => 0.6,
    "model_type" => "aphi_min"
)

#=
"reference_depth": 5000.0,
"vertical_stress": 1.0,
"min_horizontal_stress": 0.7,
"aphi_value": 0.7,
"max_horizontal_stress": null,
"model_type": "aphi_min",
"pore_pressure": 0.433,
"max_stress_azimuth": 45.0



"faults": [
        {
            "dip": 33.2,
            "northing": 19.908,
            "friction_coefficient": 0.55,
            "easting": 18.497,
            "strike": 152.5,
            "length_km": 8.6,
            "fault_id": 1.0
        }
=#

# create a fault_data dictionary
fault_data = [
    Dict(
        "fault_id" => 1.0,
        "strike" => 150.0,
        "dip" => 30.0,
        "length_km" => 8.0,
        "friction_coefficient" => 0.6,
        "easting" => 20.0,
        "northing" => 30.0
    )
]

stress_state, initial_pressure = calculate_absolute_stresses(stress_data_aphi_min, fault_data)
# print sV
println("sV: ", stress_state.principal_stresses[1])
println("sHmax: ", stress_state.principal_stresses[3])
println("sHmin: ", stress_state.principal_stresses[2])






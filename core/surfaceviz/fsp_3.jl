using JSON

# Assuming you have already imported the necessary module, GetHorFromAPhi
include("mohrs_circle.jl")  # Make sure this points to your mohrs_3D.jl file
include("get_hor_from_APhi.jl")  # Include the A-Phi module as well

# Mock data for testing
Sig0 = [8000.0, 5000.0, 3000.0]  # Initial total stresses [Sv, Shmin, SHmax]
ignoreParameter4 = 0              # Unused parameter (legacy)
p0 = 2000.0                       # Reference pore pressure
strike = [30.0, 60.0, 90.0]       # Strikes of faults in degrees
dip = [45.0, 50.0, 55.0]          # Dips of faults in degrees
SHdir = 90.0                      # Maximum horizontal stress direction
dp = [100.0, 200.0, 150.0]        # Pressure perturbations
mu = [0.6]                        # Friction coefficient for faults
biot = 1.0                        # Biot coefficient
nu = 0.25                         # Poisson's ratio
APhi = 1.5                        # A-Phi parameter (for modified A-Phi model)

# Prepare the input data for mohrs_3D function
indatacell = [Sig0, ignoreParameter4, p0, strike, dip, SHdir, dp, mu, biot, nu, APhi]

# Mock hDV data structure (simplified)
hDV = Dict(
    :plotdata => Dict(),
    :data => Dict(
        :stress => Dict(:aphi => Dict(:use => 12))  # A-Phi usage indicator
    )
)

# Call the mohrs_3D function with mock data
failout, outs, C1, C2, C3, sig_fault, tau_fault = mohrs_3D(indatacell, hDV)

# Prepare data for JSON output
mohr_circle_data = Dict(
    "mohr_circles" => Dict(
        "R1" => C1,
        "R2" => C2,
        "R3" => C3
    ),
    "fault_points" => [
        Dict("sig_fault" => sig, "tau_fault" => tau) for (sig, tau) in zip(sig_fault, tau_fault)
    ],
    "failure_envelope" => Dict("friction_coefficient" => mu[1]),
    "labels" => Dict(
        "sigma_v" => Sig0[1],
        "sigma_H" => Sig0[3],
        "sigma_h" => Sig0[2]
    )
)

# Convert data to JSON format
json_output = JSON.json(mohr_circle_data)

# Print the resulting JSON to verify
println(json_output)

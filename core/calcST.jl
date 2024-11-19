module CalcST

include("hydro_model_structs.jl")
using .HydroModelStructs

export calcST


# This module calculates storativity (S) and transmissivity (T) from the given hydro parameters
# Inputs: h: reservoir thickness in feet
#         phi: porosity fraction (given in percentage ex. 10 = 10%)
#         kap: intrinsic permeability in millidarcy
#         params::HydroModelStructs.HydroParameters: Struct containing other hydro parameters.

# Outputs: S: storativity
#          T: transmissivity
function calcST(h::Float64, phi::Float64, kap::Float64, params::HydroModelStructs.HydroParameters)
    # convert permeability (mD -> m^2)
    kap_m2 = kap * 10^(-3) * 9.9e-13 # manually entered data (-perm)

    rho = params.rho# fluid density
    g = params.g # accl gravity
    beta = params.beta # fluid compressibility
    alphav = params.alphav # vertical compressibility of aquifer
    mu = params.mu # dynamic viscosity

    # calculate storativity
    h_meters = h * 0.3048 # convert feet to meters
    S = rho * g * h_meters * (alphav + (phi * 1e-2) * beta)
    
    # calculate transmissivity
    K = kap_m2 * rho * g / mu # Saturated hydraulic conductivity
    T = K * h_meters

    return S, T
    
end



# HARDCODED DATA FOR TESTING ---------------------------------------------------------------------------------------


function test_calcST()
    # mock HydroParameters struct
    params = HydroParameters(
        100.0,  # aquifer_thickness
        10.0,   # porosity_percent
        150.0,  # permeability
        1000.0, # rho
        9.81,   # g
        0.001,  # mu
        4.5e-10,# beta
        1.0e-9  # alphav
    )

    
    h = 100.0  
    phi = 10.0 
    kap = 150.0 

    S, T = calcST(h, phi, kap, params)
    println("Storativity (S): $S")
    println("Transmissivity (T): $T")
end

end # module


using .CalcST
CalcST.test_calcST()

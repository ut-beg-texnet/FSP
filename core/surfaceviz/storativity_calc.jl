# calc_st.jl
# calculates Storativity (S) and Transmissivity (T) for a reservoir
# based on reservoir thickness, porosity, permeability, and other aquifer parameters.

module CalcST

export calcST

# Function to calculate storativity and transmissivity
function calcST(valsarr, h, phi, kap)
    # Inputs:
    #   h      : reservoir thickness (in meters)
    #   phi    : porosity fraction (dimensionless)
    #   kap    : intrinsic permeability (in millidarcies)
    #   valsarr: array of additional parameters
    
    # Extract parameters from valsarr
    rho = valsarr[5]  # fluid density (kg/m^3)
    g = valsarr[6]    # gravitational acceleration (m/s^2)
    mu = valsarr[7]   # dynamic viscosity (Pa.s)
    beta = valsarr[8] # fluid compressibility (Pa^-1)
    alphav = valsarr[9] # vertical compressibility of aquifer (Pa^-1)
    
    # Convert permeability from millidarcies (mD) to square meters (m²)
    kap_m2 = kap * 1e-3 * 9.9e-13  # mD to m² conversion

    # Calculate Storativity (S) using the provided formula
    S = rho * g * h * (alphav + phi * beta)  # Storativity (dimensionless)

    # Specific storage (optional, for comparison with MODFLOW)
    # SpecificStorage = S / h  # Specific storage = storativity / thickness

    # Calculate saturated hydraulic conductivity (K)
    K = kap_m2 * rho * g / mu  # Saturated hydraulic conductivity (m/s)

    # Calculate Transmissivity (T) based on hydraulic conductivity and thickness
    T = K * h  # Transmissivity (m²/s)

    return S, T
end

end  # module CalcST

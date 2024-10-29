# SetupData.jl (setupdata.m in MATLAB)

# Clears and sets up data structure. 
# This is everything that is saved when you save a model

module SetupData

using ..SurfaceViz  

export setupdata, get_fault_data_from_json, set_fault_data



function setupdata(self::SurfaceVizStruct)
    self.data = Dict{Symbol, Any}()

    self.data[:nwells] = 1

    # stress data
    self.data[:stress] = Dict(
        :txt => [
            "Vertical Stress Gradient [psi/ft]",
            "Min Horiz. Stress Gradient [psi/ft]",
            "Max Horiz. Stress Gradient [psi/ft]",
            "Max Hor Stress Direction [deg N CW]",
            "Reference Depth for Calculations [ft]",
            "Initial Res. Pressure Gradient [psi/ft]"
        ],
        :vals => zeros(6),
        :aphi => Dict(
            :use => 0,
            :vals => [0, 0],
            :sigvals => [0.0, 0.0]
        )
    )

    # reservoir data
    # Reservoir data setup
    self.data[:reservoir] = Dict(
        :txt => [
            "Aquifer Thickness [ft]",
            "Porosity [%]",
            "Permeability [mD]",
            "Stochastic Mag [psi]",
            "Poisson's Ratio"
        ],
        :vals => [0.0, 0.0, 0.0, 0.0, 0.5],  # Hardcode Poisson's ratio as 0.5 for now
        :loadedHeaderLines => 1,
        :importHydrology => 0  # 0 = internal hydrology calculation, 1 = external data import
    )

    # injection data
    self.data[:inject] = Dict(
        :txt => ["x [km]", "y [km]", "Inj. Rate [bbl/day]", "Start Year [yr]", "End Year [yr]"],
        :vals => zeros(100, 5)  # Placeholder, max wells = 100
    )

    # well data
    self.data[:realWellData] = Dict(
        :inputStringColumnNames => [
            "UniqueID/Name", "Easting (km)", "Northing (km)",
            "Year", "Month (1-12)", "InjectionVolume (bbl/month)"
        ],
        :stringsWellDataAdvanced => fill("", 2, 6),
        :columnIsNumber => [false, true, true, true, true, true],
        :loadedHeaderLines => 1
    )

    # fault data
    self.data[:fault] = Dict(
        :txt => [
            "Number of faults (max 500)", 
            "Friction Coefficient mu", 
            "Strike minimum [deg]", 
            "Strike maximum [deg]", 
            "Dip minimum [deg]", 
            "Dip maximum [deg]"
        ],
        :vals => [20.0, 0.6, 0.0, 0.0, 0.0, 0.0]
    )

    
    # advanced data
    self.data[:adv] = Dict(
        :txt => [
            "Min x [km]", "Max x [km]", "Min y [km]", "Max y [km]",
            "Density [kg/m^3]", "Accl Gravity [m/s^2]",
            "Dynamic Viscosity [Pa.s]", "Fluid Compressibility [Pa^-1]",
            "Rock Compressibility [Pa^-1]", "Set Random Seed?"
        ],
        :vals => [0.0, 20.0, 0.0, 20.0, 1000.0, 9.80665, 0.8e-3, 3.6e-10, 3*3.6e-10, 0.0]
    )

    return self


end

# gets the data that will update self.data[:fault]
function get_fault_data_from_json(file_path::String)
    # Load fault data from a JSON file
    fault_data = JSON.parsefile(file_path)

    # get fault IDs
    num_of_ids = Set([fault["id"] for fault in fault_data])
    num_of_faults = length(num_of_ids)

    # extract fault properties
    friction_coeff_mu = [fault["friction"] for fault in faults]
    strikes = [fault["strike"] for fault in faults]
    dips = [fault["dip"] for fault in faults]

    # not sure if this should be implemented (was used for random fault generation)
    strike_min = minimum(strikes)
    strike_max = maximum(strikes)
    dip_min = minimum(dips)
    dip_max = maximum(dips)

    # tuples are returned
    return (num_of_faults, friction_coeff_mu, strike_min, strike_max, dip_min, dip_max)
end

function set_fault_data(self::SurfaceVizStruct, num_faults, friction_coeff, strike_min, strike_max, dip_min, dip_max)
    self.data[:fault][:vals][1] = num_faults
    self.data[:fault][:vals][2] = friction_coeff
    self.data[:fault][:vals][3] = strike_min
    self.data[:fault][:vals][4] = strike_max
    self.data[:fault][:vals][5] = dip_min
    self.data[:fault][:vals][6] = dip_max
end


end

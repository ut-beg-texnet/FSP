# load_default_data.jl (callbackloaddata.m in MATLAB)

# loads the default case
# 4 random wells


module LoadDefaultData


using Dates
using Random


export loaddefaultdata

function load_default_data(hDV)

    # Set default stress values
    hDV.data[:stress][:vals] = [1.1, 0.693, 1.22, 70, 11000, 0.43]

    # Set reservoir values
    hDV.data[:reservoir][:vals] = [100, 10, 200, 0.0 * 14.5, 0.5]  # Hardcoded Poisson's Ratio to 0.5

    # Set default fault values
    hDV.data[:fault][:vals] = [20.0, 0.58, 240.0, 330.0, 45.0, 90.0]
    hDV.data[:fault][:intype] = 1
    hDV.data[:fault][:file] = nothing  # Placeholder for no file

    # Set number of wells and default injection data
    hDV.data[:nwells] = 4
    hDV.data[:inject][:vals] = [
        [17, 12, 27000, 2015, 2022],
        [3, 4, 23000, 2016, 2026],
        [14, 7, 15000, 2017, 2027],
        [9, 9, 12000, 2018, 2025]
    ]

    # Disable real well data usage
    hDV.data[:realWellData][:use] = false

    # Set advanced properties for the region and fluid
    hDV.data[:adv][:vals] = [-5, 24, -3, 20, 1000, 9.80665, 0.8e-3, 3.6e-10, 3 * 3.6e-10, 20]

    # Set probabilistic hydrology values
    hDV.data[:probHydrology][:sigvals] = [25, 3, 50, 2, 1e-6, 1e-12, 1e-11, 750]
    hDV.data[:probHydrology][:probabilistic1VsDeterministic2] = 1
    hDV.data[:probHydrology][:distributionTxt] = [
        "Aquifer Thickness [100 ft]",
        "Porosity [10 %]",
        "Perm [200 mD]",
        "Fluid Density [1000 kg/m^3]",
        "Dynamic Viscosity [0.0008 Pa.s]",
        "Fluid Compressibility [3.6e-10 Pa^-1]",
        "Rock Compressibility [1.08e-09 Pa^-1]"
    ]

    # Set datenumBarrelsPerDay for constant rate wells
    for kk in 1:hDV.data[:nwells]
        hDV.data[:inject][:datenumBarrelsPerDay][kk] = [
            (Date(hDV.data[:inject][:vals][kk][4], 1, 1), 0),
            (Date(hDV.data[:inject][:vals][kk][5], 1, 1), hDV.data[:inject][:vals][kk][3])
        ]
    end
    hDV.data[:inject][:datenumBarrelsPerDay] = hDV.data[:inject][:datenumBarrelsPerDay][1:hDV.data[:nwells]]

    # Set sigma model for Monte Carlo default
    hDV.data[:sigvals] = [0.01, 0.01, 0.01, 0.0, 0.01, 5.0, 5.0, 5.0, 0.0, 0.01, 0.0, 0.0]

    # Stress aphi values
    hDV.data[:stress][:aphi][:use] = 0
    hDV.data[:stress][:aphi][:vals] = [1.2277, 0.0]

    # Set probabilistic traffic light settings
    hDV.plotdata[:minint] = 0.0
    hDV.plotdata[:maxint] = 0.5
    # GUI-related functionality (commented out)
    # set(hDV.plotdata.pint.cmintxt, 'string', '0')
    # set(hDV.plotdata.pint.cmaxtxt, 'string', '0.5')

    # Load faults from random distribution
    nfaults = hDV.data[:fault][:vals][1]  # Number of faults
    Xmin, Xmax = hDV.data[:adv][:vals][1], hDV.data[:adv][:vals][2]
    Ymin, Ymax = hDV.data[:adv][:vals][3], hDV.data[:adv][:vals][4]
    Tmin, Tmax = hDV.data[:fault][:vals][3], hDV.data[:fault][:vals][4]
    Dmin, Dmax = hDV.data[:fault][:vals][5], hDV.data[:fault][:vals][6]

    # Set random seed for consistency
    Random.seed!(20)

    # Generate random fault positions and orientations
    hDV.data[:fault][:xf] = Xmin .+ (Xmax - Xmin) .* rand(nfaults)
    hDV.data[:fault][:yf] = Ymin .+ (Ymax - Ymin) .* rand(nfaults)
    hDV.data[:fault][:thf] = Tmin .+ (Tmax - Tmin) .* rand(nfaults)
    hDV.data[:fault][:dipf] = Dmin .+ (Dmax - Dmin) .* rand(nfaults)
    hDV.data[:fault][:lenf] = fill(2.0, nfaults)  # Set fault length to 2
    hDV.data[:fault][:muf] = fill(hDV.data[:fault][:vals][2], nfaults)  # Friction coefficient

    # Enable the calculate button (GUI-related functionality)
    # set(hDV.bCalc, 'enable', 'on')

    # Call callbackcalc to run the calculation (not implemented yet)
    # callbackcalc(hDV)
end

end # module LoadDefaultData
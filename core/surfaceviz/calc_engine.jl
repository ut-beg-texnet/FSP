# calcengine.jl
# This module runs the main calculations for different tabs in the FSP program.

module CalcEngine

export calcengine

using LinearAlgebra
using PFieldCalc  # Import the PFieldCalc module to use pfieldcalc
using CalcST   # Import the CalcST module to use calcST
using Dates
using Printf

# Dummy placeholder for mohrs_3D, which should be defined elsewhere
function mohrs_3D(inputCell, hDV)
    # Placeholder function, should be replaced with actual implementation
    return (nothing, Dict("ppfail" => rand(50, 50)), nothing, nothing, nothing, nothing, nothing)
end

# Main calcengine function
function calcengine(hDV, name::String)

  
    if name == "HYDROLOGY" || name == "INTEGRATED"
        if hDV.data[:reservoir][:importHydrology] == 0
            # Internal hydrology calculation (0 means internal calculation)
            println("Calculating internal hydrology...")

            # Initialize variables for the calculation
            xmin, xmax = hDV.data[:adv][:vals][1], hDV.data[:adv][:vals][2]
            ymin, ymax = hDV.data[:adv][:vals][3], hDV.data[:adv][:vals][4]
            npts = 50
            nwells = hDV.data[:nwells]

            # Generate mesh grid
            hDV.plotdata[:pflot][:Xgrid], hDV.plotdata[:pflot][:Ygrid] = meshgrid(linspace(xmin, xmax, npts), linspace(ymin, ymax, npts))
            r = linspace(0.5, 20, npts)

            # Get storativity (S) and transmissivity (T)
            aqthick = hDV.data[:reservoir][:vals][1] * 0.3048  # Convert ft to meters
            phi = hDV.data[:reservoir][:vals][2] / 100  # Porosity fraction
            kap = hDV.data[:reservoir][:vals][3]  # Permeability in mD
            S, T = calcST(hDV.data[:adv][:vals], aqthick, phi, kap)  # Call calcST to calculate S and T

            # Run pressure field calculations using pfieldcalc
            hDV.plotdata[:pflot][:Zgrid] = pfieldcalc(hDV)

            # Initialize pressure curve arrays
            hDV.plotdata[:pflot][:x] = zeros(nwells, length(r))
            hDV.plotdata[:pflot][:y] = zeros(nwells, length(r))

            for k in 1:hDV.data[:nwells]
                # Location of well
                if hDV.data[:realWellData][:use]
                    xwell = hDV.data[:realWellData][:XEasting][k]
                    ywell = hDV.data[:realWellData][:YNorthing][k]
                else
                    xwell = hDV.data[:inject][:vals][k, 1]
                    ywell = hDV.data[:inject][:vals][k, 2]
                end

                # Calculate radial flow model using pfieldcalc
                p = pfieldcalc(hDV, xwell + r, ywell, k)
                hDV.plotdata[:pflot][:x][k, :] = r
                hDV.plotdata[:pflot][:y][k, :] = p
            end

            # Fault failure calculation using Mohr's circle
            thf = hDV.data[:fault][:thf]  # Fault strikes
            dips = hDV.data[:fault][:dipf]  # Fault dips
            SHdir = hDV.data[:stress][:vals][4]  # Max horizontal stress direction
            muf = hDV.data[:fault][:muf]  # Fault mu

            # Get pressure from hydrology model at each fault
            nfaults = hDV.data[:fault][:vals][1]
            ppOnFault = zeros(nfaults)

            for j in 1:nfaults
                wells = 1:hDV.data[:nwells]
                x = hDV.data[:fault][:xf][j]
                y = hDV.data[:fault][:yf][j]
                ppOnFault[j] = pfieldcalc(hDV, x, y, wells)
            end

            # Update fault failure calculations
            hDV.plotdata[:pint][:ppf] = ppOnFault
            dpth = hDV.data[:stress][:vals][5]
            sig = hDV.data[:stress][:vals][1:3] * dpth
            pp0 = hDV.data[:stress][:vals][6] * dpth
            nu = hDV.data[:reservoir][:vals][5]
            biot = 1  # Biot coefficient

            if hDV.data[:stress][:aphi][:use] == 1
                inputCell = [sig', 0.00, pp0, thf, dips, SHdir, ppOnFault, muf, biot, nu, hDV.data[:stress][:aphi][:vals][1]]
            else
                inputCell = [sig', 0.00, pp0, thf, dips, SHdir, ppOnFault, muf, biot, nu]
            end

            # Run Mohr's circle calculation
            _, hDV.plotdata[:resultsH][:outs], hDV.plotdata[:resultsH][:C1], hDV.plotdata[:resultsH][:C2], hDV.plotdata[:resultsH][:C3], hDV.plotdata[:resultsH][:sig_fault], hDV.plotdata[:resultsH][:tau_fault] = mohrs_3D(inputCell, hDV)

        else
            # Importing hydrology case (1 means external hydrology import)
            println("Importing hydrology data and interpolating...")

            # Example interpolation, real logic needs griddata implementation
            ts = hDV.hdsldr[:value]
            numbersImportedHydrology = hDV.data[:reservoir][:numbersImportedHydrology]
            thisYearData = filter(row -> row[4] == ts, numbersImportedHydrology)
            thisYearX, thisYearY, thisYearPSI = thisYearData[:, 1], thisYearData[:, 2], thisYearData[:, 3]

            # Interpolate imported hydrology onto fault
            ppOnFault = griddata(thisYearX, thisYearY, thisYearPSI, hDV.data[:fault][:xf], hDV.data[:fault][:yf])
            hDV.plotdata[:pint][:ppf] = ppOnFault
        end

    elseif name == "GEOMECHANICS" || name == "PROB. GEOMECH"
        println("Calculating Geomechanics...")

        # Fault failure calculation using Mohr's circle
        thf = hDV.data[:fault][:thf]  # Fault strikes
        dips = hDV.data[:fault][:dipf]  # Fault dips
        SHdir = hDV.data[:stress][:vals][4]  # Max horizontal stress direction
        muf = hDV.data[:fault][:muf]  # Fault friction coefficients
        dpth = hDV.data[:stress][:vals][5]  # Calculation depth
        sig = hDV.data[:stress][:vals][1:3] * dpth
        pp0 = hDV.data[:stress][:vals][6] * dpth  # Native pore pressure
        nu = hDV.data[:reservoir][:vals][5]  # Poisson's Ratio
        biot = 1  # Biot coefficient

        if hDV.data[:stress][:aphi][:use] == 1
            inputCell = [sig', 0.00, pp0, thf, dips, SHdir, 0 * thf, muf, biot, nu, hDV.data[:stress][:aphi][:vals][1]]
        else
            inputCell = [sig', 0.00, pp0, thf, dips, SHdir, 0 * thf, muf, biot, nu]
        end

        # Run Mohr's circle calculation
        _, hDV.plotdata[:results][:outs], hDV.plotdata[:results][:C1], hDV.plotdata[:results][:C2], hDV.plotdata[:results][:C3], hDV.plotdata[:results][:sig_fault], hDV.plotdata[:results][:tau_fault] = mohrs_3D(inputCell, hDV)
    else
        println("No calculations required for $name")
    end
end

end  # module CalcEngine

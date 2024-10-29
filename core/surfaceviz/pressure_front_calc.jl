# PFront.jl
# calculates the pressure front generated by a set of wells with
# variable injection rates using a radially symmetric finite thickness model.

module PressureFrontCalc

export pressurefrontcalc

using LinearAlgebra
using SpecialFunctions  

#Inputs
# ts: time steps
# S: storativity
# T: transmissivity
# rho: fluid density
# g: gravitational acceleration
# thisWellDatenumBarrelsPerDay: well data


function pfront(r, yearsToCalculate, thisWellDatenumBarrelsPerDay, S, T, rho, g)
    # Default values for storativity and transmissivity if not provided
    if S === nothing
        S = 5.4e-5  # m²/sec storativity
    end
    if T === nothing
        T = 4.5e-6  # m²/sec transmissivity
    end

    # Check if yearsToCalculate is an integer, else print a warning
    if abs(round(yearsToCalculate) - yearsToCalculate) > eps()
        println("Warning: pfront expects integer years (e.g., 2019). The input year is not an integer.")
    end

    # Ensure r is a column vector, reshape it if necessary
    r_reshaped = false
    rowr, colr = size(r)
    if colr != 1
        r_reshaped = true
        r = reshape(r, :, 1)  # Reshape to a column vector
    end

    # Initialize pressure delta matrix
    delta_p_Bars = zeros(size(r))

    # Loop over the years to calculate the pressure front for each time step
    for yearCycleCount in 1:length(yearsToCalculate)
        yearOfInterest = yearsToCalculate[yearCycleCount]

        # Get well data up until the year of interest
        thisWellDataUpTilNow = filter(row -> row[1] <= Date(yearOfInterest, 1, 1), thisWellDatenumBarrelsPerDay)

        # If no data exists before the year of interest, continue with zero pressure change
        if isempty(thisWellDataUpTilNow)
            if r_reshaped
                delta_p_Bars = zeros(rowr, colr)
            else
                delta_p_Bars[1, yearCycleCount] = 0
            end
            continue
        end

        # Add a data point at the beginning of the year if necessary
        if size(thisWellDataUpTilNow, 1) < size(thisWellDatenumBarrelsPerDay, 1)
            append!(thisWellDataUpTilNow, [Date(yearOfInterest, 1, 1), thisWellDatenumBarrelsPerDay[size(thisWellDataUpTilNow, 1) + 1, 2]])
        else
            append!(thisWellDataUpTilNow, [Date(yearOfInterest, 1, 1), 0])
        end

        # Calculate elapsed time in seconds
        tdays = (thisWellDataUpTilNow[:, 1] - thisWellDataUpTilNow[1, 1]) * Day(1)
        t = tdays * 86400  # Convert days to seconds

        # Convert barrels per day to cubic meters per second
        Q = thisWellDataUpTilNow[:, 2] * 1.84013e-6

        # Initialize matrices for well function calculations
        u = zeros(length(r))
        well_func = zeros(length(r))
        tstep_sum = zeros(length(r))

        # Loop over time steps
        for ii in 2:length(t)
            u[:, 1] = ((r[:, 1].^2) * T * S) ./ (4 * T * T * (maximum(t) - t[ii - 1]))
            well_func[:, 1] = expint.(u[:, 1])
            tstep_sum .+= well_func .* (Q[ii] - Q[ii - 1])
        end

        # Final head distribution and pressure calculation
        head = tstep_sum * (1 / (4 * π * T))
        delta_p_Pascals = head * rho * g  # Pressure in Pascals

        # Convert Pascals to Bars and reshape if necessary
        if r_reshaped
            delta_p_Pascals = reshape(delta_p_Pascals, rowr, colr)
            delta_p_Bars = delta_p_Pascals * 1e-5  # Convert to Bars
        else
            delta_p_Bars[:, yearCycleCount] = delta_p_Pascals * 1e-5  # Convert to Bars
        end

        # Handle NaN and infinity values
        delta_p_Bars[isnan.(delta_p_Bars)] .= 0
        delta_p_Bars[isinf.(delta_p_Bars)] .= 0
    end

    return delta_p_Bars
end

end  # module PFront

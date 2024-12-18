module HydroCalculations


using SpecialFunctions
using Plots
using Statistics
using PrettyTables # used to print the pressure field in a beautified table

export calcST, md_to_m2, pfieldcalc, pfront

"""
Convert permeability from millidarcy (mD) to square meters (m²)
"""
function md_to_m2(permeability_md::Float64)
    return permeability_md * 10^-3 * 9.9e-13 # conversion factor from mD to m²
end




#=
Inputs: 
    h: aquifer thickness
    phi: porosity fraction
    kap: intrinsic permeability
    rho: fluid density
    mu: dynamic viscosity
    g: accl gravity
    beta: fluid compressibility
    alphav: rock compressibility

Outputs:
    S: Storativity
    T: Transmissivity
=#
function calcST(h, phi, kap, rho, mu, g, beta, alphav)

    # convert permeability from mD to m²
    kap = md_to_m2(kap)
    
    # Calculate storativity
    S = rho * g * h * (alphav + phi * beta)

    # saturated hydraulic conductivity
    K = kap * rho * g / mu

    # Calculate transmissivity
    T = K * h

    # print all arguments and results (verified)
    #=
    println("Aquifer thickness: ", h)
    println("Porosity: ", phi)
    println("Permeability: ", kap)
    println("Fluid density: ", rho)
    println("Dynamic viscosity: ", mu)
    println("Gravity: ", g)
    println("Fluid compressibility: ", beta)
    println("Rock compressibility: ", alphav)
    println("Storativity: ", S)
    println("Transmissivity: ", T)
    println("Hydraulic conductivity: ", K)
    =#
    
    return S, T
end


#=
radially symmetric finite thickness pressure front calculation
# Storativity (S) and Transmissivity (T) are calculated externally
# Calculates the pressure change (delta_p_Bars) in an aquifer over time and space due to fluid injection from one or more wells. 
This pressure change is modeled using a radially symmetric representation of the aquifer and incorporates the effects of:
    1) Fluid flow
    2) Aquifer properties (storativity and transmissivity)
    3) Variable injection rates over time
# Inputs:
    T: transmissivity (assume it's the same for X and Y directions)
    rho: fluid density in SI (1000)
    g: acceleration due to gravity in SI (9.81)
    r: radial distance from well in meters
    year_range: range of years to calculate pressure field
    r: distances from wells in meters matrix
=#
function pfront(
    r::Matrix{Float64}, 
    year_range::StepRange{Int64, Int64},
    year_of_interest::Int64,
    datenum_barrels_per_day::Matrix{Float64}, 
    S::Float64 = 5.4e-5, 
    T::Float64 = 4.5e-6, 
    rho::Float64 = 1000.0, 
    g::Float64 = 9.81)

    #println("Inside pfront function...")
    # convert r to column vector
    r_reshaped = false
    if size(r, 2) != 1
        r_reshaped = true
        r = reshape(r, :, 1)
    end

    # initialize output
    # a matrix to store the calculated pressure changes at each radial distance r for each year in year_range
    delta_p_bars = zeros(Float64, size(r,1), length(year_range))

    # loop over years
    for (year_idx, year) in enumerate(year_range)
        #year = year_range[yearidx]
        # filter data for years up to the current year
        # include only the rows where the year is less that or equal to 'year'
        filtered_data = datenum_barrels_per_day[datenum_barrels_per_day[:, 1] .<= year_of_interest, :]

        # print the filtered data
        #println("Filtered data: ", filtered_data)
        #println("Number of columns in filtered data: ", size(filtered_data, 2))
        #println("Dimensions of filtered data: ", size(filtered_data))
        #println("Dimensions of datenum_barrels_per_day: ", size(datenum_barrels_per_day))

        #println("datenum_barrels_per_day: ", datenum_barrels_per_day)
        #println("Number of columns in datenum_barrels_per_day: ", size(datenum_barrels_per_day, 2))

        # if no data exists for this year, set pressure to 0
        if isempty(filtered_data)
            println("No data for year ", year)
            return zeros(Float64, size(r, 1), 1)
        end

        # Add a datapoint at the start of the year if not present
        if size(filtered_data, 1) < size(datenum_barrels_per_day, 1)
            # update the filtered_data matrix by appending the next available value from datenum_barrels_per_day
            next_value = datenum_barrels_per_day[size(filtered_data, 1) + 1, 2]
            new_row = hcat(fill(year_of_interest, 1, size(filtered_data, 2) - 1), [next_value])
            if size(new_row, 2) != size(filtered_data, 2)
                error("Shape mismatch: Cannot append new row to filtered_data")
            end
            filtered_data = vcat(filtered_data, new_row)
        else
            # add a new row with a value of 0.0 if all original data has already been included
            new_row = hcat(fill(year_of_interest, 1, size(filtered_data, 2) - 1), [0.0])
            if size(new_row, 2) != size(filtered_data, 2)
                error("Shape mismatch: Cannot append new row to filtered_data")
            end
            filtered_data = vcat(filtered_data, new_row)
        end

        # convert elapsed time form days to seconds
        t_days = filtered_data[:, 1] .- datenum_barrels_per_day[1, 1]
        t_seconds = t_days .* 24 .* 3600


        #println("t_seconds: ", t_seconds)
        #println("t_days: ", t_days)

        # convert barrels per day to m^3 per second
        Q = filtered_data[:, 2] * 1.84013e-6

        # Well function --> models radial pressure propagation in an aquifer due to fluid injection
        # well function params
        u = zeros(Float64, size(r, 1)) # input equation into well-function (dimensionless)
        # Well-function matrix --> From Allen (1954) and Hastings (1955) where u varies greatly [Found in Practical Approximations of Well Function (Srivastava and Guzman-Guzman)]
        well_func = zeros(Float64, size(r, 1))
        
        tstep_sum = zeros(Float64, size(r, 1)) # placeholder matrix for each timestamp sum

        
        # loop over time in seconds
        for i in 2:length(t_seconds)
            u .= (r .^ 2 .* S) ./ (4 * T * (t_seconds[end] - t_seconds[i - 1])) # u matrix calculation at each x,y point
            well_func .= expint.(u) # well function calculation
            # sum previous tstep_sum with new well_function and injection rate step. Note: If Q has not changed, tstep_sum unchanged.
            tstep_sum .+= well_func .* (Q[i] - Q[i - 1]) 
        end
        
        #println("tstep_sum min: ", minimum(tstep_sum))
        #println("tstep_sum max: ", maximum(tstep_sum))
        #println("tstep_sum min: ", minimum(tstep_sum))

        

        #=
        final_time = t_seconds[end]
        u = (r .^ 2 .* S) ./ (4 * T)
        delta_t_seconds = final_time .- t_seconds[1:end-1]
        delta_Q = Q[2:end] .- Q[1:end-1]
        u_values = u ./ delta_t_seconds'
        well_func_values = expint.(u_values) # models the pressure decay with distance (use it for the graph)
        tstep_sum = well_func_values * delta_Q
        =#

        # debugging prints
        #println("min u: ", minimum(u))
        #println("max u: ", maximum(u))
        #println("S: ", S)
        #println("T: ", T)


        # head and pressure calculations
        head = tstep_sum .* (1/(4 * pi * T)) # head calculation
        #println("min value of head: ", minimum(head))
        #println("max value of head: ", maximum(head))
        delta_p_pascals = head .* rho * g # pressure calculation in pascals
        delta_p_bars[:, year_idx] .= delta_p_pascals .* 1e-5 # convert pressure from pascal to bars

        #println("min value of delta_p_bars", minimum(delta_p_bars))
        #println("max value of delta_p_bars", maximum(delta_p_bars))
        #println("min value of delta_p_pascals: ", minimum(delta_p_pascals))
        #println("max value of delta_p_pascals: ", maximum(delta_p_pascals))


        # reshape back to shape of r (if needed)
        if r_reshaped
            delta_p_bars = reshape(delta_p_bars, size(r, 1), size(delta_p_bars, 2))
        end

        #println("Pressure field calculated for year ", year, "= ", delta_p_bars)
        return delta_p_bars

    end
    #println("delta_p_bars min: ", minimum(delta_p_bars))
    #println("delta_p_bars max: ", maximum(delta_p_bars))


end

# using multiple dispatch in case the function gets called with two arguments
function pfront(T::Float64, S::Float64)
    println("pfront function called with two arguments")
    # m^2/sec storativity
    S = 5.4e-5
    # m^2/sec transmissivity
    T = 4.5e-6
end





#=
Computes the pressure field generated by a set of wells using
superposition from a simple radially symmetric pressure model
=#
function pdfieldcalc(injection_wells::Dict{String, Any}, xGrid::AbstractMatrix, yGrid::AbstractMatrix, num_wells::Int, aquifer_thickness, porosity, permeability, fluid_density, dynamic_viscosity, fluid_compressibility, rock_compressibility, year_range::StepRange{Int64, Int64}, datenum_barrels_per_day::Matrix{Float64})
    
    # Initialize pressure field
    pField = zeros(size(xGrid)) # a 2D grid of pressure values over the aquifer (50x50)

    # convert aquifer thickness from ft to m
    aquifer_thickness = aquifer_thickness * 0.3048
    S, T = calcST(aquifer_thickness, porosity, permeability, fluid_density, dynamic_viscosity, 9.81, fluid_compressibility, rock_compressibility)
    #println("Storativity from calcST: ", S)
    #println("Transmissivity from calcST: ", T)
    # setup well location and injection rate from variable or constant rate structures
    #// check if the word 'month' is used in the input file under 'injection_wells'
    if occursin("month", string(injection_wells["1"]))
        println("Monthly injection rates are used")
    else
        println("Constant injection rates are used")
        well_locations = zeros(Float64, num_wells, 2)
        #if the word 'month' is not used, then the injection rate is a constant rate
        for well in keys(injection_wells)
            # get xy coordinates from 'northing_km' and 'easting_km' under 'location'
            x = injection_wells[well]["location"]["easting_km"]
            y = injection_wells[well]["location"]["northing_km"]
            well_locations[parse(Int, well), :] = [x, y]
            #print("Well ", well, " at ", x, ", ", y, "\n")
        end

        # get injection rates from 'rate_bbl_day'
        injection_rates = zeros(num_wells)
        for well in keys(injection_wells)
            injection_rates[parse(Int, well)] = injection_wells[well]["injection_rate"]["rate_bbl_day"]
        end

        # Print well locations
        for i in 1:num_wells
            #println("Well ", i, " at x:", well_locations[i, 1], ", y:", well_locations[i, 2])
        end
        # print injection rates
        for i in 1:num_wells
            #println("Well ", i, " injection rate per day: ", injection_rates[i])
        end
    end

    # loop over each wel and calculate pressure field
    for i in 1:num_wells
        well = string(i)
        xwell = well_locations[i, 1]
        ywell = well_locations[i, 2]

        println("Well Coordinates: xwell = ", xwell, ", ywell = ", ywell)

        println("Grid Shape: ", size(xGrid), " x ", size(yGrid))



        # get data from well
        injection_rate = injection_rates[i]
        start_year = injection_wells[well]["injection_rate"]["start_year"]
        end_year = injection_wells[well]["injection_rate"]["end_year"]
        # grid distances to well
        # calculate radial distances from the well located at xwell, ywell, to various points (x,y) on the grid using complex numbers
        difflocs = xGrid .+ im .* yGrid .- (xwell + im * ywell)
        # CONTINUE FROM HERE
        #println("difflocs min: ", minimum(abs.(difflocs)))
        #println("difflocs max: ", maximum(abs.(difflocs)))

        
        # NOT THE SAME AS MATLAB
        R = abs.(difflocs) # radial distance
        #println("Radial distance min: ", minimum(R))
        #println("Radial distance max: ", maximum(R))
        rmeters = R * 1e3 # convert to meters
        #println(typeof(rmeters))

        # print the radial distances 50x50 matrix using PrettyTables
        println("Radial distances (sample rows):")
        pretty_table(rmeters[1:5, 1:5])  # Adjust range for testing
        println("Size of radial distances: ", size(rmeters))
        println("Min radial distance: ", minimum(rmeters))
        println("Max radial distance: ", maximum(rmeters))

        # rho and g constants
        rho = fluid_density
        g = 9.81 # this is from Rall's comments but on MATLAB debugger I see 9.8067

        # initialize pfront_result
        pfront_result = zeros(Float64, size(rmeters, 1), length(year_range))

        # call pfront function for pressure calculation
        # here, 'year' is the yearOfInterest in MATLAB
        for year in year_range
            #println("Calculating pressure field for year ", year)
            pfront_result = pfront(rmeters, year_range, year, datenum_barrels_per_day, S, T, rho, g)
        end


        for t = 1:size(pfront_result, 2)
            reshaped_pfront_result = reshape(pfront_result[:, t], size(pField))
            pField .+= reshaped_pfront_result # superpose pressure field
        end
    end


    # remove infinity values, typically occurs when well centers exactly match a grid location
    pField = replace(pField, Inf => 0.0)
    # remove NaN values, typically occurs when time at calc is less that the well start times
    pField = replace(pField, NaN => 0.0)
    #convert to Psi
    pField = pField * 14.5

    println("Pressure field calculated for all wells")

    return pField

    # create a heatmap of the pressure field 
    heatmap(pField, title="Pressure Field", xlabel="X", ylabel="Y", color=:viridis)
    # save the heatmap inside the graphs folder
    savefig("../graphs/pressure_field.png")


    display(heatmap(pField, title="Pressure Field", xlabel="X", ylabel="Y", color=:viridis))

end # pfieldcalc function
    
    
end # module

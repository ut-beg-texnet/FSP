using JSON
using ArgParse
using LinearAlgebra
using Statistics
using Distributions
using Dates
using ProgressBars


include("graphs/julia_fsp_graphs.jl")
include("core/hydrology_calculations.jl")
include("step4_hydrology.jl")
include("core/bill_pfront.jl")

using .JuliaFSPGraphs
using .HydroCalculations
using .DetHydrologyDriver
using .BillPFront


struct HydrologyParams
    aquifer_thickness::Float64
    porosity::Float64
    permeability::Float64
    fluid_density::Float64
    dynamic_viscosity::Float64
    fluid_compressibility::Float64
    rock_compressibility::Float64
    plus_minus::Dict{String, Float64}
    n_iterations::Int64
end



#=
For a Uniform distribution, we will need:
- Base value
- Plus/minus value
=#
function create_uniform_distribution(base_value::Union{Float64, Integer}, plus_minus::Union{Float64, Integer})
    # check for negative plus/minus values
    if plus_minus < 0
        throw(ArgumentError("Plus/minus value cannot be negative. Setting it to 0.0."))
        plus_minus = 0.0
    end
    
    min_value = base_value - plus_minus
    max_value = base_value + plus_minus
    return Uniform(min_value, max_value)
end

#=
For a Gaussian (Normal) distribution, we will need:
- Mean (base value)
- Standard deviation
- Optionally, bounds if you want a truncated normal
=#
function create_gaussian_distribution(mean_value::Union{Float64, Integer}, std_dev::Union{Float64, Integer})
    # check if the standard deviation is negative
    if std_dev < 0
        throw(ArgumentError("Standard deviation cannot be negative. Setting it to 0.0."))
        std_dev = 0.0
    end
    
    return Normal(mean_value, std_dev)
end

function run_monte_carlo_hydrology(params::HydrologyParams, distribution_type::String="uniform")

    if distribution_type == "uniform"
        

        # if we are missing the plus_minus value for a parameter, set it to 0.0
        for (key, value) in params.plus_minus
            if isnan(value)
                params.plus_minus[key] = 0.0
            end
            # if the n_iterations was missing (so it's 0.0), set it to the default 750 iterations
            if params.n_iterations == 0 || isnothing(params.n_iterations)
                params.n_iterations = 750
            end
        end

        # check that the plus_minus values are not greater than the base values
        for (key, value) in params.plus_minus
            if key == "aquifer_thickness" && value > params.aquifer_thickness ||
                key == "porosity" && value > params.porosity ||
                key == "permeability" && value > params.permeability ||
                key == "fluid_density" && value > params.fluid_density ||
                key == "dynamic_viscosity" && value > params.dynamic_viscosity ||
                key == "fluid_compressibility" && value > params.fluid_compressibility ||
                key == "rock_compressibility" && value > params.rock_compressibility
                throw(ArgumentError("plus_minus value for $key cannot be greater than the base value."))
            end
            
        end

        # create the Uniform distributions for each parameter
        distributions = Dict(
            "aquifer_thickness" => create_uniform_distribution(params.aquifer_thickness, params.plus_minus["aquifer_thickness"]),
            "porosity" => create_uniform_distribution(params.porosity, params.plus_minus["porosity"]),
            "permeability" => create_uniform_distribution(params.permeability, params.plus_minus["permeability"]),
            "fluid_density" => create_uniform_distribution(params.fluid_density, params.plus_minus["fluid_density"]),
            "dynamic_viscosity" => create_uniform_distribution(params.dynamic_viscosity, params.plus_minus["dynamic_viscosity"]),
            "fluid_compressibility" => create_uniform_distribution(params.fluid_compressibility, params.plus_minus["fluid_compressibility"]),
            "rock_compressibility" => create_uniform_distribution(params.rock_compressibility, params.plus_minus["rock_compressibility"])
        )

    elseif distribution_type == "gaussian"
        # create the Gaussian distributions for each parameter
        distributions = Dict(
            "aquifer_thickness" => create_gaussian_distribution(params.aquifer_thickness, params.plus_minus["aquifer_thickness"]),
            "porosity" => create_gaussian_distribution(params.porosity, params.plus_minus["porosity"]),
            "permeability" => create_gaussian_distribution(params.permeability, params.plus_minus["permeability"]),
            "fluid_density" => create_gaussian_distribution(params.fluid_density, params.plus_minus["fluid_density"]),
            "dynamic_viscosity" => create_gaussian_distribution(params.dynamic_viscosity, params.plus_minus["dynamic_viscosity"]),
            "fluid_compressibility" => create_gaussian_distribution(params.fluid_compressibility, params.plus_minus["fluid_compressibility"]),
            "rock_compressibility" => create_gaussian_distribution(params.rock_compressibility, params.plus_minus["rock_compressibility"])
        )
    else
        throw(ArgumentError("Invalid distribution type."))
    end


    # from the JSON file, get the number of faults
    input_data = JSON.parsefile("output/model_inputs.json")
    num_faults = length(input_data["faults"])

    # create the matrix to store the results of the Monte Carlo simulations
    ppOnFaultMC = zeros(params.n_iterations, num_faults) # rows are iterations, columns are faults

   

    # ask for the year of interest
    println("Enter the year of interest:")
    year_of_interest = parse(Int, readline())

    # run the Monte Carlo simulations
    println("Running MC iterations on faults...")
    for i in ProgressBar(1:params.n_iterations)
        println("Running iteration $i...")
        # sample the parameters from the distributions
        sampled_params = Dict(
            "aquifer_thickness" => rand(distributions["aquifer_thickness"]),
            "porosity" => rand(distributions["porosity"]),
            "permeability" => rand(distributions["permeability"]),
            "fluid_density" => rand(distributions["fluid_density"]),
            "dynamic_viscosity" => rand(distributions["dynamic_viscosity"]),
            "fluid_compressibility" => rand(distributions["fluid_compressibility"]),
            "rock_compressibility" => rand(distributions["rock_compressibility"])
        )

        # calculate storativity and transmissivity
        S, T = calcST(
            sampled_params["aquifer_thickness"], 
            sampled_params["porosity"], 
            sampled_params["permeability"], 
            sampled_params["fluid_density"], 
            sampled_params["dynamic_viscosity"], 
            9.81, 
            sampled_params["fluid_compressibility"], 
            sampled_params["rock_compressibility"]
        )


        STRho = (S, T, sampled_params["fluid_density"]) # used in presssureScenario function


        # extract fault data from the JSON
        local faults = get(input_data, "faults", [])
        if isempty(faults)
            throw(ArgumentError("Faults data is missing from the input JSON."))
        end
        
        # extract injection well data from the JSON
        injection_wells = get(input_data, "injection_wells", [])
        if isempty(injection_wells)
            throw(ArgumentError("Injection wells data is missing from the input JSON."))
        end

        

        # prepare radial distances
        r_km = range(0.5, stop=20.0, length=50)
        r_meters = r_km .* 1000


        
        # loop over faults
        for f in 1:length(faults)
            # loop over wells
            for w in 1:length(injection_wells)
                # get all data for this well
                well_data = injection_wells[string(w)]

                # get injection dates from wells
                inj_start_year = well_data["injection_rate"]["start_year"]
                inj_end_year = well_data["injection_rate"]["end_year"]

                # check if the year of interest is within the injection period
                if year_of_interest < inj_start_year
                    println("Year of interest is before injection period. Setting ppOnFault to 0.")
                    continue
                end
                local actual_end = min(inj_end_year, year_of_interest)
                if actual_end <= inj_start_year
                    println("Year of interest is after injection period. Setting ppOnFault to 0.")
                    continue
                end


                # fault coordinates
                x_fault_km = faults[f]["easting"]
                y_fault_km = faults[f]["northing"]

                # well coordinates
                x_well_km = well_data["location"]["easting_km"]
                y_well_km = well_data["location"]["northing_km"]

                # create timestamps for this well
                days, dpbs = prepare_well_data_for_pressureScenario(
                    well_data,
                    inj_start_year,
                    actual_end
                )

                if length(days) == 0 || length(dpbs) == 0
                    println("No injection days for the year of interest. Setting ppOnFault to 0.")
                    continue
                end

                local ppOnFault = pfieldcalc_all_rates(
                    x_fault_km,
                    y_fault_km,
                    STRho,
                    days,
                    dpbs,
                    x_well_km,
                    y_well_km,
                    Date(year_of_interest, 1, 1)
                )

                # store the results in the matrix
                ppOnFaultMC[i, f] += ppOnFault

            end # end of well loop
        end # end of fault loop
        
    end # end of MC loop

    # save the results to a JSON file
    open("output/prob_hydro_results.json", "w") do file
        JSON.print(file, Dict("monte_carlo_results" => ppOnFaultMC))
    end

    println("Finished running Monte Carlo simulations for hydrology.")
    println("Probabilistic Hydrology Process: max pressure on fault: ", maximum(ppOnFaultMC))
    println("Probabilistic Hydrology Process: min pressure on fault: ", minimum(ppOnFaultMC))


end # end of run_monte_carlo_hydrology function










function main()


    # first we need to check if the user want to run deterministic or probabilistic hydrology
    model_check = 1 # 0 for deterministic, 1 for probabilistic


    if model_check == 0
        # deterministic hydrology
        println("Probabilistic Hydrology Process: Running deterministic hydrology model.")

        # read the output JSON file from the probabilistic geomechanics step to get the CDF values of each fault
        prob_geo_results = JSON.parsefile("output/monte_carlo_results.json")
        # reformat it to a dictionary with fault_id as the key and the pressure values
        prob_geo_values = Dict{String, Dict{String, Any}}()

        for (fault_id, fault_data) in prob_geo_results["monte_carlo_results"]
            slip_pressures = [iter["slip_pressure"] for iter in fault_data["iterations"]]
            prob_geo_values[fault_id] = Dict(
                "pressures" => slip_pressures,
                "stats" => fault_data["slip_pressure"]
            )

        end

        det_hydro_values = JSON.parsefile("output/hydro_ppOnFault.json")

        

        # save the dictionary to a JSON file
        open("graphs/prob_geo_values.json", "w") do file
            JSON.print(file, prob_geo_values)
        end

        # plot the CDF 
        plot_cdf_det_hydro(prob_geo_values, det_hydro_values)


    else
        # probabilistic hydrology
        println("Probabilistic Hydrology Process: Running probabilistic hydrology model.")


        # prepare the CDF from the Probabilistic Geomechanics step to combine it with this one
        # read the output JSON file from the probabilistic geomechanics step to get the CDF values of each fault
        prob_geo_results = JSON.parsefile("output/monte_carlo_results.json")
        # reformat it to a dictionary with fault_id as the key and the pressure values
        prob_geo_values = Dict{String, Dict{String, Any}}()

        for (fault_id, fault_data) in prob_geo_results["monte_carlo_results"]
            slip_pressures = [iter["slip_pressure"] for iter in fault_data["iterations"]]
            prob_geo_values[fault_id] = Dict(
                "pressures" => slip_pressures,
                "stats" => fault_data["slip_pressure"]
            )

        end


        model_inputs_data = JSON.parsefile("output/model_inputs.json")

        aquifer_thickness = model_inputs_data["hydrology"]["aquifer_thickness"]
        porosity = model_inputs_data["hydrology"]["porosity"]
        permeability = model_inputs_data["hydrology"]["permeability_md"]
        fluid_density = model_inputs_data["model_parameters"]["fluid_density"]
        dynamic_viscosity = model_inputs_data["model_parameters"]["dynamic_viscosity"]
        fluid_compressibility = model_inputs_data["model_parameters"]["fluid_compressibility"]
        rock_compressibility = model_inputs_data["model_parameters"]["rock_compressibility"]


        # hardcode the parameters and distribution ranges for Now
        params = HydrologyParams(
            aquifer_thickness,    # aquifer_thickness
            porosity,     # porosity
            permeability,    # permeability
            fluid_density,     # fluid_density
            dynamic_viscosity,   # dynamic_viscosity
            fluid_compressibility,   # fluid_compressibility
            rock_compressibility,   # rock_compressibility
            Dict(    # plus_minus values
                "aquifer_thickness" => 25.0,
                "porosity" => 0.003,
                "permeability" => 50.0,
                "fluid_density" => 2.0,
                "dynamic_viscosity" => 1e-06,
                "fluid_compressibility" => 1e-12,
                "rock_compressibility" => 1e-11
            ),
            750     # n_iterations
        )


        run_monte_carlo_hydrology(params)

        prob_hydro_results = JSON.parsefile("output/prob_hydro_results.json")
        plot_prob_hydro_combined_cdf(prob_geo_values, prob_hydro_results["monte_carlo_results"])

        println("Probabilistic Hydrology Process: Finished running probabilistic hydrology model.")

    end




end





if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
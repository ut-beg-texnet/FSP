module FaultSlipCDFPlot

using Plots
using Statistics

# Input: matrix with results from MC simulation
# Each row --> a different fault
# Each column --> a different iteration for that fault
function plot_cdf(mc_results_matrix::Matrix{Float64}, output_path::String="prob_geo_cdf_plot.png")
    num_faults, num_iterations = size(mc_results_matrix)
    
    
    plt = plot(
        xlabel="Δ Pore Pressure to Slip [psi]",
        ylabel="Probability of Fault Slip",
        title="CDF of Δ Pore Pressure to Slip for Faults",
        xlim=(minimum(mc_results_matrix) * 0.95, maximum(mc_results_matrix) * 1.2),
        ylim=(0, 1),
        linewidth=2,
        legend=:outerright,  
        size=(900, 600),     # plot size
        left_margin=10Plots.mm, 
        bottom_margin=10Plots.mm 
    )

    # Define distinct colors for each fault using a color palette
    colors = distinguishable_colors(num_faults)

    # Loop over each fault (row of mc_results_matrix)
    for fault_id in 1:num_faults
        # Extract data for current fault
        pressures = mc_results_matrix[fault_id, :]

        
        sorted_pressures = sort(pressures) # might need provide an already sorted matrix so we don't have additional overhead

        # cumulative probabilities
        # create a vector from range 1 to number_iterations
        # normalize it by dividing by the number of iterations
        # this results in a vector of probabilities from 0 to 1
        probabilities = collect(1:num_iterations) / num_iterations

        # Add the current fault's CDF to the plot
        plot!(
            plt, 
            sorted_pressures, 
            probabilities, 
            label="Fault $fault_id", 
            lw=2, 
            alpha=0.8,              
            linestyle=:solid,       
            color=colors[fault_id]  # unique color to each fault
        )
    end

   
    savefig(plt, output_path)
    println("Plot saved to $output_path")
end

end # module

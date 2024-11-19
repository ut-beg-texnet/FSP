module FaultSlipCDFPlot

using Plots
using Statistics

function plot_cdf(mc_results_matrix::Matrix{Float64}, output_path::String="prob_geo_cdf_plot.png")
    num_faults, num_iterations = size(mc_results_matrix)
    #println("Plotting CDF for $num_faults faults with $num_iterations iterations each.")

    # Create a base plot with a distinct color palette
    plt = plot(
        xlabel="Δ Pore Pressure to Slip [psi]",
        ylabel="Probability of Fault Slip",
        title="CDF of Δ Pore Pressure to Slip for Faults",
        legend=:topright,
        xlim=(minimum(mc_results_matrix) * 0.95, maximum(mc_results_matrix) * 1.05),
        ylim=(0, 1),
        linewidth=2
    )

    # Define distinct colors for each fault using a color palette
    colors = distinguishable_colors(num_faults)

    # Loop over each fault (row of mc_results_matrix)
    for fault_id in 1:num_faults
        # Extract data for the current fault
        pressures = mc_results_matrix[fault_id, :]

        # Sort pressures for cumulative distribution calculation
        sorted_pressures = sort(pressures)

        # Generate cumulative probabilities
        probabilities = collect(1:num_iterations) / num_iterations

        # Debugging: Print summary for each fault being plotted
        #println("Plotting Fault $fault_id: min=$(minimum(pressures)), max=$(maximum(pressures)), mean=$(mean(pressures))")

        # Add the current fault's CDF to the plot with distinct colors and line styles
        plot!(
            plt, 
            sorted_pressures, 
            probabilities, 
            label="Fault $fault_id", 
            lw=2, 
            alpha=0.8,              # Reduced transparency for better visibility
            linestyle=:solid,       # Use solid lines for clarity; consider varying styles if needed
            color=colors[fault_id]  # Assign a unique color to each fault
        )
    end

    # Save the plot to the specified file
    savefig(plt, output_path)
    println("Plot saved to $output_path")
end

end # module

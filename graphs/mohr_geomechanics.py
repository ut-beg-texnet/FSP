#!/home/seiscomp/fsp_3/fsp_3/src/graph_helpers/venv/bin/python
import json
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os

# backend for matplotlib to allow for GUI interaction
matplotlib.use('Qt5Agg')

def plot_mohr_diagram(data):
    # Retrieve Mohr circles' x-axis min and max values for each hemisphere
    hemispheres = ["C1", "C2", "C3"]
    circles_x_values = [(data["mohr_circle_data"][h]["min"], data["mohr_circle_data"][h]["max"]) for h in hemispheres]
    
    sigma_fault = data["geomechanics_results"]["sig_fault"]
    tau_fault = data["geomechanics_results"]["tau_fault"]
    
    # Initialize plot
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_facecolor('lightgrey')
    # add grid
    ax.grid(True, linestyle='--', color='white', linewidth=0.7)
    
    # Plot Mohr circles based on provided min and max x-values
    for i, (x_min, x_max) in enumerate(circles_x_values):
        # Calculate circle properties from x_min and x_max
        center = (x_max + x_min) / 2
        radius = (x_max - x_min) / 2
        
        # Draw the circle
        circle = plt.Circle((center, 0), radius, color='grey', fill=False, linewidth=2)
        ax.add_patch(circle)
        
        # Annotate each circle
        label = r'$\sigma_h$' if i == 0 else (r'$\sigma_V$' if i == 1 else r'$\sigma_H$')
        ax.text(center, -radius - 100, label, color='cyan', fontsize=12, ha='center')
    
    # Plot frictional slip line based on start and end points in JSON
    if "frictional_slip_line" in data["mohr_circle_data"]:
        friction_start = data["mohr_circle_data"]["frictional_slip_line"]["start"]
        friction_end = data["mohr_circle_data"]["frictional_slip_line"]["end"]
        ax.plot([friction_start[0], friction_end[0]], [friction_start[1], friction_end[1]], 'b--', label="Frictional Slip Line")
    
    # Plot fault stress points
    sc = ax.scatter(sigma_fault, tau_fault, c=tau_fault, cmap="YlOrRd", edgecolor='k', s=100, label="Faults")
    
    plt.colorbar(sc, ax=ax, label=r'$\tau$ [psi]') # creates the color bar on the right (might remove later)
    
    # Axis labels and limits
    max_sigma = max([x_max for x_min, x_max in circles_x_values])  # Get the max x-value across circles
    max_tau = max(tau_fault + [friction_end[1]])  # Get the max y-value between faults and slip line end
    ax.set_xlabel(r'$\sigma$ effective [psi]')
    ax.set_ylabel(r'$\tau$ [psi]')
    ax.set_xlim(0, max_sigma + 1000)
    ax.set_ylim(0, max_tau + 500)
    ax.legend(loc="upper left")
    
    # Show plot
    plt.show()

def load_json_data(file_path):
    """Load JSON data from the specified file path."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data

# Example usage:
file_path = '../output/deterministic_geomechanics_results.json'  # Replace with the actual path to your JSON file

# Load data from the JSON file
data = load_json_data(file_path)

# Call the function to plot the Mohr diagram
plot_mohr_diagram(data)

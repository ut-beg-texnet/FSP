# CallbackCalc.jl
# This function handles calculations when triggered by the user (e.g., through a web request).

module LinkTabToCalculations

export linktabtocalculations

using Dates  # For working with date calculations
using Random # For random number generation

using CalcEngine  # For the core technical calculations based on the current tab

# links te current tab to the calculations that need to be run, and calls calc_engine.jl to run them
function link_tab_to_calculations(hDV)
    
    # Get the current tab (stored in hDV.currtab[:name])
    curr_tab = hDV.currtab[:name]

    # Call calcengine with the current data structure and the current tab name
    calcengine(hDV, curr_tab)

    # Post-calculation actions (could trigger frontend updates via web responses)
    println("Calculation completed for the current tab: $curr_tab")

end


end  # module LinkTabToCalculations

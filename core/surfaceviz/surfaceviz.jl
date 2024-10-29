# Define SurfaceViz object (surfaceviz.m)
module SurfaceViz

export SurfaceVizStruct, SurfaceViz

include("setup_data.jl")
using .SetupData


mutable struct SurfaceVizStruct
    data::Dict{Symbol, Any} # input data goes here
    plotdata::Dict{Symbol, Any} # data used for the visualization (D3 module can replace this)
    colors::Dict{Symbol, Any} # colors used for GUI
    cmapGYR::Array{Float64,2} # colormap for GUI [NEED TO IMPLEMENT]

end



# Constructor
function SurfaceVizStruct()
    self = SurfaceVizStruct(
        Dict{Symbol, Any}(), # data
        Dict{Symbol, Any}(), # plotdata
        Dict{Symbol, Any}(), # colors
        Array{Float64,2}(undef, 300, 3) # cmapGYR [PLACEHOLDER]
    )

    # set up colors
    self.colors[:red] = [1.0, 0.0, 0.0]
    self.colors[:green] = [0.0, 0.7, 0.0]
    self.colors[:white] = [1.0, 1.0, 1.0]
    self.colors[:black] = [0.0, 0.0, 0.0]
    self.colors[:yellow] = [1.0, 1.0, 0.0]

    # load colormap [TO DO: TO LOAD FROM EXTERNAL FILE]
    #self.cmapGYR = load_cmapGYR()

    # initialize default data
    self.data[:nwells_max] = 100 # max number of wells
    self.data[:NFAULTSMAX] = 500 # max number of faults

    # call data setup module
    setupdata(self)

    return self
end


end # module SurfaceViz
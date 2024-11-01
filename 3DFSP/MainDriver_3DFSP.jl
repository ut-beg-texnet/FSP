
################################################################################################################
dir=@__FILE__
mainDir="/" .* joinpath(split(dir,"/")[1:end-1])*"/"

additionalPath = ".//3DFSP_Public/Julia/LIB/Modules/"
push!(LOAD_PATH, additionalPath)

import FSP3D
import FSP3Dplots
####################################################################################################################################################################################################################
inputType="Faults"

dirFaultData=mainDir*"Data/Faults/FaultFramework/" ## Original

println("\n \t Directory containing fault data = ", dirFaultData)

@assert isdir(dirFaultData) "Dir not found "*string(dirFaultData)

fileNamesAll=FSP3D.GetFaultNamesInDirectory(dirFaultData)
fileNamesAll=fileNamesAll[1:2]


####
inputUnits = Dict(:stressGradient => "psi/ft")  ## Units of the stress gradients defined by the user
outputUnits = Dict(:depth => "ft")  ## m or ft                  

param=Dict(:mainDir => mainDir,
           :inputType => inputType,
           :dirInputData => dirFaultData,
           :IP3model => [],  #Not used in this case. User will define all the input stress parameters
           :probFailureAll => [1,5,10,50],           ## Define probability of failure to export to Petrel
           :unitsDict => Dict(:inputUnits => inputUnits,
                              :outputUnits => outputUnits
                             ),
           :nSamples => 10 ## Number of monte carlo simulations
           )


## Defining PDFs for the random variables
userDefinedPdfs = Dict(:SvGrad           =>   Dict(:pdfType => "Uniform", :pInterval => [0.98, 1.0]),
                       :PpGrad           =>   Dict(:pdfType => "Constant", :value => 0.44),
                       :ShminGrad        =>   Dict(:pdfType => "Uniform", :pInterval => [0.62, 0.65] ),                       
                       :μ                =>   Dict(:pdfType => "Uniform", :pInterval  => [0.6, 0.65]),
                       :SHmaxDir         =>   Dict(:pdfType => "Uniform", :pInterval => [80,100]),
                       :surfaceElevation =>   Dict(:pdfType => "Constant", :value => 0),
                       :Aphi             =>   Dict(:pdfType => "Uniform", :pInterval => [0.65, 0.75]),
                       )


#=
## Defining PDFs for the random variables
userDefinedPdfs = Dict(:SvGrad           =>   Dict(:pdfType => "Uniform", :pInterval => [0.90, 0.93]),
                       :PpGrad           =>   Dict(:pdfType => "Constant", :value => 0.45),
                       :ShminGrad        =>   Dict(:pdfType => "Uniform", :pInterval => [0.63, 0.73] ),                       
                       :μ                =>   Dict(:pdfType => "Uniform", :pInterval  => [0.6, 0.65]),
                       :SHmaxDir         =>   Dict(:pdfType => "Uniform", :pInterval => [65,85]),
                       :surfaceElevation =>   Dict(:pdfType => "Constant", :value => 0),
                       :Aphi             =>   Dict(:pdfType => "Uniform", :pInterval => [0.65, 0.75]),
                       )
=#
    
#=
## Defining PDFs for the random variables
userDefinedPdfs = Dict(:SvGrad           =>   Dict(:pdfType => "Gaussian", :std => 0.05, :μ => 0.92, :truncated => [0.90, 0.93]),
                       :PpGrad           =>   Dict(:pdfType => "Constant", :value => 0.45),
                       :ShminGrad        =>   Dict(:pdfType => "Gaussian", :std => 0.1, :μ => 0.68, :truncated => [0.63, 0.73] ),                       
                       :μ                =>   Dict(:pdfType => "Uniform", :pInterval  => [0.6, 0.65]),
                       :SHmaxDir         =>   Dict(:pdfType => "Uniform", :pInterval => [70,90]),
                       :surfaceElevation =>   Dict(:pdfType => "Constant", :value => 0),
                       :Aphi             =>   Dict(:pdfType => "Uniform", :pInterval => [0.7, 0.8])
                       )
=#

### Original values used Before the meeting with Yueming
#=
## Defining PDFs for the random variables
userDefinedPdfs = Dict(:Sv               =>   Dict(:pdfType => "Gaussian", :stdFraction => 0.05, :truncated => ["Sv_LS", "Sv_HS"]),
                       :Pp               =>   Dict(:pdfType => "Gaussian", :stdFraction => 0.05,  :truncated => ["Pp_LS", "Pp_HS"]),
                       :Shmin            =>   Dict(:pdfType => "Gaussian", :stdFraction => 0.1, :truncated => ["Shmin_LS", "Shmin_HS"] ),                       
                       :μ                =>   Dict(:pdfType => "Uniform", :pInterval  => [0.6, 0.65]),
                       :SHmaxDir         =>   Dict(:pdfType => "Uniform", :pInterval => [80,110]),
                       :surfaceElevation =>   Dict(:pdfType => "Constant", :value => 0),
                       :Aphi             =>   Dict(:pdfType => "Uniform", :pInterval => [0.7, 0.8])
                       )
=#

## plotting histogram of the pdf
FSP3Dplots.PlotPdfInputParameters(mainDir::String, userDefinedPdfs, param)

for fileName in fileNamesAll
    println("\n \t Working on file "*string(fileName)*" \n \n")
    FSP3D.Call3DFSPallFaults(fileName, userDefinedPdfs::Dict{Symbol, Dict{Symbol, Any}}, param::Dict{Symbol, Any})
end



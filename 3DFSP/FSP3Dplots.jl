"""
FSP3D
Author: Josimar Alves da Silva
Copyright 2024 ExxonMobil Technology and Engineering Company
"""
module FSP3Dplots

using StatsBase
using KernelDensity
using UncertaintyQuantification
using Revise
using PyCall
using PyPlot
import FSP3D
import DataFrames
import CSV
import Dierckx

const mcolors = PyNULL()
const cm = PyNULL()

function __init__()
    copy!(mcolors , pyimport("matplotlib.colors"))
    copy!(cm , pyimport("matplotlib.cm"))
end

## The set of functions below are used to create the kernel density based on the histogram
struct KDEDist{D, K} <: ContinuousUnivariateDistribution
    data::D
    kde::K
end

function SetAxisLabelsAndLimits(ax; kwargs=Dict([]))

    ax.set_xlabel(get(kwargs,:xLabel,ax.get_xlabel()))
    ax.set_ylabel(get(kwargs,:yLabel,ax.get_ylabel()))

    ax.set_title(get(kwargs,:figureTitle,ax.get_title()))

    if !haskey(kwargs,:yLim)
        ylims = ax.get_ylim()
        ax.set_ylim(top=1.1*ylims[2])
    else
        ax.set_ylim(get(kwargs,:yLim,ax.get_ylim()))
    end

    ax.set_xlim(get(kwargs,:xLim,ax.get_xlim()))

end ##SetAxisLabelsAndLimits


function PlotMultipleColoredLinesWithColorBar(x::Array{Float64},y::Array{Float64},colorparams::Union{Array{Float64}, Array{Int64}}; kwargs=Dict([]))
    #=
        This function plots multiple colored lines with a colorbar. This is useful when there are many lines
        and the use of legend is impossible.

        Reference: http://csc.ucdavis.edu/~jmahoney/matplotlib-tips.html

        Input:
        colorparams: nLines : array to construct the color of each line
        y: input x and y for each line. Expected size: nLines x nPoints?
        x: nPoints
    =#

    colormap=getproperty( cm,Symbol(get(kwargs,:colorMap,"rainbow")) ) ## if the colormapJS is not set, then use rainbow

    normalize = mcolors.Normalize(vmin=get(kwargs,:vmin,minimum(colorparams)), vmax=get(kwargs,:vmax, maximum(colorparams)))
    if get(kwargs, :logScaleColors, 0) == 1
        normalize = mcolors.LogNorm(vmin=get(kwargs,:vmin,minimum(colorparams)), vmax=get(kwargs,:vmax, maximum(colorparams)))
    end

    if !haskey(kwargs,:ax)
        fig=figure(get(kwargs,:figureLabel,"Figure = "*string(rand(Int8,1)[1])))
        clf()
        ax=gca()
    else
        ax=kwargs[:ax]
        fig=ax.figure
    end
    for i=1:size(y,1)
        color = colormap(normalize(colorparams[i]))
        if get(kwargs, :invertXYaxis, 0) == 0
            ax.plot(x[:], y[i,:], color=color, linewidth=get(kwargs,:linewidth,2), zorder=get(kwargs,:zorder,0))
        else
            ax.plot(y[i,:], x[:], color=color, linewidth=get(kwargs,:linewidth,2), zorder=get(kwargs,:zorder,0))
        end

    end

    # Colorbar setup
    s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
    s_map.set_array(colorparams)

    # If color parameters is a linspace, we can set boundaries in this way
    halfdist = 0
    if length(colorparams) > 1
        halfdist = (colorparams[2] - colorparams[1])/2.0
     end
    boundaries = range(colorparams[1] - halfdist; stop=colorparams[end] + halfdist, length=length(colorparams) + 1)


    if get(kwargs,:cbarStyle,"continuous") == "discrete"
        ## Discrete colorbar
        if get(kwargs, :logScaleColors, 0) == 0
            cbar=fig.colorbar(s_map, spacing="proportional", ticks=get(kwargs,:ticks,colorparams)[:], boundaries=boundaries, format=get(kwargs,:cbarFormat,"%2.1f")) #Integer
        else
            cbar=fig.colorbar(s_map, spacing="proportional", ticks=get(kwargs,:ticks,colorparams)[:], format=get(kwargs,:cbarFormat,"%2.1f")) #Integer
        end
    else
        # Use this to show a continuous colorbar
        #cbar = fig.colorbar(s_map, spacing="proportional", ticks=get(kwargs,:ticks,colorparams)[:], format=get(kwargs,:cbarFormat,"%2.1f"), pad=get(kwargs, :pad, 0.1)) #Float by default
        #cbar = fig.colorbar(s_map, spacing="proportional", ticks=get(kwargs,:ticks,colorparams)[:], format=get(kwargs,:cbarFormat,"%2.2f"), pad=get(kwargs, :pad, 0.05)) #Float by default
        cbar = fig.colorbar(s_map, ax=ax, spacing="proportional", ticks=get(kwargs,:ticks,[500 1500])[:], format=get(kwargs,:cbarFormat,"%2.2f"), pad=get(kwargs, :pad, 0.05)) #Float by default
        #cbar = fig.colorbar(s_map, ax=ax)
    end

    cbar.set_label(get(kwargs,:cbarLabel," "))

    SetAxisLabelsAndLimits(ax; kwargs=kwargs)

    return ax

end #PlotMultipleColoredLinesWithColorBar

function PlotSurfaceOnGrid(dispModel::Matrix{Float64}, xGrid::Vector{Float64}, yGrid::Vector{Float64}; vMin=[], vMax=[], cmapJS=[], kwargs=Dict([]))
    #=
       Plot InSAR and model surface displacement
    =#

    xminJS,xmaxJS=minimum(xGrid), maximum(xGrid)
    yminJS,ymaxJS=minimum(yGrid), maximum(yGrid)

    inputDispModel = dispModel

    if isempty(vMin)
        vMin=minimum(inputDispModel)
    end
    if isempty(vMax)
        vMax=maximum(inputDispModel)
    end
    
    if isempty(cmapJS)
        cmapJS="seismic"
    end
    
    figureName = get(kwargs, :figureName, [])
    fig=figure(get(kwargs, :figureLabel, "NoName"))
    clf()
    ax=gca()
    cb=ax.pcolormesh(xGrid, yGrid, inputDispModel' , cmap=cmapJS, vmin=vMin, vmax=vMax)
    ticks = collect(range(vMin; stop=vMax, length=8))
    cbar=colorbar(cb, ticks=ticks)
    cbar.set_label(get(kwargs, :cbarLabel, "NoLabel"))
    cbar.ax.set_yticklabels(trunc.(ticks; digits=0))

    ax.set_xlim(xminJS, xmaxJS)
    ax.set_ylim(yminJS, ymaxJS)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")    

    ax.set_aspect("equal")

    if !isempty(figureName)
        println("\n \t Saving figure = ", figureName)
        savefig(figureName, format="png", bbox_inches="tight")
    end

    return ax

    
end ## PlotSurfDispInSARandModel

function PlotPdfInputParameters(mainDir::String, userDefinedPdfs::Dict{Symbol, Dict{Symbol, Any}}, param)
    #=
      Plot Pdf of input parameters for quality control
    =#

    figureName=mainDir*"Figures/PdfInputData.png"
    FSP3Dplots.CreateDir(mainDir*"Figures/")

    nSamples = 1000000

    tmp=length(userDefinedPdfs)
    tmp=Int64(ceil(sqrt(tmp) + 1))
    nRows=tmp
    nCols=tmp

    fig=figure("Pdf of input data", figsize=(20,20))
    clf()
    fig.subplots_adjust(hspace=0.5, wspace=0.4)
    subplot(nRows,nCols,1)

    countJS=1
    for key in keys(userDefinedPdfs)
        #global countJS
        if typeof(userDefinedPdfs[key]) != Vector{Float64}
            subplot(nRows,nCols,countJS)
            ax=gca()

            if userDefinedPdfs[key][:pdfType] == "Gaussian"
                if haskey(userDefinedPdfs[Symbol(key)],:std) && haskey(userDefinedPdfs[Symbol(key)],:μ)
                    std=userDefinedPdfs[key][:std]
                    μ=userDefinedPdfs[key][:μ]
                    tmp=Distributions.Normal(μ, std)

                    if haskey(userDefinedPdfs[Symbol(key)],:truncated)
                        ls=userDefinedPdfs[Symbol(key)][:truncated][1]
                        hs=userDefinedPdfs[Symbol(key)][:truncated][2]
                        if hs > ls
                            tmp = Distributions.truncated( tmp , ls , hs )
                        end
                    end

                    xVal=range(μ .* 0.5; stop=μ .* 1.5, length=1000)
                    yVal=pdf.(tmp, xVal)
                    ax.plot(xVal, yVal, color="black")
                    countJS=countJS+1

                elseif haskey(userDefinedPdfs[Symbol(key)],:stdFraction) && haskey(userDefinedPdfs[Symbol(key)],:truncated)
                    depth=(rand(param[:IP3model][Symbol(key)][:data][:,1],1))[1]
                    μ=param[:IP3model][Symbol(key)][:spl](depth)
                    std=userDefinedPdfs[Symbol(key)][:stdFraction] .* μ

                    @assert all(std .> 0 ) "Error!!! Standard deviation is negative or zero: "*string(std)*" for property "*string(key)* " Mean Val = "*string(μ)

                    tmp=Distributions.Normal(μ, std)

                    ls=userDefinedPdfs[Symbol(key)][:truncated][1]
                    hs=userDefinedPdfs[Symbol(key)][:truncated][2]

                    ls=param[:IP3model][Symbol(ls)][:spl](depth)
                    hs=param[:IP3model][Symbol(hs)][:spl](depth)                    
                    
                    if hs > ls
                        tmp = Distributions.truncated( tmp , ls , hs )
                    elseif hs == ls
                        tmp = Distributions.truncated( tmp , ls .* 0.999 , hs .* 1.001 )
                    end

                    xVal=range(ls * 0.5; stop = hs * 1.5, length=1000)
                    yVal=pdf.(tmp, xVal)
                    ax.plot(xVal, yVal, color="black")
                    ax.plot([ls, ls], [0,maximum(yVal)], color="red", linestyle="--", label="LS")
                    ax.plot([hs, hs], [0,maximum(yVal)], color="blue", linestyle="--", label="HS")                    
                    ax.legend(loc="best")
                    countJS=countJS+1

                end

            elseif userDefinedPdfs[key][:pdfType] == "Uniform"
                if haskey(userDefinedPdfs[Symbol(key)],:pInterval)
                    pInterval=userDefinedPdfs[key][:pInterval]

                    tmp=Distributions.Uniform(pInterval[1], pInterval[2])

                    xVal=range(pInterval[1] .* 0.5; stop=pInterval[2] .* 1.5, length=1000)
                    yVal=pdf.(tmp, xVal)
                    ax.plot(xVal, yVal, color="black")
                    countJS=countJS+1
                end

            elseif userDefinedPdfs[key][:pdfType] == "Constant"
                if haskey(userDefinedPdfs[Symbol(key)],:value)
                    value=userDefinedPdfs[Symbol(key)][:value]
                    ax.plot(ones(10) * value, [0 ; ones(9)] , color="black")
                    ax.set_xlim(value .* 0.5, value .* 1.5)
                    countJS=countJS+1
                end

            elseif userDefinedPdfs[key][:pdfType] == "Histogram"
                if haskey(userDefinedPdfs[Symbol(key)],:data)
                    data=userDefinedPdfs[Symbol(key)][:data]
                    ax.hist(data, bins=100, color="black")
                    countJS=countJS+1
                end

            elseif userDefinedPdfs[key][:pdfType] == "HistogramDepth"
                if haskey(userDefinedPdfs[Symbol(key)],:data)
                    tmp=userDefinedPdfs[Symbol(key)][:data]
                    idPlot=convert(Array{Float64},rand((keys(tmp)), 1))
                    ax.hist(tmp[idPlot[1]], bins=100, color="black", label="Depth = "*string(idPlot[1])*" ft")
                    ax.legend(loc="best")

                    countJS=countJS+1
                end
                
            end

            ax.grid(true)

            if string(key) == "SHmaxDir"
                ax.set_title("SHmax dir.")
                ax.set_xlabel("SHmax dir. (degrees from North)")
            elseif string(key) == "SHmax"
                ax.set_title("SHmax")
                ax.set_xlabel("SHmax (psi)")
            elseif string(key) == "Shmin"
                ax.set_title("Shmin")
                ax.set_xlabel("Shmin (psi)")
            elseif string(key) == "Sv"
                ax.set_title("Sv")
                ax.set_xlabel("Sv (psi)")
            elseif string(key) == "Aphi"
                ax.set_title(L"A_\phi")
                ax.set_xlabel("Relative stress magnitude, "*L"A_\phi")
            elseif string(key) == "SvGrad"
                ax.set_title("Sv grad.")
                ax.set_xlabel("Sv grad. (psi/ft)")
            elseif string(key) == "ShminGrad"
                ax.set_title("Shmin grad.")
                ax.set_xlabel("Shmin grad. (psi/ft)")
            elseif string(key) == "SHmaxGrad"
                ax.set_title("SHmax grad.")
                ax.set_xlabel("SHmax grad. (psi/ft)")
            elseif string(key) == "PpGrad"
                ax.set_title("Pp grad.")
                ax.set_xlabel("Pp grad. (psi/ft)")
            elseif string(key) == "μ"
                ax.set_title("Friction coeff.")
                ax.set_xlabel("Friction coeff.")
            elseif string(key) == "faultDip"
                ax.set_title("Fault dip")
                ax.set_xlabel("Fault dip (degrees from horizontal) ")
            elseif string(key) == "faultStrike"
                ax.set_title("Fault strike")
                ax.set_xlabel("Fault strike (degrees from North)")
            elseif string(key) == "surfaceElevation"
                ax.set_title("Surface elevation")
                ax.set_xlabel("Surface elevation (m)")

            end

        end
    end


    savefig(figureName, format="png", bbox_inches="tight")


end #PlotPdfInputParameters

function PlotCummulativeDistributionPcriticalAllDepthColorLinesByDepth(PcriticalFinal::Array{Float64}, zGrid::Array{Float64}; kwargs=Dict([]))
    #=
     Plot the cummulative distribution of hte critical pore pressure to failure vs depth
    =#

    xPlot = collect(range(0; stop=maximum(zGrid), length=1000))
    yPlot = zeros(size(PcriticalFinal,1) * size(PcriticalFinal,2), length(xPlot))

    depthAll = zeros(size(zGrid,1)*size(zGrid,2))
    stepSize=2
    countIndex=1
    indexToKeep=zeros(Int64, size(zGrid,1)*size(zGrid,2))

    for i=1:stepSize:size(PcriticalFinal,1)
        #global k, yPlot, indexToKeep, countIndex
        for j=1:stepSize:size(PcriticalFinal,2)
            data = PcriticalFinal[i,j,:]
            if sum(data[:]) .> 0

                yPlot[countIndex,:] = FSP3D.ComputeCDF(data, xPlot)
                depthAll[countIndex] = zGrid[i,j]
                
                indexToKeep[countIndex] = countIndex
                countIndex=countIndex + 1
            end
        end
    end

    indexToKeep=indexToKeep[1:countIndex-1]
    yPlot=yPlot[1:countIndex - 1,:]
    depthAll=depthAll[1:countIndex - 1]

    indperm = sortperm(depthAll)
    depthAll=depthAll[indperm]
    yPlot=yPlot[indperm,:]

    minValue = minimum(depthAll)
    ticks=Dict(:ticks => range(minValue; stop=maximum(depthAll), length=10))
    
    figure(get(kwargs, :figureLabel, "None"))
    clf()
    ax=gca()
    kwargsDefault=Dict(:ax => ax, :cbarLabel => "Depth (ft)", :xLabel => L"\Delta"*"P for fault slip (psi)")
    merge!(kwargs,ticks)
    merge!(kwargs, kwargsDefault)
    
    ax=FSP3Dplots.PlotMultipleColoredLinesWithColorBar(xPlot,yPlot,depthAll; kwargs=kwargs)

    ax.grid(true)

    figureName = get(kwargs, :figureName, [])
    if !isempty(figureName)
        savefig(figureName, format="png", bbox_inches="tight")
    end

    return ax
    
end #PlotCummulativeDistributionPcriticalAllDepthColorLinesByDepth

function PlotDistributionPcritical(PcriticalFinal::Array{Float64}; kwargs=Dict([]))
    #=
      Plot Pcritical distribution for a single location
    =#
    
    xbin=Int64(round(size(PcriticalFinal,1) / 2; digits=0))
    ybin=Int64(round(size(PcriticalFinal,2) / 2; digits=0))
    
    figureName=get(kwargs, :figureName, [])
    figure("Distribution of Pcritical")
    clf()
    ax=gca()
    ax.hist(PcriticalFinal[xbin,ybin,:], bins=100, color="gray", histtype="stepfilled")
    ax.set_xlabel(L"\Delta"*"P for fault slip (psi)")
    if !isempty(figureName)
        savefig(figureName, format="png", bbox_inches="tight")
    end

end #PlotDistributionPcritical()

function PlotDistributionPcriticalAllDepthColorLineByDepth(PcriticalFinal::Array{Float64}, zGrid::Matrix{Float64}; kwargs=Dict([]))

    xPlot = collect(range(0; stop=maximum(zGrid), length=1000))
    yPlot = zeros(size(PcriticalFinal,1) * size(PcriticalFinal,2), length(xPlot))

    depthAll = zeros(size(zGrid,1)*size(zGrid,2))
    stepSize=2
    countIndex=1
    indexToKeep=zeros(Int64, size(zGrid,1)*size(zGrid,2))
    
    for i=1:stepSize:size(PcriticalFinal,1)
        for j=1:stepSize:size(PcriticalFinal,2)
            data = PcriticalFinal[i,j,:]
            if sum(data[:]) .> 0
                kde = FSP3Dplots.KDEDist(data)
                yPlot[countIndex,:] = Distributions.pdf(kde.kde,xPlot)
                depthAll[countIndex] = zGrid[i,j]
                indexToKeep[countIndex] = countIndex
                countIndex=countIndex + 1
            end
        end
    end
    

    indexToKeep=indexToKeep[1:countIndex-1]
    yPlot=yPlot[1:countIndex-1,:]
    depthAll=depthAll[1:countIndex-1]

    indperm = sortperm(depthAll)
    depthAll=depthAll[indperm]
    yPlot=yPlot[indperm,:]

    minValue = minimum(depthAll)
    ticks=Dict(:ticks => range(minValue; stop=maximum(depthAll), length=10))
    
    figure(get(kwargs, :figureLabel, "None"))
    clf()
    ax=gca()
    kwargsDefault=Dict(:ax => ax, :cbarLabel => "Depth (ft)", :xLabel => L"\Delta"*"P for fault slip (psi)")
    merge!(kwargs,ticks)
    merge!(kwargs, kwargsDefault)
    
    ax=FSP3Dplots.PlotMultipleColoredLinesWithColorBar(xPlot,yPlot,depthAll; kwargs=kwargs)

    ax.grid(true)

    figureName = get(kwargs, :figureName, [])
    if !isempty(figureName)
        savefig(figureName, format="png", bbox_inches="tight")
    end

    return ax


end #PlotDistributionPcriticalAllDepthColorLineByDepth()

function KDEDist(data::Vector{Float64})
    KDE_fit = KernelDensity.kde(data)
    ik = KernelDensity.InterpKDE(KDE_fit)
    return KDEDist(data,ik)
end #KDEDist

#=
function ComputeCDF(inputData::Array{Float64}, x::Union{Array{Float64}, Float64, Int64, Array{Int64}})
    #=
       Compute the CDF for the critical pore pressure to failure
    input:
      inputData: data that will be used to compute the ecdf
      x: value to compute the probability 
    =#

    tmpCDF = StatsBase.ecdf(inputData)
    return tmpCDF.(x)
    
end ##
=#
    
function PlotProbabilityOfFailureUserDefinedPorePressureIncrease(PcriticalFinal::Array{Float64}, grid, PpIncrease::Array{Float64}; kwargs=Dict([]))
    #=
      Plot the probability of failure for a user defined pore pressure increase
    =#

    probFailure=FSP3Dplots.ComputeProbOfFailureUserInputPressureIncrease(PcriticalFinal::Array{Float64}, PpIncrease::Array{Float64})

    inputProperty=probFailure
    tmp=inputProperty[abs.(inputProperty) .< 99]
    ax=FSP3Dplots.PlotSurfaceOnGrid(inputProperty , grid[:xGrid], grid[:yGrid]; vMin=minimum(tmp), vMax=maximum(tmp), cmapJS="jet_r", kwargs = kwargs)

    return Dict(:ax => ax, :probFailure => probFailure)

end #PlotProbabilityOfFailureUserDefinedPorePressureIncrease


function ComputeProbOfFailureUserInputPressureIncrease(PcriticalFinal::Array{Float64}, PpIncrease::Array{Float64})
    #=
      Compute the probability of failure for a user defined pressure increase
    =#

    if size(PcriticalFinal,2) != size(PpIncrease,2)
        throw(ErrorException("Error!!"))
    end

    probFailure = zeros(size(PcriticalFinal,1), size(PcriticalFinal,2))

    for i=1:size(PcriticalFinal,1)
        for j = 1 : size(PcriticalFinal,2)
            if sum(abs.(PcriticalFinal[i,j,:])) > 0 && PpIncrease[i,j] > 0
                probFailure[i,j] = FSP3D.ComputeCDF(PcriticalFinal[i,j,:],PpIncrease[i,j]) .* 100
            end
        end
    end

    return probFailure

end # ComputeProbOfFailureUserInputPressureIncrease


function PlotMaxPressureIncreaseForUserDefinedProbabilityOfFailure(mainDir::String, PcriticalFinal::Array{Float64}, grid, probFailureInput::Union{Int64, Float64}; kwargs=Dict([]))
    #=
        Plot the maximum pore pressure increase acceptable for a given user defined probability of failure
    =#

    zGrid = grid[:zGrid]
    
    ppFailureFinal=zeros(size(zGrid))
    PpAll = collect(range(0; stop=maximum(zGrid), length=1000))

    for i=1:size(zGrid,1)
        for j=1:size(zGrid,2)
            tmp = FSP3D.ComputeCDF(PcriticalFinal[i,j,:],PpAll) .* 100
            ppFailure = findall(tmp .<= probFailureInput)
            if !isempty(ppFailure)
                ppFailureFinal[i,j] = PpAll[ppFailure[end]]
            end
        end
    end

    ## plot The Z Coordinates
    inputProperty=ppFailureFinal
    tmp=inputProperty[abs.(inputProperty) .> 0]
    if !isempty(tmp)
        ax=FSP3Dplots.PlotSurfaceOnGrid(inputProperty , grid[:xGrid], grid[:yGrid]; vMin=minimum(tmp), vMax=maximum(tmp), cmapJS="jet_r", kwargs = kwargs)

        DictFixedProbFailure = Dict(:ax => ax, :ppFailure => ppFailureFinal)

        figureName = get(kwargs, :figureName, [])
        if !isempty(figureName)
            println("\n \t Saving figure = ", figureName)
            savefig(figureName, format="png", bbox_inches="tight")
        end
        return DictFixedProbFailure
    else
        return Dict([])
    end

end #PlotMaxPressureIncreaseForUserDefinedProbabilityOfFailure()

function CreateDir(DirToSave::String)
    ### Check if directory exist and if not then crates a new one

    if isdir(DirToSave) == 0
        println("\n \t Creating directory = ", DirToSave)
        mkpath(DirToSave)
    end

end  #CreateDir


function PlotPDFhistogramsInputTo3DFSP(mainDir::String, finalResults::Dict{Symbol, Array}; kwargs=Dict([]))
    #=
      Plot PDF of inputs used for 3DFSP. It also plot the histogram of the Pcritical
    =#

    inputKey=get(kwargs, :inputKey, keys(finalResults))

    tmp=collect(keys(finalResults))
    
    #iLoc=rand(1:size(finalResults[tmp[1]],1))
    iLoc=get(kwargs, :iLoc, rand(1:size(finalResults[tmp[1]],1)) )

    nKeys=ceil(Int64,sqrt(length( inputKey )))

    figureName=mainDir*"Figures/Histogram_PDFs_inputData.png"
    FSP3D.CreateDir(mainDir*"Figures/")
    
    fig, ax = plt.subplots(nKeys,nKeys, num="Histogram data ", clear = true, figsize=(15,15))
    iKey=1
    fig.subplots_adjust(wspace=0.6, hspace=0.6)

    for key in inputKey
        #global iKey
        if iKey <= length(finalResults)
            println("Plotting on key "*string(key))
            ax[iKey].hist(finalResults[key][iLoc,:], histtype="step", linewidth=3, density=true, bins=get(kwargs, :bins, 15))
            ax[iKey].set_ylabel("Density (-)")
            ax[iKey].set_xlabel(string(key))

            if get(kwargs, :plotDepth, false)
                ax[iKey].set_title( "Depth = "*string(finalResults[:depth][iLoc]) )
            end
            iKey = iKey + 1
        end
    end
    savefig(figureName, format="png", bbox_inches="tight")

    return Dict(:ax => ax, :depth => finalResults[:depth][iLoc])
    
end #PlotPDFhistogramsInputTo3DFSP()

function QCIP3model(iP3results)
    #=
      QC the results of the IP3 interpolation
    =#

    maxDepth=maximum(abs.(iP3results[:Sv][:data][:,1]))
    depthToPlot=-sort(collect(range(0; stop=maxDepth, length=20000)))
    fig, ax = plt.subplots(1,2,num="IP3 model results", clear=true)
    fig.subplots_adjust(wspace=0.5)
    for key in keys(iP3results)
        if key != :units && !occursin("BE", String(key)) && !occursin("HS", String(key)) && !occursin("LS", String(key)) 
            depth=iP3results[key][:data][:,1]
            prop=iP3results[key][:data][:,2]
            
            ax[1].plot(abs.(iP3results[key][:spl].(depthToPlot) ), depthToPlot, linewidth=2, label=string(key)*" interpolated")
            ax[1].scatter(abs.(prop), depth, s=5, marker="o",  label=string(key)*" raw data")

            ax[2].plot(abs.(iP3results[key][:spl].(depthToPlot) ./ depthToPlot), depthToPlot, linewidth=2, label=string(key) * " interpolated")
            ax[2].scatter(abs.(prop ./ depth), depth, s=5, marker="o", label=string(key)*" raw data")
        end
    end
    stressUnits=iP3results[:units][:Sv]
    depthUnits=iP3results[:units][:depth]

    ax[1].grid(true)
    ax[2].grid(true)    
    
    ax[1].legend(loc="best")
    ax[2].legend(loc="best")
    
    ax[1].set_xlabel("Stress ("*string(stressUnits)*")")
    ax[1].set_ylabel("Depth "*string(depthUnits)*")")

    ax[2].set_xlabel("Stress gradient ("*string(stressUnits)*"/"*string(depthUnits)*")")
    ax[2].set_ylabel("Depth "*string(depthUnits)*")")

    ax[2].set_xlim(0,1)

end #QCIP2model

function CreateBlockyCurvesForPlotting(x::AbstractArray,  y1::AbstractArray)
    #=
      Create blocky curves for plotting
    =#

    finalData=[0 0 ]
    for iRow=2:length(x)
        tmp = [x[iRow-1] y1[iRow] ; x[iRow] y1[iRow]]
        finalData = [finalData; tmp]
    end
    finalData = finalData[2:end,:]

    return finalData

end #CreateBlockyCurvesForPlotting

#=
function CreateBlockyCurvesForPlotting(x::AbstractArray,  y1::AbstractArray)
    #=
      Create blocky curves for plotting
    =#

    finalData=[0 0 ]
    for iRow=1:length(x)-1
        tmp = [x[iRow] y1[iRow] ; x[iRow+1] y1[iRow]]
        finalData = [finalData; tmp]
    end
    finalData = finalData[2:end,:]

    return finalData

end #CreateBlockyCurvesForPlotting
=#

function GetIndexIntoSamplesAndDepthToPlotVariationWithDepth(finalResults::Dict{Symbol, Array}; nSamplesToPlot=1000, step=10)

    nSamples=size(finalResults[:Sv],2)
    #nSamplesToPlot=1000

    if nSamplesToPlot > nSamples
        nSamplesToPlot = nSamples
    end
    
    #indSamples=rand(1:nSamples, nSamplesToPlot)
    indSamples=StatsBase.sample(1:nSamples, nSamplesToPlot; replace=false)

    tmp=finalResults[:depth][:,1]
    depthSpace=collect(range(minimum(tmp); stop=maximum(tmp), step=step))
    indMin=zeros(Int64,length(depthSpace))
    for i=1:length(depthSpace)
        minv, indMin[i] = findmin(abs.(depthSpace[i] .- tmp))
    end
    indDepth=indMin

    return indDepth, indSamples
    
end #GetIndexIntoSamplesAndDepthToPlotVariationWithDepth

function PlotIP3modelWith3DFSPresultsVsDepth(mainDir::String, finalResults::Dict{Symbol, Array}, kwargsInterpProp::Dict{Symbol, Dict{Symbol, Any}}; yLim=[], xLim=[], pValues=[5,95], kwargs=Dict([]), wells=[])
    #=
      Plot the IP3 model  and the results of the 3DFSP model
    =#

    ###############
    ## Making sure that the plot selects only a few samples
    indDepth, indSamples = GetIndexIntoSamplesAndDepthToPlotVariationWithDepth(finalResults::Dict{Symbol, Array})
    #depth=finalResults[:depth][indDepth, 1][:]
    
    ###############

    #### Automatized plot here
    keys=["Sv", "SHmax", "Shmin", "Pp"]
    colorJS=["red", "black", "green", "cyan", "gray", "yellow", "brown"]

    fig, ax=plt.subplots(1,2,num="Final QC Stress", clear=true, figsize=get(kwargs, :figSize, (5,5)))
    fig.subplots_adjust(wspace=0.5)

    countJS=1
    for key in keys

        depth=finalResults[:depth][indDepth, 1][:]
        tmpProp=finalResults[Symbol(key)][indDepth,:]
        
        pVal=zeros(size(tmpProp,1),length(pValues))
        for i=1:size(tmpProp,1)
            for j=1:length(pValues)
                pVal[i,j] = StatsBase.percentile(tmpProp[i,:], pValues[j])
            end
        end

        minShmin = pVal[:,1]
        maxShmin = pVal[:,end]

        indPerm=sortperm(depth)
        depth=depth[indPerm]
        minShmin=minShmin[indPerm]
        maxShmin=maxShmin[indPerm]

        #global countJS
        #minShmin = minimum(finalResults[Symbol(key)][indDepth, :]; dims=2)[:]
        #maxShmin = maximum(finalResults[Symbol(key)][indDepth, :]; dims=2)[:]

        ax[1].fill_betweenx(-depth, minShmin, maxShmin, alpha=0.4, color=colorJS[countJS], label=string(key))
        countJS=countJS+1
    end

    if !isempty(kwargsInterpProp)
        countJS=1
        for key in ["Shmin_LS", "Shmin_HS", "Shmin_BE", "Sv_BE", "Pp_BE", "Pp_LS", "Pp_HS"]
            if haskey(kwargsInterpProp, key)
                #global countJS
                if countJS==1
                    ax[1].scatter(abs.(kwargsInterpProp[Symbol(key)][:data][:,2]), kwargsInterpProp[Symbol(key)][:data][:,1], s=5, marker="*", color="black", label="IP3 data")
                else
                    ax[1].scatter(abs.(kwargsInterpProp[Symbol(key)][:data][:,2]), kwargsInterpProp[Symbol(key)][:data][:,1], s=5, marker="*", color="black")
                end
                countJS=countJS+1
            end
        end
    end
    

    ax[1].legend(loc="best")
    ax[1].grid(true)
    ax[1].set_xlabel("Stress (psi)")
    ax[1].set_ylabel("Depth (ft)")
    if !isempty(yLim)
        ax[1].set_ylim(yLim[1], yLim[2])
    end

    if !isempty(wells) 
        for wellName in Base.keys(wells)
            if haskey(wells[wellName], :wellTops)
                FSP3Dplots.PlotWellTop(ax[1], wells[wellName][:wellTops])
                break
            end
        end
    end
    


    #keys=["SvGrad", "SHmaxGrad", "ShminGrad", "PpGrad"]
    keys=["Sv", "SHmax", "Shmin", "Pp"]
    colorJS=["red", "black", "blue", "green", "cyan"]

    #fig, ax[]=plt.subplots(1,1,num="Final QC Stress stress gradient", clear=true)

    countJS=1
    for key in keys
        #global countJS
        depth=finalResults[:depth][indDepth, 1][:]
        tmpProp=finalResults[Symbol(key)][indDepth,:] ./ depth
        
        pVal=zeros(size(tmpProp,1),length(pValues))
        for i=1:size(tmpProp,1)
            for j=1:length(pValues)
                pVal[i,j] = StatsBase.percentile(tmpProp[i,:], pValues[j])
            end
        end

        minShmin = pVal[:,1]
        maxShmin = pVal[:,end]

        indPerm=sortperm(depth)
        depth=depth[indPerm]
        minShmin=minShmin[indPerm]
        maxShmin=maxShmin[indPerm]

        
        #minShmin = minimum(finalResults[Symbol(key)][indDepth, :]; dims=2)[:]
        #maxShmin = maximum(finalResults[Symbol(key)][indDepth, :]; dims=2)[:]

        ax[2].fill_betweenx(-depth[:,1], minShmin , maxShmin , alpha=0.4, color=colorJS[countJS], label=string(key))
        #ax[2].fill_betweenx(-depth[:,1], minShmin ./ depth[:,1], maxShmin ./ depth[:,1], alpha=0.4, color=colorJS[countJS], label=string(key))
        countJS=countJS+1
    end


    if !isempty(kwargsInterpProp)
        countJS=1
        for key in ["Shmin_LS", "Shmin_HS", "Shmin_BE", "Sv_BE", "Pp_BE", "Pp_LS", "Pp_HS"]
            #global countJS
            if haskey(kwargsInterpProp, key)
                tmpDepthToPlot=kwargsInterpProp[Symbol(key)][:data][:,1]
                propToPlot=abs.(kwargsInterpProp[Symbol(key)][:data][:,2])
                gradientToPlot = abs.(propToPlot ./ tmpDepthToPlot)

                if countJS==1
                    ax[2].scatter(gradientToPlot, tmpDepthToPlot, s=5, marker="*", color="black", label="IP3 data")
                else
                    ax[2].scatter(gradientToPlot, tmpDepthToPlot, s=5, marker="*", color="black")
                end
                countJS=countJS+1
            end
        end
    end

    ax[2].legend(loc="best")
    ax[2].grid(true)
    ax[2].set_xlabel("Stress gradient (psi / ft)")
    ax[2].set_ylabel("Depth (ft)")
    if !isempty(yLim)
        ax[2].set_ylim(yLim[1], yLim[2])
    end
    ax[2].set_xlim(0.5, 1.2)

    if !isempty(wells)
        wellNamesTMP=Base.keys(wells)
        for wellName in wellNamesTMP
            if haskey(wells[wellName], :wellTops)
                FSP3Dplots.PlotWellTop(ax[2], wells[wellName][:wellTops])
                break
            end
        end
    end
    


    figureName=mainDir*"Figures/IP3modelComparisonWith3DFSPresults.png"
    FSP3D.CreateDir(mainDir*"Figures/")
    savefig(figureName, format="png", bbox_inches="tight")

    return ax
    
end #PlotIP3modelWith3DFSPresultsVsDepth()

function PlotCriticalPorePressureToFailureVsDepth(mainDir::String, finalResults::Dict{Symbol, Array}; pValues=[5,95], nSamplesToPlot=10000, step=75, yLim=[], xLim=[], stepSizeSpace=100, kwargs=Dict([]), wells = Dict{String, Dict{Symbol, Any}}([]))
    #=
      Plot critical pore pressure to failure vs depth
    =#

    #indDepth, indSamples = GetIndexIntoSamplesAndDepthToPlotVariationWithDepth(finalResults::Dict{Symbol, Array}; nSamplesToPlot=nSamplesToPlot, step=step)
    #Pcritical=finalResults[:Pcritical][indDepth, indSamples]
    #depth=finalResults[:depth][indDepth,1]

    Pcritical=finalResults[:Pcritical]
    depth=finalResults[:depth][:,1]

    #depth, Pcritical = CreateBlockedAverageArrayForPlotting(depth,Pcritical; stepSize=50)
    
    colorJS=["red", "black", "green", "cyan", "gray", "yellow", "brown"]
    fig, ax=plt.subplots(1,1,num="Pore pressure to failure", clear=true, figsize=get(kwargs, :figSize, (5,5)))
    #fig.subplots_adjust(wspace=0.5)

    #pValues=[5,95]
    pVal=zeros(size(Pcritical,1),length(pValues))
    for i=1:size(Pcritical,1)
        for j=1:length(pValues)
            @inbounds pVal[i,j] = StatsBase.percentile(Pcritical[i,:], pValues[j])
        end
    end

    countJS=1
    i=1
    for i=1:size(pVal,2)
        #global countJS
        tmpDepth, tmpProp = FSP3Dplots.CreateBlockedAverageArrayForPlotting(depth,pVal[:,i]; stepSize=stepSizeSpace)
        finalData=FSP3Dplots.CreateBlockyCurvesForPlotting(tmpProp, tmpDepth)
        #finalData=FSP3Dplots.CreateBlockyCurvesForPlotting(pVal[:,i], depth)
        ax.plot(finalData[:,1], -finalData[:,2], linestyle="--", color=colorJS[countJS], label="P"*string(pValues[i]))
        countJS=countJS+1
    end

    #minShmin = minimum(Pcritical; dims=2)[:]
    #maxShmin = maximum(Pcritical; dims=2)[:]

    minShmin = pVal[:,1]
    maxShmin = pVal[:,end]

    #indPerm=sortperm(depth)
    #depth=depth[indPerm]
    #minShmin=minShmin[indPerm]
    #maxShmin=maxShmin[indPerm]

    tmpDepth, tmpProp = CreateBlockedAverageArrayForPlotting(depth,minShmin; stepSize=stepSizeSpace)
    finalData=FSP3Dplots.CreateBlockyCurvesForPlotting(tmpProp, tmpDepth)
    minShmin = finalData[:,1]

    tmpDepth, tmpProp = CreateBlockedAverageArrayForPlotting(depth,maxShmin; stepSize=stepSizeSpace)
    finalData=FSP3Dplots.CreateBlockyCurvesForPlotting(tmpProp, tmpDepth)
    maxShmin = finalData[:,1]

    depth = finalData[:,2]
    
    #finalData=FSP3Dplots.CreateBlockyCurvesForPlotting(maxShmin, depth)
    #maxShmin = finalData[:,1]
    #depth = finalData[:,2]

    ax.fill_betweenx(-depth, minShmin, maxShmin, alpha=0.4)
    #ax.set_xticks(collect(range(0; stop=5000, step=500)))
    ax.grid(true)
    ax.set_ylabel("Depth (ft)")
    ax.set_xlabel("Pore pressure required for fault failure (psi)")
    if !isempty(yLim)
        ax.set_ylim(yLim[1], yLim[2])
    end
    ax.set_xlim(0,5000)
    ax.legend(loc="best")

    if !isempty(wells)
        for wellName in Base.keys(wells)
            if haskey(wells[wellName], :wellTops)
                FSP3Dplots.PlotWellTop(ax, wells[wellName][:wellTops])
            end
        end
    end


    figureName=mainDir*"Figures/PorePressureToFailureVsDepth.png"
    FSP3D.CreateDir(mainDir*"Figures/")
    savefig(figureName, format="png", bbox_inches="tight")


end #PlotCriticalPorePressureToFailureVsDepth()

function PlotHistogramOfFaultDipAndStrike(mainDir::String, faultInfo::Dict{Symbol, Vector{Real}})
    
    figureName=mainDir*"Figures/HistogramFaultDipAndStrike.png"
    FSP3D.CreateDir(mainDir*"Figures")
    faultDip=faultInfo[:faultDip]
    faultStrike=faultInfo[:faultStrike]
    fig, ax = plt.subplots(1,2,num="Histogram fault dip and strike", clear=true, figsize=(15,5))
    ax[1].hist(faultDip, bins=100, density=true, color="blue")
    ax[1].set_xlabel("Fault dip (degrees from horizontal)")
    ax[2].hist(faultStrike, bins=100, density=true, color="blue")
    ax[2].set_xlabel("Fault strike (degrees from North)")
    savefig(figureName, format="png", bbox_inches="tight")
    
end

function CreateBlockedAverageArrayForPlotting(inputDepth::AbstractArray,inputProperty::AbstractArray; stepSize=500)
    #=
      Decimates an array and plot the average value at fixed steps
    =#

    indPerm=sortperm(inputDepth)
    depth=inputDepth[indPerm]

    if size(inputProperty,2) > 1
        tmpInputProperty=inputProperty[indPerm,:]
        nCols=size(Pcritical,2)
    else
        tmpInputProperty=inputProperty[indPerm]
        nCols=1
    end
    

    tmpDepth=collect(range(minimum(depth); stop=maximum(depth), step=stepSize))

    finalProp=zeros(length(tmpDepth)-1, nCols)
    for i=1:size(tmpDepth,1)-1
        ind=findall(depth .>= tmpDepth[i] .&& depth .<= tmpDepth[i+1] )
        finalProp[i,:]=mean(tmpInputProperty[ind, :];dims=1)
    end

    tmpDepth=tmpDepth[1:end-1]
    if size(finalProp,2) == 1
        finalProp=finalProp[:]
    end
    
    return tmpDepth, finalProp 
end #CreateBlockedAverageArrayForPlotting(inputDepth,inputProperty; stepSize=500)

function SetAxisLimitsWellLogs(ax, inputBoundsX::Dict{Symbol, Any}, inputBoundsY::Dict{Symbol, Any})
    #=
     Set bounds for axis plots
    =#

    xlim = inputBoundsX[:wellLogRange]
    xLabel = inputBoundsX[:wellLogName] *" ("* inputBoundsX[:wellLogUnits]*")"
    ax.set_xlim(xlim)
    ax.set_xlabel(xLabel)

    yLim = inputBoundsY[:wellLogRange]
    ax.set_ylim(yLim)
    yLabel = inputBoundsY[:wellLogName] *" ("* inputBoundsY[:wellLogUnits]*")"
    ax.set_ylabel(yLabel)

    
end #SetAxisLimitsWellLogs(bounds::Dict{Symbol, Dict{Symbol, Any}})

function PlotComparisonLogsMultipleWells(mainDir::String, kwargs)

    wells=kwargs[:data]
    wellLogsToCompare=kwargs[:wellLogsToCompare]
    wellNames=keys(wells)
    bounds=kwargs[:bounds]

    figureName=mainDir*"Figures/Compare_ElasticProperties_vs_Depth.png"
    FSP3D.CreateDir(mainDir*"Figures/")
    fig, ax=plt.subplots(1,length(wellLogsToCompare), num=get(kwargs, :figureName, "NoNme"), clear=true, figsize=get(kwargs, :figSize, (17,17)), squeeze=false, sharey=true)

    countJS=1
    for wellLogName in Symbol.(wellLogsToCompare)
        for iWell in wellNames
            df=wells[iWell][:df]
            wellColor=wells[iWell][:color]
            ax[countJS].plot(df[!, Symbol(wellLogName)], -df[!, :TVD], color=wellColor, label=iWell)

            xlim = bounds[wellLogName][:wellLogRange]
            xLabel = bounds[wellLogName][:wellLogName] *" ("* bounds[wellLogName][:wellLogUnits]*")"
            ax[countJS].set_xlim(xlim)
            ax[countJS].set_xlabel(xLabel)

            if haskey(wells[iWell][:srtData], :df)
                if wellLogName == :Shmin
                    PlotSRTdata(ax[countJS], wells[iWell][:srtData][:df])
                elseif wellLogName == :ShminGrad
                    PlotSRTdata(ax[countJS], wells[iWell][:srtData][:df]; plotGradient=true)
                end
            end

            if haskey(wells[iWell], :wellTops) 
                FSP3Dplots.PlotWellTop(ax[countJS], wells[iWell][:wellTops])
            end
            
        end
        ax[countJS].legend(loc="best")
        ax[countJS].grid(true)

        yLim = bounds[:Depth][:wellLogRange]
        ax[countJS].set_ylim(yLim)
        if countJS==1
            yLabel = bounds[wellLogName][:wellLogName] *" ("* bounds[wellLogName][:wellLogUnits]*")"
            ax[countJS].set_ylabel(yLabel)
        end
        
        countJS=countJS + 1
    end

    savefig(figureName, format="png", bbox_inches="tight")

    return ax
    
end #PlotComparisonLogsMultipleWells(kwargs)

function PlotWellLogs(df::DataFrames.DataFrame; nVals = [], kwargs=Dict([]))

    if isempty(nVals)
        nVals=names(df)
    end

    bounds=get(kwargs, :bounds, [])
    
    nPlots=length(nVals)
    fig, ax = plt.subplots(1, nPlots, num=get(kwargs, :figureName, "Well logs vs Depth"), clear=true, squeeze=false, sharey=true)
    countJS=1
    for iName in nVals
        #global countJS
        if !occursin("TVD", iName) && !occursin("DEPT", iName) && !occursin("Column", iName)
            println("Plotting "*string(iName))
            ax[countJS].plot(df[!, iName], -df[!, :TVD], color="black")
            ax[countJS].set_xlabel(bounds[Symbol(iName)][:wellLogName]*" ("*string(bounds[Symbol(iName)][:wellLogUnits])*")")
            if countJS > 1
                #ax[countJS].set_yticks([])
                #ax[countJS].set_yticks(color="w")

            end
            if countJS == 1
                #ax[countJS].set_ylabel("Depth (ft)")
                ax[countJS].set_ylabel(bounds[:Depth][:wellLogName]*" ("*string(bounds[:Depth][:wellLogUnits])*")")
            end
            if !isempty(bounds)
                if haskey(bounds, Symbol(iName))
                    #ax[countJS].set_xlim(bounds[Symbol(iName)][1], bounds[Symbol(iName)][2] )
                    ax[countJS].set_xlim(bounds[Symbol(iName)][:wellLogRange])
                end
            end

            if !isempty(get(kwargs, :wellTops, [])) 
                FSP3Dplots.PlotWellTop(ax[countJS], kwargs[:wellTops])
            end

            ax[countJS].grid(true)
            countJS=countJS+1
        end
    end

    return ax

end #PlotWellLogs(df::DataFrames.DataFrame)


function AverageWellLogsOnInterval(depth, inputVals, intervalSize::Number; pVals=[5,50,95])
    #=
     Takes in a well log as input and then computes the average, P5, P95 and so on a interval.
     This is very useful if one wants to plot well logs with average values on a interval

    It also outputs the histogram of the values within the interval

    depth:: array containing the well log depths
    inputVals:: array containing the property to average on an interval
    intervalSize: size of the interval to average the property
    =#

    tmpDepth=collect(range(minimum(depth); stop=maximum(depth), step=intervalSize))

    finalProperty=zeros(length(tmpDepth)-1, length(pVals))
    finalDepth=zeros(length(tmpDepth)-1)
    histData=Dict([])

    for i=1:length(tmpDepth)-1
        #global finalProperty, finalDepth
        ind=findall(depth .> tmpDepth[i] .&& depth .< tmpDepth[i+1] )
        if !isempty(ind)
            tmpInd = findall(isnan.(inputVals[ind]) .== 0)
            tmpShmin = inputVals[ind[tmpInd]]

            for iP = 1 : length(pVals)
                finalProperty[i,iP] = StatsBase.percentile(tmpShmin, pVals[iP])
            end


            finalDepth[i] = mean([tmpDepth[i] tmpDepth[i+1]])

            histData[finalDepth[i]] = tmpShmin
            
        end
    end

    tmpOutput = Dict(:depth => Dict(:input => depth, :output => finalDepth), :property => Dict(:input => deepcopy(inputVals), :output => deepcopy(finalProperty)), :intervalSize => intervalSize, :pVals => pVals, :histInterval => deepcopy(histData))
    return tmpOutput
    
end #AverageWellLogsOnInterval()

function CreateBulkDensityAndVerticalStress(mainDir::String, wells, bounds; a=0.345, b=1550, c=60)
    #=
      Create bulk density and vertical stresses
    =#
    
    depth=collect(range(0; stop=maximum(abs.(bounds[:Depth][:wellLogRange])), length=1000))
    RHOBfit=RhobVsDepthAnalyticalSolution.(;a=a, b=b, c=c)

    fitFun = Dict(:RHOB => RHOBfit,
                  :Depth => depth,
                  )


    ## Plotting comparison between well logs for different wells
    kwargs=Dict(:data => wells,
                :bounds => bounds,
                :wellLogsToCompare => ["RHOB"]  ,
                :figureName => "RHOB Well log comparison",
                :figSize => (5,5),
                )

    #kwargs[:figureName] = "RHOB Well log comparison"
    #kwargs[:wellLogsToCompare] = ["RHOB"]
    ax=PlotComparisonLogsMultipleWells(mainDir::String, kwargs)

    #ax[1].plot(fitFun[:RHOB](depth), -depth, color="black", label=L"\rho=60z"*L"^{0.345} + 1550")
    ax[1].plot(fitFun[:RHOB](depth), -depth, color="black", label=L"\rho="*string(c)*"z^"*string(a)*" + "*string(b))
    ax[1].legend(loc="best")

    wellNamesTMP=keys(wells)
    for wellName in wellNamesTMP
        if haskey(wells[wellName], :wellTops)
            FSP3Dplots.PlotWellTop(ax[1], wells[wellName][:wellTops])
            break
        end
    end

    figureName=mainDir*"Figures/Density_Interpolated_vs_Depth.png"
    savefig(figureName, format="png", bbox_inches="tight")
    #####

    ### Now Computing the vertical stress gradient and plotting it
    ## Creaging vertical stress gradient
    tmp=fitFun[:RHOB](depth[1:end-1])
    tmp2=diff(depth .* 0.3048)
    Sv = cumsum(tmp .* tmp2 .* 9.8)

    #Sv = cumsum(fitFun[:RHOB](depth[1:end-1]) .* diff(depth ./ 0.3048) .* 9.8)
    Sv = [Sv; Sv[end]]
    spl=FSP3D.CreateInterpolatedFunctionsFromIP3model(depth,Sv)

    fitFun[:Sv] = spl

    #SvGrad = (Sv ./ FSP3D.PsiToPa()) ./ depth
    SvGrad = (fitFun[:Sv](depth) ./ FSP3D.PsiToPa()) ./ depth
    SvGradTmp = 1

    figureName=mainDir*"Figures/Sv_vs_Depth.png"
    fig, ax=plt.subplots(1,2, num="Sv vs depth", clear=true, figsize=(5,5))
    fig.subplots_adjust(wspace=0.4)
    ax[1].plot(fitFun[:Sv](depth), -depth, color="black", label="Interpolated density")
    ax[1].plot(SvGradTmp .* depth .* FSP3D.PsiToPa() , -depth, color="red", label=string(SvGradTmp)*" psi / ft")

    wellLogName=:Sv
    xlim = bounds[wellLogName][:wellLogRange]
    xLabel = bounds[wellLogName][:wellLogName] *" ("* bounds[wellLogName][:wellLogUnits]*")"
    ax[1].set_xlim(xlim)
    ax[1].set_xlabel(xLabel)
    ax[1].grid(true)
    ax[1].legend(loc="best")

    wellNamesTMP=keys(wells)
    for wellName in wellNamesTMP
        if haskey(wells[wellName], :wellTops)
            FSP3Dplots.PlotWellTop(ax[1], wells[wellName][:wellTops])
            break
        end
    end


    ax[2].plot(SvGrad, -depth, color="black", label="Draco")
    wellLogName=:SvGrad
    xlim = bounds[wellLogName][:wellLogRange]
    xLabel = bounds[wellLogName][:wellLogName] *" ("* bounds[wellLogName][:wellLogUnits]*")"
    ax[2].set_xlim(xlim)
    ax[2].set_xlabel(xLabel)
    ax[2].grid(true)

    wellNamesTMP=keys(wells)
    for wellName in wellNamesTMP
        if haskey(wells[wellName], :wellTops)
            FSP3Dplots.PlotWellTop(ax[2], wells[wellName][:wellTops])
            break
        end
    end

    savefig(figureName, format="png", bbox_inches="tight")

    return fitFun

end #CreateBulkDensityAndVerticalStress()

function FindElasticPropertiesAtSRTdepthLocation(df::DataFrames.DataFrame,srtData::DataFrames.DataFrame)
    #=
      Find elastic properties at SRT depth location
    =#
    
    #tmpDepth=wells["Draco"][:df][!, :TVD]
    tmpDepth = df[!,:TVD]

    indAll_ElasticProp=[]
    for i=1:length(srtData.TVD)
        #global indAll_ElasticProp
        depthSRT=srtData.TVD[i]
        minv, indmin = findmin(abs.(depthSRT .- tmpDepth))
        indAll_ElasticProp=[indAll_ElasticProp; indmin]
    end
    
    E=mean(df[indAll_ElasticProp, :E])
    PR=mean(df[indAll_ElasticProp, :PR])

    return E, PR

end #FindElasticPropertiesAtSRTdepthLocation()

function QCUpscaledWellLogsForMonteCarloSampling(mainDir::String, wells::Dict{String, Dict{Symbol, Any}}, bounds; wellNames=[], kwargs=Dict([]))
    #=
      QC the upscaled well logs that will be used in the Monte Carlo sampling
    =#
    #iWell="Draco"
    if isempty(wellNames)
        wellNames=String.(keys(wells))
    end
    
    wellNameTmp=String.(keys(wells))
    nProperties=length(keys(wells[wellNameTmp[1]][:logAverage]))

    fig, ax = plt.subplots(1,nProperties,num="QC Average well logs  ", clear=true, squeeze=false, sharey=true, figsize=get(kwargs, :figSize, (10,10)))

    nWells=1
    for wellName in wellNames
        countJS=1
        for prop in keys(wells[wellName][:logAverage])
            ax[countJS].plot(wells[wellName][:df][!,Symbol(prop)], -wells[wellName][:df][!, :TVD], color="black", alpha=0.2)
            nVals=length(wells[wellName][:logAverage][Symbol(prop)][:pVals])
            for i=1:nVals
                label=wells[wellName][:logAverage][Symbol(prop)][:pVals][i]
                
                if nWells==1 && get(kwargs, :plotPercentile, false)

                    x=wells[wellName][:logAverage][Symbol(prop)][:property][:output][:,i]
                    y=wells[wellName][:logAverage][Symbol(prop)][:depth][:output]
                    finalData=CreateBlockyCurvesForPlotting(x,  y)

                    ax[countJS].plot(finalData[:,1], -finalData[:,2], color=rand(3), label="P"*string(label))
                    
                    #ax[countJS].plot(wells[wellName][:logAverage][Symbol(prop)][:property][:output][:,i], -wells[wellName][:logAverage][Symbol(prop)][:depth][:output], color=rand(3), label="P"*string(label))
                    
                end
                
            end

            if haskey(wells[wellName], :srtData)
                if  haskey(wells[wellName][:srtData], :df) && "Shmin" == String(prop)
                    PlotSRTdata(ax[countJS], wells[wellName][:srtData][:df])
                elseif haskey(wells[wellName][:srtData], :df) && "ShminGrad" == String(prop)
                    PlotSRTdata(ax[countJS], wells[wellName][:srtData][:df]; plotGradient=true)
                end
            end

            SetAxisLimitsWellLogs(ax[countJS], bounds[Symbol(prop)], bounds[:Depth])

            if haskey(wells[wellName], :wellTops)
                FSP3Dplots.PlotWellTop(ax[countJS], wells[wellName][:wellTops])
            end

            ax[countJS].grid(true)
            ax[countJS].legend(loc="best")

            countJS=countJS+1
        end
        nWells=nWells+1
    end

    figureName=mainDir*"Figures/QC_Average_Wells_Logs_WellName.png"
    savefig(figureName, format="png", bbox_inches="tight")

end #QCUpscaledWellLogsForMonteCarloSampling()


function PlotSRTdata(ax, srtData::DataFrames.DataFrame; plotGradient=false)
    x=srtData[!,:Shmin] .* FSP3D.PsiToPa()
    y=srtData[!,:TVD]

    if plotGradient
        x = srtData[!,:ShminGrad]
    end
    
    ax.scatter(x,-abs.(y), s=150, color="cyan", marker="o", label="SRT data", zorder=100)
end #PlotSRTdata(ax, srtData::DataFrames.DataFrame)

function FindStrainsForWellLogToSRTfit(wells, rhobInterp, SHmaxGrad, PpGrad, α)
    #=
      This function first the well log to the SRT data by using the poroelastic model. It first finds the strains that fits the average stress gradient provided by the user, using the elastic properties based on the well log. It then computes the stresses. The output well log stresses should be well calibrated with the SRT data.
    
    =#

    ## Now fitting the poroelastic model to the SRT data
    #SHmaxGrad=0.8
    #PpGrad=0.44
    #α=0.8

    strain_1_All = []
    strain_2_All = []

    wellName = keys(wells)
    E = [] ; PR = []
    for iWell in wellName
        #global E, PR, strain_1_All, strain_2_All
        if haskey(wells[iWell][:srtData], :df)
            if !isempty(wells[iWell][:srtData][:df])
                E, PR = FindElasticPropertiesAtSRTdepthLocation(wells[iWell][:df],wells[iWell][:srtData][:df])

                depthTmp = wells[iWell][:srtData][:df][!,:TVD]
                SvTotal = (rhobInterp[:Sv].(depthTmp))
                SHmax = (SHmaxGrad .* depthTmp .* FSP3D.PsiToPa())
                Shmin = (wells[iWell][:srtData][:df][!, :Shmin] .* FSP3D.PsiToPa())
                Pp = (PpGrad * depthTmp .* FSP3D.PsiToPa())

                strain_1, strain_2 = ComputeStrainMagnitudeUsingStressValue(PR, α, E, Shmin, SHmax, SvTotal, Pp)
                strain_1_All = [strain_1_All ; deepcopy(mean(strain_1))]
                strain_2_All = [strain_2_All ; deepcopy(mean(strain_2))]
            end

        end
    end

    if isempty(strain_1_All)
        throw(ErrorException("Error!! No SRT data found"))
    end

    strain_1 = mean(strain_1_All)
    strain_2 = mean(strain_2_All)

    for iWell in keys(wells)

        depth=wells[iWell][:df][!, :TVD]
        SvTotal=rhobInterp[:Sv].(depth)
        Pp = PpGrad .* depth .* FSP3D.PsiToPa()

        Shmin, SHmax = ComputeShminAndSHmaxMagnitudeIncludingTectonicStrain(wells[iWell][:df][!,:PR], α, wells[iWell][:df][!, :E], mean(strain_1), mean(strain_2), SvTotal, Pp)

        wells[iWell][:df][!, :Sv] = deepcopy(SvTotal)
        wells[iWell][:df][!, :Shmin] = deepcopy(Shmin)
        wells[iWell][:df][!, :SHmax] = deepcopy(SHmax)

        wells[iWell][:df][!, :SvGrad] = (wells[iWell][:df][!, :Sv]  ./ FSP3D.PsiToPa()) ./ depth
        wells[iWell][:df][!, :ShminGrad] = (wells[iWell][:df][!, :Shmin] ./ FSP3D.PsiToPa()) ./ depth
        wells[iWell][:df][!, :SHmaxGrad] = (wells[iWell][:df][!, :SHmax] ./ FSP3D.PsiToPa()) ./ depth

        wells[iWell][:units][!,:Sv] .= "Pa"
        wells[iWell][:units][!,:Shmin] .= "Pa"
        wells[iWell][:units][!,:SHmax] .= "Pa"
        wells[iWell][:units][!,:SvGrad] .= "psi/ft"
        wells[iWell][:units][!,:ShminGrad] .= "psi/ft"
        wells[iWell][:units][!,:SHmaxGrad] .= "psi/ft"
        
    end
    
end #FindStrainsForWellLogToSRTfit()

function QCwellLogToSRTfit_IndividualWells(wells, bounds)
    #=
      QC the fit of the well logs to the SRT data for each well
    =#
    
    ### Plotting the computed Shmin and SHmax
    wellName = keys(wells)
    for iWell in wellName
        kwargs=Dict(:figureName => "Shmin comparison "*string(iWell), :bounds => bounds, :wellTops => wells[iWell][:wellTops])
        ax=PlotWellLogs(wells[iWell][:df]; nVals = ["Shmin"; "SHmax"; "ShminGrad"; "SHmaxGrad"], kwargs=kwargs)
        if haskey(wells[iWell][:srtData], :df)
            PlotSRTdata(ax[1], wells[iWell][:srtData][:df])
            PlotSRTdata(ax[3], wells[iWell][:srtData][:df]; plotGradient=true)
        end
    end
    
end #QCwellLogToSRTfit_IndividualWells()

function ComputeWellLogUpscaling(wells, intervalSize::Number; pVals=[5,50,95], wellNames=[])
    #=
      Upscale the well logs by taking the average value within each interval
    =#
    if isempty(wellNames)
        wellNames=keys(wells)
    else
        wellNames=wellNames
    end

    #intervalSize=100
    logNames=["ShminGrad", "SHmaxGrad", "Shmin", "SHmax"]
    finalLogAverageAllProperty=Dict{Symbol, Dict{Symbol, Any}}([])

    for iLog in logNames
        inputPropertyAll=[]
        depthAll=[]
        for iWell in wellNames
            inputPropertyAll = [ inputPropertyAll ; wells[iWell][:df][!,Symbol(iLog)] ]
            depthAll = [depthAll ; wells[iWell][:df][!, :TVD] ]
        end
        logAverage=AverageWellLogsOnInterval(depthAll, inputPropertyAll, intervalSize; pVals=pVals)
        finalLogAverageAllProperty[Symbol(iLog)] = deepcopy(logAverage)
        for iWell in keys(wells)
            wells[iWell][:logAverage] = deepcopy(finalLogAverageAllProperty)
        end
    end

end #ComputeWellLogUpscaling()

function SaveUpscaledWellLogsToFile(mainDir::String, wells, wellName::String)
    #=
      Save  the upscaled well logs to file
    =#

    logNames=keys(wells[wellName][:logAverage])
    tmp=[]
    for iLog in logNames
        tmpStr=[]
        tmp=[]
        
        dirToSave=mainDir*"Julia/Data/Wells/"*wellName*"/Upscaled/"
        FSP3D.CreateDir(dirToSave)
        fileName=dirToSave*string(iLog)*".csv"

        println("Saving file ", fileName)
    
        tmp1=wells[wellName][:logAverage][iLog][:property][:output]
        tmp2=wells[wellName][:logAverage][iLog][:depth][:output]

        tmp=[tmp2 tmp1]
    
        pVals=wells[wellName][:logAverage][iLog][:pVals]
        tmpStr=[]
        tmpIP3keys=["LS", "BE", "HS"]
        for i=1:length(pVals)
            #tmpStr=[tmpStr; String(iLog)*"_P"*string(pVals[i])]
            tmpStr=[ tmpStr; String(iLog)*"_"*String(tmpIP3keys[i]) ]
        end

        tmpStr=["TVD" ; tmpStr]
        df=DataFrames.DataFrame([ ["ft" "Pa" "Pa" "Pa"] ; tmp], tmpStr)
        #df=DataFrames.DataFrame(tmp, tmpStr)
        CSV.write(fileName, df)

    end
        
    
end #SaveUpscaledWellLogsToFile()



function ReadWellDataForStressCalculation(dirWells::String)
    #=
      Read the well data required for the stress calculation for 3DFSP
    =#

    wells=Dict{String, Dict{Symbol, Any}}([])

    #dirName=mainDir*"Julia/Data/Wells/"
    wellNames=readdir(dirWells)

    for wellName in wellNames
        #global wells

        #dirTmp = mainDir*"Julia/Data/Wells/"*string(wellName)*"/WellLogs/"
        dirTmp = dirWells*string(wellName)*"/WellLogs/"
        if !isdir(dirTmp)
            throw(ErrorException("Error!! Directory not found "*string(dirTmp)))
        end

        fileName=readdir(dirTmp)[1]
        fileName = dirWells*string(wellName)*"/WellLogs/"*string(fileName)

        if !isfile(fileName)
            throw(ErrorException("Error!! File not found!!"))
        end

        #=
        df=CSV.File(fileName; header=1, skipto=3) |> DataFrames.DataFrame
        units=CSV.File(fileName; header=1, skipto=2, limit=1) |> DataFrames.DataFrame

        if units[1,:DTCO] != "microsecods/ft"
            throw(ErrorException("Error!! Well DTCO units should be in microseconds/ft"))
        elseif units[1,:DTSM] != "microsecods/ft"
            throw(ErrorException("Error!! Well DTSM units should be in microseconds/ft"))
        elseif units[1,:RHOB] != "kg/m3"
            throw(ErrorException("Error!! Well RHOB units should be in kg/m3"))
        elseif units[1,:TVD] != "ft"
            throw(ErrorException("Error!! Well TVD units should be in ft"))
        end
        =#

        df, units = ReadWellData(fileName::String)
        
        srtDataDict=Dict{Symbol, Any}([])
        dirTmp = dirWells*string(wellName)*"/SRTdata/"
        if isdir(dirTmp)
            fileName=readdir(dirTmp)[1]
            if !isempty(fileName)
                fileName = dirTmp .* fileName
                srtData = CSV.File(fileName; header=1, skipto=3) |> DataFrames.DataFrame
                srtData[!,:ShminGrad] = srtData[!, :Shmin] ./ srtData[!, :TVD]
                unitsSRT=CSV.File(fileName; header=1, skipto=2, limit=1) |> DataFrames.DataFrame
                srtDataDict = Dict(:df => srtData, :units => unitsSRT, :color => "cyan")
            end
        end

        fileName=dirWells*string(wellName)*"/WellTops/WellTops.csv"
        dfWellTops, unitsWellTops = ReadWellTops(fileName::String)
        
        CreateElasticProperties(df; units=units)
        wells[wellName]=Dict(:df => deepcopy(df), :color => rand(3), :srtData => deepcopy(srtDataDict), :units => units, :wellTops => Dict(:df => deepcopy(dfWellTops), :units => deepcopy(unitsWellTops)))
        
    end

    return wells
    
end #ReadWellDataForStressCalculation()

function CreateElasticProperties(df::DataFrames.DataFrame; units=Dict([]))
    #=

      Create Elastic properties based on the available well logs

    =#

    df[!,:Vp]=(1 ./ df.DTCO) .* (0.3048 ./ 1e-6)
    units[!,:Vp] .= "m/s"

    df[!,:Vs]=(1 ./ df.DTSM) .* (0.3048 ./ 1e-6)
    units[!,:Vs] .= "m/s"    

    #df[!,:RHOB] = df[!, :RHOB] .* 1000

    ν, E, K, G = ComputeElasticProperties(df[!, :Vp], df[!, :Vs], df[!, :RHOB])
    df[!, :PR] = ν
    df[!, :E] = E
    df[!, :K] = K
    df[!, :G] = G

    units[!,:PR] .= ""
    units[!,:E] .= "Pa"
    units[!,:K] .= "Pa"
    units[!,:G] .= "Pa"        
    
end #CreateElasticProperties(df::DataFrames.DataFrame)

function ComputeElasticProperties(Vp, Vs, bulkDensity)
    #=
        Computes elastic properties based on input Vp, Vs and Rhob
    =#

    ## Compute elastic properties from Vp, Vs and bulk density
    G = bulkDensity .* Vs.^2  # shear modulus
    #λ = bulkDensity .* Vp.^2 .- 2 .* G  #Lame parameter
    λ = bulkDensity .* (Vp.^2 .- 2 .* Vs .^ 2)  #Lame parameter    
    K = λ .+ (2/3) .* G # drained bulkd modulus
    E = 9 .* K .* G ./ (3*K .+ G)
    ν = λ ./ ( 2*(λ .+ G) )

    return ν, E, K, G

end #ComputeElasticProperties(Vp, Vs, density)

function ComputeShminAndSHmaxMagnitudeIncludingTectonicStrain(ν, α, E, ϵ1, ϵ2, SvTotal, Pp)
    #=
      Compute the Shmin accounting for tectonic stresses
    =#
    SvTotal = abs.(SvTotal)
    Pp = abs.(Pp)

    Shmin = (ν ./ (1 .- ν)) .* (SvTotal .- α .* Pp) .+ α .* Pp .+ E .* (1 ./ (1 .- ν .^ 2)) .* ϵ1 .+ E .* (ν ./ (1 .- ν .^ 2)) .* ϵ2
    SHmax = (ν ./ (1 .- ν)) .* (SvTotal .- α .* Pp) .+ α .* Pp .+ E .* (1 ./ (1 .- ν .^ 2)) .* ϵ2 .+ E .* (ν ./ (1 .- ν .^ 2)) .* ϵ1    

    
    return Shmin, SHmax
    
end #ComputeShminAndSHmaxMagnitudeIncludingTectonicStrain(ν, α, E, ϵ1, ϵ2, SvTotal, Pp)

function ComputeStrainMagnitudeUsingStressValue(ν, α, E, Shmin, SHmax, SvTotal, Pp)
    #=
      For a given stress value, it computes the amount of strain needed to apply the abaqus boundary conditions
    =#

    #Shmin = ShminGradient .* depth
    #SHmax = SHmaxGradient .* depth

    Shmin = abs.(Shmin)
    SHmax = abs.(SHmax)

    SvTotal = abs.(SvTotal)
    Pp = abs.(Pp)

    A = (ν ./ (1 .- ν)) .* (SvTotal .- α .* Pp) .+ α .* Pp
    B = E .* (1 ./ (1 .- ν .^ 2))
    C = E .* (ν ./ (1 .- ν .^ 2))

    D = ( (Shmin .- A) ./ B ) .- (C ./ (B .^ 2)) .* SHmax .+ C .* A ./ (B .^ 2)
    
    strain_1 = D ./ (1 .- (C .^ 2) ./ (B .^ 2))
    strain_2 = (SHmax .- A .- C .* strain_1) ./ B

    return strain_1, strain_2
    
end #ComputeStrainMagnitudeUsingStressValue()

function RhobVsDepthAnalyticalSolution(;a=0.345, b=1550, c=60)
    #=
      Create density variation with depth for Gulf of Mexico sediments
      depth: depth values should always be in ft
    =#
    rhob(depth) = c .* (depth .* 0.3048) .^ (a) .+ b

    return rhob
    
end #RhobVsDepthAnalyticalSolution()

function SaveInterpolatedVerticalStress(dirWells::String, wells, rhobInterp)
    #=
      Save the interpolatd vertical stress
    =#

    if !isdir(dirWells)
        FSP3Dplots.CreateDir(dirWells)
    end
    
    wellNames=keys(wells)

    depth=range(0; stop=30000, length=1000)

    df=[]
    for wellName in wellNames
        #global df
        dirToSave = dirWells*wellName*"/Upscaled/"
        fileName = dirToSave * "Sv.csv"
        println("Saving file "*string(fileName))
        Sv = rhobInterp[:Sv].(depth)
        df = DataFrames.DataFrame([ ["ft" "Pa" "Pa" "Pa"] ; [depth Sv Sv Sv] ], ["TVD", "Sv_LS", "Sv_BE", "Sv_HS"])
        CSV.write(fileName, df)
    end

end #SaveInterpolatedVerticalStress(dirWells::String, wells)

function ExportWellLogsToFile(dirToSave::String, wells, wellName::String )
    #=
      Export well logs to file
    =#

    #dirToSave=mainDir*"Julia/Data/Wells/"*wellName*"/ComputedStressesFromWellLogs/"
    FSP3D.CreateDir(dirToSave)
    fileNameToSave=dirToSave*wellName*".csv"
    CSV.write(fileNameToSave, [wells[wellName][:units] ; wells[wellName][:df]])

end #ExportComputedStressesUsingWellLogs()

function ReadWellData(dirToRead::Union{String, Array{String}}, wellNames::Union{String, Array{String}})

    wells=Dict{String, Dict{Symbol, Any}}([])
    
    #dirToRead=mainDir*"Julia/Data/Wells/" .* wellNames .* "/ComputedStressesFromWellLogs/"

    fileNames= dirToRead .* wellNames .* ".csv"

    tmp = FSP3Dplots.ReadWellData.(fileNames)

    countJS=1
    for wellName in wellNames
        #global countJS
        
        fileNameWellTops=dirToRead[1] * string(wellName) * "/WellTops/WellTops.csv"
        println("File = ", fileNameWellTops)
        if isfile(fileNameWellTops)
            dfWellTops, unitsWellTops = ReadWellTops(fileNameWellTops::String)
        else
            dfWellTops=Dict(:df => DataFrames.DataFrame([]), :units => DataFrames.DataFrame([]))
        end 

        wells[wellName]=Dict(:df => deepcopy(tmp[countJS][1]), 
                             :color => rand(3), 
                             :units => tmp[countJS][2],
                             :wellTops => dfWellTops
                             )
        countJS=countJS+1
    end

    return wells
    
end #ReadWellData(dirToRead::Union{String, Array{String}}, wellNames::Union{String, Array{String}})

function ReadWellData(fileName::String)

    df=CSV.File(fileName; header=1, skipto=3) |> DataFrames.DataFrame
    units=CSV.File(fileName; header=1, skipto=2, limit=1) |> DataFrames.DataFrame

    if units[1,:DTCO] != "microsecods/ft"
        throw(ErrorException("Error!! Well DTCO units should be in microseconds/ft"))
    elseif units[1,:DTSM] != "microsecods/ft"
        throw(ErrorException("Error!! Well DTSM units should be in microseconds/ft"))
    elseif units[1,:RHOB] != "kg/m3"
        throw(ErrorException("Error!! Well RHOB units should be in kg/m3"))
    elseif units[1,:TVD] != "ft"
        throw(ErrorException("Error!! Well TVD units should be in ft"))
    end

    return df, units
    
end #ReadWellData(fileName::String)

function DefineBoundsForPlotting()
    #=
      define standard bounds for plotting
    =#

    Einfo=Dict(         :wellLogName => "E",        :wellLogUnits => "Pa",     :wellLogRange => (0,50e9))
    PRinfo=Dict(        :wellLogName => "PR",       :wellLogUnits => "",       :wellLogRange => (0,0.5))
    RHOBinfo=Dict(      :wellLogName => "RHOB",     :wellLogUnits => "kg/m3",  :wellLogRange => (0,4000))
    VpInfo=Dict(        :wellLogName => "Vp",       :wellLogUnits => "m/s",    :wellLogRange => (0,5000))
    VsInfo=Dict(        :wellLogName => "Vs",       :wellLogUnits => "m/s",    :wellLogRange => (0,5000))
    ShminInfo=Dict(     :wellLogName => "Shmin",    :wellLogUnits => "Pa",     :wellLogRange => (0,10e7))
    SHmaxInfo=Dict(     :wellLogName => "SHmax",    :wellLogUnits => "Pa",     :wellLogRange => (0,10e7))
    SvInfo=Dict(        :wellLogName => "Sv",       :wellLogUnits => "Pa",     :wellLogRange => (0,10e7))
    ShminGradInfo=Dict( :wellLogName => "ShminGrad",:wellLogUnits => "psi/ft", :wellLogRange => (0,1.5))
    SHmaxGradInfo=Dict( :wellLogName => "SHmaxGrad",:wellLogUnits => "psi/ft", :wellLogRange => (0,1.5))
    SvGradInfo=Dict(    :wellLogName => "SvGrad",   :wellLogUnits => "psi/ft", :wellLogRange => (0,1.5))
    DepthInfo=Dict(     :wellLogName => "TVD",      :wellLogUnits => "ft",     :wellLogRange => (-13000,0))

    bounds=Dict(:RHOB => RHOBinfo,
                :PR => PRinfo,
                :E => Einfo,
                :Vp => VpInfo,
                :Vs => VsInfo,
                :Shmin => ShminInfo,
                :SHmax => SHmaxInfo,
                :Sv => SvInfo,
                :ShminGrad => ShminGradInfo,
                :SHmaxGrad => SHmaxGradInfo,
                :SvGrad => SvGradInfo,
                :Depth => DepthInfo)
    ####

    return bounds
    
end #DefineBoundsForPlotting()

function ReadWellTops(fileName::String)
    #=
     Read well tops for plotting
    =#

    df=CSV.File(fileName; header=1, skipto=3) |> DataFrames.DataFrame
    units=CSV.File(fileName; header=1, skipto=2, limit=1) |> DataFrames.DataFrame

    return df, units
    
end #ReadWellTops(fileName::String)

function ReadWellTops(dirWells::String, wells, wellName::String)
    #=
     Read well tops for plotting
    =#

    fileName = dirWells*wellName*"/WellTops/WellTops.csv"
    if isfile(fileName)
        df=CSV.File(fileName; header=1, skipto=3) |> DataFrames.DataFrame
        units=CSV.File(fileName; header=1, skipto=2, limit=1) |> DataFrames.DataFrame
        wells[wellName][:wellTops][:df] = df
        wells[wellName][:wellTops][:units] = units
    else
        println("File not found: ", fileName)
    end


    return df, units
    
end #ReadWellTops(fileName::String)

function PlotWellTop(ax, wellTops::Dict{Symbol, DataFrames.DataFrame})

    df=wellTops[:df]

    x=ax.get_xlim()
    for name in names(df)
        y=-abs.(df[!,name])
        y=[y,y]
        ax.plot(x, y, linestyle="--")
        ax.text(maximum(x)*0.8, y[1], name)
    end
    
end #PlotWellTop(wellTops::Dict{Symbol, DataFrames.DataFrame})


end ##FSP3Dplots

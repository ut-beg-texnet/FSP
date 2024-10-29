"""
FSP3D
Author: Josimar Alves da Silva
Copyright 2024 ExxonMobil Technology and Engineering Company
"""

module FSP3D

using CSV
using Mmap
using Distributions
import UncertaintyQuantification 
using Test
import DataFrames
using Printf
#import Polynomials
using Statistics
using Revise
using LinearAlgebra
using DelimitedFiles
using Dierckx
using ProgressMeter
using ScatteredInterpolation
using StatsBase
using EmpiricalDistributions
using Interpolations

#import FSP3Dplots

struct TectonicStressStruct
    Sv::Union{Array{Float64}, Int64, Float64}
    SHmax::Union{Array{Float64}, Int64, Float64}
    Shmin::Union{Array{Float64}, Int64, Float64}
    Pp::Union{Array{Float64}, Int64, Float64}
    SHmaxDir::Union{Array{Float64}, Int64, Float64}
    μ::Union{Array{Float64}, Int64, Float64}
    coord::Array{Float64}
end ##TectonicStressStruct

struct SurfaceWithAttributesPetrel 
    coord::Array{Float64}
    grid::Array{Int64}
    attribute::Array{Float64}
    attributeName::Array{AbstractString}
    header::Array{AbstractString}
    depthUnitSystem::AbstractString
    faultDip::Array{Float64}
    faultDipAzimuth::Array{Float64}
    faultStrike::Array{Float64}
end

struct SurfaceFromPetrelEarthVisionFormat
    #=
      Holds horizons surface from Petrel exported as EarthVision format
    =#
    coord::Array{Float64}
    grid::Array{Int64}
    attribute::Array{Float64}
    attributeName::String
    header::Array{AbstractString}
    depthUnitSystem::AbstractString
    fileName::String
    
end #SurfaceFromPetrelEarthVisionFormat


function UnitTest()
    #=
      Unit testing the implementation
    =#

    S = diagm([60, 45, 40])
    S_G = Compute_Sg(S, [0,0,90]);
    @test S_G ≈ [60.0 0 0; 0 40 0; 0 0 45]

    n, ns, nd = ComputeUnitVectors(0, 65);
    @test n ≈ [-0.          0.90630779 -0.42261826]
    @test ns ≈ [1. 0. 0.]
    @test nd ≈ [-0.          0.42261826  0.90630779]


    S = diagm([80, 70, 65])
    strike=30
    dip=-30
    angles=[-70,0,0]
    sigma_n, tau_s, tau_d, tau_mag = ComputeStressComponentsOnFault(S, strike, dip, angles)

    @test sigma_n ≈ 68.67461577598239
    @test tau_s ≈ -0.8550503583141679
    @test tau_d ≈ -6.364621222295629
    @test tau_mag ≈ 6.421799936042004


    ## Testing the principal stresses here
    #=
    S = [50 20 300 ; 20 40 0; 0 0 40]
    S1,S2,S3 = FaultTractions.ComputePrincipalStresses(S)
    println("S1 = "*string(S1)*"\t S2 = "*string(S2)*" \t S3 = "*string(S3))
    @test S1 ≈ 345.7074991
    @test S2 ≈ 39.999999999
    @test S3 ≈ -255.7074991
    =#
    
end ## UnitTest()


function Compute_Sg(S::Union{Array{Float64}, Array{Int64}}, angles::Union{Array{Float64}, Array{Int64}})
    #=
    Compute the stress tensor in geographic coordinates. The input are:
    S: Stress tensor in principal stress (S1, S2, S3)
    angles: angles that need to be set to rotate the principal stresses to the geographic coordinates
    =#
    alpha, beta, gamma = deg2rad.(angles)
    
    Rg = [ cos(alpha) * cos(beta)  sin(alpha) * cos(beta)  -sin(beta);
           cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma) sin(alpha) * sin(beta) * sin(gamma) + cos(alpha) * cos(gamma)  cos(beta) * sin(gamma);
           cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma) sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma) cos(beta) * cos(gamma)]
                  
    #return np.dot(Rg.T, np.dot(S,Rg))
    return Rg' * (S * Rg)

end #Compute_Sg(S::Array{Float64}, angles::Union{Array{Float64}, Array{Int64}})

function ComputeUnitVectors(strike::Real, dip::Real)
    #=
       Compute fault unit vectors
       strike: degrees from north
       dip: degrees from horizontal
    =#
    
    strikeRad = deg2rad(strike)
    dipRad = deg2rad(dip)
    
    n = [-sin(strikeRad) * sin(dipRad) cos(strikeRad) * sin(dipRad) -cos(dipRad) ]
    
    ns = [ cos(strikeRad) sin(strikeRad) 0 ]
    
    nd = [ -sin(strikeRad) * cos(dipRad) cos(strikeRad) * cos(dipRad) sin(dipRad) ]
    
    return n, ns, nd
end #ComputeUnitVectors(strike::Union{Float64, Int64}, dip::Union{Float64, Int64}):

function ComputeStressComponentsOnFault(S::Union{Array{Float64}, Array{Int64}}, strike::Union{Float64, Int64}, dip::Union{Float64, Int64}, angles::Union{Array{Float64}, Array{Int64}})

    S_G = Compute_Sg(S, angles)
    
    n, ns, nd = ComputeUnitVectors(strike, dip)

    sigma_n =  (S_G * n')' * n'
    tau_s =  (S_G * n')' * ns'
    tau_d = (S_G * n')' * nd'

    tau_mag = sqrt(tau_s .^ 2 .+ tau_d .^ 2)

    return sigma_n[1], tau_s[1], tau_d[1], tau_mag[1]

end #ComputeStressComponentsOnFault(S, strike, dip, angles)


function AssemblePrincipalStressesAndAngles(verticalStress::Union{Float64, Int64}, maximumHorStress::Union{Float64, Int64}, minimumHorStress::Union{Float64, Int64}, maximumHorizontalStressDirection::Union{Float64, Int64}, Pp::Union{Float64, Int64})
    #=
      The first step is to figure out the stress regime based on the user defined vertical stress, maximum horizontal stress and minimum horizontal stress

    Input:
         verticalStress: far-field vertical stress
         maximumHorStress: far-field maximum horizontal stress
         minimumHorStress: far-field minimum horizontal stress
         maximumHorizontalStressDirection: azimuth of the maximum horizontal stress direction. It should be measured from north and ranges from 0 to 360

    Output
        S1, S2, S3 : principal stress magnitudes
        α, β, γ: angles required to rotate the principal stresses to the geographic coordinate system

    Recall that:
    
       Normal stress regime:
        Sv > SHmax > Shmin    (S1=Sv, S2=SHmax, Shmin=S3)
    
       Reverse stress regime:
        SHmax > Shmin > Sv    (S1=SHmax, S2=Shmin, S3=Sv)

       Strike-slip stress regime:
        SHmax > Sv > Shmin    (S1=SHmax, S2=Sv, Shmin=Sv)
    =#
    #println("TEST1="*string(verticalStress)*"\t"*string(maximumHorStress)*"\t"*string(minimumHorStress))
    S1, S2, S3 = [] , [] , []
    α , β , γ = [], [], []
    stressRegime=[]
    if abs(verticalStress) > abs(maximumHorStress) && abs(maximumHorStress) > abs(minimumHorStress)
        ## normal stres regime
        S1 = verticalStress
        S2 = maximumHorStress
        S3 = minimumHorStress

        #=
        α = 270 + maximumHorizontalStressDirection    
        β = 270
        γ = 0
        =#
        
        α = 90 + maximumHorizontalStressDirection    
        β = -90
        γ = 0

        stressRegime="Normal"

    elseif abs(maximumHorStress) > abs(minimumHorStress) && abs(minimumHorStress) > abs(verticalStress)
        ## reverse stress regime
        S1 = maximumHorStress
        S2 = minimumHorStress
        S3 = verticalStress

        α = maximumHorizontalStressDirection
        β = 0
        γ = 0

        stressRegime="Reverse"

    elseif abs(maximumHorStress) > abs(verticalStress) && abs(verticalStress) > abs(minimumHorStress)
        ## strike-slip stress regime
        S1 = maximumHorStress
        S2 = verticalStress
        S3 = minimumHorStress

        α = maximumHorizontalStressDirection
        β = 0
        γ = 90

        stressRegime="StrikeSlip"
        
    end

    if verticalStress == 0 && maximumHorStress == 0 && minimumHorStress == 0
        stressParam = Dict(:S1 => 0, :S2 => 0, :S3 => 0, :α => 0, :β => 0, :γ => 0, :Pp => 0)
        return stressParam
    end

    if isempty(S1) || isempty(S2) || isempty(S3) || isempty(α) || isempty(β) || isempty(γ)
        #println("\n \n \t Stress values = Sv="*string(verticalStress)*"\t SHmax="*string(maximumHorStress)*"\t Shmin="*string(minimumHorStress))
        #throw(ErrorException("Error!!! Stress is empty !!!"))
        stressParam = Dict(:S1 => 0, :S2 => 0, :S3 => 0, :α => 0, :β => 0, :γ => 0, :Pp => 0)
        return stressParam
    end

    S1 = -abs(S1) + abs(Pp)
    S2 = -abs(S2) + abs(Pp)
    S3 = -abs(S3) + abs(Pp)

    stressParam = Dict(:S1 => S1, :S2 => S2, :S3 => S3, :α => α, :β => β, :γ => γ, :Pp => Pp, :stressRegime => stressRegime)

    return stressParam

end #AssemblePrincipalStressesAndAngles

function AuxFunParallelFaultTractions(i::Union{Integer, Int64}, verticalStress::AbstractArray, maximumHorStress::AbstractArray, minimumHorStress::AbstractArray, maximumHorizontalStressDirection::Real, Pp::AbstractArray, sigmaN::AbstractArray, tauS::AbstractArray, tauD::AbstractArray, tauMag::AbstractArray, faultDip::AbstractArray, faultStrike::AbstractArray)
    
    stressParam = AssemblePrincipalStressesAndAngles(verticalStress[i], maximumHorStress[i], minimumHorStress[i], maximumHorizontalStressDirection, abs(Pp[i]))

    S=diagm([stressParam[:S1], stressParam[:S2] , stressParam[:S3]])
    angles = [stressParam[:α], stressParam[:β], stressParam[:γ]]
        
    sigmaN[i], tauS[i], tauD[i], tauMag[i] = ComputeStressComponentsOnFault(S, faultStrike[i], faultDip[i], angles)
    
end #AuxFunParallelFaultTractions


function ComputeFaultTractions(verticalStress::Array{Real}, maximumHorStress::Array{Real}, minimumHorStress::Array{Real}, Pp::Array{Real}, maximumHorizontalStressDirection::Real, faultDip::Array{Real}, faultStrike::Array{Real})

    if length(verticalStress) != length(faultDip)
        throw(ErrorException("Error!!!"))
    end

    sigmaN=zeros(Real,length(faultDip))
    tauS=deepcopy(sigmaN)
    tauD = deepcopy(tauS)
    tauMag = deepcopy(tauS)

    @Threads.threads for i=1:length(faultDip)
        
        AuxFunParallelFaultTractions(i, verticalStress, maximumHorStress, minimumHorStress, maximumHorizontalStressDirection, Pp, sigmaN, tauS, tauD, tauMag, faultDip, faultStrike)
        
    end


    faultTractions=Dict(:sigmaN => sigmaN, :tauS => tauS, :tauD => tauD, :tauMag => tauMag)

    return faultTractions
    
end #ComputeFaultTractions()

function ComputeFaultTractions(verticalStress::Real, maximumHorStress::Real, minimumHorStress::Real, Pp::Real, maximumHorizontalStressDirection::Real, faultDip::Real, faultStrike::Real)

    if length(verticalStress) != length(faultDip)
        throw(ErrorException("Error!!!"))
    end

    sigmaN = zeros(1)
    tauS = deepcopy(sigmaN)
    tauD = deepcopy(tauS)
    tauMag = deepcopy(tauS)

    i=1
    AuxFunParallelFaultTractions(i, [verticalStress], [maximumHorStress], [minimumHorStress], maximumHorizontalStressDirection, [Pp], sigmaN, tauS, tauD, tauMag, [faultDip], [faultStrike])

    faultTractions=Dict(:sigmaN => sigmaN, :tauS => tauS, :tauD => tauD, :tauMag => tauMag)

    return faultTractions
    
end #ComputeFaultTractions()

function TestImplementationFaultTractionsWithPetrel(mainDir::String, testNumber::Integer)

    if testNumber == 1
        @testset "Testing Rotation Matrix" begin FaultTractions.UnitTest() end
    end

    ## Read fault information
    dirData = mainDir*"Julia/Data/BarrowIsland/ResultsFromPetrel/"
    if testNumber == 1
        fileName = dirData * "EM_Microseismic_Pink_SHmax_Dir_N30E_NormalStress"
    elseif testNumber == 2
        fileName = dirData * "EM_Microseismic_Pink_SHmax_Dir_N20E_ReverseStress"
    elseif testNumber == 3
        fileName = dirData * "EM_Microseismic_Pink_SHmax_Dir_N40W_StrikeSlipStress"
    else
        throw(ErrorException("Errro!!! There only 3 tests to try"))
    end
    
    fault = FSP3D.ReadPetrelPointsWithAttributes(fileName)

    faultStrike=fault.faultAzimuth .+ 90
    faultDip= fault.faultDip

    SvGrad = 1 ## psi / ft
    ShminGrad = 0.6 # psi / ft
    SHmaxGrad = 0.8 ## psi / ft
    PpGrad = 0.44 # psu / ft

    depth = 2000 * 3.28 # ft

    if testNumber == 1
        ## Normal stress regime, N30E
        maximumHorizontalStressDirection = 30 ## azimuth angle from  North
        verticalStress = SvGrad .* depth
        maximumHorStress = SHmaxGrad .* depth
        minimumHorStress = ShminGrad .* depth
        Pp = PpGrad .* depth
    end
    
    if testNumber == 2
        ## Reverse stress regime, N20E
        maximumHorizontalStressDirection = 20 ## azimuth angle from  North
        maximumHorStress = SvGrad .* depth
        minimumHorStress = SHmaxGrad .* depth
        verticalStress = ShminGrad .* depth
        Pp = PpGrad .* depth
    end


    if testNumber == 3
        ## Strike-slip stress regime, N40W
        maximumHorizontalStressDirection = 310 ## azimuth angle from  North
        maximumHorStress = SvGrad .* depth
        verticalStress  = SHmaxGrad .* depth
        minimumHorStress = ShminGrad .* depth
        Pp = PpGrad .* depth
    end


    ## Stresses should always be of the same size as the fault dip and fault strike arrays
    verticalStress = ones(length(faultDip)) * verticalStress
    maximumHorStress = ones(length(faultDip)) * maximumHorStress
    minimumHorStress = ones(length(faultDip)) * minimumHorStress
    Pp = ones(length(faultDip)) * Pp

    faultTractions = ComputeFaultTractions(verticalStress, maximumHorStress, minimumHorStress, Pp, maximumHorizontalStressDirection, faultDip, faultStrike)

    sigmaN = faultTractions[:sigmaN]
    tauMag = faultTractions[:tauMag]

    sigmaN = sigmaN ./ FaultTractions.BarToPsi()
    tauMag = tauMag ./ FaultTractions.BarToPsi()


    nBins=20
    figure("Shear and effective normal stress ", figsize=(20,15))
    clf()
    subplot(221)
    ax=gca()
    ax.hist(sigmaN,bins=nBins)
    ax.set_title("Min = "*string(minimum(sigmaN))*"; Max = "*string(maximum(sigmaN)))
    ax.set_xlabel("Normal stress (bar)")

    subplot(222)
    ax=gca()
    ax.hist(tauMag,bins=nBins)
    ax.set_xlabel("Shear stress (bar)")
    ax.set_title("Min = "*string(minimum(tauMag))*"; Max = "*string(maximum(tauMag)))

    subplot(223)
    ax=gca()
    ax.hist(tauMag ./ sigmaN,bins=nBins)
    ax.set_title("Min = "*string(minimum(tauMag ./ sigmaN))*"; Max = "*string(maximum(tauMag ./ sigmaN)))
    ax.set_xlabel("Slip tendency")

    subplot(224)
    ax=gca()
    ax.hist(faultStrike,bins=nBins)
    ax.set_title("Min = "*string(minimum(faultStrike))*"; Max = "*string(maximum(faultStrike)))
    ax.set_xlabel("Fault strike")



    tauMag
    Ts = tauMag ./ sigmaN

    figureName = mainDir*"Figures/OneToOnePlots_Comparison_Julia_Petrel.png"
    figure("1:1 plots; test number = "*string(testNumber), figsize=(20,20))
    clf()
    subplot(221)
    ax=gca()
    ax.scatter(fault.slipTendency[:], Ts[:], s=5, marker="o", color="black", facecolor="none")
    ax.plot([0 , 1], [0 , 1], color="blue", linewidth=2, zorder=0)
    ax.set_title("Slip tendency")
    ax.set_aspect("equal")
    ax.set_xlabel("Petrel results")
    ax.set_ylabel("Julia results")

    minV = minimum(tauMag) .* 0.5
    maxV = maximum(tauMag) .* 1.5

    subplot(222)
    ax=gca()
    ax.scatter(fault.shearStress[:], tauMag[:], s=5, marker="o", color="black", facecolor="none")
    ax.plot([minV , maxV ], [minV , maxV], color="blue", linewidth=2, zorder=0)
    ax.set_aspect("equal")
    ax.set_title("Shear stress (bar)")
    ax.set_xlabel("Petrel results")
    ax.set_ylabel("Julia results")

    minV = minimum(sigmaN) .* 0.5
    maxV = maximum(sigmaN) .* 1.5

    subplot(223)
    ax=gca()
    ax.scatter(fault.normalStress[:], sigmaN[:], s=5, marker="o", color="black", facecolor="none")
    ax.plot([minV , maxV ], [minV , maxV], color="blue", linewidth=2, zorder=0)
    ax.set_title("Normal stress (bar)")
    ax.set_aspect("equal")
    ax.set_xlabel("Petrel results")
    ax.set_ylabel("Julia results")
    savefig(figureName, format="png", bbox_inches="tight")

    figureName = mainDir*"Figures/ShearAndNormalStressComparison_Julia_Petrel.png"
    figure("Comparison shear and normal stress; Test Number = "*string(testNumber))
    clf()
    ax=gca()
    ax.scatter(sigmaN, tauMag, s=15, marker="o", color="black", facecolor="black", label="Julia")
    ax.scatter(fault.normalStress, fault.shearStress, s=15, marker="s", color="red", facecolor="none", label="Petrel")
    legend(loc="best")
    ax.set_xlabel("Normal stress (bar)")
    ax.set_ylabel("Shear stress (bar)")
    savefig(figureName, format="png", bbox_inches="tight")

    
end #TestImplementationFaultTractionsWithPetrel()

function ReadPetrelPointsWithAttributes(fileName::String)
    #=
        Read petrel points with attributes for the fault tractions
         X
         Y
         Z
         INT,Fault framework: fault
         STRING,Fault framework: fault
         DOUBLE,Fault sample area
         DOUBLE,Dip direction
         DOUBLE,Dip
         DOUBLE,Slip tendency
         DOUBLE,Critical pore fluid pressure change
         DOUBLE,Distance to failure
         DOUBLE,Shear stress magnitude: far-field stress
         DOUBLE,Normal stress: far-field stress

    =#

    colX = []
    colY = []
    colZ = []
    colFaultDip = []
    colFaultDipDirection = []
    colNormalStress = []
    colShearStress = []
    colDistanceToFailure = []
    colSlipTolerance = []
    colCriticalPoreFluid = []
    colSlipTendency = []
    colDilationTendency = []
    colDistanceToFailure = []
    countCol=1
    flag = 0
    countLines = 1
    open(fileName,"r") do io

        #global colX, countCol, colY, colZ, flag, colFaultDip, colFaultDipDirection, colNormalStress, colShearStress, colDistanceToFailure, colSlipTolerance, colCriticalPoreFluid, colSlipTendency, colDilationTendency, countLines

        for ln in eachline(io)

            if occursin("END", ln)
                rowToStartReading = countLines + 1
                break
            end
            countLines = countLines + 1

            if ln[1] == 'X' && length(ln) == 1
                colX = countCol
                flag = 1
            end
            if ln[1] == 'Y' && length(ln) == 1
                colY = countCol
            end

            if ln[1] == 'Z' && length(ln) == 1
                colZ = countCol
            end

            if occursin("Dip", ln) && !occursin("direction", ln)
                colFaultDip = countCol
            end

            if occursin("Dip direction", ln) 
                colFaultDipDirection = countCol
            end

            if occursin("Normal stress", ln) 
                colNormalStress = countCol
            end

            if occursin("Shear stress", ln) 
                colShearStress = countCol
            end

            if occursin("Distance to failure", ln)
                colDistanceToFailure = countCol
            end

            if occursin("Slip tolerance", ln)
                colSlipTolerance = countCol
            end

            if occursin("Critical pore fluid pressure change", ln)
                colCriticalPoreFluid = countCol
            end

            if occursin("Slip tendency", ln)
                colSlipTendency = countCol
            end

            if occursin("Dilation tendency", ln)
                colDilationTendency = countCol
            end

            if flag == 1
                countCol = countCol + 1
            end
        end

    end

    dataInput=readdlm(fileName, skipstart=countLines)
    data=DataFrames.DataFrame()

    if !isempty(colX)
        data.X = dataInput[:,colX]
    end

    if !isempty(colY)
        data.Y = dataInput[:,colY]
    end

    if !isempty(colZ)
        data.Z = dataInput[:,colZ]
    end

    if !isempty(colFaultDip)
        data.faultDip = dataInput[:,colFaultDip]
    end

    if !isempty(colFaultDipDirection)
        data.faultAzimuth = dataInput[:,colFaultDipDirection]
    end

    if !isempty(colNormalStress)
        data.normalStress = dataInput[:,colNormalStress]
    end

    if !isempty(colShearStress)
        data.shearStress = dataInput[:,colShearStress]
    end

    if !isempty(colDistanceToFailure)
        data.distanceToFailure = dataInput[:,colDistanceToFailure]
    end

    if !isempty(colSlipTendency)
        data.slipTendency = dataInput[:,colSlipTendency]
    end

    if !isempty(colCriticalPoreFluid)
        data.criticalPorePressure = dataInput[:,colCriticalPoreFluid]
    end

    return data
    
end ## ReadPetrelPointsWithAttributes()

function BarToPa()
    ##= convert bar to Pa
    return 100000
end ## BarToPa

function BarToPsi()
    ##= convert bar to Psi
    return 14.5038 
end ## BarToPa

function PsiToPa()
    ##= convert Psi to Pa
    return 6894.76
end ## BarToPa

function ComputeCriticalPorePressureForFailure(sigma_n::Union{Array{Real}, Real}, tau::Union{Array{Real}, Real}; mu = 0.6)
    #=
        Compute the critical pore pressure necessary for failure
    =#

    mu1 = abs(tau ./ sigma_n)
    mu2 = mu
    Pcritical = ((mu2 .- mu1) ./ mu2) .* abs(sigma_n)

    return Pcritical
end #ComputeCriticalPorePressureForFailure

#=
function MainComputePcritical(df::DataFrames.DataFrame, modelType::String; kwargs=Dict([]))
    #=
      Computes the critical pore pressure to failure
      Input:
          Dataframe containing the gradients, fault information etc...

      Output:
         Pore pressure to failure
    =#

    ## Fault friction coefficient
    μ = df.FrictionCoeff 
    
    if isempty(modelType)
        throw(ErrorException("Error!! Model type has not been defined. Available options are: \n 1) StressGradients and 2) CriticallyStressed."))
    end

    if modelType == "CriticallyStressed"
        depth = df.Depth
        xCoord = df.Xcoord
        yCoord = df.Ycoord        

        ## Getting values of the spatially interpolated properties
        ## Note that at this specific (x,y) location, all xCoord and yCoord are the same size. Therefore, it
        ## only requires to compute the interpolated value once
        propInterp = FSP3D.GetValueOfSpatiallyInterpolatedProperty(df, kwargs, xCoord[1], yCoord[1])

        ##If interpolated
        if haskey(propInterp, :SvGrad)
            Sv = propInterp[:SvGrad][1] .* abs.(depth)
        else
            if any(occursin.("SvGrad", names(df)))
                Sv = df.SvGrad .* abs.(depth)
            else
                throw(ErrorException("Error!! No value found for SvGrad"))
            end
        end

        if haskey(propInterp, :Aphi)
            Aphi = propInterp[:Aphi][1] .* ones(length(depth)) 
        else
            if any(occursin.("Aphi", names(df)))
                Aphi = df.Aphi
            else
                throw(ErrorException("Error!! No value found for Aphi"))
            end
        end

        if haskey(propInterp, :SHmaxDir)
            SHmaxDir = propInterp[:SHmaxDir][1] .* ones(length(depth))
        else
            if any(occursin.("SHmaxDir", names(df)))
                SHmaxDir = df.SHmaxDir
            else
                throw(ErrorException("Error!! No value found for SHmaxDir"))
            end
        end

        if haskey(propInterp, :PpGrad)
            Pp = propInterp[:PpGrad][1] .* abs.(depth)
        else
            if any(occursin.("PpGrad", names(df)))
                Pp = df.PpGrad .* abs.(depth)
            else
                throw(ErrorException("Error!! No value found for PpGrad"))
            end
        end

        #Pp = df.PpGrad .* abs.(depth)

        SHmax=zeros(length(Sv))
        Shmin=zeros(length(Sv))
        @Threads.threads for i=1:length(Sv)
            tmp = FSP3D.ComputeStressTensor_CS_Model.(Sv[i], Pp[i], μ[i], Aphi[i])
            SHmax[i] = tmp.SHmax[1]
            Shmin[i] = tmp.Shmin[1]
        end

        if all(SHmax .== 0)
            throw(ErrorException("Error!!! SHmax is all zeros!"))
        end

        if all(Shmin .== 0)
            throw(ErrorException("Error!!! Shmin is all zeros!"))
        end

    elseif modelType == "StressGradients"
        
        throw(ErrorException("Error!! Not implemented yet!!! "))
        
    else
        
        throw(ErrorException("Error!! Model type has not been defined. Available options are: \n 1) StressGradients and 2) CriticallyStressed."))
        
    end
    

    #=
    SvFunction = get(kwargs, :SvFunction, [])
    if typeof(SvFunction) == Spline1D
        verticalStress = SvFunction(abs.(df.Depth))
    else
        verticalStress = df.SvGrad .* abs.(df.Depth)
        #verticalStress = [] ## testing 
    end

    if any(occursin.("Aphi", names(samples)))  ##Critically stressed model
        
    else
        
    end
    
    SvFunction = get(kwargs, :SvFunction, [])
    if typeof(SvFunction) == Spline1D
        verticalStress = SvFunction(abs.(df.Depth))
    else
        verticalStress = df.SvGrad .* abs.(df.Depth)
        #verticalStress = [] ## testing 
    end
    maximumHorStress = df.SHmaxGrad .* abs.(df.Depth) 
    minimumHorStress = df.ShminGrad .* abs.(df.Depth) 
    Pp = df.PpGrad .* abs.(df.Depth) 
    SHmaxDir = df.SHmaxDir
    =#

    faultStrike = df.FaultStrike
    faultDip = df.FaultDip

    Pcritical = zeros(length(Sv))
    #println("TEST0 = ", typeof(Pcritical))
    @Threads.threads for i=1:length(Sv)
    #for i=1:length(verticalStress)
        FSP3D.AuxFunctionComputePcriticalSample(i, Pcritical, Sv, SHmax, Shmin, SHmaxDir, Pp, faultStrike, faultDip, μ)
    end

    #return Pcritical
    return Pcritical, mean(Sv) ./ depth[1], mean(SHmax) ./ depth[1], mean(Shmin) ./ depth[1], mean(SHmaxDir), mean(Pp) ./ depth[1], mean(μ), mean(Aphi)
    
end #MainComputePcritical
=#

#=
function MainComputePcritical(df::DataFrames.DataFrame; kwargs=Dict([]))
    #=
      Computes the critical pore pressure to failure
      Input:
          Dataframe containing the gradients, fault information etc...

      Output:
         Pore pressure to failure
    =#

    #PpGradPdf = Distributions.Uniform(0.43,0.45)    
    #SvGradPdf = Distributions.Uniform(0.97,1.03)
    #SHmaxGradPdf = Distributions.Uniform(1.2, 1.4)
    #ShminGradPdf = Distributions.Uniform(0.7, 0.9)
    #maximumHorizontalStressDirectionPdf = Distributions.Uniform(0, 90)

    #faultDipPdf = Distributions.Uniform(0,360)
    #faultStrikePdf = Distributions.Uniform(0,360)
    #frictionCoeffPdf = Distributions.Uniform(0.6,0.9)

    #faultDip = rand(faultDipPdf,1)[1]
    #faultStrike = rand(faultStrikePdf,1)[1]
    #frictionCoeff = rand(frictionCoeffPdf,1)[1]
    #SvGrad = rand(SvGradPdf,1)[1]
    #SHmaxGrad = rand(SHmaxGradPdf,1)[1]
    #ShminGrad = rand(ShminGradPdf,1)[1]
    #PpGrad = rand(PpGradPdf,1)[1]
    #maximumHorizontalStressDirection = rand(maximumHorizontalStressDirectionPdf,1)[1]

    SvFunction = get(kwargs, :SvFunction, [])
    if typeof(SvFunction) == Spline1D
        verticalStress = SvFunction(abs.(df.Depth))
    else
        verticalStress = df.SvGrad .* abs.(df.Depth)
        #verticalStress = [] ## testing 
    end
    #println("TEST0 = ", typeof(verticalStress))
    #verticalStress = df.SvGrad .* abs.(df.Depth)
    maximumHorStress = df.SHmaxGrad .* abs.(df.Depth) 
    minimumHorStress = df.ShminGrad .* abs.(df.Depth) 
    Pp = df.PpGrad .* abs.(df.Depth) 
    SHmaxDir = df.SHmaxDir
    faultStrike = df.FaultStrike
    faultDip = df.FaultDip
    frictionCoeff = df.FrictionCoeff


    Pcritical = zeros(length(verticalStress))
    #println("TEST0 = ", typeof(Pcritical))
    @Threads.threads for i=1:length(verticalStress)
    #for i=1:length(verticalStress)
        FSP3D.AuxFunctionComputePcriticalSample(i, Pcritical, verticalStress, maximumHorStress, minimumHorStress, SHmaxDir, Pp, faultStrike, faultDip, frictionCoeff)
    end

    return Pcritical
    
end #MainComputePcritical
=#

function AuxFunctionComputePcriticalSample(i::Integer, Pcritical::Array{Float64}, verticalStress::Array{Float64}, maximumHorStress::Array{Float64}, minimumHorStress::Array{Float64}, SHmaxDir::Array{Float64}, Pp::Array{Float64}, faultStrike::Array{Float64}, faultDip::Array{Float64}, frictionCoeff::Array{Float64}, sigmaN::Array{Float64}, tauMag::Array{Float64})    
    #=
       Aux fuinction to compute the critical pore pressure to failure in paralle
    =#

    stressParam = FSP3D.AssemblePrincipalStressesAndAngles(verticalStress[i], maximumHorStress[i], minimumHorStress[i], SHmaxDir[i], Pp[i])
    S=diagm([stressParam[:S1], stressParam[:S2] , stressParam[:S3]])
    angles = [stressParam[:α], stressParam[:β], stressParam[:γ]]

    sigmaN[i], tauS, tauD, tauMag[i] = FSP3D.ComputeStressComponentsOnFault(S, faultStrike[i], faultDip[i], angles)
    @inbounds Pcritical[i] = FSP3D.ComputeCriticalPorePressureForFailure(sigmaN[i], tauMag[i]; mu = frictionCoeff[i])

end #AuxFunctionComputePcriticalSample()


function ComputePcriticalToSlipAllFaultAzimuthsAndDip(stressParam::Dict{Symbol, Real}, frictionCoeff::Real; stepSize=5, stopStrike=[], stopDip=[])
    #=
      For a given stress tensor defined by the principal stresses (S1, S2, S3), computes the pore pressure to slip
    for faults of all azimuths and all dip
      Outputs results on a grid that is ready to be plotted on a polar plot

    Input:
       stressParam: dictionary containing the information regarding the stress tensor. See Example below

    #=
    k=1
    stressParam = FaultTractions.AssemblePrincipalStressesAndAngles(verticalStress[k], maximumHorStress[k], minimumHorStress[k], maximumHorizontalStressDirection, Pp[k])

    S=diagm([stressParam[:S1], stressParam[:S2] , stressParam[:S3]])
    angles = [stressParam[:α], stressParam[:β], stressParam[:γ]]
    =#

    =#

    if isempty(stopStrike)
        stopStrike=360
    end
    if isempty(stopDip)
        stopDip=360
    end

    S=diagm([stressParam[:S1], stressParam[:S2] , stressParam[:S3]])
    angles = [stressParam[:α], stressParam[:β], stressParam[:γ]]

    faultDipAll = collect(range(0; stop=stopDip, step=stepSize))
    faultStrikeAll = collect(range(0; stop=stopStrike, step=stepSize))
    
    
    grid = FSP3D.CreateMeshGrid(faultStrikeAll, faultDipAll, []; stepSizeX=stepSize, stepSizeY=stepSize)
    xGrid = grid[:xGrid]
    yGrid = grid[:yGrid]

    Pcritical = zeros(length(xGrid), length(yGrid) )
    @Threads.threads for i=1:length(xGrid)
        FSP3D.AuxFunctionParallelComputePcritical(Pcritical, xGrid, yGrid, S, angles, frictionCoeff, i)
    end

    return Dict(:Pcritical => Pcritical, :xGrid => xGrid, :yGrid => yGrid)

end ## ComputePcriticalToSlipAllFaultAzimuthsAndDip()

function AuxFunctionParallelComputePcritical(Pcritical::Array{Real}, xGrid::Array{Real}, yGrid::Array{Real}, S::Array{Real}, angles::Array{Real}, frictionCoeff::Real, i::Integer)
    
    for j=1:length(yGrid)
        faultDip = yGrid[j]
        faultStrike = xGrid[i]
        sigmaN, tauS, tauD, tauMag = FSP3D.ComputeStressComponentsOnFault(S, faultStrike, faultDip, angles)
        @inbounds Pcritical[i,j] = FSP3D.ComputeCriticalPorePressureForFailure(sigmaN, tauMag; mu = frictionCoeff)
    end
    
end #AuxFunctionParallelComputePcritical(Pcritical, xGrid, yGrid, S, faultStrike, faultDip, angles, frictionCoeff, i)

function CreateMeshGrid(xCoord::Vector{Float64}, yCoord::Vector{Float64}, stepSize::Integer; minX=[], maxX=[], minY=[], maxY=[], stepSizeX=[], stepSizeY=[], xGrid=[], yGrid=[])    
    #=
    Create grid to display InSAR data
    Create mesh grid that can be used to plot matrix plots using pcolormesh for example
    =#
    minValX = minimum(xCoord)
    if !isempty(minX)
        minValX = minX
    end
    maxValX = maximum(xCoord)
    if !isempty(maxX)
        maxValX = maxX
    end
    minValY = minimum(yCoord)
    if !isempty(minY)
        minValY = minY
    end
    maxValY = maximum(yCoord)
    if !isempty(maxY)
        maxValY = maxY
    end

    if isempty(stepSizeX)
        stepSizeX=stepSize
    end

    if isempty(stepSizeY)
        stepSizeY=stepSize
    end

    if isempty(xGrid) && isempty(yGrid)
        xGrid=collect(range(minValX; stop=maxValX, step=stepSizeX))
        yGrid=collect(range(minValY; stop=maxValY, step=stepSizeY))
    end

    xGridFinal = xGrid' .* ones(length(yGrid))
    yGridFinal = ones(length(xGrid))' .* yGrid

    grid=Dict(:xGrid => xGrid, :yGrid => yGrid, :xGridMatrix => xGridFinal, :yGridMatrix => yGridFinal)

    return grid
    
end ##CreateMeshGrid

function ComputePCriticalOnGrid(zGrid::Matrix{Float64}, inputParam::Dict{Symbol, Union{UncertaintyQuantification.RandomVariable, Integer}}; kwargs=Dict([]))
    #=
      For each location in the grid, computes the critical pore pressure to failure

    Input: zGrid which is a grid of size nRows vs nCols creatred using the ComuteInSARDipOnGri function
    inputParam: dictionary containing input parameters required to create the MCMC 

    Output: Critical pore pressure to failure
    =#

    PpGrad = inputParam[:PpGrad]
    SvGrad = inputParam[:SvGrad]
    SHmaxGrad = inputParam[:SHmaxGrad]
    ShminGrad = inputParam[:ShminGrad]
    SHmaxDir = inputParam[:SHmaxDir]
    faultDip = inputParam[:faultDip]
    faultStrike = inputParam[:faultStrike]
    frictionCoeff = inputParam[:frictionCoeff]
    nSamples = inputParam[:nSamples]
    
    PcriticalFinal = zeros(size(zGrid,1), size(zGrid,2), nSamples)

    ## Compute the critical pore pressure to slip for each location in the spatial grid
    @showprogress for i=1:size(zGrid,1)
        for j=1:size(zGrid,2)
            depth = UncertaintyQuantification.Parameter(abs.(zGrid[i,j]), :depth)
            if abs(zGrid[i,j]) > 0
                samples = UncertaintyQuantification.sample( [depth, PpGrad, SvGrad, ShminGrad, SHmaxGrad, SHmaxDir, faultDip, faultStrike, frictionCoeff] , nSamples )
                @inbounds PcriticalFinal[i,j,:] = FSP3D.MainComputePcritical(samples; kwargs=kwargs)
            end
        end
    end

    return PcriticalFinal

end #ComputePCriticalOnGrid(zGrid)

function AuxFunctionsParallel(dispTmp::Matrix{Float64}, i::Integer, xGrid::Vector{Float64}, yGrid::Vector{Float64}, disp::Vector{Float64}, xAll::Vector{Float64}, yAll::Vector{Float64})

    x1=xGrid[i]
    x2=xGrid[i+1]
    tmpDispAll = disp[xAll .>= x1 .&& xAll .< x2]
    yAllRef = yAll[xAll .>= x1 .&& xAll .< x2]
    
    if isempty(tmpDispAll) == false
        @inbounds for j=1:length(yGrid)-1
            y1=yGrid[j]
            y2=yGrid[j+1]
            #tmpDisp = disp[xAll .>= x1 .&& xAll .< x2 .&& yAll .>= y1 .&& yAll .< y2]
            tmpDisp = tmpDispAll[ yAllRef .>= y1 .&& yAllRef .< y2]            
            if isempty(tmpDisp) == false
                dispTmp[i,j] = mean(tmpDisp)
                #dispTmp[i,j] = minimum(tmpDisp)
            end
        end
    end
    
end

function ComputeS1FrictionStrength(S3::Array{Float64}, Pp::Array{Float64}, μ::Real)
    #=
      Computes the magnitude of the S1 principal stress based on the S3 value. This implementes equation 4.46 from Zoback
    =#

    factor = ( sqrt.( (μ .^ 2 .+ 1) ) + μ ) .^ 2
    S1 = (S3 .- Pp) .* factor .+ Pp

    return S1
    
end ## ComputeS1FrictionStrength()


function ComputePropertyOnGrid(xGrid::Vector{Float64}, yGrid::Vector{Float64}, disp::Vector{Float64}, xAll::Vector{Float64}, yAll::Vector{Float64})
    #=
       For the mesh grid defined by xGrid and yGrid, it computes the average of the samples in disp and then
    output it. The goal of this function is to create a grid that can be ploted using pcolormesh
    =#

    dispTmp = zeros(length(xGrid), length(yGrid))
    
    #@inbounds Threads.@threads for i=1:length(xGrid)-1
    @inbounds for i=1:length(xGrid)-1        
        FSP3D.AuxFunctionsParallel(dispTmp, i,xGrid, yGrid, disp, xAll, yAll)
    end
    return dispTmp
end

function ComputeGridBasedOnHorizonData(coord::Matrix{Float64}, spatialResolution::Real)
    #=
       Using the coordinates provided by the user input horizon, it computes the grid that can be used in the pcolormesh plot and in the calculation of the critical pore pressure to failure

    Input:
         coord: coordinates of the used defined horizon
         spatialResolution: spatial resolution of the grid

    Output
       grid: dictionary containing the grid information that can be plot    
    =#
    
    ## Creating grid
    grid = FSP3D.CreateMeshGrid(coord[:,1], coord[:,2], spatialResolution; minX=minimum(coord[:,1]), maxX=maximum(coord[:,1]), minY=minimum(coord[:,2]), maxY=maximum(coord[:,2]))

    propToGrid = -coord[:,3]
    zGrid = FSP3D.ComputePropertyOnGrid(grid[:xGrid], grid[:yGrid], propToGrid, coord[:,1], coord[:,2])
    grid[:zGrid]  = zGrid
    
    return grid

end #MainComputeGridBasedOnHorizonData()

function ComputeMohrCircle(sigma::Matrix{Float64})
    #=
         Compute Mohr circle
    =#
    
    S1 = sort(sigma[:]; rev = true)[1]
    S2 = sort(sigma[:]; rev = true)[2]
    S3 = sort(sigma[:]; rev = true)[3]

    theta = range(0; stop = 180, length = 180)

    sigmaN_12 = 0.5 .* (S1 .+ S2) .+ 0.5 .* (S1 .- S2) .* cosd.(theta)
    tau_12 = 0.5 .* (S1 .- S2) .* sind.(theta)

    sigmaN_23 = 0.5 .* (S2 .+ S3) .+ 0.5 .* (S2 .- S3) .* cosd.(theta)
    tau_23 = 0.5 .* (S2 .- S3) .* sind.(theta)

    sigmaN_13 = 0.5 .* (S1 .+ S3) .+ 0.5 .* (S1 .- S3) .* cosd.(theta)
    tau_13 = 0.5 .* (S1 .- S3) .* sind.(theta)

    mohrCircle = Dict(:sigmaN12 => sigmaN_12, :tau12 => tau_12, :sigmaN23 => sigmaN_23, :tau23 => tau_23 , :sigmaN13 => sigmaN_13 , :tau13 => tau_13, :theta => theta, :S1 => S1, :S2 => S2, :S3 => S3)

    return mohrCircle

end #ComputeMohrCircle

function ConvertFaultAzimuthToNorthCentered(inputAzimuth::Union{Array{Float64}, Array{Int64}})
    #=
      Convert fault azimuth to a north centered system. For example, if the fault azimuth is 120 degrees, from the North, then the convertion would yield -60, which corresponds to N60W. This is useful if I want to plot an histogram of the fault azimuth to fit a pdf for example.

    =#

    #=
    newAzimuth = deepcopy(inputAzimuth)

    for k=1:100
        for i=1:length(newAzimuth)

            if newAzimuth[i] .> 360
                newAzimuth[i] =  newAzimuth[i] .- 360
            end

           if newAzimuth[i] .>= 90 && newAzimuth[i] <= 180
               newAzimuth[i] =  newAzimuth[i] .- 180
           end

           if newAzimuth[i] .> 180 && newAzimuth[i] <= 270
               newAzimuth[i] = newAzimuth[i] .- 180
           end

           if newAzimuth[i] .> 270 && newAzimuth[i] <= 360
               newAzimuth[i] =  newAzimuth[i] .- 360
           end
        end
    end

    return newAzimuth
    =#

    
end #ConvertFaultAzimuthToNorthCentered()

function ComputePrincipalStresses(S::Matrix{Float64})
    #=
     This function computes the principal stresses based on the stress tensor. The equation that I am using can be found on page 28 of Jaeger's book

      Input: stress tensor (size 9 x 9)
    =#

    #=
    if S[1,2] != S[2,1] || S[1,3] != S[3,1] || S[2,3] != S[3,2]
        throw(ErrorException("Errro!!! Input matrix must be symmetric!"))
    end
    =#

    #=
    epsJS=1e-5
    if (S[1,2] - S[2,1]) > epsJS || (S[1,3] - S[3,1]) > epsJS || (S[2,3] - S[3,2]) > epsJS

        #println("S = "*string(S))
        #throw(ErrorException("Errro!!! Input matrix must be symmetric!"))
        
    end
    =#

    eg = LinearAlgebra.eigen(S)
    S1,S2,S3 = eg.values

    return eg.values, eg.vectors
    
end #ComputePrincipalStresses

function ConstructStressTensorJS(S::Matrix{Float64})
    # function to construct the stress tensor for the principal stress calculation
    #
    if length(S) > 6
        throw(ErrorException("Erro!!! Input to stress tensor needs to be defined by 6 elements only"))
    end
    
    sigma = [S[1] S[4] S[5] ; S[4] S[2] S[6]; S[5] S[6] S[3]] 
end # ConstructStressTensorJS

function Main_ComputePcritical(userDefinedPdfs::Dict{Symbol, Distributions.Distribution{Univariate, Continuous}}, nSamples::Real, gridJS::Dict{Symbol, Array{Float64}} )

    zGrid=gridJS[:zGrid]
    
    #println("\n \t Computing critical pore pressure to failure for a grid containing "*string(length(zGrid))*" locations and using "*string(nSamples)*" samples for each location. Total number of FSP evaluations is equal to "*string(length(zGrid)*nSamples)*"\n")
    @printf("\n \t Computing critical pore pressure to failure for a grid containing %G locations and using %G samples for each location. Total number of FSP evaluations is equal to %0.0E. See progress meter below for total computation time. \n \n", length(zGrid), nSamples, length(zGrid)*nSamples)    

    ShminGrad = UncertaintyQuantification.RandomVariable(userDefinedPdfs[:ShminGrad], :ShminGrad)
    SHmaxGrad = UncertaintyQuantification.RandomVariable(userDefinedPdfs[:SHmaxGrad], :SHmaxGrad)
    SHmaxDir = UncertaintyQuantification.RandomVariable(userDefinedPdfs[:SHmaxDir], :SHmaxDir)
    faultStrike = UncertaintyQuantification.RandomVariable(userDefinedPdfs[:FaultAzimuth], :FaultStrike)
    frictionCoeff = UncertaintyQuantification.RandomVariable(userDefinedPdfs[:FrictionCoeff], :FrictionCoeff)
    faultDip = UncertaintyQuantification.RandomVariable(userDefinedPdfs[:FaultDip], :FaultDip)
    PpGrad = UncertaintyQuantification.RandomVariable(userDefinedPdfs[:PpGrad], :PpGrad)
    SvGrad = UncertaintyQuantification.RandomVariable(userDefinedPdfs[:SvGrad], :SvGrad)

    kwargsPc = Dict([])

    inputParam = Dict{Symbol, Union{UncertaintyQuantification.RandomVariable, Integer}}(:PpGrad => PpGrad, :SvGrad => SvGrad, :SHmaxGrad => SHmaxGrad, :ShminGrad => ShminGrad, :SHmaxDir => SHmaxDir, :FaultDip => faultDip, :FaultStrike => faultStrike, :FrictionCoeff => frictionCoeff, :nSamples => nSamples)

    PcriticalFinal = FSP3D.ComputePCriticalOnGrid(zGrid, inputParam; kwargs=kwargsPc)

    println("\n \t Finished computing critical pore pressure to failure")
    
    return PcriticalFinal

end #Main_ComputePcritical()

function ReadFaultSurfaceWithAttributeHeader(fileName::AbstractString)
    #=
      Read header of surface file exported from Petrel. The file was exported as "Surfaces with attributes"
    =#

    if !isfile(fileName)
        throw(ErrorException("Error!! File not found!! "*string(fileName)))
    end
    
    ## Trying to read the header of the file
    header=AbstractString[]
    open(fileName, "r") do io
        for ln in eachline(io)
            #global header
            if occursin("Type: scattered data", ln)
                ## It means that this is an Horzions (EarthVision format) and not a surface
                break
            end
            header = [header ; ln]
            #println("test = ", ln)
            if occursin("END ATTRIBUTES", ln)
                break
            end
        end
    end

    return header
    
end #ReadFaultSurfaceWithAttributeHeader

function ReadFaultSurfaceWithAttributesFromPetrel(fileName::AbstractString)
    #= 
      Read  Petrel surface with attributes. The typical attributes expected here are the fault azimuth and the fault dip, but it can also read other attributes if available
    =#

    #Read file header
    header = FSP3D.ReadFaultSurfaceWithAttributeHeader(fileName)
    
    ## Writting function to read data from file
    coord=[]; gridJS=[]; attributes=[]
    countLine=1
    open(fileName, "r") do io
        for ln in eachline(io)
            #global countLine
            if occursin("END ATTRIBUTES", ln)
                break
            end
            countLine=countLine + 1
        end
    end
    tmp = readdlm(fileName; skipstart=countLine)
    coord = tmp[:,1:3]
    gridJS = Int64.(tmp[:,4:5])
    attributes = tmp[:,6:end]

    attInfo=FSP3D.GetAttributeNameFromHeader(header)
    depthUnits=GetUnitsOfDepthAxisBasedOnSurfaceAttributesHeader(header)

    if !isempty(attInfo[:faultDipIndex])
        faultDip=attributes[:,attInfo[:faultDipIndex]]
    else
        faultDip=zeros(size(coord,1))
    end
    
    if !isempty(attInfo[:faultDipAzimuthIndex])
        faultDipAzimuth=attributes[:,attInfo[:faultDipAzimuthIndex]]
        faultStrike=faultDipAzimuth .+ 90 ## Petrel outpus the fault dip angle, need to conver to fault strike

        ##Now convert the fault strike to degrees from North
        #faultStrike = FSP3D.ConvertFaultAzimuthToNorthCentered(faultStrike)
    else
        faultDipAzimuth=zeros(size(coord,1))
        faultStrike=zeros(size(coord,1))        
    end

    attributeName=attInfo[:attributeName]
    
    surface=FSP3D.SurfaceWithAttributesPetrel(coord, gridJS, attributes, attributeName, header, depthUnits, faultDip, faultDipAzimuth, faultStrike)
    
    return surface

end #ReadFaultSurfaceWithAttributesFromPetrel


function GetAttributeNameFromHeader(header::Vector{AbstractString})
    #=
      For a given surface with attributes in Petrel, it gets the attribute name anf check if fault dip and azimuth are present or not
      head: Input header from the surface with attributes
    =#
    
    ind1=findall(header .== "ATTRIBUTES")[1]
    ind2=findall(header .== "END ATTRIBUTES")[1]

    if isempty(ind1) || isempty(ind2)
        throw(ErrorException("Error!!! ATTRIBUTE not found!!"))
    end
    
    attributeName=AbstractString[]
    attributeIdNumber=[]
    faultDipIndex=Int64[]
    faultDipAzimuthIndex=Int64[]
    countJS=1
    for i=ind1+1:ind2-1
        tmpName = split(header[i],",")[end]
        tmpIdNumber = split(header[i],",")[1]

        attributeName=[attributeName ; tmpName]
        attributeIdNumber=[attributeIdNumber ; tmpIdNumber]

        if string(tmpIdNumber) == "7"
            faultDipIndex=countJS
        elseif string(tmpIdNumber) == "8"
            faultDipAzimuthIndex=countJS
        end
        countJS=countJS+1
    end

    return Dict(:attributeName => attributeName, :faultDipIndex => faultDipIndex, :faultDipAzimuthIndex => faultDipAzimuthIndex)
    
end #GetAttributeNameFromHeader(header::Vector{String})


function WriteSurfaceWithAttributesToFile(fileNameToSave::AbstractString, surface::SurfaceWithAttributesPetrel; additionalAttributes=Matrix{Real}(undef,0,0), additionalAttributeNames=Vector{String}(undef,0))
    #=
      Write surface with attributes to file. The output from this function can be read by Petrel
    =#

    header = surface.header
    
    if isempty(additionalAttributes)
        attribute = surface.attribute
        attributeName = surface.attributeName
    else
        attribute = [surface.attribute additionalAttributes]
        if isempty(additionalAttributeNames)
            throw(ErrorException("Erro!!! You need to specify names for all attributes" ))
        else
            if size(additionalAttributes,2) != length(additionalAttributeNames)
                throw(ErrorException("Erro!!! You need to specify names for all attributes" ))
            else
                attributeName = [surface.attributeName[:]; additionalAttributeNames[:]]
            end
        end
    end

    #println("\n TEST1 = "*string(size(attributeName))*"\t TEST2 = "*string(size(attribute)))
    coord = surface.coord
    grid = surface.grid
    
    open(fileNameToSave, "w") do io
    for i=1:length(header)
        write(io, strip(header[i])*"\n")
        #if header[i] == "END HEADER"
        if occursin("END HEADER", header[i] )
            if !isempty(attribute)
                write(io,"ATTRIBUTES\n")
                for k=1:size(attribute,2)
                    str="FLOAT,NN,General,"*string(attributeName[k])*"\n"
                    write(io,str)
                end
                write(io, "END ATTRIBUTES\n")
                break
            end
        end
    end

    for i=1:size(coord,1)

        str=[]
        for k=1:size(coord,2)
            if k==1
                str=string(coord[i,k])*" "
            else
                str=str * string(coord[i,k])*" "
            end
        end
        for k=1:size(grid,2)
            str=str * string(grid[i,k])*" "
        end

        for k=1:size(attribute,2)
            str=str * string(attribute[i,k])*" "
        end

        if i < size(coord,1)
            str=str*"\n"
        end

        write(io, str)
    end

end #open

end #WriteSurfaceWithAttributesToFile()

function GetUnitsOfDepthAxisBasedOnSurfaceAttributesHeader(header::Vector{AbstractString})
    #=
      Get units of the depth axis based on the header read from the surface with attributes
    =#
    finalUnits=[]
    for i=1:length(header)
        #global finalUnits
        if occursin("Unit in depth", header[i])
            finalUnits=strip(split(header[i],":")[2])
            break
        end
    end

    return finalUnits
    
end #GetUnitsOfDepthAxisBasedOnSurfaceAttributesHeader()



function BuildDictionaryOfUserDefinePdsForSampling(userDefinedPdfs)
    #=
       Build dictionary containing the user defined PDFs and constant values for sampling
    =#

    finalDict=Dict{Symbol, UncertaintyQuantification.RandomVariable}([])
    for key in keys(userDefinedPdfs)
        if typeof(userDefinedPdfs[key]) == Vector{Float64} || typeof(userDefinedPdfs[key]) == Float64 || typeof(userDefinedPdfs[key]) == Int64 || typeof(userDefinedPdfs[key]) == Matrix{Float64}
            tmp = UncertaintyQuantification.Parameter.(userDefinedPdfs[key], key)
            finalDict=merge(finalDict, Dict(key => tmp))
            
        elseif typeof(userDefinedPdfs[key]) == String
            finalDict=merge(finalDict, Dict(key => userDefinedPdfs[key]))
            
        else
            tmp = UncertaintyQuantification.RandomVariable(userDefinedPdfs[key], key)
            finalDict=merge(finalDict, Dict(key => tmp))
        end
    end

    return finalDict

end #BuildDictionaryOfUserDefinePdsForSampling(userDefinedPdfs)

function CreateDataFrameOfSamplesForMonteCarloSimulation(finalUserPdfs, i::Integer, nSamples::Integer)
    #=
      Create dataframe with random samples extracted from the user Pdfs

    i:: For this specific location, check if it is a PDF or aa constant value
    =#
    
    ## Build container to hold user defined PDFs for sampling
    tmp=Union{ UncertaintyQuantification.Parameter, UncertaintyQuantification.RandomVariable }[]
    for key in keys(finalUserPdfs)
        #global tmp
        if typeof(finalUserPdfs[key]) == Array{UncertaintyQuantification.Parameter,1}
            tmp = [tmp; finalUserPdfs[key][i]]  ## constant value assigned to this location
            
        elseif typeof(finalUserPdfs[key]) != String
            tmp = [tmp; finalUserPdfs[key]]  ## pdf assigned to this location            
        #else
            #tmp = [tmp; finalUserPdfs[key]]  ## pdf assigned to this location
        end
    end

    samples=UncertaintyQuantification.sample(tmp, nSamples)

    return samples
    
end #CreateDataFrameOfSamplesForMonteCarloSimulation

function CreateDataFrameOfSamplesForMonteCarloSimulation(df::DataFrames.DataFrame, nSamples::Integer)
    #=
      Create dataframe with random samples extracted from the user Pdfs

    i:: For this specific location, check if it is a PDF or aa constant value
    =#
    
    ## Build container to hold user defined PDFs for sampling
    tmp=Union{ UncertaintyQuantification.Parameter, UncertaintyQuantification.RandomVariable }[]
    for name in names(df)
        #global tmp
        if typeof(finalUserPdfs[key]) == Array{UncertaintyQuantification.Parameter,1}
            tmp = [tmp; finalUserPdfs[key][i]]  ## constant value assigned to this location
            
        elseif typeof(finalUserPdfs[key]) != String
            tmp = [tmp; finalUserPdfs[key]]  ## pdf assigned to this location            
        #else
            #tmp = [tmp; finalUserPdfs[key]]  ## pdf assigned to this location
        end
    end

    samples=UncertaintyQuantification.sample(tmp, nSamples)

    return samples
    
end #CreateDataFrameOfSamplesForMonteCarloSimulation


function ComputePcCriticalUserDefinedProbabilityOfFailure(PcriticalFinal::Matrix{Float64}, probFailureAll::Union{Vector{Int64}, Vector{Float64}})
    #=
     Compute critical pore pressure to failure for user defined set of probability thresholds
    =#

    ## Areas where depth == 0 have NaN stresses
    indTmp=findall( isnan.(PcriticalFinal) .== 1 )
    #println("\n Found "*string(length(indTmp))*" NaN values in Pore Pressure to failure. Replacing it with zeros.")
    PcriticalFinal[indTmp] .= 0
    
    ppFailureFinal=zeros(size(PcriticalFinal,1), length(probFailureAll))

    PpAll = collect(range(0; stop=maximum(PcriticalFinal), length=1000))

    for k=1:length(probFailureAll)
        probFailureInput=probFailureAll[k]
        for i=1:size(PcriticalFinal,1)
            tmp = FSP3D.ComputeCDF(PcriticalFinal[i,:],PpAll) .* 100
            ppFailure = findall(tmp .<= probFailureInput)
            if !isempty(ppFailure)
                ppFailureFinal[i,k] = PpAll[ppFailure[end]]
            end
        end
    end

    return ppFailureFinal

end #ComputePcCriticalUserDefinedProbabilityOfFailure()

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

#=
function ComputeCriticalPorePressureToFailure(surface::FSP3D.SurfaceWithAttributesPetrel, userDefinedPdfs, nSamples::Int64, modelType::String; kwargs=Dict([]) )
    
    finalUserPdfs = FSP3D.BuildDictionaryOfUserDefinePdsForSampling(userDefinedPdfs)
    #PcriticalFinal = zeros(size(surface.coord,1), nSamples)
    PcriticalFinal = zeros(length(userDefinedPdfs[:Depth]), nSamples)

    ##
    SvGrad = zeros(length(userDefinedPdfs[:Depth]))
    SHmaxGrad = deepcopy(SvGrad)
    ShminGrad = deepcopy(SvGrad)
    SHmaxDir = deepcopy(SvGrad)
    PpGrad = deepcopy(SvGrad)
    μ = deepcopy(SvGrad)
    Aphi = deepcopy(SvGrad)        
    ##

    ## Compute the critical pore pressure to slip for each location in the spatial grid
    ProgressMeter.@showprogress for i=1:size(PcriticalFinal,1)
        samples=FSP3D.CreateDataFrameOfSamplesForMonteCarloSimulation(finalUserPdfs, i,  nSamples)
        #@inbounds PcriticalFinal[i,:] = FSP3D.MainComputePcritical(samples, modelType; kwargs=kwargs)
        @inbounds PcriticalFinal[i,:], SvGrad[i], SHmaxGrad[i], ShminGrad[i], SHmaxDir[i], PpGrad[i], μ[i], Aphi[i] = FSP3D.MainComputePcritical(samples, modelType; kwargs=kwargs)
    end ##for

    propToQC = Dict(:SvGrad => SvGrad, :SHmaxGrad => SHmaxGrad, :ShminGrad => ShminGrad, :SHmaxDir => SHmaxDir, :PpGrad => PpGrad, :μ => μ, :Aphi => Aphi)

    return PcriticalFinal, propToQC
    
end #ComputeCriticalPorePressureToFailure(surface, userDefinedPdfs, nSamples; kwargs=Dict([]) )
=#

function ExportResultsPpToFailureOnSurfaceToPetrel(dirToSave::String,fileNamesAll::Vector{String}, surface::Union{Vector{Any},Vector{FSP3D.SurfaceWithAttributesPetrel}, Vector{FSP3D.SurfaceFromPetrelEarthVisionFormat}}, ppFailureFinal::Matrix{Float64}, ppFailureName::Vector{String}, finalResults::Dict{Symbol, Array}, df::DataFrames.DataFrame, param)    
    #=
      Export results of pore pressure to failure computed on surface to Petrel
    =#

    inputUnits=param[:unitsDict][:inputUnits]

    if haskey(param[:unitsDict], :outputUnits)
        outputUnits=param[:unitsDict][:outputUnits]
    else
        outputUnits=Dict([])
    end

    FSP3D.CreateDir(dirToSave)

    SvGrad = mean(finalResults[:Sv] ./ finalResults[:depth] ; dims=2)
    SHmaxGrad = mean(finalResults[:SHmax] ./ finalResults[:depth]; dims=2)
    ShminGrad = mean(finalResults[:Shmin] ./ finalResults[:depth]; dims=2)
    PpGrad = mean(finalResults[:Pp] ./ finalResults[:depth]; dims=2)
    μ=mean(finalResults[:μ] ; dims=2 )
    Aphi=mean(finalResults[:Aphi]; dims=2)
    faultDip=mean(finalResults[:faultDip]; dims=2)
    faultStrike=mean(finalResults[:faultStrike]; dims=2)
    SHmaxDir=mean(finalResults[:SHmaxDir]; dims=2)
    depth=mean(finalResults[:depth]; dims=2)

    tauMag=finalResults[:tauMag]
    sigmaN=finalResults[:sigmaN]
    Ts = tauMag ./ sigmaN

    probFailureTs=[100, 99, 95, 90, 50]
    probFailureTsName="Ts_P" .* string.(probFailureTs)
    slipTendencyPercentile = FSP3D.ComputePcCriticalUserDefinedProbabilityOfFailure(Ts, probFailureTs )

    attributesToSave = [ SvGrad SHmaxGrad ShminGrad PpGrad μ Aphi faultDip faultStrike SHmaxDir depth slipTendencyPercentile ppFailureFinal   ]

    attributeName=["SvGrad"; "SHmaxGrad"; "ShminGrad"; "PpGrad"; "FrictionCoeff"; "Aphi"; "faultDip"; "faultStrike" ; "SHmaxDir" ; "depth"; probFailureTsName[:]  ; ppFailureName[:]]

    faultIdUnique=unique(df.faultId)
    for faultId in faultIdUnique
        indAll=findall(faultId .== df.faultId)

        ## Fix output units here
        if haskey(outputUnits, :depth) && haskey(inputUnits, :depth)
            if outputUnits[:depth] == "m" && inputUnits[:depth] == "ft"
                surface[faultId].coord[:,end] = surface[faultId].coord[:,end] .* 0.3048
                replace!(surface[faultId].header, "Unit in depth: ft" => "Unit in depth: m")
            elseif outputUnits[:depth] == "ft" && inputUnits[:depth] == "m"
                surface[faultId].coord[:,end] = surface[faultId].coord[:,end] ./ 0.3048
                replace!(surface[faultId].header, "Unit in depth: m" => "Unit in depth: ft")
            end
        end 

        if hasproperty(surface[faultId], :faultDip) && hasproperty(surface[faultId], :faultStrike)
            #println("\n \t Surface is an horizon, there is no fault dip or azimuth")
            surfaceFinal=FSP3D.SurfaceWithAttributesPetrel(surface[faultId].coord, surface[faultId].grid, attributesToSave[indAll,:], attributeName, surface[faultId].header, surface[faultId].depthUnitSystem, surface[faultId].faultDip, surface[faultId].faultDipAzimuth, surface[faultId].faultStrike)
        else
            #println("\n \t Surface is an horizon, there is no fault dip or azimuth")
            tmp=surface[faultId].coord[:,1] .* 0
            surfaceFinal=FSP3D.SurfaceWithAttributesPetrel(surface[faultId].coord, surface[faultId].grid, attributesToSave[indAll,:], attributeName, surface[faultId].header, surface[faultId].depthUnitSystem, tmp, tmp, tmp)
        end

        fileNameToSave=dirToSave*fileNamesAll[faultId]
        println("Saving file "*fileNameToSave)
        FSP3D.WriteSurfaceWithAttributesToFile(fileNameToSave, surfaceFinal)
    end
    
end #ExportResultsPpToFailureOnSurfaceToPetrel(dirToSave::String, surface::FSP3D.SurfaceWithAttributesPetrel,faultDip::Vector{Float64}, surfaceName::String, faultStrike::Vector{Float64}, ppFailureFinal::Matrix{Float64})

function FixDepthUnits(surface::FSP3D.SurfaceWithAttributesPetrel, stressGradientUnits::String, depth::Vector{Float64})
    #=
     Ensure consistency in the depth units. This is important for stress gradient calculation
    =#
    
    if surface.depthUnitSystem == "m" && stressGradientUnits == "psi/ft"
        depth[:] = depth[:] * 3.28 ## meters to ft
    elseif surface.depthUnitSystem == "Foot_US" && stressGradientUnits == "psi/ft"
        depth[:] = depth[:] .* 1 ## meters to ft
    elseif surface.depthUnitSystem == "ft" && stressGradientUnits == "MPa/km"
        throw(ErrorException("Erro!! Conversion not implemented yet"))
    else
        throw(ErrorException("Erro!!"))
    end
    
    return depth
    
end #FixDepthUnits(surface::FSP3D.SurfaceWithAttributesPetrel, stressGradientUnits::String, depth)

function ComputePhiParameter(Aphi::Real, n::Real)
    #=
     Compute the ϕ parameter that is required by the critically stressed model
    =#

    ϕ = ( (Aphi .- n .- 0.5) ./ ( (-1) .^ n ) ) .+ 0.5

    return ϕ
    
end #ComputePhisParameters(Aphi, n)

function ComputeFrictionalStressLimitFactor(μ::Real)
    #=
     Compute ratios of the S1 / S3 assuming the frictional equibilibrium model
    =#
    
    if μ < 0
        throw(ErrorException("Error!! Friction coefficient cannot be negative"))
    end

    B = ( sqrt.( (μ .^ 2 .+ 1) ) .+ μ ) .^ 2
    
    return B
end #ComputeFrictionalStressLimitFactor()

function ComputeStressTensor_CS_Strike_Slip_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=1, Shmin=[], SHmax = [])
    #=
      Compute SHmax for the critically stressed faults for strike-slip faults
    =#

    ## Positive values means compressive stress
    Sv = abs.(Sv)
    Pp = abs.(Pp)
    
    if Aphi < 1 || Aphi > 2
        throw(ErrorException("Error!!! The Aphi parameter for strike slip faults should be between 1 and 2"))
    end

    ## Compute the critically stressed fault factor
    B = FSP3D.ComputeFrictionalStressLimitFactor(μ)

    ## Computing the phi parameters
    ϕ = FSP3D.ComputePhiParameter(Aphi, n)

    if Shmin == 0 && SHmax == 0
        tmp1 = (Sv ./ ϕ) .+ (Pp ./ (ϕ .* B)) .- (Pp ./ ϕ) .- (Pp ./ B) .+ Pp
        tmp2 = 1 .+ (1 ./ (ϕ .* B)) .- (1 ./ B)
        SHmax = tmp1 ./ tmp2
        Shmin = ( (SHmax .- Pp) ./ B ) .+ Pp
    elseif Shmin == 0 && SHmax != 0
        Shmin = ( (SHmax .- Pp) ./ B ) .+ Pp
    elseif Shmin != 0 && SHmax == 0
        SHmax = B .* Shmin .- Pp .* (B .- 1)
    end

    if any(isempty([Sv, Pp, SHmax, Shmin, Aphi, μ]))
        throw(ErrorException("Error!!! Stress is empty"))
    end

    return DataFrames.DataFrame(:Sv => Sv, :Pp => Pp, :μ => μ, :Aphi => Aphi, :SHmax => SHmax, :Shmin => Shmin)

end #ComputeStressTensor_CS_Strike_Slip_Faults

function ComputeStressTensor_CS_Normal_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=0, Shmin=[], SHmax=[])
    #=
    Compute SHmax for the critically stressed faults for normal faults
    =#

    ## Positive values means compressive stress
    Sv = abs.(Sv)
    Pp = abs.(Pp)

    if Aphi < 0 || Aphi > 1
        throw(ErrorException("Error!!! The Aphi parameter for normal faults should be between 0 and 1"))
    end

    ## Compute the critically stressed fault factor
    B = FSP3D.ComputeFrictionalStressLimitFactor(μ)

    # Compute the ϕ parameters
    ϕ = FSP3D.ComputePhiParameter(Aphi, n)

    if Shmin == 0 && SHmax == 0
        Shmin = ((Sv .- Pp) ./ B) .+ Pp
        SHmax = ϕ .* (Sv .- Shmin) .+ Shmin
    elseif Shmin == 0 && SHmax != 0
        Shmin = (SHmax .- ϕ .* Sv) ./ (1 .- ϕ)
    elseif Shmin != 0 && SHmax == 0
        SHmax = ϕ .* (Sv .- Shmin) .+ Shmin        
    end

    if any(isempty([Sv, Pp, SHmax, Shmin, Aphi, μ]))
        throw(ErrorException("Error!!! Stress is empty"))
    end

    return DataFrames.DataFrame(:Sv => Sv, :Pp => Pp, :μ => μ, :Aphi => Aphi, :SHmax => SHmax, :Shmin => Shmin)
    
end #ComputeStressTensor_CS_Normal_Faults

function ComputeStressTensor_CS_Reverse_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=2, Shmin = [], SHmax = [])
    #=
    Compute SHmax for the critically stressed faults for normal faults
    =#

    ## Positive values means compressive stress
    Sv = abs.(Sv)
    Pp = abs.(Pp)
    
    if Aphi < 2 || Aphi > 3
        throw(ErrorException("Error!!! The Aphi parameter for normal faults should be between 0 and 1"))
    end

    ## Compute the critically stressed fault factor
    B = FSP3D.ComputeFrictionalStressLimitFactor(μ)

    # Compute the ϕ parameters
    ϕ = FSP3D.ComputePhiParameter(Aphi, n)

    if Shmin == 0 && SHmax == 0
        SHmax = B .* (Sv .- Pp) .+ Pp
        Shmin = ϕ .* (SHmax .- Sv) .+ Sv
    elseif Shmin == 0 && SHmax != 0
        Shmin = ϕ .* (SHmax .- Sv) .+ Sv
    elseif Shmin != 0 && SHmax == 0
        SHmax = (1 ./ ϕ) * (Shmin .- Sv .* (1 - ϕ))
    end

    if any(isempty([Sv, Pp, SHmax, Shmin, Aphi, μ]))
        throw(ErrorException("Error!!! Stress is empty"))
    end

    return DataFrames.DataFrame(:Sv => Sv, :Pp => Pp, :μ => μ, :Aphi => Aphi, :SHmax => SHmax, :Shmin => Shmin)
    
end #ComputeStressTensor_CS_Reverse_Faults

function ComputeStressTensor_CS_Model(Sv::Real, SHmax::Real, Shmin::Real, Pp::Real, μ::Real, Aphi::Real)
    #=
       Compute stress tensor for the critically stressed model
    =#

    df=DataFrames.DataFrame([])
    
    if Aphi > 0 && Aphi <= 1 #normal stress regimen
        df=ComputeStressTensor_CS_Normal_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=0, Shmin=Shmin, SHmax=SHmax)
    elseif  Aphi > 1 && Aphi <= 2 #strike-slip stress regimen
        df=ComputeStressTensor_CS_Strike_Slip_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=1, Shmin=Shmin, SHmax=SHmax)
    elseif Aphi > 2 && Aphi <=3 #reverse stress regimen
        df=ComputeStressTensor_CS_Reverse_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=2, Shmin=Shmin, SHmax=SHmax)
    else
        throw(ErrorException("Error!!! Aphi parameter must be between 0 and 3"))
    end

    return df
end #ComputeStressTensor_CS_Model()

function GetFaultNamesInDirectory(dirData::String)
    #=
     Get fault names in directory
    =#
    
    fileNamesAll=readdir(dirData)
    indAll=findall( occursin.("xml", fileNamesAll) .== 0)
    if !isempty(indAll)
        fileNamesAll=sort(fileNamesAll[indAll])
    end

    return fileNamesAll
    
end #GetFaultNamesInDirectory()

function CreateDir(DirToSave::String)
    ### Check if directory exist and if not then crates a new one

    if isdir(DirToSave) == 0
        println("\n \t Creating directory = ", DirToSave)
        mkpath(DirToSave)
        #cmd=`mkdir -p $DirToSave`
        #run(cmd)
    end

end  #CreateDir

function GetValueIntepProperty(spl, xCoord::Union{Real, Vector{Float64}, Vector{Int64}}, yCoord::Union{Real, Vector{Float64}, Vector{Int64}})
    #=
      Get values of spatial properties using the interpolated spline functions
    =#

    outValue=[]
    #spl=inputFun[:spl]

    outValue = try
        spl(xCoord, yCoord)
    catch
        ScatteredInterpolation.evaluate(spl,[xCoord yCoord]')
    end

    return outValue
    
end #GetValueOfSpatiallyInterpolatedProperty(spl, xCoord, yCoord)

function GetValueOfSpatiallyInterpolatedProperty(df::DataFrames.DataFrame, kwargs::Dict{Symbol, Dict{Symbol, Any}}, xCoord::Union{Vector{Float64}, Real}, yCoord::Union{Vector{Float64}, Real})
    #=
      If spatially interpolated property is available, then get it
    =#
    
    propInterp = Dict{Symbol, Vector{Float64}}([])
    for key in keys(kwargs)
        if !isempty(get(kwargs, key, []))
            spl=kwargs[key][:spl]
            propInterp = merge(propInterp, Dict(key => FSP3D.GetValueIntepProperty(spl, xCoord, yCoord)) )
        else
            propInterp=merge(propInterp, df[!,key])
        end
    end

    return propInterp
    
end #GetValueOfSpatiallyInterpolatedProperty()

function ComputeStressInformationNeededForPorePressureToFailureCalculation(samples::DataFrames.DataFrame)
    #=
      For a given user defined stress input, it computes Sv, SHmax and Shmin that is required to compute the critical pore pressure to failure
    =#

    depth=abs.(samples.depth)
    
    Sv=zeros(length(depth))
    SHmax=deepcopy(Sv)
    Shmin=deepcopy(Sv)
    Pp = deepcopy(Sv)
    μ = deepcopy(Pp)
    Aphi = deepcopy(Sv)

    if any(names(samples) .== "SvGrad")
        Sv = samples.SvGrad .* depth
    elseif any(names(samples) .== "Sv")
        Sv = samples.Sv
    else
        throw(ErrorException("\n \t Error!! Information about the vertical stress Sv  must always be provided "))
    end

    if any(names(samples) .== "PpGrad")
        Pp = samples.PpGrad .* depth
    elseif any(names(samples) .== "Pp")
        Pp = samples.Pp
    else
        throw(ErrorException("\n \t Error!! Information about pore pressure Pp must always be provided "))
    end

    if any(names(samples) .== "SHmaxGrad")
        SHmax = samples.SHmaxGrad .* depth
    elseif any(names(samples) .== "SHmax")
        SHmax = samples.SHmax
    end

    if any(names(samples) .== "ShminGrad")
        Shmin = samples.ShminGrad .* depth
    elseif any(names(samples) .== "Shmin")
        Shmin = samples.Shmin
    end
    
    if any(names(samples) .== "μ")
        μ = samples.μ
    else
        throw(ErrorException("\n \t Error!! Information about fault friction coefficient must always be provided "))
    end

    if any(names(samples) .== "Aphi")
        Aphi = samples.Aphi
    end

    if all(SHmax .== 0) && all(Aphi .== 0)
        throw(ErrorException("\n \t Error!!! SHmax has not been defined, so you must specify a non-zero Aphi parameter"))
    end

    if all(Shmin .== 0) && all(Aphi .== 0)
        throw(ErrorException("\n \t Error!!! Shmin has not been defined, so you must specify a non-zero Aphi parameter"))
    end
    
    #=
    if any(isempty.([Sv, Pp, μ, Aphi])) && any(isempty.([Sv, Pp, μ, SHmax, Shmin]))
        throw(ErrorException("Error!!! Not all input values were defined"))
    end

    if any(isempty.([Sv, Pp, μ, Aphi])) && any(isempty.([Sv, Pp, μ, SHmax, Shmin]))
        throw(ErrorException("Error!!! Not all input values were defined"))
    end
    =#

    if all(SHmax .== 0) || all(Shmin .== 0)
        if !isempty(Aphi)
            @Threads.threads for i=1:length(Sv)
                tmp = FSP3D.ComputeStressTensor_CS_Model.(Sv[i], SHmax[i], Shmin[i], Pp[i], μ[i], Aphi[i])
                SHmax[i] = tmp.SHmax[1]
                Shmin[i] = tmp.Shmin[1]
            end
        else
            throw(ErrorException("\n \t Error!! Becauase either SHmax or Shmin are not defined, you must specify the Aphi parameter"))
        end
    end

    #if any(isempty.([Sv, Pp, SHmax, Shmin, μ]))
        #throw(ErrorException("Error!!! Not all input stress values were defined"))
    #end

    if all(SHmax .== 0) && any(abs.(depth) .> 0)
        throw(ErrorException("Error!!! SHmax is all zeros!"))
    end

    if all(Shmin .== 0) && any(abs.(depth) .> 0)
        throw(ErrorException("Error!!! Shmin is all zeros!"))
    end

    if all(Sv .== 0) && any(abs.(depth) .> 0)
        throw(ErrorException("Error!!! Sv is all zeros!"))
    end

    faultStrike=samples.faultStrike
    faultDip=samples.faultDip
    SHmaxDir=samples.SHmaxDir


    return Sv, SHmax, Shmin, Pp, SHmaxDir, faultStrike, faultDip, μ, depth, Aphi
    
end #ComputeStressInformationNeededForPorePressureToFailureCalculation(samples::DataFrames.DataFrame)

function CreateSamplesForMonteCarloSimulationUsingUserDefinedPDFs(inputDf::DataFrames.DataFrameRow{DataFrames.DataFrame, DataFrames.Index}, nSamples::Int64)
    #=
       The end goal of this function is to create samples for Monte Carlo simulation using the user defined PDFs
       Note that the input is just 1 row of the dataframes containing the use defined PDFs, where each row correspond to a different location
    =#
    
    ## Build container to hold user defined PDFs for sampling
    tmp=Union{ UncertaintyQuantification.Parameter, UncertaintyQuantification.RandomVariable }[]
    for name in names(inputDf)
        #global tmp
        if typeof(inputDf[name]) == Array{UncertaintyQuantification.Parameter,1}
            tmp = [tmp; inputDf[name]]  
            
        elseif typeof(inputDf[name]) == UncertaintyQuantification.Parameter
            tmp = [tmp; inputDf[name]]  ## constant value assigned to this location
            
        elseif typeof(inputDf[name]) == Int64 || typeof(inputDf[name]) == Float64
            tmp = [tmp; UncertaintyQuantification.Parameter.(inputDf[name], Symbol(name)) ]  ## constant value assigned to this location

        elseif typeof(inputDf[name]) != String 
            tmp = [tmp; inputDf[name]]  ## pdf assigned to this location            
                    
        #elseif typeof(inputDf[name]) != String && typeof(inputDf[name]) != Int64 && typeof(inputDf[name]) != Float64
            #tmp = [tmp; inputDf[name]]  ## pdf assigned to this location            
        #else
            #tmp = [tmp; finalUserPdfs[name]]  ## pdf assigned to this location
        end
    end

    samples=UncertaintyQuantification.sample(tmp, nSamples)

    return samples

end #CreateSamplesForMonteCarloSimulationUsingUserDefinedPDFs()

function UpdateUserDataFrameWithPdfsForSampling(df::DataFrames.DataFrame, userDefinedPdfs::Dict{Symbol, Dict{Symbol, Any}})
    #=
      This function takes as input the user defined dataframe with stress information and then updates it to included pdfs for the uncertainty quantification

    Note that below you might have to add the option for aditional probability density functions
    
    =#

    for name in names(df)
        if haskey(userDefinedPdfs,  Symbol(name))
            if !any( userDefinedPdfs[Symbol(name)][:pdfType] .== ["Gaussian", "Uniform", "Constant", "UserDefined", "Histogram", "HistogramDepth"] )
                throw(ErrorException("\n \t Error!!! User has requested a type of distribution that does not exist !! User requeste the following distribution = \n \t \t " *string(userDefinedPdfs[Symbol(name)][:pdfType])))
            end
        end
    end

    ### Adding additional fields to the dataframe from the userDefinedPDF.
    for key in String.(keys(userDefinedPdfs))
        if !any(key .== names(df))
            #println("Adding the following key = ", key)
            df[!,key] .= -9999
        end
    end

    for name in names(df)
        
        if haskey(userDefinedPdfs,  Symbol(name))
            
            if userDefinedPdfs[Symbol(name)][:pdfType] == "Gaussian"

                if haskey(userDefinedPdfs[Symbol(name)],  :std) && !haskey(userDefinedPdfs[Symbol(name)],  :μ)
                    std=userDefinedPdfs[Symbol(name)][:std]
                    #df[!,name] = UncertaintyQuantification.RandomVariable.(Distributions.Gaussian.(df[!,name], std), Symbol(name))
                    df[!,name] = Distributions.Gaussian.(df[!,name], std)
                elseif haskey(userDefinedPdfs[Symbol(name)],  :std) && haskey(userDefinedPdfs[Symbol(name)],  :μ)

                    std=userDefinedPdfs[Symbol(name)][:std]
                    meanVal=userDefinedPdfs[Symbol(name)][:μ]
                    df[!,name] .= Distributions.Gaussian.(meanVal, std)
                    
                elseif haskey(userDefinedPdfs[Symbol(name)],  :stdFraction)
                    stdFraction=userDefinedPdfs[Symbol(name)][:stdFraction]
                    std=df[!,name] .* stdFraction
                    #df[!,name] = UncertaintyQuantification.RandomVariable.(Distributions.Gaussian.(df[!,name], std), Symbol(name))
                    @assert all(std .> 0) "\n \t Error!! Standard devidation in for the Normal distribution must be > 0. Current value is "*string(std)*" for property "*string(name)
                    df[!,name] = Distributions.Gaussian.(df[!,name], std)
                else
                    throw(ErrorException("\n \t Error!! User did not defined all the necessary parameters to create a Gaussian distributions. Missing :stdFraction or :std parameters for the key name "*string(name)))
                end

                if haskey(userDefinedPdfs[Symbol(name)],  :truncated)
                    ls=userDefinedPdfs[Symbol(name)][:truncated][1]
                    hs=userDefinedPdfs[Symbol(name)][:truncated][2]

                    if typeof(ls) == String && typeof(hs) == String
                        df[!,name] = Distributions.truncated.( df[!,name] , df[!,ls], df[!,hs] )
                    elseif !(typeof(ls) == String) && typeof(hs) == String
                        df[!,name] = Distributions.truncated.( df[!,name] , ls, df[!,hs] )
                    elseif typeof(ls) <: Real && typeof(hs) <: Real
                        df[!,name] = Distributions.truncated.( df[!,name] , ls, hs )
                    else
                        throw("\n \t Errro!! Something went wrong in creating the truncated Gaussian!!!")
                    end
                    
                end

                df[!,name] = UncertaintyQuantification.RandomVariable.(df[!,name], Symbol(name))

            elseif userDefinedPdfs[Symbol(name)][:pdfType] == "Uniform"

                if haskey(userDefinedPdfs[Symbol(name)],  :pError)
                    pError=userDefinedPdfs[Symbol(name)][:pError]
                    p1 = df[!,name] .* (1 .- pError)
                    p2 = df[!,name] .* (1 .+ pError)
                    pFinal = [p1 p2]
                    pFinal = sort(pFinal; dims=2)
                end

                if haskey(userDefinedPdfs[Symbol(name)],  :pInterval)
                    pError=userDefinedPdfs[Symbol(name)][:pInterval]
                    p1 = df[!,name] .* 0 .+ pError[1]
                    p2 = df[!,name] .* 0 .+ pError[2]
                    pFinal = [p1 p2]
                    pFinal = sort(pFinal; dims=2)
                end

                if any(pFinal[:,1] .>= pFinal[:,2])
                    throw(ErrorException("\n \t Error in defining the uniform distribution with key "*string(name)*". Please check the definition of the bounds for the user defined PDF"))
                end
                
                df[!,name] = UncertaintyQuantification.RandomVariable.(Distributions.Uniform.( pFinal[:,1] , pFinal[:,2] ), Symbol(name))

            #elseif userDefinedPdfs[Symbol(name)][:pdfType] == "UserDefined"
                #pdf=userDefinedPdfs[Symbol(name)][:pdf]
                #df[!,name] .= UncertaintyQuantification.RandomVariable.(pdf, Symbol(name))
                
            elseif userDefinedPdfs[Symbol(name)][:pdfType] == "Constant"
                
                fixedValue=userDefinedPdfs[Symbol(name)][:value]
                df[!,name] .= UncertaintyQuantification.Parameter.(fixedValue, Symbol(name))
                
            #else
                #throw(ErrorException("Error!!! Pdf not defined for parateter "*string(name)))
                
            elseif userDefinedPdfs[Symbol(name)][:pdfType] == "Histogram"

                inputData=userDefinedPdfs[Symbol(name)][:data]
                uvdist=FSP3D.CreateEmpiricalDistributionBasedOnDataHistogram(inputData::AbstractArray)
                df[!,name] .= UncertaintyQuantification.RandomVariable.(uvdist, Symbol(name))
                
            elseif userDefinedPdfs[Symbol(name)][:pdfType] == "HistogramDepth"

                ## getting the depth associated with each interval in the user defined dictionary

                ## first deleting the column name to start from scratch
                DataFrames.select!(df, DataFrames.Not([name]))
                
                depthDict=Array{Float64}(collect(keys(userDefinedPdfs[Symbol(name)][:data])))

                for iVal=1:length(df.depth)

                    depthVal = df.depth[iVal]
                    if typeof(df.depth[iVal]) == UncertaintyQuantification.Parameter
                        depthVal = UncertaintyQuantification.sample(df.depth[iVal],1)[!, :depth][1]
                    end
                    
                    minv, indMin = findmin(abs.( abs(depthVal) .- abs.(depthDict)) )
                    ## getting the data from this interval to build the histogram
                    inputData = userDefinedPdfs[Symbol(name)][:data][depthDict[indMin]]
                    inputData=convert(AbstractArray{Float64}, inputData)

                    uvdist=FSP3D.CreateEmpiricalDistributionBasedOnDataHistogram(inputData::AbstractArray)

                    ### Now recreating the dataframe
                    if iVal == 1
                        df[!,name] .= UncertaintyQuantification.RandomVariable.(uvdist, Symbol(name))
                    end

                    df[iVal,name] = UncertaintyQuantification.RandomVariable(uvdist, Symbol(name))
                end #for
                
            end
            
        end

        if name == "depth"
            df[!,name] = UncertaintyQuantification.Parameter.(df[!,name], Symbol(name))
        end
        
    end #for

    ## Adding the user defined PDFs, if available
    for keyName in keys(userDefinedPdfs)
        if userDefinedPdfs[keyName][:pdfType] == "UserDefined"
            pdf=userDefinedPdfs[Symbol(keyName)][:pdf]
            df[!,keyName] .= UncertaintyQuantification.RandomVariable.(pdf, Symbol(keyName))
        end
    end
    

end #UpdateUserDataFrameWithPdfsForSampling(df::DataFrames.DataFrame, userDefinedPdfs::Dict{Symbol, Dict{Symbol, Any}})

function AssembleStressAndFaultInformationFromUser(faultShapes::Union{Dict{Int64, Dict{Symbol}}, Dict{Int64, Dict}}, unitsDict::Dict{Symbol, String}; kwargs=Dict([]))
    #=
    Takes as input the fault shapes, then creates dataframe holding all the fault and stress information
    It first searches for information in the shapefiles, then use the kwargs to add user defined information
    =#

    totalNumberOfSegments=0
    for iShape=1:length(faultShapes)
        #global totalNumberOfSegments
        totalNumberOfSegments=totalNumberOfSegments .+ faultShapes[iShape][:nSegments]
    end

    keyNamesShapeFile=sort(collect(keys(faultShapes[1][:segment])))
    keyNamesUser=sort(collect(keys(kwargs)))

    finalDict=Dict{Symbol, Matrix{Float64}}([])
    iKey=1

    for iKey = 1 : length(keyNamesShapeFile)
        #global finalDict
        if Symbol(keyNamesShapeFile[iKey]) == :coord
            finalDict = merge(finalDict, Dict(Symbol(keyNamesShapeFile[iKey])  => zeros(totalNumberOfSegments,2)) )
        else
            finalDict = merge(finalDict, Dict(Symbol(keyNamesShapeFile[iKey])  => zeros(totalNumberOfSegments)) )
        end
    end
    finalDict=merge(finalDict, Dict(:depth => zeros(totalNumberOfSegments), :segmentId => zeros(Int64, totalNumberOfSegments), :faultId => zeros(Int64, totalNumberOfSegments)) )

    ### Adding the information from the shapefile into the final dictionary
    for iKey=1:length(keyNamesShapeFile)
        #global finalDict
        countJS=1
        for iShape=1:length(faultShapes)
            #global coord, faultStrike, SHmaxDir, Aphi, SvGrad, faultId, segmentId, countJS
            for i=1:faultShapes[iShape][:nSegments]

                keyName=Symbol(keyNamesShapeFile[iKey])
                if keyName == :coord
                    if size(faultShapes[iShape][:segment][keyName],2) == 2
                        finalDict[keyName][countJS,1:2] =faultShapes[iShape][:segment][keyName][i,1:2]
                    elseif size(faultShapes[iShape][:segment][keyName],2) == 3
                        finalDict[keyName][countJS,1:2] =faultShapes[iShape][:segment][keyName][i,1:2]
                        finalDict[:depth][countJS] =faultShapes[iShape][:segment][keyName][i,3]
                    else
                        throw(ErrorException("Error!!"))
                    end
                else
                    if !isempty(faultShapes[iShape][:segment][keyName])
                        finalDict[keyName][countJS] = faultShapes[iShape][:segment][keyName][i]
                    end
                   #else
                        #finalDict[keyName][countJS] = []
                    #end
                end

                finalDict[:segmentId][countJS] = i
                finalDict[:faultId][countJS] = iShape
                countJS=countJS+1
            end
        end

    end

    ## Checking for user defined keyword and adding it to the dataframe if it does not exist
    for iKey=1:length(keyNamesUser)
        #global finalDict
        if keyNamesUser[iKey] != Symbol(:units)
            if !any(occursin.(String.(keyNamesUser[iKey]), String.(keyNamesShapeFile)))
                userVal=kwargs[keyNamesUser[iKey]]
                if !isempty(userVal)
                    finalDict=merge(finalDict, Dict( keyNamesUser[iKey] => ones(totalNumberOfSegments) .* userVal ) )
                end
            end
        end
    end

    ## Now creating the dataframe that will be used for the monte carlo simulations
    df=DataFrames.DataFrame([])
    finalKeys=collect(keys(finalDict))
    for iKey=1:length(finalKeys)
        #global df
        keyName=Symbol(finalKeys[iKey])
        if keyName != :coord
            df[!,keyName]= finalDict[keyName]
        else
            df[!,:X]= finalDict[keyName][:,1]
            df[!,:Y]= finalDict[keyName][:,2]        
        end
    end


    @assert haskey(unitsDict, :depth) "Error!! No depth units defined."

    if unitsDict[:depth] == "m" && unitsDict[:stressGradient] == "psi/ft"
        df.depth = (df.depth .* 3.28)  ## meter to feet
    elseif unitsDict[:depth] == "ft" && unitsDict[:stressGradient] == "psi/m"
        df.depth = (df.depth)  ./ 3.28
    elseif unitsDict[:depth] == "ft" && unitsDict[:stressGradient] == "psi/ft"
        df.depth = (df.depth)  
    else 
        throw(ErrorException("Error!! Units not defined."))
    end

    return df

end #AssembleStressAndFaultInformationFromUser(faultShapes::Dict{Int64, Dict{Symbol}}, kwargs)

#=
function AssembleStressAndFaultInformationFromUser(faultShapes::Dict{Int64, Dict{Symbol}}, kwargs)
    #=
      Takes as input the fault shapes, then creates dataframe holding all the fault and stress information
    =#
    totalNumberOfSegments=0
    for iShape=1:length(faultShapes)
        #global totalNumberOfSegments
        totalNumberOfSegments=totalNumberOfSegments .+ faultShapes[iShape][:nSegments]
    end

    depth=get(kwargs, :depth, [])
    faultDip=get(kwargs, :faultDip, [])
    
    PpGrad=get(kwargs, :PpGrad, [])
    μ=get(kwargs, :μ, [])

    coord=zeros(totalNumberOfSegments, 2)
    faultStrike=zeros(totalNumberOfSegments)
    
    Aphi=deepcopy(faultStrike)
    SHmaxDir=deepcopy(faultStrike)
    SvGrad=deepcopy(Aphi)
    faultId=zeros(Int64, length(Aphi))
    segmentId=deepcopy(faultId)

    μ = ones(length(faultStrike)) .* μ

    if !isempty(PpGrad)
        PpGrad=ones(totalNumberOfSegments) .* PpGrad
    else
        PpGrad = zeros(totalNumberOfSegments)
    end
    
    if !isempty(faultDip)
        faultDip = ones(length(faultStrike)) .* faultDip
    else
        faultDip = deepcopy(faultStrike)
    end

    if !isempty(depth)
        depth = ones(length(faultStrike)) .* depth
    else
        depth = deepcopy(faultStrike)
    end
    
    countJS=1
    for iShape=1:length(faultShapes)
        #global coord, faultStrike, SHmaxDir, Aphi, SvGrad, faultId, segmentId, countJS
        for i=1:size(faultShapes[iShape][:segment][:coord],1)
            coord[countJS,1:2]=faultShapes[iShape][:segment][:coord][i,1:2]
            faultStrike[countJS]=faultShapes[iShape][:segment][:faultStrike][i]

            if haskey(faultShapes[iShape][:segment], :SHmaxDir)
                SHmaxDir[countJS]=faultShapes[iShape][:segment][:SHmaxDir][i]
            end
            
            SvGrad[countJS]=faultShapes[iShape][:segment][:SvGrad][i]
            Aphi[countJS]=faultShapes[iShape][:segment][:Aphi][i]
            faultId[countJS]=iShape
            segmentId[countJS]=i

            if size(faultShapes[iShape][:segment][:coord],2) == 3
                depth[countJS]=abs(faultShapes[iShape][:segment][:coord][i,3])
            end

            if haskey(faultShapes[iShape][:segment], :faultDip)
                faultDip[countJS]=faultShapes[iShape][:segment][:faultDip][i]
            end
            
            countJS=countJS+1
        end
    end

    df=DataFrames.DataFrame(:μ => μ, :SHmaxDir => SHmaxDir, :faultStrike => faultStrike, :faultDip => faultDip, :Aphi => Aphi, :SvGrad => SvGrad, :faultId => faultId, :segmentId => segmentId, :X => coord[:,1], :Y => coord[:,2], :depth => depth, :PpGrad => PpGrad)

    return df
    
end ##AssempleStressAndFaultInformationFromUser
=#

function MainComputePcritical(inputDf::DataFrames.DataFrameRow{DataFrames.DataFrame, DataFrames.Index}, nSamples::Int64, finalResults::Dict{Symbol, Array},  iRow::Int64)
    
    #inputDf = deepcopy(df[1,:])
    samples=FSP3D.CreateSamplesForMonteCarloSimulationUsingUserDefinedPDFs(inputDf, nSamples)

    ## Now design function to compute Sv, SHmax, Shmin based on the user available information. For example: if the user has given Sv and Aphi, then use the critically stressed model
    Sv, SHmax, Shmin, Pp, SHmaxDir, faultStrike, faultDip, μ, depth, Aphi = FSP3D.ComputeStressInformationNeededForPorePressureToFailureCalculation(samples)

    sigmaN = zeros(length(Sv))
    tauMag = zeros(length(Sv))    
    Pcritical = zeros(length(Sv))
    @Threads.threads for i=1:length(Sv)
        #for i=1:length(verticalStress)
        FSP3D.AuxFunctionComputePcriticalSample(i, Pcritical, Sv, SHmax, Shmin, SHmaxDir, Pp, faultStrike, faultDip, μ, sigmaN, tauMag)
    end

    finalResults[:Pcritical][iRow,:] = Pcritical
    finalResults[:sigmaN][iRow,:] = sigmaN
    finalResults[:tauMag][iRow,:] = tauMag
    finalResults[:Sv][iRow,:] = Sv
    finalResults[:SHmax][iRow,:] = SHmax
    finalResults[:Shmin][iRow,:] = Shmin
    finalResults[:SHmaxDir][iRow,:] = SHmaxDir
    finalResults[:Pp][iRow,:] = Pp
    finalResults[:faultStrike][iRow,:] = faultStrike
    finalResults[:faultDip][iRow,:] = faultDip
    finalResults[:μ][iRow,:] = μ
    finalResults[:depth][iRow,:] = depth
    finalResults[:Aphi][iRow,:] = Aphi
    finalResults[:surfaceElevation][iRow,:] = samples.surfaceElevation
    

end ##MainComputePcritical()

function ComputeCriticalPorePressureToFailure(df::DataFrames.DataFrame, userDefinedPdfs::Dict{Symbol, Dict{Symbol, Any}}, nSamples::Int64)

    ## Updating dataframes with user defined pdfs
    FSP3D.UpdateUserDataFrameWithPdfsForSampling(df, userDefinedPdfs)

    nVals=length(df[!,1])
    finalResults=Dict(:Pcritical => zeros(nVals, nSamples),
                      :sigmaN => zeros(nVals, nSamples),
                      :tauMag => zeros(nVals, nSamples),
                      :Sv => zeros(nVals, nSamples),
                      :SHmax => zeros(nVals, nSamples),
                      :Shmin => zeros(nVals, nSamples),
                      :SHmaxDir => zeros(nVals, nSamples),
                      :Pp => zeros(nVals, nSamples),
                      :faultStrike => zeros(nVals, nSamples),
                      :faultDip => zeros(nVals, nSamples),
                      :depth => zeros(nVals, nSamples),
                      :Aphi => zeros(nVals, nSamples),
                      :μ => zeros(nVals, nSamples),
                      :surfaceElevation => zeros(nVals, nSamples),
                      :coord => [df.X df.Y],
                      :faultId => df.faultId,
                      :segmentId => df.segmentId
                      )

    ProgressMeter.@showprogress for iRow=1:size(df,1)
        FSP3D.MainComputePcritical(df[iRow,:], nSamples, finalResults, iRow)
    end

    finalResults[:SvGrad] = finalResults[:Sv] ./ finalResults[:depth]
    finalResults[:SHmaxGrad] = finalResults[:SHmax] ./ finalResults[:depth]
    finalResults[:ShminGrad] = finalResults[:Shmin] ./ finalResults[:depth]
    finalResults[:PpGrad] = finalResults[:Pp] ./ finalResults[:depth]
    finalResults[:SlipTendency] = abs.(finalResults[:tauMag] ./ finalResults[:sigmaN])

    return finalResults
    
end #ComputeCriticalPorePressureToFailure

function ReadFaultSurfacesExportedFromPetrel(dirFaultData::String; iFile=[], fileNamesAll=[])
    #=
      Read fault surfaces exported from Petrel and then create objects for MCMC 
    =#

    if isempty(fileNamesAll)
        fileNamesAll=FSP3D.GetFaultNamesInDirectory(dirFaultData)
        fileNamesAll=dirFaultData .* fileNamesAll
    else
        fileNamesAll = dirFaultData .* fileNamesAll
    end

    if !isempty(iFile)
        if length(iFile) == 1
            fileNamesAll=fileNamesAll[[iFile]]
        else
            fileNamesAll=fileNamesAll[iFile]
        end
    end
    
    println("Working with surface ", string(fileNamesAll))
    
    faultShapes = Dict{Int64, Dict{Symbol, Array{Float64}}}([])
    surface=FSP3D.SurfaceWithAttributesPetrel[]
    for iFault=1:length(fileNamesAll)
        #for iFault=1:1
        #global faultShapes, surface
        println("Reading file "*string(fileNamesAll[iFault]))
        surface = [surface; FSP3D.ReadSurfacesExportedFromPetrel(fileNamesAll[iFault])]
        
        finalCoord = surface[iFault].coord
        if hasproperty(surface[iFault], :faultDip)
            segment=Dict(:coord => finalCoord, :faultDip => surface[iFault].faultDip, :faultStrike => surface[iFault].faultStrike)
        else
            segment=Dict(:coord => finalCoord, :faultDip => [], :faultStrike => [])
        end
        
        faultShapes=merge(faultShapes, Dict(iFault => Dict(:nSegments => size(finalCoord,1), :coord => finalCoord, :paramSpl => [], :segment => segment, :units => Dict(:depth => surface[iFault].depthUnitSystem)) ))
    
    end

    return surface, faultShapes

end #ReadFaultSurfacesExportedFromPetrel()

function Save3DFSPresultsToBinary(dirToSave::String, finalResults::Dict{Symbol, Array})
    #=
      Save results of 3DFSP Run to binary file. This can be used later on to be QC
    =#

    #dirToSave = mainDir*"3DFSP/output/bin/"
    #Utils.CleanUpDirectory(dirToSave)  ## It can be dangerous
    FSP3D.CreateDir(dirToSave)

    for key in keys(finalResults)

        fileNameToSave=dirToSave .* string(key)
        println("Saving file ", fileNameToSave)

        data=finalResults[Symbol(key)]
        FSP3D.SaveFileToBinaryFormat(fileNameToSave, data)

        #tmp=ReadWriteDataGeneral.ReadBinaryFile(fileNameToSave)
        #sum(tmp .- data)

    end

end #Save3DFSPresultsToBinary()

function Read3DFSPresultsFromBinary(dirToRead::String)
    #=
     Reading results of 3DFSP runs from that were previously saved as binary
    =#

    fileNames=readdir(dirToRead)
    
    #finalDict=Dict{Symbol, Matrix{Float64}}([])
    finalDict=Dict{Symbol, Array}([])
    for key in fileNames
        fileNameToRead=dirToRead .* string(key)
        #println("Reading file ", fileNameToSave)
        if string(key) == "faultId" || string(key) == "segmentId"
            tmp=FSP3D.ReadBinaryFile(fileNameToRead; type="Int64")
        else
            tmp=FSP3D.ReadBinaryFile(fileNameToRead)            
        end
        finalDict=merge(finalDict, Dict(Symbol(key) => tmp))
    end

    return finalDict
    
end #Read3FFSPresultsFromBinary()

function ReadBinaryFile(fileName::String; type="Float64")
    #=
      Read files that were saved in binary format
    =#
    #println("Reading file ", fileName)
    ##
    s = open(fileName)
    m = read(s, Int)
    n = read(s, Int)
    if type == "Float64"
        A = Mmap.mmap(s, Matrix{Float64}, (m,n))
    else
        A = Mmap.mmap(s, Matrix{Int}, (m,n))
    end
    close(s)
    return A

end ## ReadBinaryFile

function SaveFileToBinaryFormat(fileNameToSave::String, data::Union{Matrix{Float64}, Array{Float64}, Array{Int64}})
    #=
      Save file to binary format. This is specially useful for very large files
    =#

    #fileNameToSave = dirToSave * fileName

    tmp = open(fileNameToSave, "w+")
    write(tmp, size(data,1))
    write(tmp, size(data,2))

    for i=1:size(data,2)
        write(tmp, data[:,i])
    end
    close(tmp)

end ## SaveFileToBinaryFormat

function SaveFaultTraceInformation(fileName::String, df::DataFrames.DataFrame)
    #=
      Save fault trace information for later on so that it can be displayed
    =#
    CSV.write(fileName, DataFrames.select(df,[:X, :Y, :faultId]))
    
end #SaveFaultTraceInformation()

function ReadSurfaceExportedFromPetrelEarthVisionFormat(fileName::String)
    #=
      Read surface exported from Petrel in earth vision format
    =#

    if !isfile(fileName)
        throw(ErrorException("Error!! File not found "*string(fileName)))
    end

    header=AbstractString[]
    open(fileName, "r") do io
        for ln in eachline(io)
            header = [header; ln]
            if occursin("Z_units", ln)
                break
            end
        end
    end

    gridSize=[]
    for i=1:length(header)
        if occursin("Grid_size:", header[i])
            gridSize=parse.(Int64,string.(split(header[i]," ")[[3,5]]))
        end
    end

    depthUnit = strip.(split(header[end],":"))[2]
    if depthUnit == "feet"
        depthUnit = "ft"
    end

    nVals=gridSize[1] .* gridSize[2]
    data=zeros(nVals, 5)
    countJS=1
    open(fileName, "r") do io
        for ln in eachline(io)
            if !occursin("#", ln)
                tmp = string.(split(ln," "))
                data[countJS,:] = parse.(Float64,tmp)
                countJS=countJS+1
            end
        end
    end
    data=data[1:countJS-1,:]

    gridJS=Int64.(data[:,4:5])
    data[:,3] = -abs.(data[:,3]) ## Depth is always negative
    surface=SurfaceFromPetrelEarthVisionFormat(data[:,1:3], gridJS, data[:,3], "depth", header, depthUnit, fileName)

    return surface
    
end #ReadSurfaceExportedFromPetrelEarthVisionFormat()


function ReadIP3inputFile(dirData::String)
    #=
      Read IP3 input file from user. It also creates the dictionary with the input data and the interpolation function
    =#
    #dirData = mainDir*"Julia/Data/IP3/"
    inputPropertyAll=[:Pp, :Sv, :Shmin, :SHmax]

    unitsIP3=Dict{Symbol, Any}([])
    iP3results=Dict{Symbol, Dict{Symbol, Any}}([])

    for inputProperty in inputPropertyAll
        #global iP3results, unitsIP3

        fileNameData=dirData * String(inputProperty)*".csv"
        if !isfile(fileNameData)
            println("\n \t File not found for "*string(fileNameData))
            continue
        end
        
        df=CSV.File(fileNameData; skipto=3) |> DataFrames.DataFrame

        @assert any(names(df) .== "TVD") raw"Error!! Column name TVD not found!"

        ## making sure the IP3 model is sorted according to the depth values
        indPerm=sortperm(abs.(df[!,:TVD]))
        df=df[indPerm,:]        

        if size(df,2) > 4
            throw(ErrorException("\n \t Error!! The input IP3 file should have only 4 columns: Depth, Values_LS, Value_BE, Value_HS \n "))
        end

        if any( diff( Matrix(df)[:,2:end] ; dims=2  )[:]  .< 0)

            tmpProb=[]
            for i=1:size(df,1)
                if any( diff( Matrix(df)[i,2:end] ) .< 0  )
                    #println("i = ", i)
                    tmpProb=[tmpProb; i]
                end
            end
            
            println("\n \t Error!!! The input IP3 file columns should be ordered as: depth, Values_LS, Values_BE and Values_HS. Check that the columns are correct. Here is the input values read from file:")

            println(show(df))

            println("\n \t Here are the rows in the input data with problems: \n \t" )
            for i=1:length(tmpProb)
                println("Row = "*string(tmpProb[i]))
            end
            println("\n \t ")
            
            throw(ErrorException("\n"))
        end

        ## Making sure that the depth is always negative
        df[!,1]=abs.(df[!,1]) .* -1

        ## Reading the units of the IP3 model
        dfUnits=CSV.File(fileNameData; skipto=2, limit=1) |> DataFrames.DataFrame
        unitsIP3 = merge(unitsIP3, Dict{Symbol, Any}(:depth => String(dfUnits[!,1][1]), inputProperty => String(dfUnits[!,2][1])))

        for k=2:4
            if k==2
                symbolName=Symbol(String(inputProperty)*"_LS")
            elseif k==3
                symbolName=Symbol(String(inputProperty)*"_BE")        
            elseif k==4
                symbolName=Symbol(String(inputProperty)*"_HS")        
            end

            depth=df[!,1]

            ### Convert units here
            if dfUnits[!,:TVD] == "m"
                depth = depth .* 3.28  ## meters to ft?
            end

            prop=df[!,k]

            ### Convert units here
            if dfUnits[1,symbolName] == "Pa"
                prop = prop ./ FSP3D.PsiToPa()
            end
            
            data=[depth prop]
            spl=FSP3D.CreateInterpolatedFunctionsFromIP3model(depth, prop)

            tmpFinal=Dict(symbolName => Dict(:data => data, :property => prop, :coordColumn => [3], :spl => spl))
            iP3results=merge(iP3results,  tmpFinal )
        end
    end

    iP3results=merge(iP3results, Dict(:units => unitsIP3))

    ## By default, the properties that will be used are the ones set as the "best estimate"
    if haskey(iP3results,:Sv_BE)
        iP3results=merge(iP3results, Dict(:Sv => iP3results[:Sv_BE]))
    end
    
    if haskey(iP3results,:Shmin_BE)
        iP3results=merge(iP3results, Dict(:Shmin => iP3results[:Shmin_BE]))
    end

    if haskey(iP3results,:SHmax_BE)
        iP3results=merge(iP3results, Dict(:SHmax => iP3results[:SHmax_BE]))
    end

    if haskey(iP3results,:Pp_BE)
        iP3results=merge(iP3results, Dict(:Pp => iP3results[:Pp_BE]))
    end
    
    return iP3results

end #ReadIP3inputFile

function CreateInterpolatedFunctionsFromIP3model(depth::Vector{Float64}, prop::Vector{Float64})
    #=
      Create interpolated functions for the IP3 model
    =#

    #=
    indperm=sortperm(depth)
    depth=depth[indperm]
    prop=prop[indperm]
    #spl=Dierckx.Spline1D(depth,prop; k=3, s=10000)
    spl=Dierckx.Spline1D(depth,prop; k=1)
    =#

    tmp=[depth prop]
    indPerm = sortperm(tmp[:,1])
    tmp=tmp[indPerm, :]

    ## making sure that the points are unique prior to the interpolation
    tmpUnique=unique(tmp[:,1])
    tmpFinal=zeros(size(tmpUnique,1), 2)
    countJS=1
    for i=1:length(tmpUnique)
        indAll=findall(tmpUnique[i] .== tmp[:,1])
        if !isempty(indAll) 
            tmpFinal[countJS, 1] = tmpUnique[i]
            tmpFinal[countJS, 2] = mean(tmp[indAll,2])
            countJS=countJS+1
        end
    end
    tmpFinal = tmpFinal[1:countJS-1,:]

    #tmp=unique(tmp; dims=1)
    
    depthData=tmpFinal[:,1]
    propData=tmpFinal[:,2]
    
    spl = Interpolations.linear_interpolation(depthData, propData, extrapolation_bc=Interpolations.Line())

    return spl
    
end #CreateInterpolatedFunctionsFromIP3model()

function GetStressRegime(Sv::Union{Float64, Int64}, SHmax::Union{Float64, Int64}, Shmin::Union{Float64, Int64})
    #=
      Get the stress regime asociated with the stress tensor
    =#

    stressInfo = FSP3D.AssemblePrincipalStressesAndAngles(Sv, SHmax, Shmin, 0, 0)

    return stressInfo[:stressRegime]
    
end #GetStressRegime()

function ComputeSHmaxDirectionUsingStressTensor(Sg::Union{Matrix{Float64}, Matrix{Int64}}, stressRegime::String)
    #=
      Compute the maximum horizontal stress direction for a given stress tensor. User has to specify the stress regime
    =#

    sVals, sVect = FSP3D.ComputePrincipalStresses(Sg)

    newSHmaxDir=[]
    if stressRegime == "Normal"
        k=2  ## Intermediate principal stress 
        x=sVect[:,k]
        newSHmaxDir = atand(x[2] ./ x[1])

    elseif stressRegime == "Reverse"
        k=3 ## maximum principal stress
        x=sVect[:,k]
        newSHmaxDir = atand(x[2] ./ x[1])

    elseif stressRegime == "StrikeSlip"
        k=3 ## maximum principal stress
        x=sVect[:,k]
        newSHmaxDir = atand(x[2] ./ x[1])
    else
        throw(ErrorException("Error!! Stress regime not defined !!"))
    end

    if isempty(newSHmaxDir)
        throw(ErrorException("Error!! Stress regime not defined !!"))
    end

    return Dict(:SHmaxDir => newSHmaxDir, :S1 => sVals[1], :S2 => sVals[2], :S3 => sVals[3])

end #ComputeSHmaxDirectionUsingStressTensor(Sg::Union{Matrix{Float64}, Matrix{Int64}}, stressRegime::String)

function ComputeGeometryOfCriticallyStressesFaults(μ::Float64, stressRegime::String, SHmaxDir::Union{Float64, Int64})
    #=
      Compute fault dip and fault strike for a critically stressed faults
    =#

    ϕ=atand(μ)
    β=45 - ϕ / 2

    faultStrike=[]
    faultDip=[]
    if stressRegime == "Normal"
        faultStrike=SHmaxDir
        faultDip=90 .- β
    
    elseif stressRegime == "Reverse"
        faultStrike=90 .+ SHmaxDir
        faultDip=β
        
    elseif stressRegime == "StrikeSlip"
        faultStrike=β .+ SHmaxDir
        faultDip=90
    
    else
        throw(ErrorException("Error!! Stress regime not defined !!"))
    end

    if isempty(faultDip) || isempty(faultStrike)
        throw(ErrorException("Error!! Fault dip and strike are empty "))
    end

    #println("For SHmax dir = "*string(SHmaxDir)*" and stress regime = "*string(stressRegime)*", the optimaly oriented fault dip = "*string(faultDip)*" and fault strike = "*string(faultStrike))
    
    return faultDip, faultStrike

end #ComputeGeometryOfCriticallyStressesFaults(μ::Float64, stressRegime::String, SHmaxDir::Union{Float64, Int64})

function GetFaultGeometryCriticallyStressedModel(Sg::Union{Matrix{Float64}, Matrix{Int64}}, stressRegime::String, μ::Float64)
    #=
      Get the fault geometry (fault dip and fault strike) for a critically stressed model
    =#

    tmpDict = FSP3D.ComputeSHmaxDirectionUsingStressTensor(Sg, stressRegime)
    newSHmaxDir=tmpDict[:SHmaxDir]
    faultDip, faultStrike = FSP3D.ComputeGeometryOfCriticallyStressesFaults(μ, stressRegime, newSHmaxDir)

    return faultDip, faultStrike
    
end ##  GetFaultGeometryCriticallyStressedModel(Sg::Union{Matrix{Float64}, Matrix{Int64}}, stressRegime::String, μ::Float64)

function ComputeInterpolatedPropertyAtFaultShapeSegments(faultShapes::Dict{Int64, Dict{Symbol}}, kwargsInterpProp)
    #=
     Compute interpolated proprety at fault shape segment
    =#

    for iShape=1:length(faultShapes)

        segment=faultShapes[iShape][:segment]
        coord=segment[:coord]
        
        for key in keys(kwargsInterpProp)
            #global segment
            if key != :units
                spl=kwargsInterpProp[key][:spl]
                coordColumn=kwargsInterpProp[key][:coordColumn]
                segment=merge(segment, Dict(key => FSP3D.EvaluateInterpolantAtAnyPoint(coord, spl, coordColumn) ))
            end
        end
        faultShapes[iShape][:segment] = segment
        
    end

end #ComputeInterpolatedPropertyAtFaultShapeSegments

function EvaluateInterpolantAtAnyPoint(coord::Matrix{Float64}, spl, coordColumn::Vector{Int64})
    #=
      Evaluate the interpolation function at matrix points
    =#
    
    Rgrid=zeros(size(coord,1))
    if occursin("Scattered", string(typeof(spl)))
        @Threads.threads for i=1:length(Rgrid)
            Rgrid[i]=ScatteredInterpolation.evaluate(spl,coord[i,coordColumn])[1]
        end
    else
        if length(coordColumn) > 1
            @Threads.threads for i=1:length(Rgrid)
                #Rgrid[i]=spl.(coord[i,1], coord[i,2])
                @assert 1 > 2 " Error !! You have to implement the Interpolations package !!! "
            end
        else
            @Threads.threads for i=1:length(Rgrid)
                Rgrid[i]=spl.(coord[i,coordColumn[1]])
            end
        end
    end

    return Rgrid

end ##EvaluateInterpolantAtMatrixPoints

function Read3DFSPresultsFromBinaryAllFaults(dirToRead::String; fileName=Vector{String}[], numberOfFilesToRead=[])
    #=
      Read binary file corresponding to the 3DFSP results for all faults exported
    =#
    if isempty(fileName)
        fileNames=readdir(dirToRead)
    else
        fileNames=fileName
    end

    if numberOfFilesToRead > length(fileNames)
        numberOfFilesToRead = length(fileNames)
    end

    if !isempty(numberOfFilesToRead)
        #fileNames=sort(unique(rand(fileNames, numberOfFilesToRead)))
        fileNames=sort(unique(sample(fileNames, numberOfFilesToRead; replace=false)))
    end

    println("\n \t Reading "*string(length(fileNames))*" faults")
    
    faults=Dict{Int64, Dict{Symbol,Array}}([])
    for i=1:length(fileNames)
        tmp=FSP3D.Read3DFSPresultsFromBinary(dirToRead*fileNames[i]*"/")
        faults=merge(faults, Dict(i => tmp))
    end

    return faults
end #Read3DFSPresultsFromBinaryAllFaults()

function ConcatenateResultsSeveralFaults(faults::Dict{Int64, Dict{Symbol, Array}}; keysToConcatenate=[:Pcritical, :Sv, :SHmax, :Shmin, :Pp, :depth])
    #=
      Concatenate the FSP results corresponding to several faults
    =#
    finalDict=Dict{Symbol, Array}([])

    nRows=[]
    nCols=[]
    for i=1:length(faults)
        #global nRows, nCols
        nRows=[nRows; size(faults[i][:Pcritical],1)]
        nCols=[nCols; size(faults[i][:Pcritical],2)]
    end
    totalNumberOfRows=sum(nRows)

    iStart=[1 nRows[1]]
    for k=2:length(nRows)
        #global iStart
        iStart=[iStart ; (iStart[k-1,2] + 1) (iStart[k-1,2] + nRows[k]) ]
    end

    #keysToConcatenate=[:Pcritical, :Sv, :SHmax, :Shmin, :Pp, :depth]

    for key in keysToConcatenate
        #global finalDict
        println("Reading on key "*string(key))

        if key != :coord
            tmpFinal=zeros(totalNumberOfRows,nCols[1])
        else
            tmpFinal=zeros(totalNumberOfRows,2)
        end

        @Threads.threads for iFault=1:length(faults)
            @inbounds tmpFinal[iStart[iFault,1]:iStart[iFault,2], :] .= faults[iFault][key]
        end
        finalDict=merge(finalDict, Dict(key => tmpFinal))
    end

    return finalDict

end ## ConcatenateResultsSeveralFaults

function Call3DFSPallFaults(fileName::String, userDefinedPdfs::Dict{Symbol, Dict{Symbol, Any}}, param::Dict{Symbol, Any})
    #=
      Call 3DFSP function
    =#

    #=
      Bug to fix: make sure that the IP3model dictionary does not contain unecessary information. Specially because the user might wnat to use the "userDefinePdfs" to define the probability and values.
    =#
    
    mainDir=param[:mainDir]
    inputType=param[:inputType]
    dirInputData=param[:dirInputData]
    if haskey(param, :IP3model)
        IP3model=param[:IP3model]
    end
    nSamples=param[:nSamples]
    probFailureAll=param[:probFailureAll]
    

    dirToSave=mainDir*"3DFSP/output/"*inputType*"/"
    FSP3D.CreateDir(dirToSave)
    
    ## Read and run the 3DFSP code
    surface, faultShapes = FSP3D.ReadFaultSurfacesExportedFromPetrel(dirInputData; fileNamesAll=[fileName])
    param[:unitsDict][:inputUnits][:depth] = surface[1].depthUnitSystem ## Assuming all surface have the same depth unit

    if haskey(param, :IP3model)
        FSP3D.ConverFaultSurfaceDepthUnitsTopDepthForIP3model(faultShapes, IP3model)
        param[:unitsDict][:inputUnits][:depth] = faultShapes[1][:units][:depth] ## Assuming all surfaces have the same depth unit

        ## Here, it assigns the interpolated properties to the fault segments
        FSP3D.ComputeInterpolatedPropertyAtFaultShapeSegments(faultShapes, IP3model)
    end

    ## Compute stress tensor based on the available data
    df=FSP3D.AssembleStressAndFaultInformationFromUser(faultShapes, param[:unitsDict][:inputUnits];)

    ## Run MCMC to compute the pore pressure to failure
    finalResults=FSP3D.ComputeCriticalPorePressureToFailure(df, userDefinedPdfs, nSamples::Int64)

    ## Compute the quantile from the pore pressure to failure
    ppFailureFinal = FSP3D.ComputePcCriticalUserDefinedProbabilityOfFailure(finalResults[:Pcritical], probFailureAll)

    ## Export results to Petrel
    ppFailureName="P" .* string.(probFailureAll)
    FSP3D.ExportResultsPpToFailureOnSurfaceToPetrel(dirToSave*"SurfacesToPetrel/", [fileName], surface, ppFailureFinal, ppFailureName, finalResults, df, param)

    FSP3D.Save3DFSPresultsToBinary(dirToSave*"bin/"*string(fileName)*"/", finalResults)

end #Call3DFSPallFaults

function ReadSurfacesExportedFromPetrel(fileName::String)
    #=
      Read surfaces exported from Petrel. It works for both faults and horizons
    =#

    header = FSP3D.ReadFaultSurfaceWithAttributeHeader(fileName)

    outputSurface=[]
    if isempty(header)
        ## It means that the user input surface is an horizon, not a fault surface
        outputSurface=FSP3D.ReadSurfaceExportedFromPetrelEarthVisionFormat(fileName::String)
        tmp=FSP3D.ReplaceHeaderInPetrelHorizonSurfaceWithFaultWithAttributesHeader(outputSurface.header)
        outputSurface.header[:] .= ""
        outputSurface.header[1:length(tmp)] = tmp
    else
        outputSurface=FSP3D.ReadFaultSurfaceWithAttributesFromPetrel(fileName)
    end

    if !hasproperty(outputSurface,:coord)
        throw(ErrorException("\n \t Error!!! The user specified fault surfaces does not exist!!! File name input = ", fileName))
    end

    return outputSurface

end #ReadSurfacesExportedFromPetrel()

function ReplaceHeaderInPetrelHorizonSurfaceWithFaultWithAttributesHeader(header::Array{AbstractString})
    #=
      To save the horizons in Petrel as surfaces with attributes, I need to replace the header that I read using the EarthVision formation with the same header that is used for the surfaces with attributes options
    =#

    #=
    newHeader="# Petrel Surface format
    # VERSION 1
    BEGIN HEADER
    # Field 1: X
    # Field 2: Y
    # Field 3: Z
    # Field 4: column
    # Field 5: row
    Unit in X and Y direction: Foot_US
    Unit in depth: ft
    Elevation depth
    Grid_size: 105 x 25
    STRUCTURED

    END HEADER"
    =#

    newHeader="# Petrel Surface format
    # VERSION 1
    BEGIN HEADER
    # Field 1: X
    # Field 2: Y
    # Field 3: Z
    # Field 4: column
    # Field 5: row
    Unit in X and Y direction: Foot_US
    Unit in depth: ft
    Elevation depth
    Grid_size: 105 x 25
    STRUCTURED

    END HEADER"

    #fileName=fileNamesAll[1]
    ## Read and run the 3DFSP code
    #surface, faultShapes = FSP3D.ReadFaultSurfacesExportedFromPetrel(dirFaultData; fileNamesAll=[fileName])
    #header=surface[1].header

    ind=findall(occursin.("Grid_size", header))

    tmp=string(header[ind][1])
    strFinal=strip(split(tmp,"#")[2])

    finalHeader=replace(newHeader, "Grid_size: 105 x 25" => strFinal)
    finalHeader=split(String(finalHeader), "\n")

    tmpFinal=String[]
    for i=1:length(finalHeader)
        #global tmpFinal
        tmpFinal=[tmpFinal ; String(finalHeader[i])]
    end
    finalHeader=tmpFinal
    
    return finalHeader

end #ReplaceHeaderInPetrelHorizonSurfaceWithFaultWithAttributesHeader()

function CreateEmpiricalDistributionBasedOnDataHistogram(inputData::AbstractArray)
    #=
      Create empirical distribution based on data histogram
    =#
    #data=faultDip
    #nSamples=1000000

    uvhist = fit(Histogram, inputData, nbins=length(inputData))
    uvdist = EmpiricalDistributions.UvBinnedDist(uvhist)

    return uvdist

end #CreateEmpiricalDistributionBasedOnDataHistogram

function GetFaultDipAndStrikeAllFaults(dirFaultData::String)
    #=
      Reads the fault data and create a distribution of fault dip and azimuth
    =#

    fileNamesAll=FSP3D.GetFaultNamesInDirectory(dirFaultData)

    faultDip=Real[]
    faultStrike=Real[]
    for fileName in fileNamesAll
        #global faultDip, faultStrike
        surface, faultShapes = FSP3D.ReadFaultSurfacesExportedFromPetrel(dirFaultData; fileNamesAll=[fileName])
        faultDip=[faultDip; surface[1].faultDip]
        faultStrike=[faultStrike ; surface[1].faultStrike]
    end

    return Dict(:faultDip => faultDip, :faultStrike => faultStrike)

end #GetFaultDipAndStrikeAllFaults

function CheckIfDirExist(dirName::String)
    #=
      Check if direcotry exist
    =#
    if !isdir(dirName)
        throw(ErrorException("\n \t Error!! Directory does not exist: "*string(dirName)))
    end
end #CheckIfDirExist(dirName::String)

function ConverFaultSurfaceDepthUnitsTopDepthForIP3model(faultShapes, IP3model)

    for iFault = 1 : length(faultShapes)
        if haskey(faultShapes[iFault], :units)
            if haskey(faultShapes[iFault][:units],:depth)
                if faultShapes[iFault][:units][:depth] == "m" && IP3model[:units][:depth] == "ft"
                    faultShapes[iFault][:coord][:,end] = faultShapes[iFault][:coord][:,end] .* 3.28 
                    
                    ## BIG BUG!! Don't change the coordinates here. Recall that this is a dictionary, so 
                    ## faultShapes[iFault][:segment][:coord] has the same memory location as faultShapes[iFault][:coord]
                    ## because I did not use "deepcopy"
                    #faultShapes[iFault][:segment][:coord][:,end] = faultShapes[iFault][:segment][:coord][:,end] .* 3.28
                    ##

                    faultShapes[iFault][:units][:depth] = "ft"
                elseif faultShapes[iFault][:units][:depth] == "ft" && IP3model[:units][:depth] == "m"
                    faultShapes[iFault][:coord][:,end] = faultShapes[iFault][:coord][:,end] ./ 3.28 
                    #faultShapes[iFault][:segment][:coord][:,end] = faultShapes[iFault][:segment][:coord][:,end] ./ 3.28
                    faultShapes[iFault][:units][:depth] = "m"
                else
                    @assert 1 > 2 "Error!! TEST1 = "*string(faultShapes[iFault][:units][:depth])*" TEST2 = "*string(IP3model[:units][:depth])
                end
            end
        end
    end

end #ConverFaultSurfaceDepthUnitsTopDepthForIP3model()

end ## module FSP3D

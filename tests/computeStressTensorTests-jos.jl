dir=@__FILE__
#mainDir="/" .* joinpath(split(dir,"/")[1:end-1])*"/"

additionalPath = "C:/Users/bakirtzisn/Desktop/FSP_dev_test/3D_FSP/"
push!(LOAD_PATH, additionalPath)

import FSP3D
#import FSP3Dplots


# test inputs for ComputeStressTensor_CS_Normal_Faults
Sv=60
Pp=10
μ=0.6
Aphi=0.5
Shmin=20
SHmax=0


df=FSP3D.ComputeStressTensor_CS_Normal_Faults(Sv::Real, Pp::Real, μ::Real, Aphi::Real; n=0, Shmin=Shmin, SHmax=SHmax)

println(df)
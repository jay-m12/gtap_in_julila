module GTAPinJulia

using JuMP, PATHSolver, CSV, DataFrames

PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")

include("macros.jl")

export @extract


include("gtap.jl")

export load_gtap_data


include("model.jl")

export GTAPModel

export JuMP

end
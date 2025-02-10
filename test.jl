include("src/GTAPinJulia.jl")

using .GTAPinJulia
using JuMP

#   Choose a database and model options
True = false
if True
    inFolder    = "Data/DBUG0"
    BaseName    = "DBUG0"
    resName     = :USA
    capName     = :cap
    nrsFact     = [:nrs]
    savfClosure = :capFix
else
    inFolder    = "Data/3x3"
    BaseName    = "3x3"
    resName     = :USA
    capName     = :CAP
    nrsFact     = [:xxx]
    savfClosure = :capFix
end
inscale  = 1e-6
popscale = 1e-3

#   Get the GTAP data, initialize the model variables and calibrate the model parameters
data = GTAPinJulia.load_gtap_data(; inFolder = inFolder, BaseName = BaseName, 
                                    inscale = inscale, popscale = popscale, resName = resName, 
                                    capName = capName, nrsFact = nrsFact)


#   Define the model
#   The capital closure options are:
#    1. RoRFlex--standard GTAP closure with foreign savings to equate expected rates of return
#    2. capFix--Foreign savings fixed (wrt to model numéraire)
#    3. capShrFix--Foreign savings fixed wrt to GDP (requires a residual region)
GTAPModel = GTAPinJulia.GTAPModel(data, resName = resName, 
                                capName = capName, savfClosure = savfClosure)

#   Fix the exogenous variables
GTAPinJulia.fix.(GTAPModel[:ptax], data[:ptax0]; force=true)
GTAPinJulia.fix.(GTAPModel[:tfe], data[:tfe0]; force=true)
GTAPinJulia.fix.(GTAPModel[:tinc],data[:tinc0]; force=true)

GTAPinJulia.fix.(GTAPModel[:tfd], data[:tfd0]; force=true)
GTAPinJulia.fix.(GTAPModel[:tfm], data[:tfm0]; force=true)
GTAPinJulia.fix.(GTAPModel[:tpd], data[:tpd0]; force=true)
GTAPinJulia.fix.(GTAPModel[:tpm], data[:tpm0]; force=true)
GTAPinJulia.fix.(GTAPModel[:tgd], data[:tgd0]; force=true)
GTAPinJulia.fix.(GTAPModel[:tgm], data[:tgm0]; force=true)
GTAPinJulia.fix.(GTAPModel[:tid], data[:tid0]; force=true)
GTAPinJulia.fix.(GTAPModel[:tim], data[:tim0]; force=true)

GTAPinJulia.fix.(GTAPModel[:txs], data[:txs0]; force=true)
GTAPinJulia.fix.(GTAPModel[:tmarg], data[:tmarg0]; force=true)
GTAPinJulia.fix.(GTAPModel[:tms], data[:tms0]; force=true)

GTAPinJulia.fix.(GTAPModel[:chi_qe], 1.0; force=true)
GTAPinJulia.fix.(GTAPModel[:pop], data[:pop0]; force=true)
GTAPinJulia.fix.(GTAPModel[:pnum], 1.0; force=true)

# Benchmark
GTAPinJulia.fix(GTAPModel[:pnum], 1; force=true)

GTAPinJulia.set_attribute(GTAPModel, "cumulative_iteration_limit", 10_000)
GTAPinJulia.optimize!(GTAPModel)

#   Include the output functions
include("src/reporting.jl")
using CSV, DataFrames

acts = data[:acts] 
result_df = name_of_operation(GTAPModel, data, 2025)

println(result_df)
display(result_df)

CSV.write("GTAP_results.csv", result_df)
println("Results saved to GTAP_results.csv")

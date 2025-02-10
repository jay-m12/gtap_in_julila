#  The GTAP model in Julia
#  The model is documented in Corong et al 2017 (https://jgea.org/ojs/index.php/jgea/article/view/47/30)
#  Contributors: Dominique van der Mensbrugghe with advice from Mitch Phillipson

#   Model pre-amble: location of files, model options, user functions

#   Options linked to a specific database
if false
    inFolder  = "v:/Julia/GTAPinJulia/Data/3x3"
    outFolder = "v:/Julia/GTAPinJulia/Data/3x3"
    BaseName = "3x3"
#   Provide name of the capital label -- N.B. Julia is case sensitive
    capName  = Symbol("CAP")
#   One region is set as the residual region, provide the label
    resName  = Symbol("USA")
#   One (or more) factors can be treated as sector-specific, provide the labels
    nrsFact  = [Symbol("NRS")]
else
    inFolder  = "v:/Julia/GTAPinJulia/Data/DBUG0"
    outFolder = "v:/Julia/GTAPinJulia/Data/DBUG0"
    BaseName = "DBUG0"
#   Provide name of the capital label -- N.B. Julia is case sensitive
    capName  = Symbol("cap")
#   One region is set as the residual region, provide the label
    resName  = Symbol("USA")
#   One (or more) factors can be treated as sector-specific, provide the labels
    nrsFact  = [Symbol("nrs")]
    
end

#   Model options
simName  = "COMP"
savfOptions = [:RoRFlex, :capFix, :capShrFix]
savfClosure = :capFix
inscale  = 1e-6
popscale = 1e-3
using JuMP, Complementarity, DataFrames, CSV
using PATHSolver
# using JLD2, FileIO
# save("test.jld2", Dict("arrays") = arrays)
# X = load("test.jld2")
PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")

function getSet(setName)
    #   Get the set definitions for the set named 'setname'
    #   Function returns the Julia DataFrame as read in from the CSV file,
    #       the set labels transformed into Symbol form, and a Dictionary that maps the
    #       set labels to an index (not sure if this will be used)
    set       = CSV.read(inFolder * "/" * BaseName * setName * ".csv", DataFrames.DataFrame)
    setDict   = Dict{String, Int}()
    setLabels = Array{Symbol}(undef, nrow(set))
    for (i,row) in enumerate(eachrow(set))
        #   Tried to use push! but got an error message, so initialized it with a fixed row dimension
        setLabels[i] = Symbol(row[1])
        merge!(setDict, Dict(String(row[1]) => i))
    end
    return setLabels, set, setDict
end

#   Overloaded getData function

function getData(parmName,dim1)
    #   Read the data from the CSV file
    csv_data = CSV.read(inFolder * "/" * BaseName * parmName * ".csv", DataFrames.DataFrame)
    #   Create the container
    x = Containers.@container([dim1], 0.0)
    #   Loop over each row and fill the array
    for row in eachrow(csv_data)
        x[Symbol(row[1])] = row[:value]
    end
    return x
end

function getData(parmName,dim1,dim2)
    #   Read the data from the CSV file
    csv_data = CSV.read(inFolder * "/" * BaseName * parmName * ".csv", DataFrames.DataFrame)
    #   Create the container
    x = Containers.@container([dim1,dim2], 0.0)
    #   Loop over each row and fill the array
    for row in eachrow(csv_data)
        x[Symbol(row[1]), Symbol(row[2])] = row[:value]
    end
    return x
end

function getData(parmName,dim1,dim2,dim3)
    #   Read the data from the CSV file
    csv_data = CSV.read(inFolder * "/" * BaseName * parmName * ".csv", DataFrames.DataFrame)
    #   Create the container
    x = Containers.@container([dim1,dim2,dim3], 0.0)
    #   Loop over each row and fill the array
    for row in eachrow(csv_data)
        x[Symbol(row[1]), Symbol(row[2]), Symbol(row[3])] = row[:value]
    end
    return x
end

function getData(parmName,dim1,dim2,dim3,dim4)
    #   Read the data from the CSV file
    csv_data = CSV.read(inFolder * "/" * BaseName * parmName * ".csv", DataFrames.DataFrame)
    #   Create the container
    x = Containers.@container([dim1,dim2,dim3,dim4], 0.0)
    #   Loop over each row and fill the array
    for row in eachrow(csv_data)
        x[Symbol(row[1]), Symbol(row[2]), Symbol(row[3]), Symbol(row[4])] = row[:value]
    end
    return x
end

#   Overloaded 'initVar' function

function initVar(ifZero, dim1)
    if ifZero
        x = Containers.@container([dim1], 0.0)
    else
        x = Containers.@container([dim1], 1.0)
    end
    return x
end

function initVar(ifZero, dim1, dim2)
    if ifZero
        x = Containers.@container([dim1,dim2], 0.0)
    else
        x = Containers.@container([dim1,dim2], 1.0)
    end
    return x
end

function initVar(ifZero, dim1, dim2, dim3)
    if ifZero
        x = Containers.@container([dim1,dim2,dim3], 0.0)
    else
        x = Containers.@container([dim1,dim2,dim3], 1.0)
    end
    return x
end

#   Read model sets and database

#   Get the set labels
#   Dropping the margins subset

acts, activities, acts_ndx = getSet("ACTS")
endw, endowments, endw_ndx = getSet("ENDW")
comm, commodities, comm_ndx = getSet("COMM")
reg, regions, reg_ndx = getSet("REG")

#   Get the GTAP data
vdfb = inscale * getData("VDFB", comm, acts, reg)
vdfp = inscale * getData("VDFP", comm, acts, reg)
vmfb = inscale * getData("VMFB", comm, acts, reg)
vmfp = inscale * getData("VMFP", comm, acts, reg)

evfb = inscale * getData("EVFB", endw, acts, reg)
evfp = inscale * getData("EVFP", endw, acts, reg)
evos = inscale * getData("EVOS", endw, acts, reg)

maks = inscale * getData("MAKS", comm, acts, reg)
makb = inscale * getData("MAKB", comm, acts, reg)
ptax = inscale * getData("PTAX", comm, acts, reg)

vdpb = inscale * getData("VDPB", comm, reg)
vdpp = inscale * getData("VDPP", comm, reg)
vmpb = inscale * getData("VMPB", comm, reg)
vmpp = inscale * getData("VMPP", comm, reg)

vdgb = inscale * getData("VDGB", comm, reg)
vdgp = inscale * getData("VDGP", comm, reg)
vmgb = inscale * getData("VMGB", comm, reg)
vmgp = inscale * getData("VMGP", comm, reg)

vdib = inscale * getData("VDIB", comm, reg)
vdip = inscale * getData("VDIP", comm, reg)
vmib = inscale * getData("VMIB", comm, reg)
vmip = inscale * getData("VMIP", comm, reg)

vxsb = inscale * getData("VXSB", comm, reg, reg)
vfob = inscale * getData("VFOB", comm, reg, reg)
vcif = inscale * getData("VCIF", comm, reg, reg)
vmsb = inscale * getData("VMSB", comm, reg, reg)

vst  = inscale * getData("VST", comm, reg)
vtwr = inscale * getData("VTWR", comm, comm, reg, reg)

vkb  = inscale * getData("VKB", reg)
vdep = inscale * getData("VDEP", reg)
save = inscale * getData("SAVE", reg)
pop0 = popscale * getData("POP", reg)

#   Get the GTAP parameters

ESUBT  = getData("ESUBT", acts, reg)
@. ifelse(ESUBT == 1.0, 1.01, ESUBT)
ESUBC  = getData("ESUBC", acts, reg)
@. ifelse(ESUBC == 1.0, 1.01, ESUBC)
ESUBVA = getData("ESUBVA", acts, reg)
@. ifelse(ESUBVA == 1.0, 1.01, ESUBVA)
ETRAQ  = -getData("ETRAQ", acts, reg)
ESUBQ  = getData("ESUBQ", comm, reg)
@. ifelse(ESUBQ == 1.0, 1.01, ESUBQ)
ESUBG  = getData("ESUBG", reg)
@. ifelse(ESUBG == 1.0, 1.01, ESUBG)
ESUBI  = getData("ESUBI", reg)
@. ifelse(ESUBI == 1.0, 1.01, ESUBI)
INCPAR = getData("INCPAR", comm, reg)
SUBPAR = getData("SUBPAR", comm, reg)
ESUBD  = getData("ESUBD", comm, reg)
@. ifelse(ESUBD == 1.0, 1.01, ESUBD)
ESUBM  = getData("ESUBM", comm, reg)
@. ifelse(ESUBM == 1.0, 1.01, ESUBM)
ETRAE  = getData("ETRAE", endw, reg)
ESUBS  = getData("ESUBS", comm)
RORFLEX = getData("RORFLEX", reg)

#   Initialize variables

#   Production

pfa0  = initVar(false, comm, acts, reg)
qfa0  = (vdfp .+ vmfp) ./ pfa0
pes0  = initVar(false, endw, acts, reg)
peb0  = initVar(false, endw, acts, reg)
pfe0  = initVar(false, endw, acts, reg)
qfe0  = evos ./ pes0
tinc0 = evfb .- evos
tfe0  = evfp .- evfb

tinc0 = @. ifelse(evfb .!= 0, tinc0/evfb,0)
peb0  = pes0 ./ (1 .- tinc0)
tfe0  = @. ifelse(evfb .!= 0, tfe0 / evfb, 0.0)
pfe0  = peb0 .* (1.0 .+ tfe0)

#   We should read this in future iterations
#   Elasticity of sector specific resources
#   Let's initialize to 0
ETAFF = initVar(true, endw, acts, reg)

pe0 = initVar(false, endw, reg)

mobFact = Containers.@container([endw], true)
for n in nrsFact
    for f in endw
        if n == f
            mobFact[f] = false
        end
    end
end

#   The expression below removes the labels
# qe0 = [(sum(pes0[f,a,r] * qfe0[f,a,r] for a in acts)) for f in endw, r in reg]
#   The labels come back here
# qe0 = qe0 ./ pe0
#   This is an alternative solution...
# qe0 = Containers.DenseAxisArray(sum.(pes0[f,acts,r] .* qfe0[f,acts,r] for f in endw, r in reg), endw, reg)
qe0 = sum.(pes0[f,acts,r] .* qfe0[f,acts,r] for f in endw, r in reg) ./ pe0

pint0 = initVar(false, acts, reg)
qint0 = sum.(pfa0[comm,a,r] .* qfa0[comm,a,r] for a in acts, r in reg) ./ pint0
pva0  = initVar(false, acts, reg)
qva0  = sum.(evfp[endw,a,r] for a in acts, r in reg) ./ pva0
po0   = initVar(false, acts, reg)
qo0   = (pint0 .* qint0) .+ (pva0 .* qva0) ./ po0

#   Make
ps0   = initVar(false, comm, acts, reg)
qca0  = maks ./ ps0
ptax0 = @. ifelse(qca0 .!= 0, (makb .- maks) ./ maks, 1)
pca0  = ps0 .* (1 .+ ptax0)
pds0  = initVar(false, comm, reg)
qc0   = sum.(pca0[i,acts,r] .* qca0[i,acts,r] for i in comm, r in reg) ./ pds0

#   Final demand
ppa0 = initVar(false, comm, reg)
qpa0 = (vdpp .+ vmpp) ./ ppa0
yp0 = Containers.DenseAxisArray(sum.(ppa0[comm,r] .* qpa0[comm,r] for r in reg), reg)
ppa0 = initVar(false, comm, reg)
conshr0 = initVar(false, comm, reg)
for i in comm, r in reg
    conshr0[i,r] = ppa0[i,r]*qpa0[i,r] / yp0[r]
end

pga0  = initVar(false, comm, reg)
qga0  = (vdgp .+ vmgp) ./ pga0
yg0   = Containers.DenseAxisArray(sum.(pga0[comm,r] .* qga0[comm,r] for r in reg), reg)
pgov0 = initVar(false, reg)
xg0   = yg0 ./ pgov0

pia0  = initVar(false, comm, reg)
qia0  = (vdip .+ vmip) ./ pia0
yi0   = Containers.DenseAxisArray(sum.(pia0[comm,r] .* qia0[comm,r] for r in reg), reg)
pinv0 = initVar(false, reg)
qinv0 = yi0 ./ pinv0

#   Armington module

tfd0 = @. ifelse(vdfb .!= 0, (vdfp .- vdfb) ./ vdfb,0)
tfm0 = @. ifelse(vmfb .!= 0, (vmfp .- vmfb) ./ vmfb,0)
tpd0 = @. ifelse(vdpb .!= 0, (vdpp .- vdpb) ./ vdpb,0)
tpm0 = @. ifelse(vmpb .!= 0, (vmpp .- vmpb) ./ vmpb,0)
tgd0 = @. ifelse(vdgb .!= 0, (vdgp .- vdgb) ./ vdgb,0)
tgm0 = @. ifelse(vmgb .!= 0, (vmgp .- vmgb) ./ vmgb,0)
tid0 = @. ifelse(vdib .!= 0, (vdip .- vdib) ./ vdib,0)
tim0 = @. ifelse(vmib .!= 0, (vmip .- vmib) ./ vmib,0)

pms0 = initVar(false, comm, reg)

#   Probably can be converted to something more succinct
pfd0 = (initVar(false, comm, acts, reg))
pfm0 = initVar(false, comm, acts, reg)
ppd0 = initVar(false, comm, reg)
ppm0 = initVar(false, comm, reg)
pgd0 = initVar(false, comm, reg)
pgm0 = initVar(false, comm, reg)
pid0 = initVar(false, comm, reg)
pim0 = initVar(false, comm, reg)
for i in comm
    for r in reg
        for a in acts
            pfd0[i,a,r] = (1 + tfd0[i,a,r]) * pds0[i,r]
            pfm0[i,a,r] = (1 + tfm0[i,a,r]) * pms0[i,r]
        end
        ppd0[i,r] = (1 + tpd0[i,r]) * pds0[i,r]
        ppm0[i,r] = (1 + tpm0[i,r]) * pms0[i,r]
        pgd0[i,r] = (1 + tgd0[i,r]) * pds0[i,r]
        pgm0[i,r] = (1 + tgm0[i,r]) * pms0[i,r]
        pid0[i,r] = (1 + tid0[i,r]) * pds0[i,r]
        pim0[i,r] = (1 + tim0[i,r]) * pms0[i,r]
    end
end

qfd0 = vdfp ./ pfd0
qfm0 = vmfp ./ pfm0
qpd0 = vdpp ./ ppd0
qpm0 = vmpp ./ ppm0
qgd0 = vdgp ./ pgd0
qgm0 = vmgp ./ pgm0
qid0 = vdip ./ pid0
qim0 = vmip ./ pim0

save0 = deepcopy(save)
savf0 = yi0 .- save0 .- vdep
#   Add any residual inconsistency to the largest capital account balance (in absolute value)
#   !!!! There is most likely a more elegant way to do this in Julia
resid = sum(savf0[r] for r in reg)
rmax  = 0
sfmax = 0
for r in reg
    if abs(savf0[r] > sfmax)
        global sfmax = abs(savf0[r])
        global rmax = r
    end
end
savf0[rmax] = savf0[rmax] .- resid

uepriv0 = Containers.DenseAxisArray(sum.(conshr0[comm,r] .* INCPAR[comm,r] for r in reg), reg)

taxrpc0  = Containers.DenseAxisArray(sum.((vdpp[comm,r] .- vdpb[comm,r]) .+ (vmpp[comm,r] .- vmpb[comm,r]) for r in reg), reg)
taxrgc0  = Containers.DenseAxisArray(sum.((vdgp[comm,r] .- vdgb[comm,r]) .+ (vmgp[comm,r] .- vmgb[comm,r]) for r in reg), reg)
taxric0  = Containers.DenseAxisArray(sum.((vdip[comm,r] .- vdib[comm,r]) .+ (vmip[comm,r] .- vmib[comm,r]) for r in reg), reg)
taxriu0  = Containers.DenseAxisArray(sum.((vdfp[comm,acts,r] .- vdfb[comm,acts,r]) .+ (vmfp[comm,acts,r] .- vmfb[comm,acts,r]) for r in reg), reg)
taxrfu0  = Containers.DenseAxisArray(sum.((evfp[endw,acts,r] .- evfb[endw,acts,r]) for r in reg), reg)
taxrout0 = Containers.DenseAxisArray(sum.((makb[comm,acts,r] .- maks[comm,acts,r]) for r in reg), reg)
taxrexp0 = Containers.DenseAxisArray(sum.((vfob[comm,r,reg] .- vxsb[comm,r,reg]) for r in reg), reg)
taxrimp0 = Containers.DenseAxisArray(sum.((vmsb[comm,reg,r] .- vcif[comm,reg,r]) for r in reg), reg)
indtax0  = taxrpc0 .+ taxrgc0 .+ taxric0 .+ taxriu0 .+ taxrfu0 .+ taxrout0 .+ taxrexp0 .+ taxrimp0

depr     = vdep ./ vkb
kb0      = deepcopy(vkb)
ke0      = (1 .- depr) .* kb0 .+ qinv0
pinv0    = initVar(false, reg)
fincome0 = Containers.DenseAxisArray(sum.((peb0[endw,acts,r] .* qfe0[endw,acts,r]) for r in reg), reg)
fincome0 = fincome0 .- depr .* pinv0 .* kb0
y0       = fincome0 .+ indtax0
gdpfc0   = deepcopy(y0)

dppriv = yp0 ./ y0
dpgov  = yg0 ./ y0
dpsave = (1.0 .- dppriv .- dpgov)
uelas0 = 1.0 ./ ((dppriv ./ uepriv0) .+ dpgov .+ dpsave)

u0  = initVar(false, reg)  
up0 = initVar(false, reg)  
ug0 = initVar(false, reg)  
us0 = initVar(false, reg)

au = ((up0 .^ dppriv) .* (ug0 .^ dpgov) .* (us0 .^ dpsave))
au = u0 ./ au
aus = us0 .* pop0 ./ save0

#   tmarg0 will be defined with respect to the base volume, not the vfob value
txs0   = @. ifelse(vxsb .!= 0, (vfob .- vxsb) ./ vxsb, 0)
tmarg0 = @. ifelse(vxsb .!= 0, (vcif .- vfob) ./ vxsb, 0)
tms0   = @. ifelse(vcif .!= 0, (vmsb .- vcif) ./ vcif, 0)

#   Set some tolerance level for tmarg0
tmarg0 = @. ifelse(abs.(tmarg0) > 1e-6, tmarg0, 0)
for i ∈ comm, s ∈ reg, d ∈ reg
    if tmarg0[i,s,d] == 0.0
        for m ∈ comm
            vtwr[m,i,s,d] = 0
        end
    end
end

pfob0 = initVar(false, comm, reg, reg)
for i in comm, s in reg, d in reg
    pfob0[i,s,d] = (1 + txs0[i,s,d]) * pds0[i,s]
end
pcif0 = pfob0 .+ tmarg0
pmds0 = pcif0 .* (1 .+ tms0)

pms0 = initVar(false, comm, reg)
qms0 = sum.(qfm0[i,acts,r] for i in comm, r in reg) .+ qpm0 .+ qgm0 .+ qim0
qxs0 = vfob ./ pfob0

ptrans0 = initVar(false, comm, reg, reg)
qst0    = vst ./ pds0
qds0    = sum.(qfd0[i,acts,r] for i in comm, r in reg) .+ qpd0 .+ qgd0 .+ qid0

pt0     = initVar(false, comm)
qtm0    = initVar(true, comm)
qtmfsd0 = deepcopy(vtwr)
for m in comm, i in comm, s in reg, d in reg
    qtmfsd0[m,i,s,d] = qtmfsd0[m,i,s,d] / pt0[m]
    qtm0[m] = qtm0[m] + qtmfsd0[m,i,s,d]
end

#   Capital account closure

rental0 = initVar(true, reg)
rental0 = sum.(pes0[capName,acts,r] .* qfe0[capName,acts,r] for r in reg) ./ kb0

rorc0 = rental0 ./ pinv0 .- depr
rore0 = rorc0 .* (kb0 ./ ke0) .^ RORFLEX
rorg0 = sum(rore0[r]*pinv0[r]*(qinv0[r] - depr[r]*kb0[r]) for r in reg) / sum(pinv0[r]*(qinv0[r] - depr[r]*kb0[r]) for r in reg)
risk  = rorg0 ./ rore0
ifRes = initVar(true, reg)
ifRes[resName] = 1
pinvLag  = initVar(false, reg)
psaveLag = initVar(false, reg)
invwgt   = initVar(true, reg)
savwgt   = initVar(true, reg)

for r in reg
    isum = sum(pinv0[d]*(qinv0[d] - depr[d]*kb0[d]) for d in reg)
    invwgt[r] = pinv0[r]*(qinv0[r] - depr[r]*kb0[r]) / isum
    ssum = sum(save[d] for d in reg)
    savwgt[r] = save[r] / ssum
end
pfact0  = initVar(false, reg)
pfactw0 = 1
Walras  = 0

#   Calibrate parameters

shr_qint = @. ifelse(qo0 .!= 0, (pint0 .* qint0) ./ (pint0 .* qint0 .+ pva0 .* qva0),0)
shr_qva  = @. ifelse(qo0 .!= 0, (pva0 .* qva0) ./ (pint0 .* qint0 .+ pva0 .* qva0),0)

shr_qfa = initVar(true, comm, acts, reg)
shr_qfe = initVar(true, endw, acts, reg)
shr_qcas = initVar(true, comm, acts, reg)
shr_qcad = initVar(true, comm, acts, reg)
shr_gov = initVar(true, comm, reg)
shr_inv = initVar(true, comm, reg)

for r in reg
    #   Production
    for a in acts
        for i in comm
            if qint0[a,r] != 0
                shr_qfa[i,a,r] = pfa0[i,a,r]*qfa0[i,a,r] / sum(pfa0[j,a,r]*qfa0[j,a,r] for j in comm)
            end
        end
        for f in endw
            if qva0[a,r] != 0
                shr_qfe[f,a,r] = pfe0[f,a,r]*qfe0[f,a,r] / sum(pfe0[e,a,r]*qfe0[e,a,r] for e in endw)
            end
        end
        for i in comm
            if qo0[a,r] != 0
                shr_qcas[i,a,r] = ps0[i,a,r]*qca0[i,a,r] / sum(ps0[j,a,r]*qca0[j,a,r] for j in comm)
            end
        end
    end
    for i in comm
        #   Make
        for a in acts
            if qc0[i,r] != 0
                shr_qcad[i,a,r] = pca0[i,a,r]*qca0[i,a,r] / sum(pca0[i,b,r]*qca0[i,b,r] for b in acts)
            end
        end
        #   Final demand
        if yg0[r] != 0
            shr_gov[i,r] = pga0[i,r]*qga0[i,r] / sum(pga0[j,r]*qga0[j,r] for j in comm)
        end
        if yi0[r] != 0
            shr_inv[i,r] = pia0[i,r]*qia0[i,r] / sum(pia0[j,r]*qia0[j,r] for j in comm)
        end
    end
end

#   Top level Armington shares

shr_qfd  = (pfd0 .* qfd0) ./ (pfd0 .* qfd0 .+ pfm0 .* qfm0)
shr_qfm  = (pfm0 .* qfm0) ./ (pfd0 .* qfd0 .+ pfm0 .* qfm0)
shr_qpd  = (ppd0 .* qpd0) ./ (ppd0 .* qpd0 .+ ppm0 .* qpm0)
shr_qpm  = (ppm0 .* qpm0) ./ (ppd0 .* qpd0 .+ ppm0 .* qpm0)
shr_qgd  = (pgd0 .* qgd0) ./ (pgd0 .* qgd0 .+ pgm0 .* qgm0)
shr_qgm  = (pgm0 .* qgm0) ./ (pgd0 .* qgd0 .+ pgm0 .* qgm0)
shr_qid  = (pid0 .* qid0) ./ (pid0 .* qid0 .+ pim0 .* qim0)
shr_qim  = (pim0 .* qim0) ./ (pid0 .* qid0 .+ pim0 .* qim0)

#   Second level Armington shares

shr_qxs  = initVar(true, comm, reg, reg)
for d in reg, i in comm
    totimp = sum(pmds0[i,s,d]*qxs0[i,s,d] for s in reg)
    if totimp != 0
        for s in reg
            shr_qxs[i,s,d] = pmds0[i,s,d] * qxs0[i,s,d] / totimp
        end
    end
end

#   Factor allocation shares

shr_qfes  = initVar(true, endw, acts, reg)
for r in reg, f in endw
    totinc = sum(pes0[f,a,r] * qfe0[f,a,r] for a in acts)
    if totinc != 0
        for a in acts
            shr_qfes[f,a,r] = pes0[f,a,r] * qfe0[f,a,r] / totinc
        end
    end
end

#   Transportation shares
shr_vtwr  = deepcopy(qtmfsd0)

for i in comm, s in reg, d in reg
    wsum = sum(vtwr[m,i,s,d] for m in comm)
    if wsum != 0
        #   There is transport for this node
        for m in comm
            shr_vtwr[m,i,s,d] = vtwr[m,i,s,d] / wsum 
        end
    else
        #   There is no transport for this node
        for m in comm
            shr_vtwr[m,i,s,d] = 0 
        end
    end
end

shr_qst  = initVar(true, comm, reg)
A_qtm    = initVar(true, comm)
for m in comm
    totmarg = sum(vst[m,r] for r in reg)
    if totmarg != 0
        for r in reg
            shr_qst[m,r] = vst[m,r] / totmarg
        end
        if ESUBS[m] == 1.0
            A_qtm[m] = pt0[m] / prod((pds0[m,r]/shr_qst[m,r])^shr_qst[m,r] for r in reg)
        else
            A_qtm[m] = 1.0
        end
    end
end

#   Private consumption

alphac = initVar(true, comm, reg)

for r in reg
    denom = sum(ifelse(SUBPAR[j,r] != 0.0, conshr0[j,r] / SUBPAR[j,r], 0.0) for j in comm)
    for i in comm
        if SUBPAR[i,r] != 0.0
            alphac[i,r] = (conshr0[i,r] / SUBPAR[i,r]) * (((yp0[r] / pop0[r]) / ppa0[i,r]) ^ SUBPAR[i,r]) * (up0[r]^(-INCPAR[i,r] * SUBPAR[i,r])) / denom
        end
    end
end

zshr0  = initVar(true, comm, reg)
for r in reg, i in comm
    zshr0[i,r] = alphac[i,r] * SUBPAR[i,r] * (up0[r] ^ (INCPAR[i,r] * SUBPAR[i,r])) * ((ppa0[i,r] / (yp0[r] / pop0[r])) ^ SUBPAR[i,r])    
end

#   Government utility
aug = pop0 ./ xg0

#   Define the model 

function GTAPModel()
    GTAPMod = JuMP.Model(PATHSolver.Optimizer)

    @variables GTAPMod begin
        qint[a = acts, r = reg], (start = 1)                                    # Aggregate intermediate demand
        qva[a = acts, r =  reg], (start = 1)                                    # Demand for aggregate value added bundle
        po[a = acts, r = reg], (start = 1)                                      # Unit cost of production
        qfa[i = comm, a = acts, r = reg], (start = 1)                           # (Armington) demand for inputs
        pint[a = acts, r = reg], (start = 1)                                    # Aggregate cost of inputs
        qfe[f = endw, a = acts, r = reg], (start = 1)                           # Factor demand
        pva[a = acts, r = reg], (start = 1)                                     # Aggregate price of value added
        qca[i = comm, a = acts, r = reg], (start = 1)                           # Commodity i supplied by a
        qo[a = acts, r = reg], (start = 1)                                      # Output by activity a
        pca[i = comm, a = acts, r = reg], (start = 1)                           # Tax-inclusive price of commodity i supplied by a
        pds[i = comm, r = reg], (start = 1)                                     # Price of supply of commodity i
        ps[i = comm, a = acts, r = reg], (start = 1)                            # Tax-exclusive price of commodity i supplied by a
        taxrout[r = reg], (start = 1)                                           # Value of output tax revenues
        taxrfu[r = reg], (start = 1)                                            # Value of factor-use tax revenues
        taxriu[r = reg], (start = 1)                                            # Value of tax revenues on intermediate demand
        taxrpc[r = reg], (start = 1)                                            # Value of tax revenues on private demand
        taxrgc[r = reg], (start = 1)                                            # Value of tax revenues on government demand
        taxric[r = reg], (start = 1)                                            # Value of tax revenues on investment demand
        taxrexp[r = reg], (start = 1)                                           # Value of tax revenues on exports
        taxrimp[r = reg], (start = 1)                                           # Value of tax revenues on imports
        indtax[r = reg], (start = 1)                                            # Total value of indirect tax revenues
        fincome[r = reg], (start = 1)                                           # Factor income at basic prices net of depreciation
        y[r = reg], (start = 1)                                                 # Regional income
        uepriv[r = reg], (start = uepriv0[r])                                   # Elasticity of cost wrt utility from private consumption
        uelas[r = reg], (start = uelas0[r])                                     # Elasticity of cost of utility wrt utility
        yp[r = reg], (start = 1)                                                # Private consumption expenditure
        yg[r = reg], (start = 1)                                                # Government consumption expenditure
        save[r = reg], (start = 1)                                              # Regional supply of nominal savings
        u[r = reg], (start = u0[r])                                             # Total utility (per capita)
        us[r = reg], (start = us0[r])                                           # Utility derived from savings
        zshr[i = comm, r = reg], (start = 1)                                    # Auxiliary share variable for private consumption
        conshr[i = comm, r = reg], (start = 1)                                  # Private budget shares
        qpa[i = comm, r = reg], (start = 1)                                     # Private consumption (at the Armington level)
        up[r = reg], (start = 1)                                                # Private utility
        qga[i = comm, r = reg], (start = 1)                                     # Government consumption (at the Armington level)
        pgov[r = reg], (start = 1)                                              # Government expenditure price deflator
        xg[r = reg], (start = 1)                                                # Aggregate volume of government expenditures
        ug[r = reg], (start = 1)                                                # Per capita utility derived from public expenditures
        qia[i = comm, r = reg], (start = 1)                                     # Investment expenditure (at the Armington level)
        pinv[r = reg], (start = 1)                                              # Investment expenditure price deflator
        qinv[r = reg], (start = 1)                                              # Aggregate volume of investment expenditures
        pfd[i = comm, a = acts, r = reg], (start = 1)                           # End user price for domestic intermediate goods
        pfm[i = comm, a = acts, r = reg], (start = 1)                           # End user price for imported intermediate goods
        ppd[i = comm, r = reg], (start = 1)                                     # End user price for domestic private goods
        ppm[i = comm, r = reg], (start = 1)                                     # End user price for imported private goods
        pgd[i = comm, r = reg], (start = 1)                                     # End user price for domestic government goods
        pgm[i = comm, r = reg], (start = 1)                                     # End user price for imported government goods
        pid[i = comm, r = reg], (start = 1)                                     # End user price for domestic investment goods
        pim[i = comm, r = reg], (start = 1)                                     # End user price for imported investment goods
        pfa[i = comm, a = acts, r = reg], (start = 1)                           # Armington price for intermediate goods
        ppa[i = comm, r = reg], (start = 1)                                     # Armington price for private goods
        pga[i = comm, r = reg], (start = 1)                                     # Armington price for government goods
        pia[i = comm, r = reg], (start = 1)                                     # Armington price for investment goods
        qfd[i = comm, a = acts, r = reg], (start = 1)                           # Demand for domestic intermediate goods
        qfm[i = comm, a = acts, r = reg], (start = 1)                           # Demand for imported intermediate goods
        qpd[i = comm, r = reg], (start = 1)                                     # Demand for domestic private goods
        qpm[i = comm, r = reg], (start = 1)                                     # Demand for imported private goods
        qgd[i = comm, r = reg], (start = 1)                                     # Demand for domestic government goods
        qgm[i = comm, r = reg], (start = 1)                                     # Demand for imported government goods
        qid[i = comm, r = reg], (start = 1)                                     # Demand for domestic investment goods
        qim[i = comm, r = reg], (start = 1)                                     # Demand for imported investment goods
        qms[i = comm, r = reg], (start = 1)                                     # Aggregate import demand
        qxs[i = comm, s = reg, d = reg], (start = 1)                            # Export/import from region s towards region d
        pms[i = comm, d = reg], (start = 1)                                     # Aggregate import price
        pfob[i = comm, s = reg, d = reg], (start = 1)                           # FOB export price
        pcif[i = comm, s = reg, d = reg], (start = 1)                           # CIF import price
        pmds[i = comm, s = reg, d = reg], (start = 1)                           # Tariff inclusive bilateral import price
        qds[i = comm, r = reg], (start = 1)                                     # Aggregate domestic demand for domestic goods (/x qst)
        qc[i = comm, r = reg], (start = 1)                                      # Market equilibrium for domestically produced goods
        pes[f = endw, a = acts, r = reg], (start = 1)                           # Factor market allocation and equilibrium
        pe[f = endw, r = reg], (start = 1)                                      # Aggregate factor price
        pfe[f = endw, a = acts, r = reg], (start = 1)                           # Factor returns at purchasers prices
        peb[f = endw, a = acts, r = reg], (start = 1)                           # Factor returns at basic prices
        qe[f = endw, r = reg], (start = 1)                                      # Total factor stock
        qtmfsd[m = comm, i = comm, s = reg, d = reg], (start = 1)               # Demand for service m to transport i from s to d
        ptrans[i = comm, s = reg, d = reg], (start = 1)                         # Average price to transport i from s to d
        qtm[m = comm], (start = 1)                                              # Total demand for service m
        qst[m = comm, r = reg], (start = 1)                                     # Regional supply of service m
        pt[m = comm], (start = 1)                                               # Aggregate world price of service m
        gdpfc[r = reg], (start = 1)                                             # Nominal GDP at factor cost
        kb[r = reg], (start = 1)                                                # Non-normalized aggregate capital stock
        ke[r = reg], (start = 1)                                                # End-of-period capital stock
        rental[r = reg], (start = 1)                                            # After-tax gross rate of return to capital
        rorc[r = reg], (start = 1)                                              # Net rate of return to capital
        rore[r = reg], (start = 1)                                              # Expected rate of return to capital
        savf[r = reg], (start = savf0[r])                                       # Default capital account closure--determination of capital flows
        rorg, (start = 1)                                                       # Default capital account closure--determination of the global rate of return
        yi[r = reg], (start = 1)                                                # Nominal gross investment
        chisave, (start = 1)                                                    # Savings adjustment factor
        psave[r = reg], (start = 1)                                             # Price index for regional savings
        pfact[r = reg], (start = 1)                                             # Average regional factor price
        pfactw, (start = 1)                                                     # Average global factor price
        Walras, (start = 0)                                                     # Check on Walras' Law

        #   (Normally) exogenous variables
        ptax[i = comm, a = acts, r = reg], (start = ptax0[i,a,r])               # Output tax
        tfe[f = endw, a = acts, r = reg], (start = tfe0[f,a,r])                 # Tax rate on factor use
        tinc[f = endw, a = acts, r = reg], (start = tinc0[f,a,r])               # Tax on factor returns
        tfd[i = comm, a = acts, r = reg], (start = tfd0[i,a,r])                 # Tax on intermediate use of domestic goods
        tfm[i = comm, a = acts, r = reg], (start = tfm0[i,a,r])                 # Tax on intermediate use of imported goods
        tpd[i = comm, r = reg], (start = tpd0[i,r])                             # Tax on private use of domestic goods
        tpm[i = comm, r = reg], (start = tpm0[i,r])                             # Tax on private use of imported goods
        tgd[i = comm, r = reg], (start = tgd0[i,r])                             # Tax on public use of domestic goods
        tgm[i = comm, r = reg], (start = tgm0[i,r])                             # Tax on public use of imported goods
        tid[i = comm, r = reg], (start = tid0[i,r])                             # Tax on investment use of domestic goods
        tim[i = comm, r = reg], (start = tim0[i,r])                             # Tax on investment use of imported goods
        txs[i = comm, s = reg, d = reg], (start = txs0[i,s,d])                  # Tax on exports by source and destination
        tmarg[i = comm, s = reg, d = reg], (start = tmarg0[i,s,d])              # Total transport margin by source and destination
        tms[i = comm, s = reg, d = reg], (start = tms0[i,s,d])                  # Tax on imports by source and destination
        pop[r = reg], (start = pop0[r])                                         # Population
        chi_qe[f = endw, r = reg], (start = 1)                                  # Shifter on aggregate factor stock
        pnum, (start = 1)                                                       # Model numeraire--fixed outside of model
    end
    #   Production module --------------------------------------------------------------------------

    #   Intermediate demand bundle
    @constraint(GTAPMod, eq_qint[a in acts, r in reg; qint0[a,r] != 0],
        qint[a,r] - qo[a,r]*(po[a,r]/pint[a,r])^ESUBT[a,r] ⟂ qint[a,r])
    for a ∈ acts, r ∈ reg
        if qint0[a,r] == 0
            fix(qint[a,r], 0, force=true)
        end
    end

    #   Value added bundle
    @constraint(GTAPMod, eq_qva[a in acts, r in reg; qva0[a,r] != 0],
        qva[a,r] - qo[a,r]*(po[a,r]/pva[a,r])^ESUBT[a,r] ⟂ qva[a,r])
    for a ∈ acts, r ∈ reg
        if qva0[a,r] == 0
            fix(qva[a,r], 0, force=true)
        end
    end
    
    #   Unit cost of production
    @constraint(GTAPMod, eq_po[a in acts, r in reg; qo0[a,r] != 0],
        po[a,r]^(1 - ESUBT[a,r]) - (shr_qint[a,r]*pint[a,r]^(1 - ESUBT[a,r]) + shr_qva[a,r]*pva[a,r]^(1 - ESUBT[a,r])) ⟂ po[a,r])
    for a ∈ acts, r ∈ reg
        if qo0[a,r] == 0
            fix(po[a,r], 1, force=true)
        end
    end

    #   Armington intermediate demand
    @constraint(GTAPMod, eq_qfa[i in comm, a in acts, r in reg; qfa0[i,a,r] != 0],
        qfa[i,a,r] - qint[a,r]*(pint[a,r]/pfa[i,a,r])^ESUBC[a,r] ⟂ qfa[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qfa0[i,a,r] == 0
            fix(qfa[i,a,r], 0, force=true)
        end
    end

     #   Aggregate price of inputs
     @constraint(GTAPMod, eq_pint[a in acts, r in reg; qint0[a,r] != 0],
        pint[a,r]^(1 - ESUBC[a,r]) - (sum(shr_qfa[i,a,r]*pfa[i,a,r]^(1 - ESUBC[a,r]) for i in comm)) ⟂ pint[a,r])
    for a ∈ acts, r ∈ reg
        if qint0[a,r] == 0
            fix(pint[a,r], 1, force=true)
        end
    end

    #   Factor demand
    @constraint(GTAPMod, eq_qfe[f in endw, a in acts, r in reg; qfe0[f,a,r] != 0],
        qfe[f,a,r] - qva[a,r]*(pva[a,r]/pfe[f,a,r])^ESUBVA[a,r] ⟂ qfe[f,a,r])
    for f ∈ endw, a ∈ acts, r ∈ reg
        if qfe0[f,a,r] == 0
            fix(qfe[f,a,r], 0, force=true)
        end
    end

     #   Aggregate price of value added
    @constraint(GTAPMod, eq_pva[a in acts, r in reg; qva0[a,r] != 0],
        pva[a,r]^(1 - ESUBVA[a,r]) - (sum(shr_qfe[f,a,r]*pfe[f,a,r]^(1 - ESUBVA[a,r]) for f in endw)) ⟂ pva[a,r])
    for a ∈ acts, r ∈ reg
        if qva0[a,r] == 0
            fix(pva[a,r], 1, force=true)
        end
    end

     #   Make module --------------------------------------------------------------------------------

    #   Supply of commodity i by activity a
    @constraint(GTAPMod, eq_qca[i in comm, a in acts, r in reg; qca0[i,a,r] != 0],
        if ETRAQ[a,r] != Inf
            qca[i,a,r] - qo[a,r]*(ps[i,a,r]/po[a,r])^ETRAQ[a,r]
        else
            ps[i,a,r] - po[a,r]
        end
        ⟂ qca[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qca0[i,a,r] == 0
            fix(qca[i,a,r], 0, force=true)
        end
    end

    #   Defines output by activity a
    @constraint(GTAPMod, eq_qo[a in acts, r in reg; qo0[a,r] != 0],
        if ETRAQ[a,r] != Inf
            po[a,r]^(1 + ETRAQ[a,r]) - (sum(shr_qcas[i, a, r]*ps[i,a,r]^(1 + ETRAQ[a,r]) for i in comm)) 
        else
            qo[a,r] - sum(shr_qcas[i, a, r]*qca[i,a,r] for i in comm)
        end
        ⟂ qo[a,r])
    for a ∈ acts, r ∈ reg
        if qo0[a,r] == 0
            fix(qo[a,r], 0, force=true)
        end
    end

    #  Demand for commodity i produced by activity a (inverse CES)
    @constraint(GTAPMod, eq_pca[i in comm, a in acts, r in reg; qca0[i,a,r] != 0],
        if ESUBQ[i,r] != Inf
            qca[i,a,r] - qc[i,r]*(pds[i,r]/pca[i,a,r])^ESUBQ[i,r] 
        else
            pca[i,a,r] - pds[i,r]
        end    
        ⟂ pca[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qca0[i,a,r] == 0
            fix(pca[i,a,r], 1, force=true)
        end
    end

    #   Defines supply price of commodity i
    @constraint(GTAPMod, eq_pds[i in comm, r in reg; qc0[i,r] != 0],
        if ESUBQ[i,r] != Inf
            pds[i,r]^(1 - ESUBQ[i,r]) - (sum(shr_qcad[i,a,r]*pca[i,a,r]^(1 - ESUBQ[i,r]) for a in acts))
        else
            qc[i,r] - sum(shr_qcad[i,a,r]*qca[i,a,r] for a in acts)    
        end
        ⟂ pds[i,r])
    for i ∈ comm, r ∈ reg
        if qc0[i,r] == 0
            fix(pds[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_ps[i in comm, a in acts, r in reg; qca0[i,a,r] != 0],
        pca[i,a,r] - ps[i,a,r]*(1 + ptax[i,a,r])*(ps0[i,a,r]/pca0[i,a,r]) ⟂ ps[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qca0[i,a,r] == 0
            fix(ps[i,a,r], 1, force=true)
        end
    end

    #   Income and its distribution ----------------------------------------------------------------

    #   Tax revenues from intermediate demand
    @constraint(GTAPMod, eq_taxriu[r in reg],
        taxriu[r]*taxriu0[r] - (sum(tfd[i,a,r]*pds[i,r]*qfd[i,a,r]*pds0[i,r]*qfd0[i,a,r]
                             +      tfm[i,a,r]*pms[i,r]*qfm[i,a,r]*pms0[i,r]*qfm0[i,a,r] for i in comm, a in acts)) ⟂ taxriu[r])

    #   Tax revenues from private demand
    @constraint(GTAPMod, eq_taxrpc[r in reg],
        taxrpc[r]*taxrpc0[r] - (sum(tpd[i,r]*pds[i,r]*qpd[i,r]*pds0[i,r]*qpd0[i,r]
                             +      tpm[i,r]*pms[i,r]*qpm[i,r]*pms0[i,r]*qpm0[i,r] for i in comm)) ⟂ taxrpc[r])

    #   Tax revenues from government demand
    @constraint(GTAPMod, eq_taxrgc[r in reg],
        taxrgc[r]*taxrgc0[r] - (sum(tgd[i,r]*pds[i,r]*qgd[i,r]*pds0[i,r]*qgd0[i,r]
                             +      tgm[i,r]*pms[i,r]*qgm[i,r]*pms0[i,r]*qgm0[i,r] for i in comm)) ⟂ taxrgc[r])

    #   Tax revenues from investment demand
    @constraint(GTAPMod, eq_taxric[r in reg],
        taxric[r]*taxric0[r] - (sum(tid[i,r]*pds[i,r]*qid[i,r]*pds0[i,r]*qid0[i,r]
                             +      tim[i,r]*pms[i,r]*qim[i,r]*pms0[i,r]*qim0[i,r] for i in comm)) ⟂ taxric[r])

    #   Output tax revenues
    @constraint(GTAPMod, eq_taxrout[r in reg],
       taxrout[r]*taxrout0[r] - (sum((ps0[i,a,r]*qca0[i,a,r])*ptax[i,a,r]*ps[i,a,r]*qca[i,a,r] for i in comm, a in acts)) ⟂ taxrout[r])

    #   Factor use tax revenues
    @constraint(GTAPMod, eq_taxrfu[r in reg],
       taxrfu[r]*taxrfu0[r] - (sum(tfe[f,a,r]*peb[f,a,r]*qfe[f,a,r]*peb0[f,a,r]*qfe0[f,a,r] for f in endw, a in acts)) ⟂ taxrfu[r])

    #   Export tax revenues
    @constraint(GTAPMod, eq_taxrexp[r in reg],
        taxrexp[r]*taxrexp0[r] - (sum(txs[i,r,d]*pds[i,r]*qxs[i,r,d]*pds0[i,r]*qxs0[i,r,d] for i in comm, d in reg)) ⟂ taxrexp[r])

    #   Import tax revenues
    @constraint(GTAPMod, eq_taxrimp[r in reg],
        taxrimp[r]*taxrimp0[r] - (sum(tms[i,s,r]*pcif[i,s,r]*qxs[i,s,r]*pcif0[i,s,r]*qxs0[i,s,r] for i in comm, s in reg)) ⟂ taxrimp[r])

    #   Total indirect taxes
    @constraint(GTAPMod, eq_indtax[r in reg],
        indtax[r] - (taxrpc[r]*(taxrpc0[r]/indtax0[r]) + taxrgc[r]*(taxrgc0[r]/indtax0[r]) + taxric[r]*(taxric0[r]/indtax0[r]) + taxriu[r]*(taxriu0[r]/indtax0[r])
            + taxrfu[r]*(taxrfu0[r]/indtax0[r]) + taxrout[r]*(taxrout0[r]/indtax0[r]) + taxrexp[r]*(taxrexp0[r]/indtax0[r]) + taxrimp[r]*(taxrimp0[r]/indtax0[r])) ⟂ indtax[r])

    @constraint(GTAPMod, eq_fincome[r in reg],
       fincome[r] - (sum(peb[f,a,r]*qfe[f,a,r]*(peb0[f,a,r]*qfe0[f,a,r]/fincome0[r]) for f in endw, a in acts) - depr[r]*pinv[r]*kb[r]*(pinv0[r]*kb0[r]/fincome0[r])) ⟂ fincome[r])

    @constraint(GTAPMod, eq_y[r in reg],
       y[r] - (fincome[r]*(fincome0[r]/y0[r]) + indtax[r]*(indtax0[r]/y0[r])) ⟂ y[r])

    #   Disposition of income and final demand -----------------------------------------------------

    #   Elasticity of expenditure wrt utility from private consumption
    @constraint(GTAPMod, eq_uepriv[r in reg],
        uepriv[r] - sum(conshr[i,r]*conshr0[i,r]*INCPAR[i,r] for i in comm) ⟂ uepriv[r])

    #   Elasticity of total expenditure wrt to utility
    @constraint(GTAPMod, eq_uelas[r in reg],
        uelas[r]*(dppriv[r]/uepriv[r] + dpgov[r] + dpsave[r]) - 1 ⟂ uelas[r])

    #   Aggregate nominal private consumption
    @constraint(GTAPMod, eq_yp[r in reg],
        yp[r] - dppriv[r]*(uelas[r]/uepriv[r])*y[r]*(y0[r]/yp0[r]) ⟂ yp[r])

    #   Aggregate nominal government consumption
    @constraint(GTAPMod, eq_yg[r in reg],
        yg[r] - dpgov[r]*uelas[r]*y[r]*(y0[r]/yg0[r]) ⟂ yg[r])

    #   Domestic supply of savings
    @constraint(GTAPMod, eq_qsave[r in reg],
        save[r]*(save0[r]/y0[r]) - dpsave[r]*uelas[r]*y[r] ⟂ save[r])

    #   Utility function
    @constraint(GTAPMod, eq_u[r in reg],
        log(u[r]) - (log(au[r]) +  dppriv[r]*log(up[r]) + dpgov[r]*log(ug[r]) + dpsave[r]*log(us[r])) ⟂ u[r])

    #   Utility from savings
    @constraint(GTAPMod, eq_us[r in reg],
        us[r] - (save0[r]*aus[r]/pop[r])*(save[r]/psave[r]) ⟂ us[r])

    #   Auxiliary consumption variable
    @constraint(GTAPMod, eq_zshr[i in comm, r in reg; qpa0[i,r] != 0],
        zshr[i,r] - (alphac[i,r]*SUBPAR[i,r]*(((ppa0[i,r]/yp0[r])^SUBPAR[i,r])/zshr0[i,r]))*(up[r]^(INCPAR[i,r]*SUBPAR[i,r]))*((pop[r]*ppa[i,r]/yp[r])^SUBPAR[i,r]) ⟂ zshr[i,r])
    for i ∈ comm, r ∈ reg
        if qpa0[i,r] == 0
            fix(zshr[i,r], 0, force=true)
        end
    end

    #   Private consumption budget shares
    @constraint(GTAPMod, eq_conshr[i in comm, r in reg; qpa0[i,r] != 0],
        conshr[i,r] - (zshr[i,r]*(zshr0[i,r]/conshr0[i,r]))/(sum(zshr[j,r]*zshr0[j,r] for j in comm)) ⟂ conshr[i,r])
    for i ∈ comm, r ∈ reg
        if qpa0[i,r] == 0
            fix(conshr[i,r], 0, force=true)
        end
    end

    #   Household demand for goods and services
    @constraint(GTAPMod, eq_qpa[i in comm, r in reg; qpa0[i,r] != 0],
        ppa[i,r]*qpa[i,r] - conshr[i,r]*yp[r]*((conshr0[i,r]*yp0[r])/(ppa0[i,r]*qpa0[i,r])) ⟂ qpa[i,r])
    for i ∈ comm, r ∈ reg
        if qpa0[i,r] == 0
            fix(qpa[i,r], 0, force=true)
        end
    end

    #   Private utility
    @constraint(GTAPMod, eq_up[r in reg],
        1 - sum(ifelse(zshr0[i,r] != 0.0, (zshr0[i,r]/SUBPAR[i,r])*zshr[i,r], 0.0) for i in comm) ⟂ up[r])

    #   Government demand for goods and services
    @constraint(GTAPMod, eq_qga[i in comm, r in reg; qga0[i,r] != 0],
        qga[i,r] - (xg[r] * (pgov[r] / pga[i,r])^ESUBG[r]) ⟂ qga[i,r])
    for i ∈ comm, r ∈ reg
        if qga0[i,r] == 0
            fix(qga[i,r], 0, force=true)
        end
    end

    #   Government expenditure price deflator
    @constraint(GTAPMod, eq_pgov[r in reg],
        pgov[r]*xg[r] - (sum(shr_gov[i,r]*pga[i,r]*qga[i,r] for i in comm)) ⟂ pgov[r])

    #   Aggregate volume of government expenditures
    @constraint(GTAPMod, eq_xg[r in reg],
        pgov[r]*xg[r] - yg[r]*(yg0[r]/(pgov0[r]*xg0[r])) ⟂ xg[r])

    @constraint(GTAPMod, eq_ug[r in reg],
        ug[r] - (aug[r]*xg0[r]/pop[r])*xg[r] ⟂ ug[r])

    #   Investment demand for goods and services
    @constraint(GTAPMod, eq_qia[i in comm, r in reg; qia0[i,r] != 0],
        qia[i,r] - (qinv[r] * (pinv[r] / pia[i,r])^ESUBI[r]) ⟂ qia[i,r])
    for i ∈ comm, r ∈ reg
        if qia0[i,r] == 0
            fix(qia[i,r], 0, force=true)
        end
    end

    #   Investment expenditure price deflator
    @constraint(GTAPMod, eq_pinv[r in reg],
        pinv[r]*qinv[r] - (sum(shr_inv[i,r]*pia[i,r]*qia[i,r] for i in comm)) ⟂ pinv[r])

    #   Aggregate volume of investment expenditures
    @constraint(GTAPMod, eq_qinv[r in reg],
        pinv[r]*qinv[r] - yi[r]*(yi0[r]/(pinv0[r]*qinv0[r])) ⟂ qinv[r])

    #   Armington module ---------------------------------------------------------------------------

    #   End use prices for domestic and imported goods
    #   Assumes all users face same basic prices: pds for domestic goods and pms for imported goods

    @constraint(GTAPMod, eq_pfd[i in comm, a in acts, r in reg; qfd0[i,a,r] != 0],
        pfd[i,a,r] - (1 + tfd[i,a,r]) * pds[i,r]*(pds0[i,r]/pfd0[i,a,r]) ⟂ pfd[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qfd0[i,a,r] == 0
            fix(pfd[i,a,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pfm[i in comm, a in acts, r in reg; qfm0[i,a,r] != 0],
        pfm[i,a,r] - (1 + tfm[i,a,r]) * pms[i,r]*(pms0[i,r]/pfm0[i,a,r]) ⟂ pfm[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qfm0[i,a,r] == 0
            fix(pfm[i,a,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_ppd[i in comm, r in reg; qpd0[i,r] != 0],
        ppd[i,r] - (1 + tpd[i,r]) * pds[i,r]*(pds0[i,r]/ppd0[i,r]) ⟂ ppd[i,r])
    for i ∈ comm, r ∈ reg
        if qpd0[i,r] == 0
            fix(ppd[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_ppm[i in comm, r in reg; qpm0[i,r] != 0],
        ppm[i,r] - (1 + tpm[i,r]) * pms[i,r]*(pms0[i,r]/ppm0[i,r]) ⟂ ppm[i,r])
    for i ∈ comm, r ∈ reg
        if qpm0[i,r] == 0
            fix(ppm[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pgd[i in comm, r in reg; qgd0[i,r] != 0],
        pgd[i,r] - (1 + tgd[i,r]) * pds[i,r]*(pds0[i,r]/pgd0[i,r]) ⟂ pgd[i,r])
    for i ∈ comm, r ∈ reg
        if qgd0[i,r] == 0
            fix(pgd[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pgm[i in comm, r in reg; qgm0[i,r] != 0],
        pgm[i,r] - (1 + tgm[i,r]) * pms[i,r]*(pms0[i,r]/pgm0[i,r]) ⟂ pgm[i,r])
    for i ∈ comm, r ∈ reg
        if qgm0[i,r] == 0
            fix(pgm[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pid[i in comm, r in reg; qid0[i,r] != 0],
        pid[i,r] - (1 + tid[i,r]) * pds[i,r]*(pds0[i,r]/pid0[i,r]) ⟂ pid[i,r])
    for i ∈ comm, r ∈ reg
        if qid0[i,r] == 0
            fix(pid[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pim[i in comm, r in reg; qim0[i,r] != 0],
        pim[i,r] - (1 + tim[i,r]) * pms[i,r]*(pms0[i,r]/pim0[i,r]) ⟂ pim[i,r])
    for i ∈ comm, r ∈ reg
        if qim0[i,r] == 0
            fix(pim[i,r], 1, force=true)
        end
    end

    #   Armington price

    @constraint(GTAPMod, eq_pfa[i in comm, a in acts, r in reg; qfa0[i,a,r] != 0],
        pfa[i,a,r]^(1 - ESUBD[i,r]) - (shr_qfd[i,a,r]*pfd[i,a,r]^(1 - ESUBD[i,r]) + shr_qfm[i,a,r]*pfm[i,a,r]^(1 - ESUBD[i,r])) ⟂ pfa[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qfa0[i,a,r] == 0
            fix(pfa[i,a,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_ppa[i in comm, r in reg; qpa0[i,r] != 0],
        ppa[i,r]^(1 - ESUBD[i,r]) - (shr_qpd[i,r]*ppd[i,r]^(1 - ESUBD[i,r]) + shr_qpm[i,r]*ppm[i,r]^(1 - ESUBD[i,r])) ⟂ ppa[i,r])
    for i ∈ comm, r ∈ reg
        if qpa0[i,r] == 0
            fix(ppa[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pga[i in comm, r in reg; qga0[i,r] != 0],
        pga[i,r]^(1 - ESUBD[i,r]) - (shr_qgd[i,r]*pgd[i,r]^(1 - ESUBD[i,r]) + shr_qgm[i,r]*pgm[i,r]^(1 - ESUBD[i,r])) ⟂ pga[i,r])
    for i ∈ comm, r ∈ reg
        if qga0[i,r] == 0
            fix(pga[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pia[i in comm, r in reg; qia0[i,r] != 0],
        pia[i,r]^(1 - ESUBD[i,r]) - (shr_qid[i,r]*pid[i,r]^(1 - ESUBD[i,r]) + shr_qim[i,r]*pim[i,r]^(1 - ESUBD[i,r])) ⟂ pia[i,r])
    for i ∈ comm, r ∈ reg
        if qia0[i,r] == 0
            fix(pia[i,r], 1, force=true)
        end
    end

    #   Top level sourcing
    @constraint(GTAPMod, eq_qfd[i in comm, a in acts, r in reg; qfd0[i,a,r] != 0],
        qfd[i,a,r] - (qfa[i,a,r]  * (pfa[i,a,r] / pfd[i,a,r])^ESUBD[i,r]) ⟂ qfd[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qfd0[i,a,r] == 0
            fix(qfd[i,a,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qfm[i in comm, a in acts, r in reg; qfm0[i,a,r] != 0],
        qfm[i,a,r] - (qfa[i,a,r]  * (pfa[i,a,r] / pfm[i,a,r])^ESUBD[i,r]) ⟂ qfm[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qfm0[i,a,r] == 0
            fix(qfm[i,a,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qpd[i in comm, r in reg; qpd0[i,r] != 0],
        qpd[i,r] - (qpa[i,r]  * (ppa[i,r] / ppd[i,r])^ESUBD[i,r]) ⟂ qpd[i,r])
    for i ∈ comm, r ∈ reg
        if qpd0[i,r] == 0
            fix(qpd[i,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qpm[i in comm, r in reg; qpm0[i,r] != 0],
        qpm[i,r] - (qpa[i,r]  * (ppa[i,r] / ppm[i,r])^ESUBD[i,r]) ⟂ qpm[i,r])
    for i ∈ comm, r ∈ reg
        if qpm0[i,r] == 0
            fix(qpm[i,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qgd[i in comm, r in reg; qgd0[i,r] != 0],
        qgd[i,r] - (qga[i,r]  * (pga[i,r] / pgd[i,r])^ESUBD[i,r]) ⟂ qgd[i,r])
    for i ∈ comm, r ∈ reg
        if qgd0[i,r] == 0
            fix(qgd[i,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qgm[i in comm, r in reg; qgm0[i,r] != 0],
        qgm[i,r] - (qga[i,r]  * (pga[i,r] / pgm[i,r])^ESUBD[i,r]) ⟂ qgm[i,r])
    for i ∈ comm, r ∈ reg
        if qgm0[i,r] == 0
            fix(qgm[i,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qid[i in comm, r in reg; qid0[i,r] != 0],
        qid[i,r] - (qia[i,r]  * (pia[i,r] / pid[i,r])^ESUBD[i,r]) ⟂ qid[i,r])
    for i ∈ comm, r ∈ reg
        if qid0[i,r] == 0
            fix(qid[i,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qim[i in comm, r in reg; qim0[i,r] != 0],
        qim[i,r] - (qia[i,r]  * (pia[i,r] / pim[i,r])^ESUBD[i,r]) ⟂ qim[i,r])
    for i ∈ comm, r ∈ reg
        if qim0[i,r] == 0
            fix(qim[i,r], 0, force=true)
        end
    end

    #   Second level Armington -- sourcing by region
    #   Aggregate import demand
    @constraint(GTAPMod, eq_qms[i in comm, r in reg; qms0[i,r] != 0],
        qms[i,r] - (sum(qfm[i,a,r]*(qfm0[i,a,r]/qms0[i,r]) for a in acts) + qpm[i,r]*(qpm0[i,r]/qms0[i,r]) 
                        + qgm[i,r]*(qgm0[i,r]/qms0[i,r]) + qim[i,r]*(qim0[i,r]/qms0[i,r])) ⟂ qms[i,r])
    for i ∈ comm, r ∈ reg
        if qms0[i,r] == 0
            fix(qms[i,r], 0, force=true)
        end
    end

    #   Demand for imports by region d sourced from region s
    @constraint(GTAPMod, eq_qxs[i in comm, s in reg, d in reg; qxs0[i,s,d] != 0],
        qxs[i,s,d] - (qms[i,d] * (pms[i,d] / pmds[i,s,d])^ESUBM[i,d]) ⟂ qxs[i,s,d])
    for i ∈ comm, s ∈ reg, d ∈ reg
        if qxs0[i,s,d] == 0
            fix(qxs[i,s,d], 0, force=true)
        end
    end

    #   Aggregate price of imports
    @constraint(GTAPMod, eq_pms[i in comm, d in reg; qms0[i,d] != 0],
        pms[i,d]^(1 - ESUBM[i,d]) - (sum(shr_qxs[i,s,d]*pmds[i,s,d]^(1 - ESUBM[i,d]) for s in reg)) ⟂ pms[i,d])
    for i ∈ comm, d ∈ reg
        if qms0[i,d] == 0
            fix(pms[i,d], 1, force=true)
        end
    end

    #   Bilateral prices ---------------------------------------------------------------------------

    #   Export price at FOB
    @constraint(GTAPMod, eq_pfob[i in comm, s in reg, d in reg; qxs0[i,s,d] != 0],
        pfob[i,s,d] - ((1 + txs[i,s,d])*pds[i,s]*(pds0[i,s]/pfob0[i,s,d])) ⟂ pfob[i,s,d])
    for i ∈ comm, s ∈ reg, d ∈ reg
        if qxs0[i,s,d] == 0
            fix(pfob[i,s,d], 1, force=true)
        end
    end

    #   Import price at CIF
    @constraint(GTAPMod, eq_pcif[i in comm, s in reg, d in reg; qxs0[i,s,d] != 0],
        pcif[i,s,d] - (pfob[i,s,d]*(pfob0[i,s,d]/pcif0[i,s,d]) + ptrans[i,s,d]*tmarg[i,s,d]*(ptrans0[i,s,d]/pcif0[i,s,d])) ⟂ pcif[i,s,d])
    for i ∈ comm, s ∈ reg, d ∈ reg
        if qxs0[i,s,d] == 0
            fix(pcif[i,s,d], 1, force=true)
        end
    end

    #   Import price tariff inclusive
    @constraint(GTAPMod, eq_pmds[i in comm, s in reg, d in reg; qxs0[i,s,d] != 0],
        pmds[i,s,d] - ((1 + tms[i,s,d])*pcif[i,s,d]*(pcif0[i,s,d]/pmds0[i,s,d])) ⟂ pmds[i,s,d])
    for i ∈ comm, s ∈ reg, d ∈ reg
        if qxs0[i,s,d] == 0
            fix(pmds[i,s,d], 1, force=true)
        end
    end

    #   Global transport services ------------------------------------------------------------------

    @constraint(GTAPMod, eq_qtmfsd[m in comm, i in comm, s in reg, d in reg; qtmfsd0[m,i,s,d] !=0],
        qtmfsd[m,i,s,d] - tmarg[i,s,d]*qxs[i,s,d]*(qxs0[i,s,d]/qtmfsd0[m,i,s,d]) ⟂ qtmfsd[m,i,s,d])
    for m ∈ comm, i ∈ comm, s ∈ reg, d ∈ reg
        if qtmfsd0[m,i,s,d] == 0
            fix(qtmfsd[m,i,s,d], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_ptrans[i in comm, s in reg, d in reg; tmarg0[i,s,d] != 0],
        ptrans[i,s,d] - sum(shr_vtwr[m,i,s,d]*pt[m] for m in comm) ⟂ ptrans[i,s,d])
    for i ∈ comm, s ∈ reg, d ∈ reg
        if tmarg0[i,s,d] == 0
            fix(ptrans[i,s,d], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_qtm[m in comm; qtm0[m] != 0],
        qtm[m] - sum(qtmfsd[m,i,s,d]*(qtmfsd0[m,i,s,d]/qtm0[m]) for i in comm, s in reg, d in reg) ⟂ qtm[m])
    for m ∈ comm
        if qtm0[m] == 0
            fix(qtm[m], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qst[m in comm, r in reg; qst0[m,r] != 0],
        qst[m,r] - (qtm[m] * (pt[m] / pds[m,r])^ESUBS[m]) ⟂ qst[m,r])
    for m ∈ comm, r ∈ reg
        if qst0[m,r] == 0
            fix(qst[m,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_pt[m in comm; qtm0[m] != 0],
        if ESUBS[m] == 1.0
            (pt[m] - (A_qtm[m]/pt0[m]) * prod((pds[m,r]*(pds0[m,r]/shr_qst[m,r]))^shr_qst[m,r] for r in reg))
        else
            (pt[m]^(1-ESUBS[m]) - sum(shr_qst[m,r]*pds[m,r]^(1-ESUBS[m]) for r in reg))  
        end
        ⟂ pt[m])
    for m ∈ comm
        if qtm0[m] == 0
            fix(pt[m], 1, force=true)
        end
    end

    #   Domestic goods market equilibrium ----------------------------------------------------------

    #   Total domestic demand for domestic goods
    @constraint(GTAPMod, eq_qds[i in comm, r in reg; qds0[i,r] != 0],
        qds[i,r] - (sum(qfd[i,a,r]*(qfd0[i,a,r]/qds0[i,r]) for a in acts) + qpd[i,r]*(qpd0[i,r]/qds0[i,r]) 
                    + qgd[i,r]*(qgd0[i,r]/qds0[i,r]) + qid[i,r]*(qid0[i,r])/qds0[i,r]) ⟂ qds[i,r])
    for i ∈ comm, r ∈ reg
        if qds0[i,r] == 0
            fix(qds[i,r], 0, force=true)
        end
    end

    #   Market clearing for domestically produced goods
    @constraint(GTAPMod, eq_qc[i in comm, r in reg; qc0[i,r] != 0],
        qc[i,r] - (qds[i,r]*(qds0[i,r]/qc0[i,r]) + sum(qxs[i,r,d]*(qxs0[i,r,d]/qc0[i,r]) for d in reg) + qst[i,r]*(qst0[i,r]/qc0[i,r])) ⟂ qc[i,r])
    for i ∈ comm, r ∈ reg
        if qc0[i,r] == 0
            fix(qc[i,r], 0, force=true)
        end
    end

    #   Factor market equilibrium ------------------------------------------------------------------

    #   Factor supply and equilibrium
    @constraint(GTAPMod, eq_pes[f in endw, a in acts, r in reg; qfe0[f,a,r] != 0],
        if mobFact[f]
            if ETRAE[f,r] != Inf
                qfe[f,a,r] - (qe[f,r] * (pes[f,a,r]/pe[f,r])^ETRAE[f,r])
            else
                pes[f,a,r] - pe[f,r]
            end
        else
            qfe[f,a,r] - (pes[f,a,r] / pfact[r])^ETAFF[f,a,r]
        end
        ⟂ pes[f,a,r])
    for f ∈ endw, a ∈ acts, r ∈ reg
        if qfe0[f,a,r] == 0
            fix(pes[f,a,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pe[f in endw, r in reg; qe0[f,r] != 0],
        if mobFact[f]
            if ETRAE[f,r] != Inf
                pe[f,r]^(1 + ETRAE[f,r]) - sum(shr_qfes[f,a,r]*pes[f,a,r]^(1 + ETRAE[f,r]) for a in acts)
            else
                qe[f,r] - sum(shr_qfes[f,a,r]*qfe[f,a,r] for a in acts)
            end
        else
            pe[f,r] - sum(shr_qfes[f,a,r]*pes[f,a,r] for a in acts)
        end
        ⟂ pe[f,r])
    for f ∈ endw, r ∈ reg
        if qe0[f,r] == 0
            fix(pe[f,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_qe[f in endw, r in reg; qe0[f,r] != 0],
    if mobFact[f]
        #   By default, aggregate factor stocks are fixed
        qe[f,r] - chi_qe[f,r]
    else
        #   For sector specific factor, this is purely definitional
        pe[f,r]*qe[f,r] - sum(shr_qfes[f,a,r]*pes[f,a,r]*qfe[f,a,r] for a in acts)
    end
    ⟂ qe[f,r])
    for f ∈ endw, r ∈ reg
        if qe0[f,r] == 0
            fix(qe[f,r], 0, force=true)
        end
    end

    #   Factor returns at purchasers prices
    @constraint(GTAPMod, eq_pfe[f in endw, a in acts, r in reg; qfe0[f,a,r] != 0],
        pfe[f,a,r] - (1 + tfe[f,a,r])*peb[f,a,r]*(peb0[f,a,r]/pfe0[f,a,r]) ⟂ pfe[f,a,r])
    for f ∈ endw, a ∈ acts, r ∈ reg
        if qfe0[f,a,r] == 0
            fix(pfe[f,a,r], 1, force=true)
        end
    end

    #   Factor returns at basic prices
    @constraint(GTAPMod, eq_peb[f in endw, a in acts, r in reg; qfe0[f,a,r] != 0],
        pes[f,a,r] - (1 - tinc[f,a,r])*peb[f,a,r]*(peb0[f,a,r]/pes0[f,a,r]) ⟂ peb[f,a,r])
    for f ∈ endw, a ∈ acts, r ∈ reg
        if qfe0[f,a,r] == 0
            fix(peb[f,a,r], 1, force=true)
        end
    end

    #   Macro aggregates
    @constraint(GTAPMod, eq_gdpfc[r in reg],
        gdpfc[r] - y[r] ⟂ gdpfc[r])

    #   Capital account closure --------------------------------------------------------------------

    #   Non-normalized capital stock--change is equal to change in normalized capital stock
    #   !!!! In the GEMPACK code, isn't VES(ENDWC,r) / GROSSCAP(r) always equal to 1⟂
    @constraint(GTAPMod, eq_kb[r in reg],
        kb[r] - qe[capName,r] ⟂ kb[r])

    #   End-of-period capital stock
    @constraint(GTAPMod, eq_ke[r in reg],
        ke[r] - ((1 - depr[r])*kb[r]*(kb0[r]/ke0[r]) + qinv[r]*(qinv0[r]/ke0[r])) ⟂ ke[r])

    #   After-tax gross rate of return to capital
    @constraint(GTAPMod, eq_rental[r in reg],
        rental[r]*kb[r] - (sum(pes[capName,a,r]*qfe[capName,a,r]*(pes0[capName,a,r]*qfe0[capName,a,r]/(rental0[r]*kb0[r])) for a in acts)) ⟂ rental[r])

    #   Net rate of return to capital
    @constraint(GTAPMod, eq_rorc[r in reg],
        rorc[r] - ((rental0[r]/(rorc0[r]*pinv0[r]))*rental[r]/pinv[r] - depr[r]/rorc0[r]) ⟂ rorc[r])

    #   Expected rate of return to capital
    @constraint(GTAPMod, eq_rore[r in reg],
        rore[r] - (((rorc0[r]/rore0[r])*(kb0[r]/ke0[r])^RORFLEX[r])*rorc[r]*(kb[r]/ke[r])^RORFLEX[r]) ⟂ rore[r])

    #   Default capital account closure--determination of capital flows
    
    if savfClosure == :RoRFlex
        #   Equal rates of return (for all regions)
        @constraint(GTAPMod, eq_savf[r in reg],
            rorg - risk[r]*rore[r]*(rore0[r]/rorg0) ⟂ savf[r])
    elseif savfClosure == :capFix
        #   Fixed flows in real terms for all regions save residual
        @constraint(GTAPMod, eq_savf[r in reg; r != resName],
            savf[r] - pnum*savf0[r] ⟂ savf[r])
    elseif savfClosure == :capShrFix
        #   Fixed flows wrt GDP for all regions save residual
        @constraint(GTAPMod, eq_savf[r in reg; r != resName],
            savf[r] - savf0[r]*gdpfc[r] ⟂ savf[r])
    end

    #   Foreign saving must add to 0
    if savfClosure == :RoRFlex
        #   Determines the global rate of return
        @constraint(GTAPMod, eq_rorg,
            sum(savf[r] for r in reg) ⟂ rorg)
    else
        #   Determines foreign saving for the residual region
        @constraint(GTAPMod, eq_savfRes[r in reg; r == resName],
            sum(savf[s] for s in reg) ⟂ savf[r])
        #   rorg is the weighted average of all the rore's
        @constraint(GTAPMod, eq_rorg,
            rorg*rorg0*sum((pinv[r]*qinv[r]*pinv0[r]*qinv0[r] - depr[r]*kb[r]*kb0[r]) for r in reg) 
                - sum(rore[r]*rore0[r]*(pinv[r]*qinv[r]*pinv0[r]*qinv0[r] - depr[r]*kb[r]*kb0[r]) for r in reg) ⟂ rorg)
    end

    #   Savings = investment
    @constraint(GTAPMod, eq_yi[r in reg],
        yi[r] - (save[r]*(save0[r]/yi0[r]) + savf[r]/yi0[r] + depr[r]*pinv[r]*kb[r]*(pinv0[r]*kb0[r]/yi0[r])) + ifRes[r]*Walras/yi0[r] ⟂ yi[r])

    @constraint(GTAPMod, eq_chisave,
        chisave - sum(invwgt[r]*pinv[r]/pinvLag[r] for r in reg) / sum(savwgt[r]*psave[r]/psaveLag[r] for r in reg) ⟂ chisave)

    @constraint(GTAPMod, eq_psave[r in reg],
        psave[r]/psaveLag[r] - chisave*(pinv[r]/pinvLag[r]) ⟂ psave[r])

    @constraint(GTAPMod, eq_pfact[r in reg],
        pfact[r] - pfact0[r]*sqrt((
            (sum(peb0[f,a,r]*qfe0[f,a,r]*peb[f,a,r] for f in endw, a in acts))/
            (sum(peb0[f,a,r]*qfe0[f,a,r] for f in endw, a in acts)))
            *(
            (sum(peb0[f,a,r]*qfe0[f,a,r]*peb[f,a,r]*qfe[f,a,r] for f in endw, a in acts))/
            (sum(peb0[f,a,r]*qfe0[f,a,r]*qfe[f,a,r] for f in endw, a in acts)))) ⟂ pfact[r])

    @constraint(GTAPMod, eq_pfactw,
            pfactw - pfactw0*sqrt((
                (sum(peb0[f,a,r]*qfe0[f,a,r]*peb[f,a,r] for f in endw, a in acts, r in reg))/
                (sum(peb0[f,a,r]*qfe0[f,a,r] for f in endw, a in acts, r in reg))
                )*(
                (sum(peb0[f,a,r]*qfe0[f,a,r]*peb[f,a,r]*qfe[f,a,r] for f in endw, a in acts, r in reg))/
                (sum(peb0[f,a,r]*qfe0[f,a,r]*qfe[f,a,r] for f in endw, a in acts, r in reg))
                )) ⟂ pfactw)

        #   Numeraire definition
    @constraint(GTAPMod, eq_Walras,
        pnum - pfactw ⟂ Walras)

    return GTAPMod
end

#   Function to save SAM as a CSV file

function saveSAM(fileName, outscale, append, results, year)
    #   Full file name
    #   Scale factor for the output SAM
    #   Boolean: true means add to the existing CSV file, false means create new CSV file
    #   Results from a model run
    #   'Year' of the model run
    if append
        file = open(fileName, "a")
    else
        file = open(fileName, "w")
        write(file,"Sim,Reg,RLab,CLab,Year,Value\n")
    end
   
    for r ∈ reg
        itax = 0.0
        vtax = 0.0
        ptax = 0.0
        mtax = 0.0
        etax = 0.0
        dtax = 0.0
        for a ∈ acts
            for i ∈ comm
                x = outscale*(value(results[:qfd][i,a,r])*value(results[:pds][i,r])*qfd0[i,a,r]*pds0[i,r]+value(results[:qfm][i,a,r])*value(results[:pms][i,r])*qfm0[i,a,r]*pms0[i,r])
                str = simName * "," * String(r) * "," * String(i) * "," * string(a) * ","  * string(year) * ","  * string(x) * "\n"
                write(file,str)
            end
            for f in endw
                x = outscale*value(GTAPMod[:qfe][f,a,r])*value(GTAPMod[:peb][f,a,r])*qfe0[f,a,r]*peb0[f,a,r]
                str = simName * "," * String(r) * "," * String(f) * "," * string(a) * ","  * string(year) * ","  * string(x) * "\n"
                write(file,str)
            end
            x = outscale*(sum(value(GTAPMod[:qfd][i,a,r])*value(GTAPMod[:pds][i,r])*qfd0[i,a,r]*pds0[i,r]*value(GTAPMod[:tfd][i,a,r]) 
                + value(GTAPMod[:qfm][i,a,r])*value(GTAPMod[:pms][i,r])*qfm0[i,a,r]*pms0[i,r]*value(GTAPMod[:tfm][i,a,r]) for i in comm))
            str = simName * "," * String(r) * "," * "itax" * "," * string(a) * ","  * string(year) * ","  * string(x) * "\n"
            write(file,str)
            itax += x
            x = outscale*(sum(value(GTAPMod[:qfe][f,a,r])*value(GTAPMod[:peb][f,a,r])*qfe0[f,a,r]*peb0[f,a,r]*value(GTAPMod[:tfe][f,a,r]) for f in endw))
            str = simName * "," * String(r) * "," * "vtax" * "," * String(a) * "," * string(year) * "," * string(x) * "\n"
            vtax += x
            write(file,str)

            #   Make
            for i in comm
                x = outscale*(value(GTAPMod[:qca][i,a,r])*value(GTAPMod[:ps][i,a,r])*qca0[i,a,r]*ps0[i,a,r])
                str = simName * "," * String(r) * "," * String(a) * "," * String(i) * "," * string(year) * "," * string(x) * "\n"
                write(file,str)
            end
        end
        for i in comm
            x = outscale*(sum(value(GTAPMod[:qca][i,a,r])*value(GTAPMod[:ps][i,a,r])*qca0[i,a,r]*ps0[i,a,r]*value(GTAPMod[:ptax][i,a,r]) for a in acts))
            str = simName * "," * String(r) * "," * "ptax" * "," * String(i) * "," * string(year) * "," * string(x) * "\n"
            write(file,str)
            ptax += x
        end
        #   Imports/exports
        for i in comm
            for s in reg
                x = outscale*(value(GTAPMod[:qxs][i,s,r])*value(GTAPMod[:pcif][i,s,r])*qxs0[i,s,r]*pcif0[i,s,r])
                str = simName * "," * String(r) * "," * String(s) * "," * String(i) * "," * string(year) * "," * string(x) * "\n"
                write(file,str)
            end
            # x = outscale*(sum(value(GTAPMod[:qxs][i,s,r])*(value(GTAPMod[:pcif][i,s,r])*pcif0[i,s,r] - value(GTAPMod[:pfob][i,s,r])*pfob0[i,s,r])*qxs0[i,s,r] for s in reg))
            # str = simName * "," * String(r) * "," * "tmg" * "," * String(i) * "," * string(year) * "," * string(x) * "\n"
            # write(file,str)
            x = outscale*(sum(value(GTAPMod[:qxs][i,s,r])*value(GTAPMod[:pcif][i,s,r])*qxs0[i,s,r]*pcif0[i,s,r]*value(GTAPMod[:tms][i,s,r]) for s in reg))
            str = simName * "," * String(r) * "," * "mtax" * "," * String(i) * "," * string(year) * "," * string(x) * "\n"
            write(file,str)
            mtax += x
            for d in reg
                x = outscale*(value(GTAPMod[:qxs][i,r,d])*value(GTAPMod[:pfob][i,r,d])*qxs0[i,r,d]*pfob0[i,r,d])
                str = simName * "," * String(r) * "," * String(i) * "," * String(d) * "," * string(year) * "," * string(x) * "\n"
                write(file,str)
            end
            x = outscale*(sum(value(GTAPMod[:qxs][i,r,d])*value(GTAPMod[:pds][i,r])*qxs0[i,r,d]*pds0[i,r]*value(GTAPMod[:txs][i,r,d]) for d in reg))
            str = simName * "," * String(r) * "," * "etax" * "," * String(i) * "," * string(year) * "," * string(x) * "\n"
            write(file,str)
            etax += x
        end
        for d in reg
            x = outscale*(sum(value(GTAPMod[:qxs][i,r,d])*value(GTAPMod[:pfob][i,r,d])*qxs0[i,r,d]*pfob0[i,r,d] for i in comm))
            str = simName * "," * String(r) * "," * String(d) * "," * "BoP" * "," * string(year) * "," * string(x) * "\n"
            write(file,str)
        end
        for s in reg
            x = outscale*(sum(value(GTAPMod[:qxs][i,s,r])*value(GTAPMod[:pcif][i,s,r])*qxs0[i,s,r]*pcif0[i,s,r] for i in comm))
            str = simName * "," * String(r) * "," * "BoP" * "," * String(s) * "," * string(year) * "," * string(x) * "\n"
            write(file,str)
        end
        #   Income distribution
        for f in endw
            x = outscale*(sum(value(GTAPMod[:qfe][f,a,r])*value(GTAPMod[:pes][f,a,r])*qfe0[f,a,r]*pes0[f,a,r] for a in acts))
            str = simName * "," * String(r) * "," * "RegY" * "," * String(f) * "," * string(year) * "," * string(x) * "\n"
            write(file,str)
            x = outscale*(sum(value(GTAPMod[:qfe][f,a,r])*value(GTAPMod[:peb][f,a,r])*qfe0[f,a,r]*peb0[f,a,r]*value(GTAPMod[:tinc][f,a,r]) for a in acts))
            str = simName * "," * String(r) * "," * "dtax" * "," * String(f) * "," * string(year) * "," * string(x) * "\n"
            write(file,str)
            dtax += x
        end
        x = outscale*(value(GTAPMod[:yp][r])*yp0[r])
        str = simName * "," * String(r) * "," * "hhd" * "," * "RegY" * "," * string(year) * "," * string(x) * "\n"
        write(file,str)
        x = outscale*(value(GTAPMod[:save][r])*save0[r])
        str = simName * "," * String(r) * "," * "inv" * "," * "RegY" * "," * string(year) * "," * string(x) * "\n"
        write(file,str)
        x = outscale*(value(GTAPMod[:yg][r])*yg0[r])
        str = simName * "," * String(r) * "," * "gov" * "," * "RegY" * "," * string(year) * "," * string(x) * "\n"
        write(file,str)
        x = outscale*(value(GTAPMod[:kb][r])*kb0[r]*value(GTAPMod[:pinv][r])*pinv0[r]*depr[r])
        str = simName * "," * String(r) * "," * "deprY" * "," * "RegY" * "," * string(year) * "," * string(x) * "\n"
        write(file,str)
        str = simName * "," * String(r) * "," * "inv" * "," * "deprY" * "," * string(year) * "," * string(x) * "\n"
        write(file,str)
        #   Domestic final demand
        for i in comm
            x = outscale*(value(GTAPMod[:qpd][i,r])*value(GTAPMod[:pds][i,r])*qpd0[i,r]*pds0[i,r]+value(GTAPMod[:qpm][i,r])*value(GTAPMod[:pms][i,r])*qpm0[i,r]*pms0[i,r])
            str = simName * "," * String(r) * "," * String(i) * "," * "hhd" * "," * string(year) * "," * string(x) * "\n"
            write(file,str)
            x = outscale*(value(GTAPMod[:qgd][i,r])*value(GTAPMod[:pds][i,r])*qgd0[i,r]*pds0[i,r]+value(GTAPMod[:qgm][i,r])*value(GTAPMod[:pms][i,r])*qgm0[i,r]*pms0[i,r])
            str = simName * "," * String(r) * "," * String(i) * "," * "gov" * "," * string(year) * "," * string(x) * "\n"
            write(file,str)
            x = outscale*(value(GTAPMod[:qid][i,r])*value(GTAPMod[:pds][i,r])*qid0[i,r]*pds0[i,r]+value(GTAPMod[:qim][i,r])*value(GTAPMod[:pms][i,r])*qim0[i,r]*pms0[i,r])
            str = simName * "," * String(r) * "," * String(i) * "," * "inv" * "," * string(year) * "," * string(x) * "\n"
            write(file,str)
        end
        x = outscale*(sum(value(GTAPMod[:qpd][i,r])*value(GTAPMod[:pds][i,r])*qpd0[i,r]*pds0[i,r]*value(GTAPMod[:tpd][i,r]) 
            + value(GTAPMod[:qpm][i,r])*value(GTAPMod[:pms][i,r])*qpm0[i,r]*pms0[i,r]*value(GTAPMod[:tpm][i,r]) for i in comm))
        str = simName * "," * String(r) * "," * "itax" * "," * "hhd" * "," * string(year) * "," * string(x) * "\n"
        write(file,str)
        itax += x
        x = outscale*(sum(value(GTAPMod[:qgd][i,r])*value(GTAPMod[:pds][i,r])*qgd0[i,r]*pds0[i,r]*value(GTAPMod[:tgd][i,r])
            + value(GTAPMod[:qgm][i,r])*value(GTAPMod[:pms][i,r])*qgm0[i,r]*pms0[i,r]*value(GTAPMod[:tgm][i,r]) for i in comm))
        str = simName * "," * String(r) * "," * "itax" * "," * "gov" * "," * string(year) * "," * string(x) * "\n"
        write(file,str)
        itax += x
        x = outscale*(sum(value(GTAPMod[:qid][i,r])*value(GTAPMod[:pds][i,r])*qid0[i,r]*pds0[i,r]*value(GTAPMod[:tid][i,r])
            + value(GTAPMod[:qim][i,r])*value(GTAPMod[:pms][i,r])*qim0[i,r]*pms0[i,r]*value(GTAPMod[:tim][i,r]) for i in comm))
        str = simName * "," * String(r) * "," * "itax" * "," * "inv" * "," * string(year) * "," * string(x) * "\n"
        write(file,str)
        itax += x
        #   Margins and BoP closure
        for m in comm
            x = outscale*(value(GTAPMod[:qst][m,r])*value(GTAPMod[:pds][m,r])*qst0[m,r]*pds0[m,r])
            str = simName * "," * String(r) * "," * String(m) * "," * "tmg" * "," * string(year) * "," * string(x) * "\n"
            write(file,str)
        end
        x = outscale*(sum(value(GTAPMod[:qst][m,r])*value(GTAPMod[:pds][m,r])*qst0[m,r]*pds0[m,r] for m in comm))
        str = simName * "," * String(r) * "," * "tmg" * "," * "BoP" * "," * string(year) * "," * string(x) * "\n"
        write(file,str)
        x = outscale*(value(GTAPMod[:savf][r]))
        str = simName * "," * String(r) * "," * "inv" * "," * "BoP" * "," * string(year) * "," * string(x) * "\n"
        write(file,str)
        #   Tax revenue closure
        str = simName * "," * String(r) * "," * "RegY" * "," * "itax" * "," * string(year) * "," * string(itax) * "\n"
        write(file,str)
        str = simName * "," * String(r) * "," * "RegY" * "," * "vtax" * "," * string(year) * "," * string(vtax) * "\n"
        write(file,str)
        str = simName * "," * String(r) * "," * "RegY" * "," * "ptax" * "," * string(year) * "," * string(ptax) * "\n"
        write(file,str)
        str = simName * "," * String(r) * "," * "RegY" * "," * "mtax" * "," * string(year) * "," * string(mtax) * "\n"
        write(file,str)
        str = simName * "," * String(r) * "," * "RegY" * "," * "etax" * "," * string(year) * "," * string(etax) * "\n"
        write(file,str)
        str = simName * "," * String(r) * "," * "RegY" * "," * "dtax" * "," * string(year) * "," * string(dtax) * "\n"
        write(file,str)
      end              

    close(file)
end

#   'Compile the model
GTAPMod = GTAPModel()
# set_attribute(GTAPMod, "cumulative_iteration_limit",0)

#   Fix the exogenous variables
fix.(GTAPMod[:ptax], ptax0; force=true)
fix.(GTAPMod[:tfe], tfe0; force=true)
fix.(GTAPMod[:tinc], tinc0; force=true)

fix.(GTAPMod[:tfd], tfd0; force=true)
fix.(GTAPMod[:tfm], tfm0; force=true)
fix.(GTAPMod[:tpd], tpd0; force=true)
fix.(GTAPMod[:tpm], tpm0; force=true)
fix.(GTAPMod[:tgd], tgd0; force=true)
fix.(GTAPMod[:tgm], tgm0; force=true)
fix.(GTAPMod[:tid], tid0; force=true)
fix.(GTAPMod[:tim], tim0; force=true)

fix.(GTAPMod[:txs], txs0; force=true)
fix.(GTAPMod[:tmarg], tmarg0; force=true)
fix.(GTAPMod[:tms], tms0; force=true)

fix.(GTAPMod[:chi_qe], 1.0; force=true)
fix.(GTAPMod[:pop], pop0; force=true)
fix.(GTAPMod[:pnum], 1.0; force=true)

#   Set up the number of shocks, typically ':0' is the benchmark

years = [:0, :1]

variables = [:Walras, :pe, :qe]
for t ∈ years
    if t == :1
        fix(GTAPMod[:pnum], 1.2; force=true)
        #fix.(GTAPMod[:chi_qe][:LAB, reg], 1.04, force=true)  
        #fix.(GTAPMod[:ptax], 0.95*ptax0, force=true)
    else
        fix(GTAPMod[:pnum], 1.0; force=true)
        #fix.(GTAPMod[:chi_qe][:LAB, reg], 1.0, force=true)
        #fix.(GTAPMod[:ptax], ptax0, force=true)         
    end

    if false && t == :0
        open("stdout.txt", "w") do io
            redirect_stdout(io) do
                println(GTAPMod)
            end
        end
    end
    if true && t == :0
        open("stdout.txt","w") do f
            for c in vcat([all_constraints(GTAPMod, t...) for t in list_of_constraint_types(GTAPMod)]...)
               println(f , c)
            end
        end
    end

    # set_silent(GTAPMod)
    optimize!(GTAPMod)

    println(termination_status(GTAPMod))
    println(raw_status(GTAPMod))
    println(solution_summary)

    if true
        for var in variables
            println("$var:")
            println(value.(GTAPMod[var]))
        end
    end

    if(t == :0)
        saveSAM(outFolder * "/" * simName * "SAM.csv", 1.0/inscale, false, GTAPMod, t)
    else
        saveSAM(outFolder * "/" * simName * "SAM.csv", 1.0/inscale, true, GTAPMod, t)
    end
end
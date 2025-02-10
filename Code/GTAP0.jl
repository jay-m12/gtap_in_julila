#  The GTAP model in Julia

inFolder = "V:/Julia/GTAPinJulia/Data/3x3"
outFolder = "V:/Julia/GTAPinJulia/Data/3x3"
BaseName = "3x3"
capName  = "CAP"
resName  = "USA"
simName  = "COMP"
inscale  = 1e-6
popscale = 1e-3
using JuMP, Complementarity, DataFrames, CSV
using PATHSolver
PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")

function getSet(setName)
    set = CSV.read(inFolder * "/" * BaseName * setName * ".csv", DataFrames.DataFrame)
    ndx = Dict{String, Int}()
    set1 = set[:,"Label"]
    for i in eachindex(set1)
        merge!(ndx, Dict(set1[i] => i))
    end
    return set, ndx
end

function getData1(parmName,dim1,ndx1)
    csv_data = CSV.read(inFolder * "/" * BaseName * parmName * ".csv", DataFrames.DataFrame)
    x = zeros((length(ndx1)))
    for row in eachrow(csv_data)
        x[get!(ndx1,row[dim1], 1)] = row[:value]
    end
    return x
end

function getData2(parmName,dim1,dim2,ndx1,ndx2)
    csv_data = CSV.read(inFolder * "/" * BaseName * parmName * ".csv", DataFrames.DataFrame)
    x = zeros((length(ndx1),length(ndx2)))
    for row in eachrow(csv_data)
        x[get!(ndx1,row[dim1], 1), get!(ndx2,row[dim2], 1)] = row[:value]
    end
    return x
end

function getData3(parmName,dim1,dim2,dim3,ndx1,ndx2,ndx3)
    csv_data = CSV.read(inFolder * "/" * BaseName * parmName * ".csv", DataFrames.DataFrame)
    x = zeros((length(ndx1),length(ndx2),length(ndx3)))
    for row in eachrow(csv_data)
        x[get!(ndx1,row[dim1], 1), get!(ndx2,row[dim2], 1), get!(ndx3,row[dim3], 1)] = row[:value]
    end
    return x
end

function getData4(parmName,dim1,dim2,dim3,dim4,ndx1,ndx2,ndx3,ndx4)
    csv_data = CSV.read(inFolder * "/" * BaseName * parmName * ".csv", DataFrames.DataFrame)
    x = zeros((length(ndx1),length(ndx2),length(ndx3),length(ndx4)))
    for row in eachrow(csv_data)
        x[get!(ndx1,row[dim1], 1), get!(ndx2,row[dim2], 1), get!(ndx3,row[dim3], 1), get!(ndx4,row[dim4], 1)] = row[:value]
    end
    return x
end

function initVar1(ifZero, ndx1)
    if ifZero
        x = zeros(length(ndx1))
    else
        x = ones(length(ndx1))
    end
    return x
end

function initVar2(ifZero, ndx1, ndx2)
    if ifZero
        x = zeros(length(ndx1),length(ndx2))
    else
        x = ones(length(ndx1),length(ndx2))
    end
    return x
end

function initVar3(ifZero, ndx1, ndx2, ndx3)
    if ifZero
        x = zeros(length(ndx1),length(ndx2),length(ndx3))
    else
        x = ones(length(ndx1),length(ndx2),length(ndx3))
    end
    return x
end

function setFlag(x)
    for i in eachindex(x)
        x[i] = (x[i] > 0.0) ? 1 : 0
    end
end

#   There are probably more efficient ways of inverting the dictionary
#   d is a dictionary and n is the index for which we are looking for the label (i.e., key)
function get_key(d,n)
    for item in d
        if item.second == n
            return(item.first)
        end
    end
end

#   Get the set labels

acts, acts_ndx = getSet("ACTS")
comm, comm_ndx = getSet("COMM")
marg, marg_ndx = getSet("MARG")
endw, endw_ndx = getSet("ENDW")
reg, reg_ndx   = getSet("REG")

activities  = [i for i=1:length(acts_ndx)]
commodities = [i for i=1:length(comm_ndx)]
margins     = [i for i=1:length(marg_ndx)]
factors     = [i for i=1:length(endw_ndx)]
regions     = [i for i=1:length(reg_ndx)]

#   Get the index for the capital factor ("CAP")
cap = get!(endw_ndx, capName, 1)
resReg = get!(reg_ndx, resName, 1)

#   Get the GTAP data
vdfb = inscale * getData3("VDFB", "COMM", "ACTS", "REG", comm_ndx, acts_ndx, reg_ndx)
vdfp = inscale * getData3("VDFP", "COMM", "ACTS", "REG", comm_ndx, acts_ndx, reg_ndx)
vmfb = inscale * getData3("VMFB", "COMM", "ACTS", "REG", comm_ndx, acts_ndx, reg_ndx)
vmfp = inscale * getData3("VMFP", "COMM", "ACTS", "REG", comm_ndx, acts_ndx, reg_ndx)

evfb = inscale * getData3("EVFB", "ENDW", "ACTS", "REG", endw_ndx, acts_ndx, reg_ndx)
evfp = inscale * getData3("EVFP", "ENDW", "ACTS", "REG", endw_ndx, acts_ndx, reg_ndx)
evos = inscale * getData3("EVOS", "ENDW", "ACTS", "REG", endw_ndx, acts_ndx, reg_ndx)

maks = inscale * getData3("MAKS", "COMM", "ACTS", "REG", comm_ndx, acts_ndx, reg_ndx)
makb = inscale * getData3("MAKB", "COMM", "ACTS", "REG", comm_ndx, acts_ndx, reg_ndx)
ptax = inscale * getData3("PTAX", "COMM", "ACTS", "REG", comm_ndx, acts_ndx, reg_ndx)

vdpb = inscale * getData2("VDPB", "COMM", "REG", comm_ndx, reg_ndx)
vdpp = inscale * getData2("VDPP", "COMM", "REG", comm_ndx, reg_ndx)
vmpb = inscale * getData2("VMPB", "COMM", "REG", comm_ndx, reg_ndx)
vmpp = inscale * getData2("VMPP", "COMM", "REG", comm_ndx, reg_ndx)

vdgb = inscale * getData2("VDGB", "COMM", "REG", comm_ndx, reg_ndx)
vdgp = inscale * getData2("VDGP", "COMM", "REG", comm_ndx, reg_ndx)
vmgb = inscale * getData2("VMGB", "COMM", "REG", comm_ndx, reg_ndx)
vmgp = inscale * getData2("VMGP", "COMM", "REG", comm_ndx, reg_ndx)

vdib = inscale * getData2("VDIB", "COMM", "REG", comm_ndx, reg_ndx)
vdip = inscale * getData2("VDIP", "COMM", "REG", comm_ndx, reg_ndx)
vmib = inscale * getData2("VMIB", "COMM", "REG", comm_ndx, reg_ndx)
vmip = inscale * getData2("VMIP", "COMM", "REG", comm_ndx, reg_ndx)

vxsb = inscale * getData3("VXSB", "COMM", "SRC", "DST", comm_ndx, reg_ndx, reg_ndx)
vfob = inscale * getData3("VFOB", "COMM", "SRC", "DST", comm_ndx, reg_ndx, reg_ndx)
vcif = inscale * getData3("VCIF", "COMM", "SRC", "DST", comm_ndx, reg_ndx, reg_ndx)
vmsb = inscale * getData3("VMSB", "COMM", "SRC", "DST", comm_ndx, reg_ndx, reg_ndx)

vst  = inscale * getData2("VST", "COMM", "REG", comm_ndx, reg_ndx)
vtwr = inscale * getData4("VTWR", "MARG", "COMM", "SRC", "DST", comm_ndx, comm_ndx, reg_ndx, reg_ndx)

vkb  = inscale * getData1("VKB", "REG", reg_ndx)
vdep = inscale * getData1("VDEP", "REG", reg_ndx)
save = inscale * getData1("SAVE", "REG", reg_ndx)
pop0 = popscale * getData1("POP", "REG", reg_ndx)

#   Get the GTAP parameters

ESUBT  = getData2("ESUBT", "ACTS", "REG", acts_ndx, reg_ndx)
ESUBC  = getData2("ESUBC", "ACTS", "REG", acts_ndx, reg_ndx)
ESUBVA = getData2("ESUBVA", "ACTS", "REG", acts_ndx, reg_ndx)
ETRAQ  = -getData2("ETRAQ", "ACTS", "REG", acts_ndx, reg_ndx)
ESUBQ  = getData2("ESUBQ", "COMM", "REG", comm_ndx, reg_ndx)
ESUBG  = getData1("ESUBG", "REG", reg_ndx)
ESUBI  = getData1("ESUBI", "REG", reg_ndx)
INCPAR = getData2("INCPAR", "COMM", "REG", comm_ndx, reg_ndx)
SUBPAR = getData2("SUBPAR", "COMM", "REG", comm_ndx, reg_ndx)
ESUBD  = getData2("ESUBD", "COMM", "REG", comm_ndx, reg_ndx)
ESUBM  = getData2("ESUBM", "COMM", "REG", comm_ndx, reg_ndx)
ETRAE  = getData2("ETRAE", "ENDW", "REG", endw_ndx, reg_ndx)
ESUBS  = getData1("ESUBS", "MARG", comm_ndx)
RORFLEX = getData1("RORFLEX", "REG", reg_ndx)

#   Initialize variables

#   Production

pfa0  = initVar3(false, comm_ndx, acts_ndx, reg_ndx)
qfa0  = (vdfp + vmfp) ./ pfa0
pes0  = initVar3(false, endw_ndx, acts_ndx, reg_ndx)
peb0  = initVar3(false, endw_ndx, acts_ndx, reg_ndx)
pfe0  = initVar3(false, endw_ndx, acts_ndx, reg_ndx)
qfe0  = evos ./ pes0
tinc0 = evfb .- evos
tfe0  = evfp .- evfb
for i in eachindex(qfe0)
    if (evfb[i] > 0)
        tinc0[i] = (tinc0[i] / evfb[i])
        peb0[i]  = pes0[i] / (1 - tinc0[i])
        tfe0[i]  = (tfe0[i] / evfb[i])
        pfe0[i]  = peb0[i] * (1 + tfe0[i])
    end
end
pe0 = initVar2(false, endw_ndx, reg_ndx)
qe0 = dropdims(sum(pes0 .* qfe0, dims=2), dims=2) ./ pe0

ifPerf = initVar2(false, endw_ndx, reg_ndx)
for r in regions
    for f in factors
        if ETRAE[f,r] == Inf
            ifPerf[f,r] = 1
        else
            ifPerf[f,r] = 0
        end
    end
end

pint0 = initVar2(false, acts_ndx, reg_ndx)
qint0 = dropdims(sum(qfa0, dims=1), dims=1) ./ pint0
pva0  = initVar2(false, acts_ndx, reg_ndx)
qva0  = dropdims(sum(evfp, dims=1), dims=1) ./ pva0
po0   = initVar2(false, acts_ndx, reg_ndx)
qo0   = (pint0 .* qint0) .+ (pva0 .* qva0) ./ po0

#   Make
ps0   = initVar3(false, comm_ndx, acts_ndx, reg_ndx)
qca0  = maks ./ ps0
ptax0 = makb .- maks
for i in eachindex(qca0)
    if(qca0[i] > 0)
        ptax0[i] = (ptax0[i] ./ (ps0[i] .* qca0[i]))
    end
end
pca0 = ps0 .* (1 .+ ptax0)
pds0  = initVar2(false, comm_ndx, reg_ndx)
qc0   = dropdims(sum(pca0 .* qca0, dims=2), dims=2) ./ pds0

#   Final demand
ppa0    = initVar2(false, comm_ndx, reg_ndx)
qpa0    = (vdpp + vmpp) ./ ppa0
yp0     = dropdims(sum(ppa0 .* qpa0, dims=1), dims=1)
conshr0 = ppa0 .* qpa0 ./ sum(ppa0 .* qpa0, dims=1)

pga0  = initVar2(false, comm_ndx, reg_ndx)
qga0  = (vdgp + vmgp) ./ pga0
yg0   = dropdims(sum(pga0 .* qga0, dims=1), dims=1)
pgov0 = initVar1(false, reg_ndx)
xg0   = yg0 ./ pgov0

pia0  = initVar2(false, comm_ndx, reg_ndx)
qia0  = (vdip + vmip) ./ pia0
yi0   = dropdims(sum(pia0 .* qia0, dims=1), dims=1)
pinv0 = initVar1(false, reg_ndx)
qinv0 = yi0 ./ pinv0

#   Armington module

tfd0 = (vdfp .- vdfb) ./ vdfb
tfm0 = (vmfp .- vmfb) ./ vmfb
tpd0 = (vdpp .- vdpb) ./ vdpb
tpm0 = (vmpp .- vmpb) ./ vmpb
tgd0 = (vdgp .- vdgb) ./ vdgb
tgm0 = (vmgp .- vmgb) ./ vmgb
tid0 = (vdip .- vdib) ./ vdib
tim0 = (vmip .- vmib) ./ vmib

pms0 = initVar2(false, comm_ndx, reg_ndx)

#   Probably can be converted to something more succinct
pfd0 = initVar3(false, comm_ndx, acts_ndx, reg_ndx)
pfm0 = initVar3(false, comm_ndx, acts_ndx, reg_ndx)
ppd0 = initVar2(false, comm_ndx, reg_ndx)
ppm0 = initVar2(false, comm_ndx, reg_ndx)
pgd0 = initVar2(false, comm_ndx, reg_ndx)
pgm0 = initVar2(false, comm_ndx, reg_ndx)
pid0 = initVar2(false, comm_ndx, reg_ndx)
pim0 = initVar2(false, comm_ndx, reg_ndx)
for i in commodities
    for r in regions
        for a in activities
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
savf0 = dropdims(sum(pia0 .* qia0, dims=1), dims=1) - save0 - vdep
#   Add any residual inconsistency to the largest capital account balance (in absolute value)
#   !!!! There is most likely a more elegant way to do this in Julia
resid = dropdims(sum(savf0, dims=1), dims=1)
rmax = 0
sfmax = 0
for r in eachindex(savf0)
    if abs(savf0[r] > sfmax)
        global sfmax = abs(savf0[r])
        global rmax = r
    end
end
savf0[rmax] = savf0[rmax] .- resid

uepriv0 = dropdims(sum(conshr0 .* INCPAR, dims=1), dims=1)

taxrpc0  = dropdims(sum((vdpp .- vdpb) + (vmpp .- vmpb), dims=1), dims=1)
taxrgc0  = dropdims(sum((vdgp .- vdgb) + (vmgp .- vmgb), dims=1), dims=1)
taxric0  = dropdims(sum((vdip .- vdib) + (vmip .- vmib), dims=1), dims=1)
taxriu0  = dropdims(sum(dropdims(sum((vdfp .- vdfb) + (vmfp .- vmfb), dims=1), dims=1), dims=1), dims=1)
taxrfu0  = dropdims(sum(dropdims(sum((evfp .- evfb), dims=1), dims=1), dims=1), dims=1)
taxrout0 = dropdims(sum(dropdims(sum((makb .- maks), dims=1), dims=1), dims=1), dims=1)
taxrexp0 = dropdims(sum(dropdims(sum((vfob .- vxsb), dims=3), dims=3), dims=1), dims=1)
taxrimp0 = dropdims(sum(dropdims(sum((vmsb .- vcif), dims=2), dims=2), dims=1), dims=1)
indtax0  = taxrpc0 .+ taxrgc0 .+ taxric0 .+ taxriu0 .+ taxrfu0 .+ taxrout0 .+ taxrexp0 .+ taxrimp0

depr     = vdep ./ vkb
kb0      = deepcopy(vkb)
ke0      = (1 .- depr) .* kb0 .+ qinv0
pinv0    = initVar1(false, reg_ndx)
fincome0 = dropdims(sum(dropdims(sum(peb0 .* qfe0, dims=1), dims=1),dims=1),dims=1) - depr .* pinv0 .* kb0
y0       = fincome0 .+ indtax0

dppriv = yp0 ./ y0
dpgov  = yg0 ./ y0
dpsave = (1 .- dppriv .- dpgov)
uelas0 = 1 ./ (dppriv ./ uepriv0 .+ dpgov .+ dpsave)

u0  = initVar1(false, reg_ndx)  
up0 = initVar1(false, reg_ndx)  
ug0 = initVar1(false, reg_ndx)  
us0 = initVar1(false, reg_ndx)
#au  = u0 .* ((up0 .^ (-dppriv)) .* (ug0 .^ (-dpgov)) .* (us0 .^ (-dpsave)))
au = ((up0 .^ dppriv) .* (ug0 .^ dpgov) .* (us0 .^ dpsave))
au = u0 ./ au
aus = us0 .* pop0 ./ save0

#   tmarg0 will be defined with respect to the base volume, not the vfob value
txs0   = @. ifelse(vxsb .!= 0, (vfob .- vxsb) ./ vxsb, 0)
tmarg0 = @. ifelse(vxsb .!= 0, (vcif .- vfob) ./ vxsb, 0)
tms0   = @. ifelse(vcif .!= 0, (vmsb .- vcif) ./ vcif, 0)

#   Set some tolerance level for tmarg0
tmarg0 = @. ifelse(abs.(tmarg0) > 1e-6, tmarg0, 0)

pfob0 = initVar3(false, comm_ndx, reg_ndx, reg_ndx)
for i in commodities
    for s in regions
        for d in regions
            pfob0[i,s,d] = (1 + txs0[i,s,d]) * pds0[i,s]
        end
    end
end
pcif0 = pfob0 .+ tmarg0
pmds0 = pcif0 .* (1 .+ tms0)

pms0 = initVar2(false, comm_ndx, reg_ndx)
qms0 = dropdims(sum(qfm0, dims=2), dims=2) .+ qpm0 .+ qgm0 .+ qim0
qxs0 = vfob ./ pfob0

ptrans0 = initVar3(false, comm_ndx, reg_ndx, reg_ndx)
qst0    = vst ./ pds0
qds0    = dropdims(sum(qfd0, dims=2), dims=2) .+ qpd0 .+ qgd0 .+ qid0

pt0     = initVar1(false, comm_ndx)
qtmfsd0 = vtwr ./ pt0
qtm0    = initVar1(true, comm_ndx)
for m in commodities
    for i in commodities
        for s in regions
            for d in regions
                qtm0[m] = qtm0[m] + qtmfsd0[m,i,s,d]
            end
        end
    end
end

#   Flags
qfaFlag = deepcopy(qfa0)
qfeFlag = deepcopy(qfe0)
qcaFlag = deepcopy(qca0)
qgaFlag = deepcopy(qga0)
qiaFlag = deepcopy(qia0)
qtmfsdFlag = deepcopy(qtmfsd0)
qxsFlag    = deepcopy(qxs0)
qtmFlag    = deepcopy(qtm0)
qstFlag    = deepcopy(qst0)
# qfeFlag .= (qfeFlag > 0.0) ? 1 : 0
for i in eachindex(qfeFlag)
    qfeFlag[i] = (qfeFlag[i] > 0.0) ? 1 : 0
end
for i in eachindex(qfaFlag)
    qfaFlag[i] = (qfaFlag[i] > 0.0) ? 1 : 0
end
for i in eachindex(qcaFlag)
    qcaFlag[i] = (qcaFlag[i] > 0.0) ? 1 : 0
end
for i in eachindex(qgaFlag)
    qgaFlag[i] = (qgaFlag[i] > 0.0) ? 1 : 0
end
for i in eachindex(qiaFlag)
    qiaFlag[i] = (qiaFlag[i] > 0.0) ? 1 : 0
end

setFlag(qtmfsdFlag)
setFlag(qxsFlag)
setFlag(qtmFlag)
setFlag(qstFlag)

#   Capital account closure

rental0 = initVar1(true, reg_ndx)
for r in regions
    rental0[r] = sum(pes0[cap,a,r]*qfe0[cap,a,r] for a in activities) / kb0[r]
end
rorc0 = rental0 ./ pinv0 .- depr
rore0 = rorc0 .* (kb0 ./ ke0) .^ RORFLEX
rorg0 = sum(rore0[r]*pinv0[r]*(qinv0[r] - depr[r]*kb0[r]) for r in regions) / sum(pinv0[r]*(qinv0[r] - depr[r]*kb0[r]) for r in regions)
risk  = rorg0 ./ rore0
ifRes = initVar1(true, reg_ndx)
ifRes[resReg] = 1
pinvLag = initVar1(false, reg_ndx)
psaveLag = initVar1(false, reg_ndx)
invwgt = initVar1(true, reg_ndx)
savwgt = initVar1(true, reg_ndx)

for r in regions
    isum = sum(pinv0[d]*(qinv0[d] - depr[d]*kb0[d]) for d in regions)
    invwgt[r] = pinv0[r]*(qinv0[r] - depr[r]*kb0[r]) / isum
    ssum = sum(save[d] for d in regions)
    savwgt[r] = save[r] / ssum
end
pfact0  = initVar1(false, reg_ndx)
pfactw0 = 1
Walras  = 0

#   Calibrate parameters

#   Production

#   !!!! Need to control for division by zero

shr_qint = (pint0 .* qint0) ./ (po0 .* qo0)
shr_qva  = (pva0 .* qva0) ./ (po0 .* qo0)
shr_qfa  = (pfa0 .* qfa0) ./ sum(pfa0 .* qfa0, dims=1)
shr_qfe  = (pfe0 .* qfe0) ./ sum(pfe0 .* qfe0, dims=1)
shr_qcas = (ps0 .* qca0) ./ sum(ps0 .* qca0, dims=1)
shr_qcad = (pca0 .* qca0) ./ sum(pca0 .* qca0, dims=2)
shr_gov  = (pga0 .* qga0) ./ sum(pga0 .* qga0, dims=1)
shr_inv  = (pia0 .* qia0) ./ sum(pia0 .* qia0, dims=1)
shr_qfd  = (pfd0 .* qfd0) ./ (pfd0 .* qfd0 .+ pfm0 .* qfm0)
shr_qfm  = (pfm0 .* qfm0) ./ (pfd0 .* qfd0 .+ pfm0 .* qfm0)
shr_qpd  = (ppd0 .* qpd0) ./ (ppd0 .* qpd0 .+ ppm0 .* qpm0)
shr_qpm  = (ppm0 .* qpm0) ./ (ppd0 .* qpd0 .+ ppm0 .* qpm0)
shr_qgd  = (pgd0 .* qgd0) ./ (pgd0 .* qgd0 .+ pgm0 .* qgm0)
shr_qgm  = (pgm0 .* qgm0) ./ (pgd0 .* qgd0 .+ pgm0 .* qgm0)
shr_qid  = (pid0 .* qid0) ./ (pid0 .* qid0 .+ pim0 .* qim0)
shr_qim  = (pim0 .* qim0) ./ (pid0 .* qid0 .+ pim0 .* qim0)
shr_qxs  = (pmds0 .* qxs0) ./ sum(pmds0 .* qxs0, dims=2)
shr_qfes = (pes0 .* qfe0) ./ sum(pes0 .* qfe0, dims=2)
shr_vtwr = deepcopy(qtmfsd0)
transFlag = initVar3(true, comm_ndx, reg_ndx, reg_ndx)

for i in commodities
    for s in regions
        for d in regions
            wsum = 0
            for m in commodities
                wsum = wsum + vtwr[m,i,s,d]
            end
            if wsum > 0
                #   There is transport for this node
                for m in commodities
                    shr_vtwr[m,i,s,d] = vtwr[m,i,s,d] / wsum 
                end
                transFlag[i,s,d] = 1
            else
                #   There is no transport for this node
                for m in commodities
                    shr_vtwr[m,i,s,d] = 0 
                end
                transFlag[i,s,d] = 0
            end
        end
    end
end

shr_qst  = vst ./ sum(vst, dims=2)
#for i in eachindex(shr_vtwr)
#    shr_vtwr[i] = (isnan(shr_vtwr[i]) ? 0 : shr_vtwr[i])
#end       
for i in eachindex(shr_qst)
    shr_qst[i] = (isnan(shr_qst[i]) ? 0 : shr_qst[i])
end       

#   Private consumption

alphac = initVar2(true, comm_ndx, reg_ndx)
for r in regions
    for i in commodities
        alphac[i,r] = (conshr0[i,r] / SUBPAR[i,r]) * (((yp0[r] / pop0[r]) / ppa0[i,r]) ^ SUBPAR[i,r]) * (up0[r]^(-INCPAR[i,r] * SUBPAR[i,r])) / sum(conshr0[j,r] / SUBPAR[j,r] for j in commodities)    
    end
end
# Can we replace the dual loops above and below with a more succinct representation?
# alphac = (conshr0 ./ SUBPAR) .* (((yp0 ./ pop0) ./ ppa0) .^ SUBPAR) .* (up0 .^ (-INCPAR .* SUBPAR)) ./ sum(conshr0 ./ SUBPAR, dims=1)
zshr0  = initVar2(true, comm_ndx, reg_ndx)
for r in regions
    for i in commodities
        zshr0[i,r] = alphac[i,r] * SUBPAR[i,r] * (up0[r] ^ (INCPAR[i,r] * SUBPAR[i,r])) * ((ppa0[i,r] / (yp0[r] / pop0[r])) ^ SUBPAR[i,r])    
    end
end
# zshr0  = alphac .* SUBPAR .* (up0 .^ (INCPAR .* SUBPAR)) .* ((ppa0 ./ (yp0 ./ pop0)) .^ SUBPAR)

#   Government utility
aug = pop0 ./ xg0

function solve_GTAP()
    GTAPMod = Model(PATHSolver.Optimizer)
   
    #   Normally exogenous variables
    ptax  = deepcopy(ptax0)
    tfe   = deepcopy(tfe0)
    tinc  = deepcopy(tinc0)
    tfd   = deepcopy(tfd0)
    tfm   = deepcopy(tfm0)
    tpd   = deepcopy(tpd0)
    tpm   = deepcopy(tpm0)
    tgd   = deepcopy(tgd0)
    tgm   = deepcopy(tgm0)
    tid   = deepcopy(tid0)
    tim   = deepcopy(tim0)
    txs   = deepcopy(txs0)
    tmarg = deepcopy(tmarg0)
    tms   = deepcopy(tms0)
    qe    = ones(length(factors), length(regions))
    pop   = deepcopy(pop0)
    pnum  = 1.2
    # ptax  = 0.5 .* ptax0
    # tmarg = 0.05 .* tmarg0
        
    @variables GTAPMod begin
        qint[a = activities, r = regions], (start = 1)                              # Aggregate intermediate demand
        qva[a = activities, r =  regions], (start = 1)                              # Demand for aggregate value added bundle
        po[a = activities, r = regions], (start = 1)                                # Unit cost of production
        qfa[i = commodities, a = activities, r = regions], (start = 1)              # (Armington) demand for inputs
        pint[a = activities, r = regions], (start = 1)                              # Aggregate cost of inputs
        qfe[f = factors, a = activities, r = regions], (start = 1)                  # Factor demand
        pva[a = activities, r = regions], (start = 1)                               # Aggregate price of value added
        qca[i = commodities, a = activities, r = regions], (start = qcaFlag[i,a,r]) # Commodity i supplied by a
        qo[a = activities, r = regions], (start = 1)                                # Output by activity a
        pca[i = commodities, a = activities, r = regions], (start = 1)              # Tax-inclusive price of commodity i supplied by a
        pds[i = commodities, r = regions], (start = 1)                              # Price of supply of commodity i
        ps[i = commodities, a = activities, r = regions], (start = 1)               # Tax-exclusive price of commodity i supplied by a
        taxrout[r = regions], (start = 1)                                           # Value of output tax revenues
        taxrfu[r = regions], (start = 1)                                            # Value of factor-use tax revenues
        taxriu[r = regions], (start = 1)                                            # Value of tax revenues on intermediate demand
        taxrpc[r = regions], (start = 1)                                            # Value of tax revenues on private demand
        taxrgc[r = regions], (start = 1)                                            # Value of tax revenues on government demand
        taxric[r = regions], (start = 1)                                            # Value of tax revenues on investment demand
        taxrexp[r = regions], (start = 1)                                           # Value of tax revenues on exports
        taxrimp[r = regions], (start = 1)                                           # Value of tax revenues on imports
        indtax[r = regions], (start = 1)                                            # Total value of indirect tax revenues
        fincome[r = regions], (start = 1)                                           # Factor income at basic prices net of depreciation
        y[r = regions], (start = 1)                                                 # Regional income
        uepriv[r = regions], (start = uepriv0[r])                                   # Elasticity of cost wrt utility from private consumption
        uelas[r = regions], (start = uelas0[r])                                     # Elasticity of cost of utility wrt utility
        yp[r = regions], (start = 1)                                                # Private consumption expenditure
        yg[r = regions], (start = 1)                                                # Government consumption expenditure
        save[r = regions], (start = 1)                                              # Regional supply of nominal savings
        u[r = regions], (start = u0[r])                                             # Total utility (per capita)
        us[r = regions], (start = us0[r])                                           # Utility derived from savings
        zshr[i = commodities, r = regions], (start = 1)                             # Auxiliary share variable for private consumption
        conshr[i = commodities, r = regions], (start = 1)                           # Private budget shares
        qpa[i = commodities, r = regions], (start = 1)                              # Private consumption (at the Armington level)
        up[r = regions], (start = 1)                                                # Private utility
        qga[i = commodities, r = regions], (start = 1)                              # Government consumption (at the Armington level)
        pgov[r = regions], (start = 1)                                              # Government expenditure price deflator
        xg[r = regions], (start = 1)                                                # Aggregate volume of government expenditures
        ug[r = regions], (start = 1)                                                # Per capita utility derived from public expenditures
        qia[i = commodities, r = regions], (start = 1)                              # Investment expenditure (at the Armington level)
        pinv[r = regions], (start = 1)                                              # Investment expenditure price deflator
        qinv[r = regions], (start = 1)                                              # Aggregate volume of investment expenditures
        pfd[i = commodities, a = activities, r = regions], (start = 1)              # End user price for domestic intermediate goods
        pfm[i = commodities, a = activities, r = regions], (start = 1)              # End user price for imported intermediate goods
        ppd[i = commodities, r = regions], (start = 1)                              # End user price for domestic private goods
        ppm[i = commodities, r = regions], (start = 1)                              # End user price for imported private goods
        pgd[i = commodities, r = regions], (start = 1)                              # End user price for domestic government goods
        pgm[i = commodities, r = regions], (start = 1)                              # End user price for imported government goods
        pid[i = commodities, r = regions], (start = 1)                              # End user price for domestic investment goods
        pim[i = commodities, r = regions], (start = 1)                              # End user price for imported investment goods
        pfa[i = commodities, a = activities, r = regions], (start = 1)              # Armington price for intermediate goods
        ppa[i = commodities, r = regions], (start = 1)                              # Armington price for private goods
        pga[i = commodities, r = regions], (start = 1)                              # Armington price for government goods
        pia[i = commodities, r = regions], (start = 1)                              # Armington price for investment goods
        qfd[i = commodities, a = activities, r = regions], (start = 1)              # Demand for domestic intermediate goods
        qfm[i = commodities, a = activities, r = regions], (start = 1)              # Demand for imported intermediate goods
        qpd[i = commodities, r = regions], (start = 1)                              # Demand for domestic private goods
        qpm[i = commodities, r = regions], (start = 1)                              # Demand for imported private goods
        qgd[i = commodities, r = regions], (start = 1)                              # Demand for domestic government goods
        qgm[i = commodities, r = regions], (start = 1)                              # Demand for imported government goods
        qid[i = commodities, r = regions], (start = 1)                              # Demand for domestic investment goods
        qim[i = commodities, r = regions], (start = 1)                              # Demand for imported investment goods
        qms[i = commodities, r = regions], (start = 1)                              # Aggregate import demand
        qxs[i = commodities, s = regions, d = regions], (start = 1)                 # Export/import from region s towards region d
        pms[i = commodities, d = regions], (start = 1)                              # Aggregate import price
        pfob[i = commodities, s = regions, d = regions], (start = 1)                # FOB export price
        pcif[i = commodities, s = regions, d = regions], (start = 1)                # CIF import price
        pmds[i = commodities, s = regions, d = regions], (start = 1)                # Tariff inclusive bilateral import price
        qds[i = commodities, r = regions], (start = 1)                              # Aggregate domestic demand for domestic goods (/x qst)
        qc[i = commodities, r = regions], (start = 1)                               # Market equilibrium for domestically produced goods
        pes[f = factors, a = activities, r = regions], (start = 1)                  # Factor market allocation and equilibrium
        pe[f = factors, r = regions], (start = 1)                                   # Aggregate factor price
        pfe[f = factors, a = activities, r = regions], (start = 1)                  # Factor returns at purchasers prices
        peb[f = factors, a = activities, r = regions], (start = 1)                  # Factor returns at basic prices
        qtmfsd[m = commodities, i = commodities, s = regions, d = regions], (start = qtmfsdFlag[m,i,s,d]) # Demand for service m to transport i from s to d
        ptrans[i = commodities, s = regions, d = regions], (start = 1)              # Average price to transport i from s to d
        qtm[m = commodities], (start = qtmFlag[m])                                  # Total demand for service m
        qst[m = commodities, r = regions], (start = qstFlag[m,r])                   # Regional supply of service m
        pt[m = commodities], (start = 1)                                            # Aggregate world price of service m
        kb[r = regions], (start = 1)                                                # Non-normalized aggregate capital stock
        ke[r = regions], (start = 1)                                                # End-of-period capital stock
        rental[r = regions], (start = 1)                                            # After-tax gross rate of return to capital
        rorc[r = regions], (start = 1)                                              # Net rate of return to capital
        rore[r = regions], (start = 1)                                              # Expected rate of return to capital
        savf[r = regions], (start = savf0[r])                                       # Default capital account closure--determination of capital flows
        # rorg, (start = 1)                                                           # Default capital account closure--determination of the global rate of return 
        yi[r = regions], (start = 1)                                                # Nominal gross investment
        chisave, (start = 1)                                                        # Savings adjustment factor
        psave[r = regions], (start = 1)                                             # Price index for regional savings
        pfact[r = regions], (start = 1)                                             # Average regional factor price
        pfactw, (start = 1)                                                         # Average global factor price
        Walras, (start = 0)                                                         # Check on Walras' Law
    end
    #   Production module --------------------------------------------------------------------------

    #   Intermediate demand bundle
    @constraint(GTAPMod, eq_qint[a in activities, r in regions], 
        qint[a,r] - qo[a,r]*(po[a,r]/pint[a,r])^ESUBT[a,r] ⟂ qint[a,r])
    
    #   Value added bundle
    @constraint(GTAPMod, eq_qva[a in activities, r in regions], 
        qva[a,r] - qo[a,r]*(po[a,r]/pva[a,r])^ESUBT[a,r] ⟂ qva[a,r])
    
    #   Unit cost of production
    @constraint(GTAPMod, eq_po[a in activities, r in regions], 
        po[a,r]^(1 - ESUBT[a,r]) - (shr_qint[a,r]*pint[a,r]^(1 - ESUBT[a,r]) + shr_qva[a,r]*pva[a,r]^(1 - ESUBT[a,r])) ⟂ po[a,r])
        # po[a,r]*qo[a,r] - (shr_qint[a,r]*pint[a,r]*qint[a,r] + shr_qva[a,r]*pva[a,r]*qva[a,r]) ⟂ po[a,r])
        # po[a,r]*qo[a,r]*po0[a,r]*qo0[a,r] - (sum(pfa[i,a,r]*qfa[i,a,r]*pfa0[i,a,r]*qfa0[i,a,r] for i in commodities) 
        #    + sum(pfe[f,a,r]*qfe[f,a,r]*pfe0[f,a,r]*qfe0[f,a,r] for f in factors)) ⟂ po[a,r])

    #   Armington input demand
    @constraint(GTAPMod, eq_qfa[i in commodities, a in activities, r in regions], 
        qfaFlag[i,a,r]*(qfa[i,a,r] - qint[a,r]*(pint[a,r]/pfa[i,a,r])^ESUBC[a,r]) + (1 - qfaFlag[i,a,r])*(qfa[i,a,r]) ⟂ qfa[i,a,r])

     #   Aggregate price of inputs
     @constraint(GTAPMod, eq_pint[a in activities, r in regions], 
        pint[a,r]^(1 - ESUBC[a,r]) - (sum(shr_qfa[i,a,r]*pfa[i,a,r]^(1 - ESUBC[a,r]) for i in commodities)) ⟂ pint[a,r])
        # pint[a,r]*qint[a,r] - (sum(shr_qfa[i,a,r]*pfa[i,a,r]*qfa[i,a,r] for i in commodities)) ⟂ pint[a,r])
    
    #   Factor demand
    @constraint(GTAPMod, eq_qfe[f in factors, a in activities, r in regions], 
        qfeFlag[f,a,r]*(qfe[f,a,r] - qva[a,r]*(pva[a,r]/pfe[f,a,r])^ESUBVA[a,r]) + (1 - qfeFlag[f,a,r])*(qfe[f,a,r]) ⟂ qfe[f,a,r])
    
     #   Aggregate price of value added
    @constraint(GTAPMod, eq_pva[a in activities, r in regions], 
        pva[a,r]^(1 - ESUBVA[a,r]) - (sum(shr_qfe[f,a,r]*pfe[f,a,r]^(1 - ESUBVA[a,r]) for f in factors)) ⟂ pva[a,r])
        # pva[a,r]*qva[a,r] - (sum(shr_qfe[f,a,r]*pfe[f,a,r]*qfe[f,a,r] for f in factors)) ⟂ pva[a,r])
    
    #   Make module --------------------------------------------------------------------------------

    #   Supply of commodity i by activity a
    @constraint(GTAPMod, eq_qca[i in commodities, a in activities, r in regions], 
        (qcaFlag[i,a,r])*(qca[i,a,r] - qo[a,r]*(ps[i,a,r]/po[a,r])^ETRAQ[a,r]) + (1-qcaFlag[i,a,r])*(qca[i,a,r]) ⟂ qca[i,a,r])
    
    #   Defines output by activity a
    @constraint(GTAPMod, eq_qo[a in activities, r in regions], 
        po[a,r]^(1 + ETRAQ[a,r]) - (sum(shr_qcas[i, a, r]*ps[i,a,r]^(1 + ETRAQ[a,r]) for i in commodities)) ⟂ qo[a,r])
        # po[a,r]*qo[a,r] - (sum(shr_qcas[i,a,r]*ps[i,a,r]*qca[i,a,r] for i in commodities)) ⟂ qo[a,r])

    #  Demand for commodity i produced by activity a (inverse CES)
    @constraint(GTAPMod, eq_pca[i in commodities, a in activities, r in regions], 
        (qcaFlag[i,a,r])*(pca[i,a,r] - pds[i,r]*(qc[i,r]/qca[i,a,r])^ESUBQ[i,r]) + (1-qcaFlag[i,a,r])*(pca[i,a,r] - 1) ⟂ pca[i,a,r])
    
    @constraint(GTAPMod, eq_ps[i in commodities, a in activities, r in regions], 
        (qcaFlag[i,a,r])*(pca[i,a,r] - ps[i,a,r]*(1 + ptax[i,a,r])*(ps0[i,a,r]/pca0[i,a,r])) + (1-qcaFlag[i,a,r])*(ps[i,a,r] - 1) ⟂ ps[i,a,r])

    #   Defines supply price of commodity i
    @constraint(GTAPMod, eq_pds[i in commodities, r in regions], 
        pds[i,r]*qc[i,r] - (sum(shr_qcad[i,a,r]*pca[i,a,r]*qca[i,a,r] for a in activities)) ⟂ pds[i,r])

    #   Income and its distribution ----------------------------------------------------------------
 
    #   Tax revenues from intermediate demand
    @constraint(GTAPMod, eq_taxriu[r in regions],
        taxriu[r]*taxriu0[r] - (sum(tfd[i,a,r]*pds[i,r]*qfd[i,a,r]*pds0[i,r]*qfd0[i,a,r]
                             +      tfm[i,a,r]*pms[i,r]*qfm[i,a,r]*pms0[i,r]*qfm0[i,a,r] for i in commodities, a in activities)) ⟂ taxriu[r])

    #   Tax revenues from private demand
    @constraint(GTAPMod, eq_taxrpc[r in regions],
        taxrpc[r]*taxrpc0[r] - (sum(tpd[i,r]*pds[i,r]*qpd[i,r]*pds0[i,r]*qpd0[i,r] 
                             +      tpm[i,r]*pms[i,r]*qpm[i,r]*pms0[i,r]*qpm0[i,r] for i in commodities)) ⟂ taxrpc[r])

    #   Tax revenues from government demand
    @constraint(GTAPMod, eq_taxrgc[r in regions],
        taxrgc[r]*taxrgc0[r] - (sum(tgd[i,r]*pds[i,r]*qgd[i,r]*pds0[i,r]*qgd0[i,r] 
                             +      tgm[i,r]*pms[i,r]*qgm[i,r]*pms0[i,r]*qgm0[i,r] for i in commodities)) ⟂ taxrgc[r])

    #   Tax revenues from investment demand
    @constraint(GTAPMod, eq_taxric[r in regions],
        taxric[r]*taxric0[r] - (sum(tid[i,r]*pds[i,r]*qid[i,r]*pds0[i,r]*qid0[i,r] 
                             +      tim[i,r]*pms[i,r]*qim[i,r]*pms0[i,r]*qim0[i,r] for i in commodities)) ⟂ taxric[r])

    #   Output tax revenues
    @constraint(GTAPMod, eq_taxrout[r in regions],
       taxrout[r] - (sum((ps0[i,a,r]*qca0[i,a,r]/taxrout0[r])*ptax[i,a,r]*ps[i,a,r]*qca[i,a,r]*qcaFlag[i,a,r] for i in commodities, a in activities)) ⟂ taxrout[r])

    #   Factor use tax revenues
    @constraint(GTAPMod, eq_taxrfu[r in regions], 
       taxrfu[r]*taxrfu0[r] - (sum(tfe[f,a,r]*peb[f,a,r]*qfe[f,a,r]*peb0[f,a,r]*qfe0[f,a,r] for f in factors, a in activities)) ⟂ taxrfu[r])

    #   Export tax revenues
    @constraint(GTAPMod, eq_taxrexp[r in regions], 
        taxrexp[r]*taxrexp0[r] - (sum(txs[i,r,d]*pds[i,r]*qxs[i,r,d]*pds0[i,r]*qxs0[i,r,d] for i in commodities, d in regions)) ⟂ taxrexp[r])

    #   Import tax revenues
    @constraint(GTAPMod, eq_taxrimp[r in regions], 
        taxrimp[r]*taxrimp0[r] - (sum(tms[i,s,r]*pcif[i,s,r]*qxs[i,s,r]*pcif0[i,s,r]*qxs0[i,s,r] for i in commodities, s in regions)) ⟂ taxrimp[r])

    #   Total indirect taxes
    @constraint(GTAPMod, eq_indtax[r in regions], 
        indtax[r]*indtax0[r] - (taxrpc0[r]*taxrpc[r] + taxrgc0[r]*taxrgc[r] + taxric0[r]*taxric[r] + taxriu0[r]*taxriu[r] 
            + taxrfu0[r]*taxrfu[r] + taxrout0[r]*taxrout[r] + taxrexp0[r]*taxrexp[r] + taxrimp0[r]*taxrimp[r]) ⟂ indtax[r])

    @constraint(GTAPMod, eq_fincome[r in regions], 
       fincome[r]*fincome0[r] - (sum(peb[f,a,r]*qfe[f,a,r]*peb0[f,a,r]*qfe0[f,a,r] for f in factors, a in activities) - depr[r]*pinv[r]*pinv0[r]*kb[r]*kb0[r]) ⟂ fincome[r])  

    @constraint(GTAPMod, eq_y[r in regions], 
       y[r]*y0[r] - (fincome[r]*fincome0[r] + indtax[r]*indtax0[r]) ⟂ y[r])

    #   Disposition of income and final demand -----------------------------------------------------
 
    #   Elasticity of expenditure wrt utility from private consumption
    @constraint(GTAPMod, eq_uepriv[r in regions],
        uepriv[r] - sum(conshr[i,r]*conshr0[i,r]*INCPAR[i,r] for i in commodities) ⟂ uepriv[r])

    #   Elasticity of total expenditure wrt to utility
    @constraint(GTAPMod, eq_uelas[r in regions],
        uelas[r]*(dppriv[r]/uepriv[r] + dpgov[r] + dpsave[r]) - 1 ⟂ uelas[r])

    #   Aggregate nominal private consumption
    @constraint(GTAPMod, eq_yp[r in regions],
        yp[r]*yp0[r] - dppriv[r]*(uelas[r]/uepriv[r])*y[r]*y0[r] ⟂ yp[r])

    #   Aggregate nominal government consumption
    @constraint(GTAPMod, eq_yg[r in regions],
        yg[r]*yg0[r] - dpgov[r]*uelas[r]*y[r]*y0[r] ⟂ yg[r])

    #   Domestic supply of savings
    @constraint(GTAPMod, eq_qsave[r in regions],
        save[r]*save0[r] - dpsave[r]*uelas[r]*y[r]*y0[r] ⟂ save[r])

    #   Utility function
    @constraint(GTAPMod, eq_u[r in regions],
        log(u[r]) - (log(au[r]) +  dppriv[r]*log(up[r]) + dpgov[r]*log(ug[r]) + dpsave[r]*log(us[r])) ⟂ u[r])
    
    #   Utility from savings
    @constraint(GTAPMod, eq_us[r in regions],
        us[r] - (save0[r]*aus[r]/pop[r])*(save[r]/psave[r]) ⟂ us[r])

    #   Auxiliary consumption variable
    @constraint(GTAPMod, eq_zshr[i in commodities, r in regions],
        zshr[i,r]*zshr0[i,r] - (alphac[i,r]*SUBPAR[i,r]*((pop[r]*ppa0[i,r]/yp0[r])^SUBPAR[i,r])*(up[r]^(INCPAR[i,r]*SUBPAR[i,r]))*((ppa[i,r]/yp[r])^SUBPAR[i,r])) ⟂ zshr[i,r])
        
    #   Private consumption budget shares
    @constraint(GTAPMod, eq_conshr[i in commodities, r in regions],   
        conshr[i,r]*conshr0[i,r]*(sum(zshr[j,r]*zshr0[j,r] for j in commodities)) - (zshr[i,r]*zshr0[i,r]) ⟂ conshr[i,r])
    
    #   Household demand for goods and services
    @constraint(GTAPMod, eq_qpa[i in commodities, r in regions],
        ppa[i,r]*qpa[i,r]*ppa0[i,r]*qpa0[i,r] - (conshr0[i,r]*conshr[i,r]*yp0[r]*yp[r]) ⟂ qpa[i,r])
    
    #   Private utility
    @constraint(GTAPMod, eq_up[r in regions],
        1 - sum(zshr0[i,r]*zshr[i,r]/SUBPAR[i,r] for i in commodities) ⟂ up[r])

    #   Government demand for goods and services
    @constraint(GTAPMod, eq_qga[i in commodities, r in regions],   
        (qgaFlag[i,r]*(qga[i,r] - (xg[r] * (pgov[r] / pga[i,r])^ESUBG[r]))) + ((1-qgaFlag[i,r])*qga[i,r]) ⟂ qga[i,r])
        #qga[i,r] - (xg[r] * (pgov[r] / pga[i,r])^ESUBG[r]) ⟂ qga[i,r])
    
    #   Government expenditure price deflator
    @constraint(GTAPMod, eq_pgov[r in regions],
        pgov[r]*xg[r] - (sum(shr_gov[i,r]*pga[i,r]*qga[i,r] for i in commodities)) ⟂ pgov[r])
    
    #   Aggregate volume of government expenditures
    @constraint(GTAPMod, eq_xg[r in regions],
        pgov[r]*xg[r]*pgov0[r]*xg0[r] - yg[r]*yg0[r] ⟂ xg[r])

    @constraint(GTAPMod, eq_ug[r in regions],
        ug[r] - (aug[r]*xg0[r]/pop[r])*xg[r] ⟂ ug[r])

    #   Investment demand for goods and services
    @constraint(GTAPMod, eq_qia[i in commodities, r in regions],   
        (qiaFlag[i,r]*(qia[i,r] - (qinv[r] * (pinv[r] / pia[i,r])^ESUBI[r]))) + ((1-qiaFlag[i,r])*qia[i,r]) ⟂ qia[i,r])

    #   Investment expenditure price deflator
    @constraint(GTAPMod, eq_pinv[r in regions],
        pinv[r]*qinv[r] - (sum(shr_inv[i,r]*pia[i,r]*qia[i,r] for i in commodities)) ⟂ pinv[r])

    #   Aggregate volume of investment expenditures
    @constraint(GTAPMod, eq_qinv[r in regions],
        pinv[r]*qinv[r]*pinv0[r]*qinv0[r] - yi[r]*yi0[r] ⟂ qinv[r])

    #   Armington module ---------------------------------------------------------------------------

    #   End use prices for domestic and imported goods
    #   Assumes all users face same basic prices: pds for domestic goods and pms for imported goods
    
    @constraint(GTAPMod, eq_pfd[i in commodities, a in activities, r in regions],
        pfd[i,a,r]*pfd0[i,a,r] - (1 + tfd[i,a,r]) * pds[i,r]*pds0[i,r] ⟂ pfd[i,a,r])

    @constraint(GTAPMod, eq_pfm[i in commodities, a in activities, r in regions],
        pfm[i,a,r]*pfm0[i,a,r] - (1 + tfm[i,a,r]) * pms[i,r]*pms0[i,r] ⟂ pfm[i,a,r])

    @constraint(GTAPMod, eq_ppd[i in commodities, r in regions],
        ppd[i,r]*ppd0[i,r] - (1 + tpd[i,r]) * pds[i,r]*pds0[i,r] ⟂ ppd[i,r])

    @constraint(GTAPMod, eq_ppm[i in commodities, r in regions],
        ppm[i,r]*ppm0[i,r] - (1 + tpm[i,r]) * pms[i,r]*pms0[i,r] ⟂ ppm[i,r])

    @constraint(GTAPMod, eq_pgd[i in commodities, r in regions],
        pgd[i,r]*pgd0[i,r] - (1 + tgd[i,r]) * pds[i,r]*pds0[i,r] ⟂ pgd[i,r])

    @constraint(GTAPMod, eq_pgm[i in commodities, r in regions],
        pgm[i,r]*pgm0[i,r] - (1 + tgm[i,r]) * pms[i,r]*pms0[i,r] ⟂ pgm[i,r])

    @constraint(GTAPMod, eq_pid[i in commodities, r in regions],
        pid[i,r]*pid0[i,r] - (1 + tid[i,r]) * pds[i,r]*pds0[i,r] ⟂ pid[i,r])

    @constraint(GTAPMod, eq_pim[i in commodities, r in regions],
        pim[i,r]*pim0[i,r] - (1 + tim[i,r]) * pms[i,r]*pms0[i,r] ⟂ pim[i,r])

    #   Armington price

    @constraint(GTAPMod, eq_pfa[i in commodities, a in activities, r in regions],
        pfa[i,a,r]*qfa[i,a,r] - (shr_qfd[i,a,r]*pfd[i,a,r]*qfd[i,a,r] + shr_qfm[i,a,r]*pfm[i,a,r]*qfm[i,a,r]) ⟂ pfa[i,a,r])

    @constraint(GTAPMod, eq_ppa[i in commodities, r in regions],
        ppa[i,r]*qpa[i,r] - (shr_qpd[i,r]*ppd[i,r]*qpd[i,r] + shr_qpm[i,r]*ppm[i,r]*qpm[i,r]) ⟂ ppa[i,r])

    @constraint(GTAPMod, eq_pga[i in commodities, r in regions],
        pga[i,r]*qga[i,r] - (shr_qgd[i,r]*pgd[i,r]*qgd[i,r] + shr_qgm[i,r]*pgm[i,r]*qgm[i,r]) ⟂ pga[i,r])

    @constraint(GTAPMod, eq_pia[i in commodities, r in regions],
        pia[i,r]*qia[i,r] - (shr_qid[i,r]*pid[i,r]*qid[i,r] + shr_qim[i,r]*pim[i,r]*qim[i,r]) ⟂ pia[i,r])

    #   Top level sourcing
    @constraint(GTAPMod, eq_qfd[i in commodities, a in activities, r in regions],
        qfd[i,a,r] - (qfa[i,a,r]  * (pfa[i,a,r] / pfd[i,a,r])^ESUBD[i,r]) ⟂ qfd[i,a,r])

    @constraint(GTAPMod, eq_qfm[i in commodities, a in activities, r in regions],
        qfm[i,a,r] - (qfa[i,a,r]  * (pfa[i,a,r] / pfm[i,a,r])^ESUBD[i,r]) ⟂ qfm[i,a,r])

    @constraint(GTAPMod, eq_qpd[i in commodities, r in regions],
        qpd[i,r] - (qpa[i,r]  * (ppa[i,r] / ppd[i,r])^ESUBD[i,r]) ⟂ qpd[i,r])

    @constraint(GTAPMod, eq_qpm[i in commodities, r in regions],
        qpm[i,r] - (qpa[i,r]  * (ppa[i,r] / ppm[i,r])^ESUBD[i,r]) ⟂ qpm[i,r])

    @constraint(GTAPMod, eq_qgd[i in commodities, r in regions],
        qgd[i,r] - (qga[i,r]  * (pga[i,r] / pgd[i,r])^ESUBD[i,r]) ⟂ qgd[i,r])

    @constraint(GTAPMod, eq_qgm[i in commodities, r in regions],
        qgm[i,r] - (qga[i,r]  * (pga[i,r] / pgm[i,r])^ESUBD[i,r]) ⟂ qgm[i,r])

    @constraint(GTAPMod, eq_qid[i in commodities, r in regions],
        qid[i,r] - (qia[i,r]  * (pia[i,r] / pid[i,r])^ESUBD[i,r]) ⟂ qid[i,r])

    @constraint(GTAPMod, eq_qim[i in commodities, r in regions],
        qim[i,r] - (qia[i,r]  * (pia[i,r] / pim[i,r])^ESUBD[i,r]) ⟂ qim[i,r])

    #   Second level Armington -- sourcing by region
    #   Aggregate import demand
    @constraint(GTAPMod, eq_qms[i in commodities, r in regions],
        qms[i,r]*qms0[i,r] - (sum(qfm[i,a,r]*qfm0[i,a,r] for a in activities) + qpm[i,r]*qpm0[i,r] + qgm[i,r]*qgm0[i,r] + qim[i,r]*qim0[i,r]) ⟂ qms[i,r])        

    #   Demand for imports by region d sourced from region s
    @constraint(GTAPMod, eq_qxs[i in commodities, s in regions, d in regions],
        qxs[i,s,d] - (qms[i,d] * (pms[i,d] / pmds[i,s,d])^ESUBM[i,d]) ⟂ qxs[i,s,d])
 
    #   Aggregate price of imports
    @constraint(GTAPMod, eq_pms[i in commodities, d in regions],
        pms[i,d]^(1 - ESUBM[i,d]) - (sum(shr_qxs[i,s,d]*pmds[i,s,d]^(1 - ESUBM[i,d]) for s in regions)) ⟂ pms[i,d])      
        # pms[i,d]*qms[i,d] - (sum(shr_qxs[i,s,d]*pmds[i,s,d]*qxs[i,s,d] for s in regions)) ⟂ pms[i,d])      

    #   Bilateral prices ---------------------------------------------------------------------------

    #   Export price at FOB
    @constraint(GTAPMod, eq_pfob[i in commodities, s in regions, d in regions],
        pfob[i,s,d]*pfob0[i,s,d] - ((1 + txs[i,s,d])*pds[i,s]*pds0[i,s]) ⟂ pfob[i,s,d])      

    #   Import price at CIF
    @constraint(GTAPMod, eq_pcif[i in commodities, s in regions, d in regions],
        pcif[i,s,d]*pcif0[i,s,d] - (pfob[i,s,d]*pfob0[i,s,d] + ptrans0[i,s,d]*ptrans[i,s,d]*tmarg[i,s,d]) ⟂ pcif[i,s,d])      

    #   Import price tariff inclusive
    @constraint(GTAPMod, eq_pmds[i in commodities, s in regions, d in regions],
        pmds[i,s,d]*pmds0[i,s,d] - ((1 + tms[i,s,d])*pcif[i,s,d]*pcif0[i,s,d]) ⟂ pmds[i,s,d])      

    #   Global transport services ------------------------------------------------------------------

    @constraint(GTAPMod, eq_qtmfsd[m in commodities, i in commodities, s in regions, d in regions],
        qtmfsdFlag[m,i,s,d]*(qtmfsd[m,i,s,d]*qtmfsd0[m,i,s,d] - tmarg[i,s,d]*qxs[i,s,d]*qxs0[i,s,d]) + (1 - qtmfsdFlag[m,i,s,d])*qtmfsd[m,i,s,d] ⟂ qtmfsd[m,i,s,d])      

    @constraint(GTAPMod, eq_ptrans[i in commodities, s in regions, d in regions],
        transFlag[i,s,d] * (ptrans[i,s,d] - sum(shr_vtwr[m,i,s,d]*pt[m] for m in commodities)) + (1 - transFlag[i,s,d])*(ptrans[i,s,d] - 1) ⟂ ptrans[i,s,d])      

    @constraint(GTAPMod, eq_qtm[m in commodities],
        qtmFlag[m]*(qtm[m]*qtm0[m] - sum(qtmfsd[m,i,s,d]*qtmfsd0[m,i,s,d] for i in commodities, s in regions, d in regions)) + (1 - qtmFlag[m])*qtm[m] ⟂ qtm[m])      
    
    @constraint(GTAPMod, eq_qst[m in commodities, r in regions],
        qstFlag[m,r]*(qst[m,r] - (qtm[m] * (pt[m] / pds[m,r])^ESUBS[m])) + (1 - qstFlag[m,r])*qst[m,r] ⟂ qst[m,r])      

    @constraint(GTAPMod, eq_pt[m in commodities],
        qtmFlag[m]*(pt[m]*qtm[m] - (sum(shr_qst[m,r]*pds[m,r]*qst[m,r] for r in regions))) + (1 - qtmFlag[m])*(pt[m] - 1) ⟂ pt[m])      

    #   Domestic goods market equilibrium ----------------------------------------------------------

    #   Total domestic demand for domestic goods
    @constraint(GTAPMod, eq_qds[i in commodities, r in regions],
        qds[i,r]*qds0[i,r] - (sum(qfd[i,a,r]*qfd0[i,a,r] for a in activities) + qpd[i,r]*qpd0[i,r] + qgd[i,r]*qgd0[i,r] + qid[i,r]*qid0[i,r]) ⟂ qds[i,r])        

    #   Market clearing for domestically produced goods
    @constraint(GTAPMod, eq_qc[i in commodities, r in regions],
        qc[i,r]*qc0[i,r] - (qds[i,r]*qds0[i,r] + sum(qxs[i,r,d]*qxs0[i,r,d] for d in regions) + qst[i,r]*qst0[i,r]) ⟂ qc[i,r])        

    #   Factor market equilibrium ------------------------------------------------------------------

    #   Factor supply and equilibrium
    @constraint(GTAPMod, eq_pes[f in factors, a in activities, r in regions], 
        (qfa[f,a,r] - (qe[f,r] * (pes[f,a,r]/pe[f,r])^ETRAE[f,r]))*(1 - ifPerf[f,r]) 
        + (pes[f,a,r]*pes0[f,a,r] - pe[f,r]*pe0[f,r])*ifPerf[f,r] ⟂ pes[f,a,r]) 
 
    @constraint(GTAPMod, eq_pe[f in factors, r in regions], 
        pe[f,r]*qe[f,r] - (sum(shr_qfes[f,a,r]*pes[f,a,r]*qfe[f,a,r] for a in activities)) ⟂ pe[f,r]) 
 
    #   Factor returns at purchasers prices
    @constraint(GTAPMod, eq_pfe[f in factors, a in activities, r in regions], 
        pfe[f,a,r]*pfe0[f,a,r] - (1 + tfe[f,a,r])*peb[f,a,r]*peb0[f,a,r] ⟂ pfe[f,a,r]) 

    #   Factor returns at basic prices
    @constraint(GTAPMod, eq_peb[f in factors, a in activities, r in regions], 
        pes[f,a,r]*pes0[f,a,r] - (1 - tinc[f,a,r])*peb[f,a,r]*peb0[f,a,r] ⟂ peb[f,a,r]) 

    #   Capital account closure --------------------------------------------------------------------

    #   Non-normalized capital stock--change is equal to change in normalized capital stock
    #   !!!! In the GEMPACK code, isn't VES(ENDWC,r) / GROSSCAP(r) always equal to 1?
    @constraint(GTAPMod, eq_kb[r in regions],
        kb[r] - qe[cap,r] ⟂ kb[r])

    #   End-of-period capital stock
    @constraint(GTAPMod, eq_ke[r in regions],
        ke[r]*ke0[r] - ((1 - depr[r])*kb[r]*kb0[r] + qinv[r]*qinv0[r]) ⟂ ke[r])

    #   After-tax gross rate of return to capital
    @constraint(GTAPMod, eq_rental[r in regions],
        rental[r]*kb[r]*rental0[r]*kb0[r] - (sum(pes[cap,a,r]*qfe[cap,a,r]*pes0[cap,a,r]*qfe0[cap,a,r] for a in activities)) ⟂ rental[r])

    #   Net rate of return to capital
    @constraint(GTAPMod, eq_rorc[r in regions],
        rorc[r]*rorc0[r] - ((rental0[r]/pinv0[r])*rental[r]/pinv[r] - depr[r]) ⟂ rorc[r])

    #   Expected rate of return to capital
    @constraint(GTAPMod, eq_rore[r in regions],
        rore[r] - (((rorc0[r]/rore0[r])*(kb0[r]/ke0[r])^RORFLEX[r])*rorc[r]*(kb[r]/ke[r])^RORFLEX[r]) ⟂ rore[r])

    #   Default capital account closure--determination of capital flows
    #@constraint(GTAPMod, eq_savf[r in regions],
    #    risk[r]*rore[r]*rore0[r] - rorg*rorg0 ⟂ savf[r])

    #   Default capital account closure--determination of the global rate of return
    #@constraint(GTAPMod, eq_rorg,
    #    sum(savf[r] for r in regions) ⟂ rorg)

    @constraint(GTAPMod, eq_savf[r in regions],
        savf[r] - pnum*savf0[r] ⟂ savf[r])

    #   Savings = investment
    @constraint(GTAPMod, eq_yi[r in regions],
        yi[r]*yi0[r] - (save[r]*save0[r] + savf[r] + depr[r]*pinv[r]*kb[r]*pinv0[r]*kb0[r]) + ifRes[r]*Walras ⟂ yi[r])
    
    @constraint(GTAPMod, eq_chisave,
        chisave - sum(invwgt[r]*pinv[r]/pinvLag[r] for r in regions) / sum(savwgt[r]*psave[r]/psaveLag[r] for r in regions) ⟂ chisave)
 
    @constraint(GTAPMod, eq_psave[r in regions],
        psave[r]/psaveLag[r] - chisave*(pinv[r]/pinvLag[r]) ⟂ psave[r])

    @constraint(GTAPMod, eq_pfact[r in regions],
        pfact[r] - pfact0[r]*sqrt((
            (sum(peb0[f,a,r]*qfe0[f,a,r]*peb[f,a,r] for f in factors, a in activities))/
            (sum(peb0[f,a,r]*qfe0[f,a,r] for f in factors, a in activities)))
            *(
            (sum(peb0[f,a,r]*qfe0[f,a,r]*peb[f,a,r]*qfe[f,a,r] for f in factors, a in activities))/
            (sum(peb0[f,a,r]*qfe0[f,a,r]*qfe[f,a,r] for f in factors, a in activities)))) ⟂ pfact[r])
     
    @constraint(GTAPMod, eq_pfactw,
            pfactw - pfactw0*sqrt((
                (sum(peb0[f,a,r]*qfe0[f,a,r]*peb[f,a,r] for f in factors, a in activities, r in regions))/
                (sum(peb0[f,a,r]*qfe0[f,a,r] for f in factors, a in activities, r in regions))
                )*(
                (sum(peb0[f,a,r]*qfe0[f,a,r]*peb[f,a,r]*qfe[f,a,r] for f in factors, a in activities, r in regions))/
                (sum(peb0[f,a,r]*qfe0[f,a,r]*qfe[f,a,r] for f in factors, a in activities, r in regions))
                )) ⟂ pfactw)
         
        #   Numeraire definition
    @constraint(GTAPMod, eq_Walras,
        pnum - pfactw ⟂ Walras)
    
    if true
        open("stdout.txt", "w") do io
            redirect_stdout(io) do
                println(GTAPMod)
            end
        end
    end
    
    optimize!(GTAPMod)
    println(termination_status(GTAPMod))
    println(raw_status(GTAPMod))
    # @assert is_solved_and_feasible(GTAPMod)
    println(solution_summary)
    println("Factor prices across activities (pfe):")
    println(JuMP.value.(pfe))
    println("Aggregate factor prices (pe):")
    println(JuMP.value.(pe))
    println("Walras:")
    println(JuMP.value(Walras))
    println("Numeraire:")
    println(JuMP.value(pfactw))

    outscale = 1 / inscale
    open(outFolder * "/" * simName * "SAM.csv", "w") do file
        write(file, "Sim,Reg,RLab,CLab,Value\n")
        #   Cost structure
        for r in regions
            reglab = get_key(reg_ndx,r)
            for a in activities
                actlab = get_key(acts_ndx,a)
                #   Cost structure
                for i in commodities
                    comlab = get_key(comm_ndx,i)
                    x = outscale*(value(GTAPMod[:qfd][i,a,r])*value(GTAPMod[:pds][i,r])*qfd0[i,a,r]*pds0[i,r]+value(GTAPMod[:qfm][i,a,r])*value(GTAPMod[:pms][i,r])*qfm0[i,a,r]*pms0[i,r])
                    str = simName * "," * reglab * "," * comlab * "," * actlab * ","  * string(x) * "\n"
                    write(file,str)
                end
                for f in factors
                    factlab = get_key(endw_ndx,f)
                    x = outscale*value(GTAPMod[:qfe][f,a,r])*value(GTAPMod[:peb][f,a,r])*qfe0[f,a,r]*peb0[f,a,r]
                    str = simName * "," * reglab * "," * factlab * "," * actlab * ","  * string(x) * "\n"
                    write(file,str)
                end
                x = outscale*(sum(value(GTAPMod[:qfd][i,a,r])*value(GTAPMod[:pds][i,r])*qfd0[i,a,r]*pds0[i,r]*tfd[i,a,r] + value(GTAPMod[:qfm][i,a,r])*value(GTAPMod[:pms][i,r])*qfm0[i,a,r]*pms0[i,r]*tfm[i,a,r] for i in commodities))
                str = simName * "," * reglab * "," * "itax" * "," * actlab * ","  * string(x) * "\n"
                write(file,str)
                x = outscale*(sum(value(GTAPMod[:qfe][f,a,r])*value(GTAPMod[:peb][f,a,r])*qfe0[f,a,r]*peb0[f,a,r]*tfe[f,a,r] for f in factors))
                str = simName * "," * reglab * "," * "vtax" * "," * actlab * ","  * string(x) * "\n"
                write(file,str)
                #   Make
                for i in commodities
                    comlab = get_key(comm_ndx,i)
                    x = outscale*(value(GTAPMod[:qca][i,a,r])*value(GTAPMod[:ps][i,a,r])*qca0[i,a,r]*ps0[i,a,r])
                    str = simName * "," * reglab * "," * actlab * "," * comlab * ","  * string(x) * "\n"
                    write(file,str)
                end
                
            end
            for i in commodities
                comlab = get_key(comm_ndx,i)
                x = outscale*(sum(value(GTAPMod[:qca][i,a,r])*value(GTAPMod[:ps][i,a,r])*qca0[i,a,r]*ps0[i,a,r]*ptax[i,a,r] for a in activities))
                str = simName * "," * reglab * "," * "ptax" * "," * comlab * ","  * string(x) * "\n"
                write(file,str)
            end
            #   Imports/exports
            for i in commodities
                comlab = get_key(comm_ndx,i)
                for s in regions
                    srclab = get_key(reg_ndx,s)
                    x = outscale*(value(GTAPMod[:qxs][i,s,r])*value(GTAPMod[:pcif][i,s,r])*qxs0[i,s,r]*pcif0[i,s,r])
                    # x = outscale*(value(GTAPMod[:qxs][i,s,r])*value(GTAPMod[:pfob][i,s,r])*qxs0[i,s,r]*pfob0[i,s,r])
                    str = simName * "," * reglab * "," * srclab * "," * comlab * ","  * string(x) * "\n"
                    write(file,str)
                end
                # x = outscale*(sum(value(GTAPMod[:qxs][i,s,r])*(value(GTAPMod[:pcif][i,s,r])*pcif0[i,s,r] - value(GTAPMod[:pfob][i,s,r])*pfob0[i,s,r])*qxs0[i,s,r] for s in regions))
                # str = simName * "," * reglab * "," * "tmg" * "," * comlab * ","  * string(x) * "\n"
                # write(file,str)
                x = outscale*(sum(value(GTAPMod[:qxs][i,s,r])*value(GTAPMod[:pcif][i,s,r])*qxs0[i,s,r]*pcif0[i,s,r]*tms[i,s,r] for s in regions))
                str = simName * "," * reglab * "," * "mtax" * "," * comlab * ","  * string(x) * "\n"
                write(file,str)
                for d in regions
                    dstlab = get_key(reg_ndx,d)
                    x = outscale*(value(GTAPMod[:qxs][i,r,d])*value(GTAPMod[:pfob][i,r,d])*qxs0[i,r,d]*pfob0[i,r,d])
                    str = simName * "," * reglab * "," * comlab * "," * dstlab * ","  * string(x) * "\n"
                    write(file,str)
                end
                x = outscale*(sum(value(GTAPMod[:qxs][i,r,d])*value(GTAPMod[:pds][i,r])*qxs0[i,r,d]*pds0[i,r]*txs[i,r,d] for d in regions))
                str = simName * "," * reglab * "," * "etax" * "," * comlab * ","  * string(x) * "\n"
                write(file,str)
            end
            for d in regions
                dstlab = get_key(reg_ndx, d)
                x = outscale*(sum(value(GTAPMod[:qxs][i,r,d])*value(GTAPMod[:pfob][i,r,d])*qxs0[i,r,d]*pfob0[i,r,d] for i in commodities))
                str = simName * "," * reglab * "," * dstlab * "," * "BoP" * ","  * string(x) * "\n"
                write(file,str)
            end
            for s in regions
                srclab = get_key(reg_ndx, s)
                x = outscale*(sum(value(GTAPMod[:qxs][i,s,r])*value(GTAPMod[:pcif][i,s,r])*qxs0[i,s,r]*pcif0[i,s,r] for i in commodities))
                str = simName * "," * reglab * "," * "BoP" * "," * srclab * ","  * string(x) * "\n"
                write(file,str)
            end
            #   Income distribution
            for f in factors
                factlab  = get_key(endw_ndx,f)
                x = outscale*(sum(value(GTAPMod[:qfe][f,a,r])*value(GTAPMod[:pes][f,a,r])*qfe0[f,a,r]*pes0[f,a,r] for a in activities))
                str = simName * "," * reglab * "," * "RegY" * "," * factlab * ","  * string(x) * "\n"
                write(file,str)
                x = outscale*(sum(value(GTAPMod[:qfe][f,a,r])*value(GTAPMod[:peb][f,a,r])*qfe0[f,a,r]*peb0[f,a,r]*tinc[f,a,r] for a in activities))
                str = simName * "," * reglab * "," * "dtax" * "," * factlab * ","  * string(x) * "\n"
                write(file,str)
            end
            x = outscale*(value(GTAPMod[:yp][r])*yp0[r])
            str = simName * "," * reglab * "," * "hhd" * "," * "RegY" * ","  * string(x) * "\n"
            write(file,str)
            x = outscale*(value(GTAPMod[:save][r])*save0[r])
            str = simName * "," * reglab * "," * "inv" * "," * "RegY" * ","  * string(x) * "\n"
            write(file,str)
            x = outscale*(value(GTAPMod[:yg][r])*yg0[r])
            str = simName * "," * reglab * "," * "gov" * "," * "RegY" * ","  * string(x) * "\n"
            write(file,str)
            x = outscale*(value(GTAPMod[:kb][r])*kb0[r]*value(GTAPMod[:pinv][r])*pinv0[r]*depr[r])
            str = simName * "," * reglab * "," * "deprY" * "," * "RegY" * ","  * string(x) * "\n"
            write(file,str)
            str = simName * "," * reglab * "," * "inv" * "," * "deprY" * ","  * string(x) * "\n"
            write(file,str)
            #   Domestic final demand
            for i in commodities
                comlab = get_key(comm_ndx,i)
                x = outscale*(value(GTAPMod[:qpd][i,r])*value(GTAPMod[:pds][i,r])*qpd0[i,r]*pds0[i,r]+value(GTAPMod[:qpm][i,r])*value(GTAPMod[:pms][i,r])*qpm0[i,r]*pms0[i,r])
                str = simName * "," * reglab * "," * comlab * "," * "hhd" * ","  * string(x) * "\n"
                write(file,str)
                x = outscale*(value(GTAPMod[:qgd][i,r])*value(GTAPMod[:pds][i,r])*qgd0[i,r]*pds0[i,r]+value(GTAPMod[:qgm][i,r])*value(GTAPMod[:pms][i,r])*qgm0[i,r]*pms0[i,r])
                str = simName * "," * reglab * "," * comlab * "," * "gov" * ","  * string(x) * "\n"
                write(file,str)
                x = outscale*(value(GTAPMod[:qid][i,r])*value(GTAPMod[:pds][i,r])*qid0[i,r]*pds0[i,r]+value(GTAPMod[:qim][i,r])*value(GTAPMod[:pms][i,r])*qim0[i,r]*pms0[i,r])
                str = simName * "," * reglab * "," * comlab * "," * "inv" * ","  * string(x) * "\n"
                write(file,str)
            end            
            x = outscale*(sum(value(GTAPMod[:qpd][i,r])*value(GTAPMod[:pds][i,r])*qpd0[i,r]*pds0[i,r]*tpd[i,r] + value(GTAPMod[:qpm][i,r])*value(GTAPMod[:pms][i,r])*qpm0[i,r]*pms0[i,r]*tpm[i,r] for i in commodities))
            str = simName * "," * reglab * "," * "itax" * "," * "hhd" * ","  * string(x) * "\n"
            write(file,str)
            x = outscale*(sum(value(GTAPMod[:qgd][i,r])*value(GTAPMod[:pds][i,r])*qgd0[i,r]*pds0[i,r]*tgd[i,r] + value(GTAPMod[:qgm][i,r])*value(GTAPMod[:pms][i,r])*qgm0[i,r]*pms0[i,r]*tgm[i,r] for i in commodities))
            str = simName * "," * reglab * "," * "itax" * "," * "gov" * ","  * string(x) * "\n"
            write(file,str)
            x = outscale*(sum(value(GTAPMod[:qid][i,r])*value(GTAPMod[:pds][i,r])*qid0[i,r]*pds0[i,r]*tid[i,r] + value(GTAPMod[:qim][i,r])*value(GTAPMod[:pms][i,r])*qim0[i,r]*pms0[i,r]*tim[i,r] for i in commodities))
            str = simName * "," * reglab * "," * "itax" * "," * "inv" * ","  * string(x) * "\n"
            write(file,str)
            for i in commodities
                comlab = get_key(comm_ndx,i)
                x = outscale*(value(GTAPMod[:qst][i,r])*value(GTAPMod[:pds][i,r])*qst0[i,r]*pds0[i,r])
                str = simName * "," * reglab * "," * comlab * "," * "tmg" * ","  * string(x) * "\n"
                write(file,str)
            end
            x = outscale*(sum(value(GTAPMod[:qst][i,r])*value(GTAPMod[:pds][i,r])*qst0[i,r]*pds0[i,r] for i in commodities))
            str = simName * "," * reglab * "," * "tmg" * "," * "BoP" * ","  * string(x) * "\n"
            write(file,str)            
            x = outscale*(value(GTAPMod[:savf][r]))
            str = simName * "," * reglab * "," * "inv" * "," * "BoP" * ","  * string(x) * "\n"
            write(file,str)            
       end
    end

    open(outFolder * "/" * simName * ".csv", "w") do file
        write(file, "Sim,Var,Reg,RLab,CLab,Year,Value\n")
        #   Cost structure
        for r in regions
            reglab = get_key(reg_ndx,r)
            for a in activities
                actlab = get_key(acts_ndx,a)
                # Output
                x = outscale*(value(GTAPMod[:qo][a,r])*qo0[a,r])
                str = simName * "," * "QO," * reglab * "," * actlab * "," * ",Comp,"  * string(x) * "\n"
                write(file,str)
                x = outscale*(qo0[a,r])
                str = simName * "," * "QO," * reglab * "," * actlab * "," * ",Base,"  * string(x) * "\n"
                write(file,str)
                x = value(GTAPMod[:po][a,r]*po0[a,r])
                str = simName * "," * "PO," * reglab * "," * actlab * "," * ",Comp,"  * string(x) * "\n"
                write(file,str)
                x = value(po0[a,r])
                str = simName * "," * "PO," * reglab * "," * actlab * "," * ",Base,"  * string(x) * "\n"
                write(file,str)
                for i in commodities
                    comlab = get_key(comm_ndx,i)
                    x = outscale*(value(GTAPMod[:qca][i,a,r])*qca0[i,a,r])
                    str = simName * "," * "QCA," * reglab * "," * comlab * "," * actlab * ",Comp,"  * string(x) * "\n"
                    write(file,str)
                    x = outscale*(qca0[i,a,r])
                    str = simName * "," * "QCA," * reglab * "," * comlab * "," * actlab * ",Base,"  * string(x) * "\n"
                    write(file,str)
                    x = (value(GTAPMod[:ps][i,a,r])*ps0[i,a,r])
                    str = simName * "," * "PS," * reglab * "," * comlab * "," * actlab * ",Comp,"  * string(x) * "\n"
                    write(file,str)
                    x = (ps0[i,a,r])
                    str = simName * "," * "PS," * reglab * "," * comlab * "," * actlab * ",Base,"  * string(x) * "\n"
                    write(file,str)
                    x = (value(GTAPMod[:pca][i,a,r])*pca0[i,a,r])
                    str = simName * "," * "PCA," * reglab * "," * comlab * "," * actlab * ",Comp,"  * string(x) * "\n"
                    write(file,str)
                    x = (pca0[i,a,r])
                    str = simName * "," * "PCA," * reglab * "," * comlab * "," * actlab * ",Base,"  * string(x) * "\n"
                    write(file,str)
                end
            end
            for i in commodities
                comlab = get_key(comm_ndx,i)
                # Supply
                x = outscale*(value(GTAPMod[:qc][i,r])*qc0[i,r])
                str = simName * "," * "QC," * reglab * "," * comlab * "," * ",Comp,"  * string(x) * "\n"
                write(file,str)
                x = outscale*(qc0[i,r])
                str = simName * "," * "QC," * reglab * "," * comlab * "," * ",Base,"  * string(x) * "\n"
                write(file,str)
                x = value(GTAPMod[:pds][i,r]*pds0[i,r])
                str = simName * "," * "PDS," * reglab * "," * comlab * "," * ",Comp,"  * string(x) * "\n"
                write(file,str)
                x = value(pds0[i,r])
                str = simName * "," * "PDS," * reglab * "," * comlab * "," * ",Base,"  * string(x) * "\n"
                write(file,str)
            end
        end
    end
end

solve_GTAP()

stdLabels = Dict("itax" => 1, "vtax" => 2, "ptax" => 3, "mtax" => 4, "etax" => 5, "dtax" => 6, 
                "regY" => 7, "hhd" => 8, "gov" => 9, "inv" => 10, "deprY" => 11, "tmg" => 12)
bopLabel  = Dict("bop" => 1)

samLabels = merge(acts_ndx, comm_ndx, endw_ndx, stdLabels, reg_ndx, bopLabel)
# println(samLabels)


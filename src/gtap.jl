
"""
    getSet(setName::String)

Get the set definitions for the set named 'setname'
Function returns the Julia DataFrame as read in from the CSV file,
the set labels transformed into Symbol form, and a Dictionary that maps the
set labels to an index (not sure if this will be used)
"""
function getSet(setName::String; inFolder = "Data/DBUG0",BaseName = "DBUG0")
    set       = CSV.read(inFolder * "/" * BaseName * setName * ".csv", DataFrames.DataFrame)
    setLabels = Symbol.(set[!,1])
    setDict = Dict(label => i for (i,label) in enumerate(setLabels))
    return setLabels#, set, setDict
end


"""
    getData(parmName, dims...; scale = 1, inFolder = "Data/DBUG0", BaseName = "DBUG0")

Load the parameter `parmName` with dimensions `dims`.
"""
function getData(parmName, dims...; scale = 1, inFolder = "Data/DBUG0", BaseName = "DBUG0")
    #   Read the data from the CSV file
    #csv_data = CSV.read(inFolder * "/" * BaseName * parmName * ".csv", DataFrames.DataFrame)
    csv_data = CSV.File(inFolder * "/" * BaseName * parmName * ".csv")
    #   Create the container
    x = Containers.DenseAxisArray(zeros([length(d) for d∈dims]...), dims...)
    #   Loop over each row and fill the array
    for row in csv_data
        x[[Symbol(row[i]) for i∈1:length(dims)]...] = row[:value]*scale
    end
    return x
end

"""
    initVar(dims...; ifZero = true)

Create a dense array of dimensions dims with either 0's or 1's
"""
function initVar(dims...; ifZero = true)
    if ifZero
        data = zeros([length(d) for d∈dims]...)
    else
        data = ones([length(d) for d∈dims]...)
    end
    return Containers.DenseAxisArray(data, dims...)
end

"""
    initialize_sets(; inFolder = "Data/DBUG0",BaseName = "DBUG0")

Initialize the model sets--for the moment only 4 of them
The set definitions are read from the relevant CSV file, e.g., '[infolder/BaseName]ACTS.csv'
"""

function initialize_sets(; inFolder = "Data/DBUG0",BaseName = "DBUG0")

    sets = [:acts => "ACTS", :endw => "ENDW", :comm => "COMM", :reg => "REG"]

    data = Dict()

    for (set, setName) in sets
        data[set] = getSet(setName; inFolder = inFolder, BaseName = BaseName)
    end

    return data
end

"""
    read_GTAPdata!(data; inscale = 1e-6, popscale = 1e-3, inFolder = "Data/DBUG0",BaseName = "DBUG0")

Read the GTAP data:
    1. The individual components are each stored in a CSV file, e.g., '[infolder/BaseName]VDFB.csv'
       These are scaled by inscalie
    2. Population is also in a CSV file, but scaled with a different factor
    3. The GTAP elasticities are each in a separate CSV file, these are not scaled
       [For the moment] CES elasticities = 1 are converted to 1.01
"""

function read_GTAPdata!(data; inscale = 1e-6, popscale = 1e-3, inFolder = "Data/DBUG0",BaseName = "DBUG0")

    inscale_data = Dict(
        :vdfb => ("VDFB", [:comm, :acts, :reg]),
        :vdfp => ("VDFP", [:comm, :acts, :reg]),
        :vmfb => ("VMFB", [:comm, :acts, :reg]),
        :vmfp => ("VMFP", [:comm, :acts, :reg]),
        :evfb => ("EVFB", [:endw, :acts, :reg]),
        :evfp => ("EVFP", [:endw, :acts, :reg]),
        :evos => ("EVOS", [:endw, :acts, :reg]),
        :maks => ("MAKS", [:comm, :acts, :reg]),
        :makb => ("MAKB", [:comm, :acts, :reg]),
        :ptax => ("PTAX", [:comm, :acts, :reg]),
        :vdpb => ("VDPB", [:comm, :reg]),
        :vdpp => ("VDPP", [:comm, :reg]),
        :vmpb => ("VMPB", [:comm, :reg]),
        :vmpp => ("VMPP", [:comm, :reg]),
        :vdgb => ("VDGB", [:comm, :reg]),
        :vdgp => ("VDGP", [:comm, :reg]),
        :vmgb => ("VMGB", [:comm, :reg]),
        :vmgp => ("VMGP", [:comm, :reg]),
        :vdib => ("VDIB", [:comm, :reg]),
        :vdip => ("VDIP", [:comm, :reg]),
        :vmib => ("VMIB", [:comm, :reg]),
        :vmip => ("VMIP", [:comm, :reg]),
        :vxsb => ("VXSB", [:comm, :reg, :reg]),
        :vfob => ("VFOB", [:comm, :reg, :reg]),
        :vcif => ("VCIF", [:comm, :reg, :reg]),
        :vmsb => ("VMSB", [:comm, :reg, :reg]),
        :vst => ("VST",   [:comm, :reg]),
        :vtwr => ("VTWR", [:comm, :comm, :reg, :reg]),
        :vkb => ("VKB", [:reg]),
        :vdep => ("VDEP", [:reg]),
        :save => ("SAVE", [:reg]),
    )

    #   Retrieve the 'SAM' components of the GTAP database and scale using the user-input 'inscale'
    for (var, (parmName, dims)) in inscale_data
        # print("$var => $dims")
        data[var] = getData(parmName, [data[d] for d in dims]...; scale = inscale, inFolder = inFolder,BaseName = BaseName)
    end

    #   Retrieve the GTAP-based population and scale using the user-input 'popscale'
    data[:pop0] = getData("POP", data[:reg]; scale = popscale, inFolder = inFolder,BaseName = BaseName)

    no_scale_data = Dict(
        :ESUBT => ("ESUBT",   [:acts, :reg]),
        :ESUBC => ("ESUBC",   [:acts, :reg]),
        :ESUBVA => ("ESUBVA", [:acts, :reg]),
        :ETRAQ => ("ETRAQ",   [:acts, :reg]),
        :ESUBQ => ("ESUBQ",   [:comm, :reg]),
        :ESUBG => ("ESUBG",   [:reg]),
        :ESUBI => ("ESUBI",   [:reg]),
        :INCPAR => ("INCPAR", [:comm, :reg]),
        :SUBPAR => ("SUBPAR", [:comm, :reg]),
        :ESUBD => ("ESUBD",   [:comm, :reg]),
        :ESUBM => ("ESUBM",   [:comm, :reg]),
        :ETRAE => ("ETRAE",   [:endw, :reg]),
        :ESUBS => ("ESUBS",   [:comm]),
        :RORFLEX => ("RORFLEX", [:reg]),
    )

    #   Retrieve the GTAP-based elasticities
    for (var, (parmName, dims)) in no_scale_data
        data[var] = getData(parmName, [data[d] for d in dims]...; scale = 1, inFolder = inFolder,BaseName = BaseName)
    end

    #   Adjust CES substitution parameters to avoid Cobb-Douglas
    parameters_to_adjust = [
        :ESUBT, :ESUBC, :ESUBVA, :ESUBQ, :ESUBG, :ESUBI, :ESUBD, :ESUBM
    ]

    for parm in parameters_to_adjust
        data[parm] = @. ifelse(data[parm] == 1.0, 1.01, data[parm])
    end

    return data
end


"""
    initialize_production!(data; nrsFact = [:nrs])

Function to initialize the production structure of the standard GTAP model
    1. The function requires that the user input the set label(s) for sector-specific
       factors. In the standard database, there is a single factor that is sector
       specific, i.e., natural resources. The set label(s) is (are) user-determined and
       part of the aggregation facility. In the standard GTAP database it is labeled
       'NatlRes'. The default used here is 'nrs'. N.B. The user must input the
       label(s) as an array of Symbols, i.e., within brackets, even if there is
       only one set label.
"""
function initialize_production!(data; nrsFact = [:nrs])

    #   Extract the needed data for this function

    @extract(data,
        comm,
        acts,
        reg,
        endw,
        vdfb,
        vdfp,
        vmfp,
        evos,
        evfb,
        evfp,
    )

    data[:pfa0]  = initVar(comm, acts, reg; ifZero = false)
    data[:qfa0]  = (vdfp .+ vmfp) ./ data[:pfa0]
    data[:pes0]  = initVar(endw, acts, reg; ifZero = false)
    data[:peb0]  = initVar(endw, acts, reg; ifZero = false)
    data[:pfe0]  = initVar(endw, acts, reg; ifZero = false)
    data[:qfe0]  = evos ./ data[:pes0]
    data[:tinc0] = evfb .- evos
    data[:tfe0]  = evfp .- evfb
    
    data[:tinc0] = @. ifelse(evfb .!= 0, data[:tinc0]/evfb,0)
    data[:peb0]  = data[:pes0] ./ (1 .- data[:tinc0])
    data[:tfe0]  = @. ifelse(evfb .!= 0, data[:tfe0] / evfb, 0.0)
    data[:pfe0]  = data[:peb0] .* (1.0 .+ data[:tfe0])
    
    #   We should read this in future iterations
    #   Elasticity of sector specific resources
    #   Let's initialize to 0
    data[:ETAFF] = initVar(endw, acts, reg; ifZero = true)
    
    data[:pe0] = initVar(endw, reg; ifZero = false)
    
    data[:mobFact] = Containers.@container([endw], true)
    
    for n in nrsFact
        for f in endw
            if n == f
                data[:mobFact][f] = false
            end
        end
    end
    
    data[:qe0] = sum.(data[:pes0][f,acts,r] .* data[:qfe0][f,acts,r] for f in endw, r in reg) ./ data[:pe0]
    
    #   Initialize the intermediate CES bundles--qint and qva
    data[:pint0] = initVar(acts, reg; ifZero = false)
    data[:qint0] = sum.(data[:pfa0][comm,a,r] .* data[:qfa0][comm,a,r] for a in acts, r in reg) ./ data[:pint0]
    data[:pva0]  = initVar(acts, reg; ifZero = false)
    data[:qva0]  = sum.(evfp[endw,a,r] for a in acts, r in reg) ./ data[:pva0]
    data[:po0]   = initVar(acts, reg; ifZero = false)
    data[:qo0]   = (data[:pint0] .* data[:qint0]) .+ (data[:pva0] .* data[:qva0]) ./ data[:po0]

    return data
end


function initialize_make!(data)

    @extract(data,
        comm,
        acts,
        reg,
        maks,
        makb,
    )

    data[:ps0]   = initVar(comm, acts, reg; ifZero = false)
    data[:qca0]  = maks ./ data[:ps0]
    data[:ptax0] = @. ifelse(data[:qca0] .!= 0, (makb .- maks) ./ maks, 1)
    data[:pca0]  = data[:ps0] .* (1 .+ data[:ptax0])
    data[:pds0]  = initVar(comm, reg; ifZero = false)
    data[:qc0]   = sum.(data[:pca0][i,acts,r] .* data[:qca0][i,acts,r] for i in comm, r in reg) ./ data[:pds0]

    return data
end

function initialize_domestic_demand!(data)

    @extract(data,
        comm,
        acts,
        reg,
        vdpp,
        vmpp,
        vmpb,
        vdgp,
        vmgp,
        vdip,
        vmip,
    )

    data[:ppa0] = initVar(comm, reg; ifZero = false)
    data[:qpa0] = (vdpp .+ vmpp) ./ data[:ppa0]
    data[:yp0] = Containers.DenseAxisArray(sum.(data[:ppa0][comm,r] .* data[:qpa0][comm,r] for r in reg), reg)
    data[:ppa0] = initVar(comm, reg; ifZero = false)
    data[:conshr0] = initVar(comm, reg; ifZero = false)
    for i in comm, r in reg
        data[:conshr0][i,r] = data[:ppa0][i,r]*data[:qpa0][i,r] / data[:yp0][r]
    end

    data[:pga0]  = initVar(comm, reg; ifZero = false)
    data[:qga0]  = (vdgp .+ vmgp) ./ data[:pga0]
    data[:yg0]   = Containers.DenseAxisArray(sum.(data[:pga0][comm,r] .* data[:qga0][comm,r] for r in reg), reg)
    data[:pgov0] = initVar(reg; ifZero = false)
    data[:xg0]   = data[:yg0] ./ data[:pgov0]

    data[:pia0]  = initVar(comm, reg; ifZero = false)
    data[:qia0]  = (vdip .+ vmip) ./ data[:pia0]
    data[:yi0]   = Containers.DenseAxisArray(sum.(data[:pia0][comm,r] .* data[:qia0][comm,r] for r in reg), reg)
    data[:pinv0] = initVar(reg; ifZero = false)
    data[:qinv0] = data[:yi0] ./ data[:pinv0]

    return data

end

function initialize_trade!(data)

    @extract(data,
            comm,
            reg,
            acts,
            endw,
            vdfp,
            vdfb,
            vmfp,
            vmfb,
            vdpp,
            vdpb,
            vmpp,
            vmpb,
            vmgp,
            vmgb,
            vdgp,
            vdgb,
            vmip,
            vdib,
            vmip,
            vmib,
            pds0,
            yi0,
            INCPAR,
            conshr0,
            vdip,
            vdib,
            vfob,
            evfp,
            evfb,
            makb,
            maks,
            vxsb,
            vmsb,
            vcif,
            vdep,
            vkb,
            qinv0,
            peb0,
            qfe0,
            yp0,
            yg0,
            pop0,
            vtwr,
            vst,
            save
    )

    data[:tfd0] = @. ifelse(vdfb .!= 0, (vdfp .- vdfb) ./ vdfb,0)
    data[:tfm0] = @. ifelse(vmfb .!= 0, (vmfp .- vmfb) ./ vmfb,0)
    data[:tpd0] = @. ifelse(vdpb .!= 0, (vdpp .- vdpb) ./ vdpb,0)
    data[:tpm0] = @. ifelse(vmpb .!= 0, (vmpp .- vmpb) ./ vmpb,0)
    data[:tgd0] = @. ifelse(vdgb .!= 0, (vdgp .- vdgb) ./ vdgb,0)
    data[:tgm0] = @. ifelse(vmgb .!= 0, (vmgp .- vmgb) ./ vmgb,0)
    data[:tid0] = @. ifelse(vdib .!= 0, (vdip .- vdib) ./ vdib,0)
    data[:tim0] = @. ifelse(vmib .!= 0, (vmip .- vmib) ./ vmib,0)

    data[:pms0] = initVar(comm, reg; ifZero = false)

    #   Probably can be converted to something more succinct
    data[:pfd0] = (initVar(comm, acts, reg; ifZero = false))
    data[:pfm0] = initVar(comm, acts, reg; ifZero = false)
    data[:ppd0] = initVar(comm, reg; ifZero = false)
    data[:ppm0] = initVar(comm, reg; ifZero = false)
    data[:pgd0] = initVar(comm, reg; ifZero = false)
    data[:pgm0] = initVar(comm, reg; ifZero = false)
    data[:pid0] = initVar(comm, reg; ifZero = false)
    data[:pim0] = initVar(comm, reg; ifZero = false)
    for i in comm
        for r in reg
            for a in acts
                data[:pfd0][i,a,r] = (1 + data[:tfd0][i,a,r]) * pds0[i,r]       
                data[:pfm0][i,a,r] = (1 + data[:tfm0][i,a,r]) * data[:pms0][i,r]
            end
            data[:ppd0][i,r] = (1 + data[:tpd0][i,r]) * pds0[i,r]
            data[:ppm0][i,r] = (1 + data[:tpm0][i,r]) * data[:pms0][i,r]        
            data[:pgd0][i,r] = (1 + data[:tgd0][i,r]) * pds0[i,r]
            data[:pgm0][i,r] = (1 + data[:tgm0][i,r]) * data[:pms0][i,r]        
            data[:pid0][i,r] = (1 + data[:tid0][i,r]) * pds0[i,r]
            data[:pim0][i,r] = (1 + data[:tim0][i,r]) * data[:pms0][i,r]        
        end
    end

    data[:qfd0] = vdfp ./ data[:pfd0]
    data[:qfm0] = vmfp ./ data[:pfm0]
    data[:qpd0] = vdpp ./ data[:ppd0]
    data[:qpm0] = vmpp ./ data[:ppm0]
    data[:qgd0] = vdgp ./ data[:pgd0]
    data[:qgm0] = vmgp ./ data[:pgm0]
    data[:qid0] = vdip ./ data[:pid0]
    data[:qim0] = vmip ./ data[:pim0]

    data[:save0] = deepcopy(save)
    data[:savf0] = yi0 .- data[:save0] .- vdep
    #   Add any residual inconsistency to the largest capital account balance (in absolute value)
    #   !!!! There is most likely a more elegant way to do this in Julia
    resid = sum(data[:savf0][r] for r in reg)
    rmax  = 0
    sfmax = 0
    for r in reg
        if abs(data[:savf0][r] > sfmax)
            sfmax = abs(data[:savf0][r])
            rmax = r
        end
    end
    data[:savf0][rmax] = data[:savf0][rmax] .- resid

    data[:uepriv0] = Containers.DenseAxisArray(sum.(conshr0[comm,r] .* INCPAR[comm,r] for r in reg), reg)

    data[:taxrpc0]  = Containers.DenseAxisArray(sum.((vdpp[comm,r] .- vdpb[comm,r]) .+ (vmpp[comm,r] .- vmpb[comm,r]) for r in reg), reg)
    data[:taxrgc0]  = Containers.DenseAxisArray(sum.((vdgp[comm,r] .- vdgb[comm,r]) .+ (vmgp[comm,r] .- vmgb[comm,r]) for r in reg), reg)
    data[:taxric0]  = Containers.DenseAxisArray(sum.((vdip[comm,r] .- vdib[comm,r]) .+ (vmip[comm,r] .- vmib[comm,r]) for r in reg), reg)
    data[:taxriu0]  = Containers.DenseAxisArray(sum.((vdfp[comm,acts,r] .- vdfb[comm,acts,r]) .+ (vmfp[comm,acts,r] .- vmfb[comm,acts,r]) for r in reg), reg)
    data[:taxrfu0]  = Containers.DenseAxisArray(sum.((evfp[endw,acts,r] .- evfb[endw,acts,r]) for r in reg), reg)
    data[:taxrout0] = Containers.DenseAxisArray(sum.((makb[comm,acts,r] .- maks[comm,acts,r]) for r in reg), reg)
    data[:taxrexp0] = Containers.DenseAxisArray(sum.((vfob[comm,r,reg] .- vxsb[comm,r,reg]) for r in reg), reg)
    data[:taxrimp0] = Containers.DenseAxisArray(sum.((vmsb[comm,reg,r] .- vcif[comm,reg,r]) for r in reg), reg)
    data[:indtax0]  = data[:taxrpc0] .+ data[:taxrgc0] .+ data[:taxric0] .+ data[:taxriu0] .+ data[:taxrfu0] .+ data[:taxrout0] .+ data[:taxrexp0] .+ data[:taxrimp0]

    data[:depr]     = vdep ./ vkb
    data[:kb0]      = deepcopy(vkb)
    data[:ke0]      = (1 .- data[:depr]) .* data[:kb0] .+ qinv0
    data[:pinv0]    = initVar(reg; ifZero = false)
    data[:fincome0] = Containers.DenseAxisArray(sum.((peb0[endw,acts,r] .* qfe0[endw,acts,r]) for r in reg), reg)
    data[:fincome0] = data[:fincome0] .- data[:depr] .* data[:pinv0] .* data[:kb0]
    data[:y0]       = data[:fincome0] .+ data[:indtax0]
    data[:gdpfc0]   = deepcopy(data[:y0])

    data[:dppriv] = yp0 ./ data[:y0]
    data[:dpgov]  = yg0 ./ data[:y0]
    data[:dpsave] = (1.0 .- data[:dppriv] .- data[:dpgov])
    data[:uelas0] = 1.0 ./ ((data[:dppriv] ./ data[:uepriv0]) .+ data[:dpgov] .+ data[:dpsave])

    data[:u0]  = initVar(reg; ifZero = false)
    data[:up0] = initVar(reg; ifZero = false)
    data[:ug0] = initVar(reg; ifZero = false)
    data[:us0] = initVar(reg; ifZero = false)

    data[:au] = ((data[:up0] .^ data[:dppriv]) .* (data[:ug0] .^ data[:dpgov]) .* (data[:us0] .^ data[:dpsave]))
    data[:au] = data[:u0] ./ data[:au]
    data[:aus] = data[:us0] .* pop0 ./ data[:save0]

    #   data[:tmarg0] will be defined with respect to the base volume, not the vfob value
    data[:txs0]   = @. ifelse(vxsb .!= 0, (vfob .- vxsb) ./ vxsb, 0)
    data[:tmarg0] = @. ifelse(vxsb .!= 0, (vcif .- vfob) ./ vxsb, 0)
    data[:tms0]   = @. ifelse(vcif .!= 0, (vmsb .- vcif) ./ vcif, 0)

    #   Set some tolerance level for data[:tmarg0]
    data[:tmarg0] = @. ifelse(abs.(data[:tmarg0]) > 1e-6, data[:tmarg0], 0)
    for i ∈ comm, s ∈ reg, d ∈ reg
        if data[:tmarg0][i,s,d] == 0.0
            for m ∈ comm
                vtwr[m,i,s,d] = 0
            end
        end
    end

    data[:pfob0] = initVar(comm, reg, reg; ifZero = false)
    for i in comm, s in reg, d in reg
        data[:pfob0][i,s,d] = (1 + data[:txs0][i,s,d]) * pds0[i,s]
    end
    data[:pcif0] = data[:pfob0] .+ data[:tmarg0]
    data[:pmds0] = data[:pcif0] .* (1 .+ data[:tms0])

    data[:pms0] = initVar(comm, reg; ifZero = false)
    data[:qms0] = sum.(data[:qfm0][i,acts,r] for i in comm, r in reg) .+ data[:qpm0] .+ data[:qgm0] .+ data[:qim0]
    data[:qxs0] = vfob ./ data[:pfob0]

    data[:ptrans0] = initVar(comm, reg, reg; ifZero = false)
    data[:qst0]    = vst ./ pds0
    data[:qds0]    = sum.(data[:qfd0][i,acts,r] for i in comm, r in reg) .+ data[:qpd0] .+ data[:qgd0] .+ data[:qid0]

    data[:pt0]     = initVar(comm; ifZero = false)
    data[:qtm0]    = initVar(comm; ifZero = true)
    data[:qtmfsd0] = deepcopy(vtwr)
    for m in comm, i in comm, s in reg, d in reg
        data[:qtmfsd0][m,i,s,d] = data[:qtmfsd0][m,i,s,d] / data[:pt0][m]
        data[:qtm0][m] = data[:qtm0][m] + data[:qtmfsd0][m,i,s,d]
    end

    return data

end

"""
    initialize_capital_closure!(data; resName = :USA, capName = :cap)

Function to initialize the model's capital closure
    1. The function requires that the user input the set label for the
       residual region. This is a single set label and can be any of the
       model's regions (from the aggregation).
       The function requires that the user input the set label for the
        capital factor. This is a single set label and is user-determined
        and is derived from the aggregation facility. In the standard
        GTAP database it is labeled 'Capital'.
"""
function initialize_capital_closure!(data; resName = :USA, capName = :cap)

    @extract(data,
        comm,
        acts,
        reg,
        pes0,
        qfe0,
        kb0,
        ke0,
        pinv0,
        depr,
        RORFLEX,
        qinv0,
        save
    )
    
    data[:rental0] = initVar(reg; ifZero = true)
    data[:rental0] = sum.(pes0[capName,acts,r] .* qfe0[capName,acts,r] for r in reg) ./ kb0
    
    data[:rorc0] = data[:rental0] ./ pinv0 .- depr
    data[:rore0] = data[:rorc0] .* (kb0 ./ ke0) .^ RORFLEX
    data[:rorg0] = sum(data[:rore0][r]*pinv0[r]*(qinv0[r] - depr[r]*kb0[r]) for r in reg) / sum(pinv0[r]*(qinv0[r] - depr[r]*kb0[r]) for r in reg)
    data[:risk]  = data[:rorg0] ./ data[:rore0]
    data[:ifRes] = initVar(reg; ifZero = true)
    data[:ifRes][resName] = 1
    data[:pinvLag]  = initVar(reg; ifZero = false)
    data[:psaveLag] = initVar(reg; ifZero = false)
    data[:invwgt]   = initVar(reg; ifZero = true)
    data[:savwgt]   = initVar(reg; ifZero = true)
    
    for r in reg
        isum = sum(pinv0[d]*(qinv0[d] - depr[d]*kb0[d]) for d in reg)
        data[:invwgt][r] = pinv0[r]*(qinv0[r] - depr[r]*kb0[r]) / isum
        ssum = sum(save[d] for d in reg)
        data[:savwgt][r] = save[r] / ssum
    end
    data[:pfact0]  = initVar(reg; ifZero = false)
    data[:pfactw0] = 1

    return data
end

function calibrate_parameters!(data)

    @extract(data,
        comm,
        reg,
        acts,
        endw,
        qo0,
        pint0,
        qint0,
        pva0,
        qva0,
        pva0,
        qva0,
        pfa0,
        qfa0,
        pfa0,
        qfa0,
        pfe0,
        qfe0,
        ps0,
        qca0,
        pca0,
        qca0,
        pga0,
        qga0,
        pia0,
        qia0,
        qc0,
        yg0,
        yi0
    )



    data[:shr_qint] = @. ifelse(qo0 .!= 0, (pint0 .* qint0) ./ (pint0 .* qint0 .+ pva0 .* qva0),0)
    data[:shr_qva]  = @. ifelse(qo0 .!= 0, (pva0 .* qva0) ./ (pint0 .* qint0 .+ pva0 .* qva0),0)
    
    data[:shr_qfa] = initVar(comm, acts, reg; ifZero = true)
    data[:shr_qfe] = initVar(endw, acts, reg; ifZero = true)
    data[:shr_qcas] = initVar(comm, acts, reg; ifZero = true)
    data[:shr_qcad] = initVar(comm, acts, reg; ifZero = true)
    data[:shr_gov] = initVar(comm, reg; ifZero = true)
    data[:shr_inv] = initVar(comm, reg; ifZero = true)
    
    for r in reg
        #   Production
        for a in acts
            for i in comm
                if qint0[a,r] != 0
                    data[:shr_qfa][i,a,r] = pfa0[i,a,r]*qfa0[i,a,r] / sum(pfa0[j,a,r]*qfa0[j,a,r] for j in comm)
                end
            end
            for f in endw
                if qva0[a,r] != 0
                    data[:shr_qfe][f,a,r] = pfe0[f,a,r]*qfe0[f,a,r] / sum(pfe0[e,a,r]*qfe0[e,a,r] for e in endw)
                end
            end
            for i in comm
                if qo0[a,r] != 0
                    data[:shr_qcas][i,a,r] = ps0[i,a,r]*qca0[i,a,r] / sum(ps0[j,a,r]*qca0[j,a,r] for j in comm)
                end
            end
        end
        for i in comm
            #   Make
            for a in acts
                if qc0[i,r] != 0
                    data[:shr_qcad][i,a,r] = pca0[i,a,r]*qca0[i,a,r] / sum(pca0[i,b,r]*qca0[i,b,r] for b in acts)
                end
            end
            #   Final demand
            if yg0[r] != 0
                data[:shr_gov][i,r] = pga0[i,r]*qga0[i,r] / sum(pga0[j,r]*qga0[j,r] for j in comm)
            end
            if yi0[r] != 0
                data[:shr_inv][i,r] = pia0[i,r]*qia0[i,r] / sum(pia0[j,r]*qia0[j,r] for j in comm)
            end
        end
    end

    return data
end



function calibrate_trade_shares!(data)

    @extract(data,
        comm,
        reg,
        pfd0,
        qfd0,
        pfm0,
        qfm0,
        ppd0,
        qpd0,
        ppm0,
        qpm0,
        pgd0,
        qgd0,
        pgm0,
        qgm0,
        pid0,
        qid0,
        pim0,
        qim0,
        pmds0,
        qxs0,
        qtmfsd0,
        vst,
        vtwr,
        pt0,
        pds0,
        ESUBS
    )
    
    data[:shr_qfd]  = (pfd0 .* qfd0) ./ (pfd0 .* qfd0 .+ pfm0 .* qfm0)
    data[:shr_qfm]  = (pfm0 .* qfm0) ./ (pfd0 .* qfd0 .+ pfm0 .* qfm0)
    data[:shr_qpd]  = (ppd0 .* qpd0) ./ (ppd0 .* qpd0 .+ ppm0 .* qpm0)
    data[:shr_qpm]  = (ppm0 .* qpm0) ./ (ppd0 .* qpd0 .+ ppm0 .* qpm0)
    data[:shr_qgd]  = (pgd0 .* qgd0) ./ (pgd0 .* qgd0 .+ pgm0 .* qgm0)
    data[:shr_qgm]  = (pgm0 .* qgm0) ./ (pgd0 .* qgd0 .+ pgm0 .* qgm0)
    data[:shr_qid]  = (pid0 .* qid0) ./ (pid0 .* qid0 .+ pim0 .* qim0)
    data[:shr_qim]  = (pim0 .* qim0) ./ (pid0 .* qid0 .+ pim0 .* qim0)


    data[:shr_qxs]  = initVar(comm, reg, reg; ifZero = true)
    for d in reg, i in comm
        totimp = sum(pmds0[i,s,d]*qxs0[i,s,d] for s in reg)
        if totimp != 0
            for s in reg
                data[:shr_qxs][i,s,d] = pmds0[i,s,d] * qxs0[i,s,d] / totimp
            end
        end
    end

    #   Transportation shares
    data[:shr_vtwr]  = deepcopy(qtmfsd0)

    for i in comm, s in reg, d in reg
        wsum = sum(vtwr[m,i,s,d] for m in comm)
        if wsum != 0
            #   There is transport for this node
            for m in comm
                data[:shr_vtwr][m,i,s,d] = vtwr[m,i,s,d] / wsum 
            end
        else
            #   There is no transport for this node
            for m in comm
                data[:shr_vtwr][m,i,s,d] = 0 
            end
        end
    end

    data[:shr_qst]  = initVar(comm, reg; ifZero = true)
    data[:A_qtm]    = initVar(comm; ifZero = true)
    for m in comm
        totmarg = sum(vst[m,r] for r in reg)
        if totmarg != 0
            for r in reg
                data[:shr_qst][m,r] = vst[m,r] / totmarg
            end
            if ESUBS[m] == 1.0
                data[:A_qtm][m] = pt0[m] / prod((pds0[m,r]/data[:shr_qst][m,r])^data[:shr_qst][m,r] for r in reg)
            else
                data[:A_qtm][m] = 1.0
            end
        end
    end

    return data
end


function calibrate_factor_supply_shares!(data)


    @extract(data,
        reg,
        acts,
        endw,
        pes0,
        qfe0        
    )

#   Factor allocation shares

    data[:shr_qfes]  = initVar(endw, acts, reg; ifZero = true)
    for r in reg, f in endw
        totinc = sum(pes0[f,a,r] * qfe0[f,a,r] for a in acts)
        if totinc != 0
            for a in acts
                data[:shr_qfes][f,a,r] = pes0[f,a,r] * qfe0[f,a,r] / totinc
            end
        end
    end

    return data
end

function calibrate_consumption!(data)

    @extract(data,
        comm,
        reg,
        yp0,
        pop0,
        ppa0,
        up0,
        xg0,
        conshr0,
        SUBPAR,
        INCPAR
    )

    #   Private consumption

    data[:alphac] = initVar(comm, reg; ifZero = true)

    for r in reg
        denom = sum(ifelse(SUBPAR[j,r] != 0.0, conshr0[j,r] / SUBPAR[j,r], 0.0) for j in comm)
        for i in comm
            if SUBPAR[i,r] != 0.0
                data[:alphac][i,r] = (conshr0[i,r] / SUBPAR[i,r]) * (((yp0[r] / pop0[r]) / ppa0[i,r]) ^ SUBPAR[i,r]) * (up0[r]^(-INCPAR[i,r] * SUBPAR[i,r])) / denom
            end
        end
    end

    data[:zshr0]  = initVar(comm, reg; ifZero = true)
    for r in reg, i in comm
        data[:zshr0][i,r] = data[:alphac][i,r] * SUBPAR[i,r] * (up0[r] ^ (INCPAR[i,r] * SUBPAR[i,r])) * ((ppa0[i,r] / (yp0[r] / pop0[r])) ^ SUBPAR[i,r])    
    end

    #   Government utility
    data[:aug] = pop0 ./ xg0

    return data
end


function load_gtap_data(; inFolder = "Data/DBUG0", BaseName = "DBUG0", inscale = 1e-6, popscale = 1e-3, resName = :USA, capName = :cap, nrsFact = [:nrs])

    #   Initialize the set definitions
    data = initialize_sets(inFolder = inFolder, BaseName = BaseName)
    
    #   Read the GTAP database
    read_GTAPdata!(data; inscale = inscale, popscale = popscale, inFolder = inFolder, BaseName = BaseName, )
   
    #   Initialize model variables
    initialize_production!(data; nrsFact = nrsFact)
    initialize_make!(data)
    initialize_domestic_demand!(data)
    initialize_trade!(data)
    initialize_capital_closure!(data; resName = resName, capName = capName)
    
    #   Calibrate share and other parameters
    calibrate_parameters!(data) 
    calibrate_trade_shares!(data)
    calibrate_factor_supply_shares!(data)
    calibrate_consumption!(data)
    
    return data
end
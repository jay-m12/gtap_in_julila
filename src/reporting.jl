using JuMP
using DataFrames

function convert_denseaxis_to_dataframe(data; dim_names::Vector{Symbol} = Vector{Symbol}(), value_col::Symbol=:value)
    if isa(data, Float64)  
        println("Warning: values is a Float64, converting to a single-row DataFrame")
        return DataFrame(:value => [data], :reg => ["Unknown"], :comm => ["Unknown"], :acts => ["Unknown"])
    elseif isa(data, JuMP.Containers.DenseAxisArray)
        if length(dim_names) == 0
            dim_names = [Symbol("dim$i") for i in 1:length(data.axes)]
        end
        if length(dim_names) != length(data.axes)
            throw(ArgumentError("Length of given name list does not fit the number of variable dimensions"))
        end
        tup_dim = (dim_names...,)
        ind = reshape([collect(k[i] for i in 1:length(dim_names)) for k in Base.Iterators.product(data.axes...)],:,1)
        return DataFrame(merge(NamedTuple{tup_dim}(i), NamedTuple{(value_col,)}(data[(i...,)...])) for i in ind)
    else
        throw(ArgumentError("Unsupported data type: $(typeof(data))"))
    end
end

function name_of_operation(GTAPMod, data, year; outscale = 1/1e-6)
    out = []
    acts = data[:acts]
    comms = data[:comm]
    regs = data[:reg]

    qfd_val = JuMP.value.(GTAPMod[:qfd])
    pds_val = JuMP.value.(GTAPMod[:pds])
    qfd0_val = JuMP.value.(data[:qfd0])
    pds0_val = JuMP.value.(data[:pds0])

    qfm_val = JuMP.value.(GTAPMod[:qfm])
    pms_val = JuMP.value.(GTAPMod[:pms])
    qfm0_val = JuMP.value.(data[:qfm0])
    pms0_val = JuMP.value.(data[:pms0])

    for a ∈ acts
        for c ∈ comms
            for r ∈ regs
                values = @. outscale * (
                    qfd_val[c,a,r] * pds_val[c,r] * qfd0_val[c,a,r] * pds0_val[c,r] +
                    qfm_val[c,a,r] * pms_val[c,r] * qfm0_val[c,a,r] * pms0_val[c,r]
                )

                println("Type of values: ", typeof(values))

                if isa(values, Float64)
                    println("Warning: values is a Float64, converting to a single-row DataFrame")
                    df = DataFrame(:value => [values], :reg => [r], :comm => [c], :acts => [a])
                else
                    df = convert_denseaxis_to_dataframe(values, dim_names = [:reg, :comm, :acts])
                end

                if :reg ∉ names(df)
                    df[!, :reg] .= "Unknown"
                end
                if :comm ∉ names(df)
                    df[!, :comm] .= "Unknown"
                end
                if :acts ∉ names(df)
                    df[!, :acts] .= "Unknown"
                end

                df = df |>
                    x -> transform(x, :reg => ByRow(y -> a) => :acts, :reg => ByRow(y -> year) => :year) |>
                    x -> select(x, [:reg, :comm, :acts, :year, :value])

                push!(out, df)
            end
        end
    end

    return vcat(out...)
end

"""

Credit to [this answer on discourse](https://discourse.julialang.org/t/extracting-jump-results/51429/6)
"""
function convert_denseaxis_to_dataframe(
    data::JuMP.Containers.DenseAxisArray;
    dim_names::Vector{Symbol} = Vector{Symbol}(), 
    value_col::Symbol=:value)

    if isempty(data)
        return DataFrame()
    end

    if length(dim_names) == 0
        dim_names = [Symbol("dim$i") for i in 1:length(data.axes)]
    end

    if length(dim_names) != length(data.axes)
        throw(ArgumentError("Length of given name list does not fit the number of variable dimensions"))
    end

    tup_dim = (dim_names...,)
    ind = reshape([collect(k[i] for i in 1:length(dim_names)) for k in Base.Iterators.product(data.axes...)],:,1)

    df = DataFrame(merge(NamedTuple{tup_dim}(i), NamedTuple{(value_col,)}(data[(i...,)...])) for i in ind)

    return df
end






"""
    name_of_operation(GTAPMod, data, year; outscale = 1/1e-6)

I don't know what this should be called, but it is the first step of the `saveSAM` function.
This should be done for each of the for loops. 

Returns a dataframe with the columns `:reg`, `:comm`, `:acts`, `:year`, and `:value`.
"""
function name_of_operation(GTAPMod, data, year; outscale = 1/1e-6)

    out = []

    for a∈acts
        values = @. outscale*value(GTAPMod[:qfd][comm,a,reg]*GTAPMod[:pds][comm,reg]*data[:qfd0][comm,a,reg]*data[:pds0][comm,reg]+GTAPMod[:qfm][comm,a,reg]*GTAPMod[:pms][comm,reg]*data[:qfm0[comm,a,reg]*pms0[comm,reg]])
        df = convert_denseaxis_to_dataframe(values, dim_names = [:reg, :comm]) |>
            x -> transform(x,
                :reg => ByRow(y -> a) => :acts,
                :reg => ByRow(y -> year) => :year
            ) |>
            x -> select(x, [:reg,:comm, :acts, :year, :value])
        push!(out, df)
    end


    return vcat(out...)
end


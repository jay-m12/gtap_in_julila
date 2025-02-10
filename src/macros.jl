macro extract(data, vars...)
    code = quote end

    for val in vars
        push!(code.args, :($(esc(val)) = $(esc(data))[$(QuoteNode(val))])) 
    end

    return code

end
export rfftsize
export rfftfreq

#
# Real FFT Frequency
#
@inline function rfftsize(dims::Number)::Tuple
    return tuple(dims ÷ 2 + 1)
end

@inline function rfftsize(dims::Number...)::Tuple
    if length(dims) > 1
        return (dims[1] ÷ 2 + 1, dims[2:end]...)
    end
end

@inline function AbstractFFTs.rfftfreq(dims)::Tuple
    return Tuple((
        if i == 1
            rfftfreq(dims[i], 1)
        else
            fftfreq(dims[i], 1)
        end
    ) for i in 1:length(dims))
end

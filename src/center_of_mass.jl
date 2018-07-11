
function centerofmass(pos::Vector{T},
                         lo::T,
                         L::T,
                         weight::Vector{Float64} = ones(Float64,length(pos))
                        ) where T
    if length(pos)==1 ; return pos[1] ; end
    
    length(weight) == length(pos) || error("Weights and positions mismatch")
    
    invtotweight = 1.0/sum(weight)
    i2pi = 0.5/pi
    
    a = zero(eltype(pos))
    b = zero(eltype(pos))

    for i in eachindex(pos)
        t = (pos[i]-lo)./L*2pi 
        a -= cos.(t) * weight[i]
        b -= sin.(t) * weight[i]
    end
    a = a*invtotweight
    b = b*invtotweight

    tcom = atan.(b,a)+pi
    com = tcom.*i2pi.*L + lo
end



function center_of_mass(
    pos::Vector{T},
    box::SimulationBox{T},
    weight = ones(Float64,length(pos))
) where T
    if length(pos)==1 ; return pos[1] ; end

    L = box.lengths
    
    length(weight) == length(pos) || error("Weights and positions mismatch")
    
    invtotweight = 1.0/sum(weight)
    i2pi = 0.5/pi
    
    a = zero(eltype(pos))
    b = zero(eltype(pos))

    for i in eachindex(pos)
        t = pos[i]./L*2pi 
        a -= cos.(t) * weight[i]
        b -= sin.(t) * weight[i]
    end
    a = a*invtotweight
    b = b*invtotweight

    tcom = atan.(b,a)+pi
    com = tcom.*i2pi.*L
end


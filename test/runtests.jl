include("../src/SimulationBoxes.jl")

module TestSimulationBoxes

using FactCheck
using SimulationBoxes
import FixedSizeArrays
noerror = not(NaN)

#a 3D vector for testing
typealias Vec FixedSizeArrays.Vec{3,Float64}

typealias BoxV Box{Vec,3,(true,true,true)}
typealias BoxA Box{Vector{Float64},3,(true,true,true)}

facts("Test Box constructor") do
    origin = Vec(0,1,2)
    vectors  = ( Vec(3,0,0),
                 Vec(0,4,0),
                 Vec(0,0,5))
    avectors = ([3.,0.,0.],
                [0.,4.,0.],
                [0.,0.,5.])
    avectors2 = ([3.,0.],
                 [0.,4.])
    ppp = (true,true,true)
    pp = (true,true)
    @fact Box{Vec,3,ppp}(vectors, origin) --> noerror
    box = BoxV(vectors, origin) 
    @fact_throws MethodError Box{Vec,2,ppp}(vectors, origin)
    @fact_throws ErrorException Box{Vec,3,pp}(vectors, origin)

    @fact Box{Vector{Float64},3,ppp}(avectors, Float64[0,1,2]) --> noerror
    @fact Box{Vector{Float64},2,pp}(avectors2, Float64[0,1]) --> noerror
    @fact_throws MethodError Box{Vector{Float64},2,pp}(
        avectors, Float64[0,1])
    @fact_throws DimensionMismatch Box{Vector{Float64},2,pp}(
        avectors2, Float64[0,1,2])
    @fact BoxV(vectors, origin) --> Box(Vec(3,4,5), origin)
end

facts("Test getters") do
    p = (true,false,true)
    vectors  = ( Vec(3,0,0),
                 Vec(0,4,0),
                 Vec(0,0,5))
    origin = Vec(0,1,2)
    box = Box{Vec,3,p}(vectors, origin)
    @fact isperiodic(box) --> p
    @fact getorigin(box) --> origin
    @fact getvectors(box) --> vectors
    @fact getdiagonal(box) --> convert(Vec, collect(vectors[i][i] for i in 1:3) )
end

for (makeV,V) in zip([  (x,y,z) -> Vec(x,y,z), (x,y,z) -> Float64[x,y,z]  ],
                     [  Vec, Vector{Float64}       ]
                     )

    facts("Test functions with V=$V") do
        origin, sides = makeV(0,1,2), makeV(3,4,5)
        hi = origin + sides

        box = Box(sides, origin, (true,true,true))
        boxpfp = Box(sides,origin, (true,false,true))

        v = makeV(-1,3,8)
        v1 = origin+makeV(0.5,1,1)
        v2 = origin+sides+makeV(-0.5,-1,-1)
        
        context("test dimensions/coordinates access functions") do
            @fact getorigin(box) --> origin
            @fact getdiagonal(box) --> sides
        end
        
        
        context("Test wrapping and unwrapping") do
            @fact wrap(v,box) --> (makeV(2,3,3),(-1,0,1))
            @fact wrap(v1,box) --> (v1,(0,0,0))
            @fact wrap(origin,box) --> (origin,(0,0,0))
            @fact wrap(origin+sides,box) --> (origin,(1,1,1))
            for img in ((0,2,4), [0,2,4], makeV(0,2,4))
                @fact unwrap(v1, img, box) --> makeV(0.5, 10, 23)
                @fact unwrap(v1, img, boxpfp) --> makeV(0.5, 10, 23)
            end
        end
        
        context("Test separation") do
            @fact separation(v1,v2,box) --> makeV(1,-2,2)
            @pending separation(need,more,tests) --> makeV(x,x,x)
        end
        
#       context("Test isinside") do
#           @fact isinside(v1, box) --> true
#           @fact isinside(makeV(-10,-10,-10),box) --> false
#           @fact isinside(makeV(1,-10,-10),box) --> false
#           @fact isinside(lo, box) --> true
#           @fact isinside(hi, box) --> false
#       end
        
    end
end

end


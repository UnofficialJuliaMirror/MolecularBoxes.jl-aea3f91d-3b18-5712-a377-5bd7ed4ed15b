
using SimulationBoxes

using Test
using StaticArrays

#a 3D vector for testing
const Vec = SVector{3,Float64}

const BoxV =  Box{Vec,3,(true,true,true)}

@testset "Test Box constructor" begin
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
    @test Box{Vec,3,ppp}(vectors, origin) != nothing
    box = BoxV(vectors, origin) 
    @test_throws MethodError Box{Vec,2,ppp}(vectors, origin)
    @test_throws ErrorException Box{Vec,3,pp}(vectors, origin)

    @test Box{Vector{Float64},3,ppp}(avectors, Float64[0,1,2]) != nothing
    @test Box{Vector{Float64},2,pp}(avectors2, Float64[0,1]) != nothing
    @test_throws MethodError Box{Vector{Float64},2,pp}(avectors, [0.,1.])
    @test_throws DimensionMismatch Box{Vector{Float64},2,pp}(avectors2, [0.,1.,2.])
    @test BoxV(vectors, origin) == Box(Vec(3,4,5), origin)
end

@testset "Test getters" begin
    p = (true,false,true)
    vectors  = (
        Vec(3,0,0),
        Vec(0,4,0),
        Vec(0,0,5),
    )
    origin = Vec(0,1,2)
    box = Box{Vec,3,p}(vectors, origin)
    @test isperiodic(box) == p
    @test getorigin(box) == origin
    @test getvectors(box) == vectors
    @test getdiagonal(box) == convert(Vec, collect(vectors[i][i] for i in 1:3) )
end


origin, sides = Vec(0,1,2), Vec(3,4,5)
hi = origin + sides

box = Box(sides, origin, (true,true,true))
boxpfp = Box(sides,origin, (true,false,true))

v = Vec(-1,3,8)
v1 = origin+Vec(0.5,1,1)
v2 = origin+sides+Vec(-0.5,-1,-1)

@testset "test dimensions/coordinates access functions" begin
    @test getorigin(box) == origin
    @test getdiagonal(box) == sides
end

@testset "Test wrapping and unwrapping" begin
    @test wrap(v,box) == (Vec(2,3,3),(-1,0,1))
    @test wrap(v1,box) == (v1,(0,0,0))
    @test wrap(origin,box) == (origin,(0,0,0))
    @test wrap(origin+sides,box) == (origin,(1,1,1))
    for img in ((0,2,4), [0,2,4], Vec(0,2,4))
        @test unwrap(v1, img, box) == Vec(0.5, 10, 23)
        @test unwrap(v1, img, boxpfp) == Vec(0.5, 10, 23)
    end
end

@testset "Test separation" begin
    @test separation(v1,v2,box) == Vec(1,-2,2)
    #@pending separation(need,more,tests) --> Vec(x,x,x)
end
   


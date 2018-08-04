
using SimulationBoxes

using Test
using StaticArrays

#a 3D vector for testing
const Vec = SVector{3,Float64}

const BoxV =  Box{Vec,3,(true,true,true)}

@testset "Test Box constructor" begin
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
    @test Box{Vec,3,ppp}(vectors) != nothing
    box = BoxV(vectors) 
    @test_throws MethodError Box{Vec,2,ppp}(vectors)
    @test_throws ErrorException Box{Vec,3,pp}(vectors)

    @test BoxV(vectors) == Box(Vec(3,4,5))
end

@testset "Test getters" begin
    p = (true,false,true)
    vectors  = (
        Vec(3,0,0),
        Vec(0,4,0),
        Vec(0,0,5),
    )
    box = Box{Vec,3,p}(vectors)
    @test isperiodic(box) == p
    @test box.vectors == vectors
    @test box.lengths == convert(Vec, collect(vectors[i][i] for i in 1:3) )
end


lengths = Vec(3,4,5)
hi = lengths

box = Box(lengths, (true,true,true))
boxpfp = Box(lengths, (true,false,true))

v = Vec(-1,3,8)
v1 = Vec(0.5,1,1)
v2 = lengths+Vec(-0.5,-1,-1)

@testset "test dimensions/coordinates access functions" begin
    @test box.lengths == lengths
end

@testset "Test wrapping and unwrapping" begin
    @test wrap(v,box) == (Vec(2,3,3),(-1,0,1))
    @test wrap(v1,box) == (v1,(0,0,0))
    @test wrap(zero(Vec),box) == (zero(Vec),(0,0,0))
    @test wrap(lengths,box) == (zero(Vec),(1,1,1))
    for img in ((0,2,4), [0,2,4], Vec(0,2,4))
        @test unwrap(v1, img, box) == Vec(0.5, 9, 21)
        @test unwrap(v1, img, boxpfp) == Vec(0.5, 9, 21)
    end
end

@testset "Test separation" begin
    @test separation(v1,v2,box) == Vec(1,-2,2)
    #@pending separation(need,more,tests) --> Vec(x,x,x)
end
   


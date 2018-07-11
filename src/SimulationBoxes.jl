module SimulationBoxes

export SimulationBox, Box
export getorigin, getvectors, getdiagonal
export setorigin
export isperiodic
export wrap, wrap!, unwrap, separation
export centerofmass

"""
SimulationBox has three type parameters:
    * V: An immutable vector-like type with a Base.length(::V) method
        and a dimfields method returning the names of
        the fields associated with each dimension
    * N: The number of dimensions
    * P: A tuple of booleans indicating which dimensions are periodic

Concrete types inheriting from SimulationBox
must implement 3 functions: origin(::Box), vectors(::Box),
 si(::Box)
"""
abstract type SimulationBox{V,N,P} end

include("center_of_mass.jl")

@inline isperiodic(b::SimulationBox{V,N,P}) where {V,N,P} = P
@inline isperiodic(b::SimulationBox{V,N,P}, dim) where {V,N,P} = P[dim]

@inline Base.eltype(::SimulationBox{V}) where V = V

struct Box{V,N,P} <: SimulationBox{V,N,P}
    origin::V
    vectors::NTuple{N,V}
    diagonal::V
    function Box{V,N,P}(vectors::NTuple{N,V}, origin = zero(V)) where {N,V,P}
        if !isa(N,Integer)
            error("N parameter of a SimulationBox must be an Integer")
        end
        if !isa(P,NTuple{N,Bool})
            error(string("P=$P parameter of a `SimulationBox` must be a tuple ",
                         "of `Bool`s, with one `Bool` for each dimension"))
        end
        for vec in vectors
            if !(length(origin)==length(vec)==N)
                msg = "origin and box vector arguments do not have equal lengths"
                throw(DimensionMismatch(msg))
            end
        end
        for i in 1:N
            if vectors[i][i] <=0 
                error("Diagonal box matrix elements must be >= 0")
            end
            for j in 1:i-1
                if vectors[j][i] != 0 
                    error("box vector element [$(j)][$(i)] must be zero")
                end
            end
            for j in i+1:N
                if vectors[j][i] != 0 
                    warn("Triclinic box support not fully tested ($j,$i)")
                end
            end
        end
        diagonal = convert(V, collect(vectors[i][i] for i in 1:N) )
        new{V,N,P}(origin, vectors, diagonal)
    end
end

function Box(sides::V, origin::V=zero(V),
                  periodic=( (true for x in 1:length(sides))..., )) where V
    N = length(sides)
    vectors = map( (i,L)-> begin 
        v = zeros(eltype(V), N)
        v[i] = L
        convert(V,v)
    end, 1:N, sides)
    Box{V,N,(periodic...,)}((vectors...,), origin)
end

function setorigin(box::Box{V,N,P}, origin::V) where {V,N,P}
    Box{V,N,P}(getvectors(box), origin)
end

@inline getorigin(box::Box) = box.origin
@inline getvectors(box::Box) = box.vectors
@inline function getdiagonal(box::SimulationBox{V,N}) :: V where {V,N}
    v = vectors(box)
    collect( v[i][i] for i in 1:N )
end
@inline getdiagonal(box::Box) = box.diagonal

@inline _wrap(
        x::V,
        origin::V,
        vector::V,
        dim::Integer,
        p::Type{Val{false}},
        ) where V = x, 0

@inline function _wrap(
        x::V,
        origin::V,
        vector::V,
        dim::Integer,
        p::Type{Val{true}},
        ) ::Tuple{V,Int} where V
    image, x0 = divrem(x[dim]-origin[dim], vector[dim])
    image -= ifelse(x0<0, 1, 0)
    x - image*vector, convert(Int, image)
end

Base.convert(::Type{Vector{T}}, x::NTuple{N,T}) where {N,T} = collect(x)

@generated function wrap(x::V, box::SimulationBox{V,N,P}) where {V,N,P}
    lines = Expr[]
    images = [ Symbol("img$i") for i in 1:N ]
    for i in 1:N
        line = quote
            x, $(images[i]) = _wrap(x, origin, vectors[$i], $i, Val{$(P[i])})
        end
        push!(lines, line)
    end
    quote
        :(Expr(:meta, :inline))
        origin, vectors = getorigin(box), getvectors(box)
        $(lines...)
        x,( $(images...), )
    end
end

function wrap!(x::Vector{V}, box::SimulationBox{V}) where V
    for i in eachindex(x)
        x[i],_ = wrap(x[i], box)
    end
    x
end

function unwrap(x::V, image, box::SimulationBox{V,N}) :: V where {V,N}
    if length(image) != N
        throw(DimensionMismatch("image argument wrong length"))
    end
    if length(x) != N
        throw(DimensionMismatch("x argument wrong length"))
    end
    vectors = getvectors(box)
    x + sum( vectors[i] * image[i] for i in 1:N )
end

"""
Unwraps periodic boundaries given the coordinate in image 0.
Does not check for periodic boundaries, which should have an image flag of zero.
"""
unwrap

@inline _separation(x1::Real,
                   x2::Real,
                   origin::Real,
                   diagonal::Real,
                   periodic::Type{Val{false}}) = x2-x1

@inline function _separation(x1::Real,
                            x2::Real,
                            origin::Real,
                            diagonal::Real,
                            periodic::Type{Val{true}})
    r = x1-x2
    hdiag = 0.5diagonal
    r + ifelse(r>hdiag, -diagonal,
                 ifelse(r<-hdiag, diagonal, zero(diagonal)))
end

@generated function separation(x1::V, x2::V,
                               box::SimulationBox{V,N,P}) where {V,N,P}
    args = [:(_separation(x1[$i],
                          x2[$i],
                          origin[$i],
                          diagonal[$i],
                          Val{$(P[i])}))
            for i in 1:N ]
    meta = Expr(:meta, :inline)
    quote
        $meta
        origin, diagonal = getorigin(box), getdiagonal(box)
        convert(V, ( ($(args...),) ))
    end
end

"""
separation{V,N::Int,P::Tuple{Vararg{Bool}}}(v1::V, v2::V, box::SimulationBox{V,N,P})
      -> returns a vector of type V corresponding to v1-v2
         taking into account periodic boundaries using the
         nearest image convention, assuming v1 and v2 are in the central box.
"""
separation

"""
Returns wraps the given coordinate
into the central box
according to any periodic boundaries
"""
wrap

end #module

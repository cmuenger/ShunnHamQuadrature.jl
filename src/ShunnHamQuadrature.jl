module ShunnHamQuadrature

using LinearAlgebra
using StaticArrays

include("Quadrature/Triangle_ref.jl")
include("Quadrature/Triangle.jl")

include("Quadrature/Tetrahedron_ref.jl")
include("Quadrature/Tetrahedron.jl")

include("Quadrature/Pentatope_ref.jl")
include("Quadrature/Pentatope.jl")

include("Quadrature/Hexateron.jl")

export(QScheme)

"""
    QScheme
A quadrature scheme (points and weights). Schemes are created by calling
`shunnham*D` gives the points in bariycentric coordinates. With 'get_pts_wts' 
they can be converted in cartisian coordinates for a given simplex.
"""
struct QScheme{N,T}
    weights::Vector{T}
    points::Vector{SVector{N,T}}
end


@inbounds function calc_vol(vertices::SMatrix{D,N,T}) where {D,N,T}
    @assert N == D + 1
    X = SMatrix{D,D,T}(vertices[i, j + 1] - vertices[i, 1]
                       for i in 1:D, j in 1:D)
    vol = det(X)
    return vol
end

export get_pts_wts

function get_pts_wts(scheme::QScheme{N,T},  vertices::SMatrix{D,N,U}) where {N,T,D,U}
    @assert N > 0
    @assert N >= D + 1

    ws = scheme.weights
    ps = scheme.points
    @assert length(ws) == length(ps)

    xs  = Vector{SVector{D,T}}()
    for i in 1:length(ws)

        p = ps[i]
       
        push!(xs,vertices * p)
    
    end

    vol = calc_vol(vertices) / factorial(N - 1)

    return xs, ws*vol
end

function get_pts_wts(scheme::QScheme{N,T}, vertices::SVector{N,SVector{D,U}}) where {N,T,D,U}
    return get_pts_wts( scheme,SMatrix{D,N,U}(vertices[n][d] for d in 1:D, n in 1:N))
end

function get_pts_wts( scheme::QScheme{N}, vertices::AbstractVector) where {N}
    @assert length(vertices) == N
    @assert N > 0
    D = length(vertices[1])
    @assert N >= D + 1
    vertices′ = SVector{N}(map(SVector{D}, vertices))
    return get_pts_wts(scheme, vertices′)
end



export shunnham5D

export shunnham4D
export shunnham4D_ref 
export shunnham3D
export shunnham3D_ref 
export shunnham2D
export shunnham2D_ref 


"""
    shunnhamXD(n)
Returns the n-th quadrature rule on the X-dimensional simplex. Returns a Qudarature Scheme with
the quadrature points in barycentric coordiantes and the normalized quadrature weights.
"""

function shunnham5D(n)
    @assert 1 <= n <= length(hexateronAn)
    p = hexateronAn[n]
    w = hexateronWn[n]

    return QScheme(w, p)
end


function shunnham4D_ref(n)
    @assert 1 <= n <= length(pentatopeA)
    p = pentatopeA[n]
    w = pentatopeW[n]

    return QScheme(w, p)
end

function shunnham4D(n)
    @assert 1 <= n <= length(pentatopeAn)
    p = pentatopeAn[n]
    w = pentatopeWn[n]

    return QScheme(w, p)
end


function shunnham3D_ref(n)
    @assert 1 <= n <= length(tetrahedronA)
    p = tetrahedronA[n]
    w = tetrahedronW[n]

    return QScheme(w, p)
end

function shunnham3D(n)
    @assert 1 <= n <= length(tetrahedronAn)
    p = tetrahedronAn[n]
    w = tetrahedronWn[n]

    return QScheme(w, p)
end

function shunnham2D_ref(n)
    @assert 1 <= n <= length(triangleA)
    p = triangleA[n]
    w = triangleW[n]

    return QScheme(w, p)
end


function shunnham2D(n)
    @assert 1 <= n <= length(triangleAn)
    p = triangleAn[n]
    w = triangleWn[n]

    return QScheme(w, p)
end

end







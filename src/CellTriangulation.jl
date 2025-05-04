# A CellTriangulation triangulates a face on a reference cell. If non-triangular,
# we start by splitting it into triangles. Triangles are then recursively refined.
# TODO: This should rather be `FaceTriangulation`, 
# potentially contained inside `FaceTriangulations`. 
# Might also make sense then to have `FaceValues` in a similar fashion to `FacetValues`. 
"""
    CellTriangulation{RefShape, T}(refinement_levels::Int)

A `CellTriangulation` defines how a standard cell can be represented by a triangulation of linear
triangles, suitable for visualization. Essientially, each vertex of the triangulation of the cell 
consists of points, ξ, in the reference coordinate system of the cell. 
This can also be used as a QuadratureRule in a
CellValues object, in order to evaluate a function for each vertex of its triangles. However, the
integration weights will be NaN, and should not be used.
"""
CellTriangulation

struct CellTriangulation{RefShape, QR <: QuadratureRule{RefShape}}
    qr::QR
    triangles::Matrix{Int}  # (3, n), q_points in each of the n triangles
    edge::Vector{Int}       # q_point::Int on the edge
    function CellTriangulation{RefShape}(points::Vector{Vec{rdim, T}}, triangles::Matrix{Int}, edge::Vector{Int}) where {RefShape, T, rdim}
        weights = [convert(T, NaN) for _ in 1:length(points)]
        qr = QuadratureRule{RefShape}(weights, points)
        return new{RefShape, typeof(qr)}(qr, triangles, edge)
    end
end

function CellTriangulation{RefShape}(refinement_levels::Int) where {RefShape}
    points, triangles, edge = triangulate_cell(RefShape, refinement_levels)
    return CellTriangulation{RefShape}(points, triangles, edge)
end

Ferrite.getpoints(ct::CellTriangulation) = Ferrite.getpoints(ct.qr)
gettriangles(ct::CellTriangulation) = ct.triangles
getedge(ct::CellTriangulation) = ct.edge
Ferrite.CellValues(ct::CellTriangulation, args...) = Ferrite.CellValues(ct.qr, args...)

function triangulate_cell(::Type{RefQuadrilateral}, refinement_levels::Int)
    n = (refinement_levels + 1)         # Number of quad corner subnodes in each direction
    n_cells = 4 * refinement_levels^2   # Total number of subcells
    n_nodes = n^2 + refinement_levels^2 # Total number of subnodes
    T = Vec{2, Float64}
    points = Vector{Vec{2, Float64}}(undef, n_nodes)
    fill!(points, zero(Vec{2}) * NaN)
    coord_fun(index) = -1 + 2 * (index - 1) / refinement_levels
    for j in 1:n
        j_offset = n * (j - 1) + (n - 1) * (j - 1)
        for i in 1:n
            points[j_offset + i] = T((coord_fun(i), coord_fun(j)))
            if i < n && j < n
                points[j_offset + n + i] = T((coord_fun(i) + coord_fun(i + 1), coord_fun(j) + coord_fun(j + 1))) / 2
            end
        end
    end
    edge = Vector{Int}(undef, 4 * (n - 1) + 1)
    edge[1:n] .= 1:n
    edge[n:(2n-1)] .= n:(n + refinement_levels):n_nodes
    edge[(2n-1):(3n-2)] .= n_nodes:-1:(n_nodes-n+1)
    edge[(3n-2):(4n-4)] .= (n_nodes-n+1):(-(n + refinement_levels)):2 # 2 instead of 1 to skip last
    edge[end] = 1

    triangles = Matrix{Int}(undef, 3, n_cells)
    offset = 0
    for j in 1:refinement_levels
        for i in 1:refinement_levels
            bot = (2n - 1) * (j - 1) + i .+ (0, 1)  # Two nodes at bottom of quad
            top = bot .+ (2n - 1)                   # Two nodes at top of quad
            mid = bot[1] + n                        # Inserted middle node
            triangles[:, offset + 1] .= (bot[1], bot[2], mid)
            triangles[:, offset + 2] .= (bot[2], top[2], mid)
            triangles[:, offset + 3] .= (top[2], top[1], mid)
            triangles[:, offset + 4] .= (top[1], bot[1], mid)
            offset += 4
        end
    end
    return points, triangles, edge
end

function triangulate_cell(::Type{RefTriangle}, refinement_levels::Int)
    n_segments = 2^(refinement_levels - 1)
    vertices = Ferrite.reference_coordinates(Lagrange{RefTriangle, 1}())
    v_base = vertices[3] # Lower left corner
    r_base = vertices[1] # Lower right corner
    dx_left = (vertices[2] - v_base) / n_segments
    dx_horizontal = (r_base - v_base) / n_segments
    points = Vector{Vec{2, Float64}}(undef, (n_segments + 2) * (n_segments + 1) ÷ 2)
    fill!(points, NaN * zero(Vec{2})) # Temporary to check
    edge = Vector{Int}(undef, n_segments * 3 + 1)
    fill!(edge, -1)
    cnt = 1
    for j in 0:n_segments
        v0 = v_base + dx_left * j
        edge[3 * n_segments + 1 - j] = cnt
        for i in 0:(n_segments - j)
            points[cnt] = v0 + i * dx_horizontal
            j == 0 && (edge[i+1] = cnt)
            cnt += 1
        end
        edge[n_segments + 1 + j] = cnt - 1
    end
    n_triangles = (n_segments + 1) * n_segments ÷ 2 + n_segments * (n_segments - 1) ÷ 2
    triangles = Matrix{Int}(undef, 3, n_triangles)
    fill!(triangles, -1) # Temporary to check
    num_added = 0
    for j in 1:n_segments
        n_remaining = n_segments - j + 1
        i0_bot = length(points) - (n_remaining + 2) * (n_remaining + 1) ÷ 2 + 1
        i0_top = i0_bot + (n_segments - j + 2)
        for i in 1:(n_segments - j + 1)
            triangles[:, num_added + 2 * (i - 1) + 1] .= (i0_bot + i - 1, i0_bot + i, i0_top + i - 1)
        end
        for i in 1:(n_segments - j)
            triangles[:, num_added + 2 * (i - 1) + 2] .= (i0_bot + i, i0_top + i, i0_top + i - 1)
        end
        num_added += 2 * (n_segments - j) + 1
    end
    return points, triangles, edge
end

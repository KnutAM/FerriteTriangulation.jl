using FerriteTriangulation: Triangulation, SubTriangulation, get_data_indicies

function _create_data!(f, data, grid, a, cvs, subtria::SubTriangulation)
    c1 = first(subtria.faces)[1]
    x = copy(getcoordinates(grid, c1))
    dofs = copy(celldofs(dh, c1))
    ae = zeros(eltype(a), length(dofs))
    for (i, (cellnr, facenr)) in enumerate(subtria.faces)
        cv = cvs[facenr]
        getcoordinates!(x, grid, cellnr)
        reinit!(cv, getcells(grid, cellnr), x)
        celldofs!(dofs, dh, cellnr)
        copyto!(ae, view(a, dofs))
        for (i_point, i_data) in enumerate(get_data_indicies(subtria, i))
            data[i_data] = f(function_value(cv, i_point, ae))
        end
    end
end

"""
    create_data(tr::Triangulation, grid::AbstractGrid, a::Vector{<:Number}, ::NTuple{N, <:Interpolation}; 
        f = identity)

Create scalar data by evaluating `f(function_value(...))` at each triangulation node in the `grid`.
"""
function create_data(tr::Triangulation, grid, a, ips; f = identity)
    data = zeros(length(tr.nodes))
    for (ip, subtria) in zip(ips, tr.sub_triangulation)
        cvs = [CellValues(cr, ip, geometric_interpolation(getcelltype(subtria.sdh))) for cr in subtria.rules]
        _create_data!(f, data, grid, a, cvs, subtria)
    end
    return data
end
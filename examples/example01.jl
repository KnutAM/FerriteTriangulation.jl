using Ferrite, Tensors
using FerriteTriangulation: Triangulation

include(joinpath(@__DIR__, "example_utils.jl")) # create_data(tr, grid, a, ips)

CT = QuadraticTriangle #QuadraticQuadrilateral #Triangle #Quadrilateral
RS = Ferrite.getrefshape(CT)

ip = Lagrange{RS, 2}()
grid = generate_grid(CT, (2000, 1000))
transform_coordinates!(grid, x -> x + Vec((0, 1)) * cospi(x[1]/2)/2)
dh = close!(add!(DofHandler(grid), :u, ip))
a = zeros(ndofs(dh))
apply_analytical!(a, dh, :u, x -> x ⋅ x)

tr = Triangulation(dh, 3)

data = create_data(tr, grid, a, (ip,))

import GeometryBasics as GB
import CairoMakie as Plt

fig = Plt.Figure()
ax = Plt.Axis(fig[1,1]; aspect=Plt.DataAspect())

nodes = [GB.Point(x.data) for x in tr.nodes]
tri = [GB.GLTriangleFace(idx) for idx in eachcol(tr.triangles)];

m = Plt.mesh!(ax, nodes, tri; color=data)
Plt.Colorbar(fig[1,2], m)
#for i in 2:length(tr.tri_edges)
#    Plt.lines!(view(nodes, view(tr.edges, tr.tri_edges[i-1]:(tr.tri_edges[i]-1))); color=:black)
#end
#fig

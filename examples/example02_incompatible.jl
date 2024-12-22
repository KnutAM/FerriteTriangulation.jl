using Ferrite, Tensors
using FerriteTriangulation: Triangulation, SubTriangulation
include(joinpath(@__DIR__, "example_utils.jl")) # create_data(tr, grid, a, ips)

CT = Triangle
RS = Ferrite.getrefshape(CT)

ip1 = Lagrange{RS, 1}()
ip2 = Lagrange{RS, 2}()

grid = generate_grid(CT, (2, 2))
addcellset!(grid, "left", x -> x[1] < 1e-3)
addcellset!(grid, "right", setdiff(1:getncells(grid), getcellset(grid, "left")))
dh = DofHandler(grid)
sdh1 = SubDofHandler(dh, getcellset(grid, "left"))
add!(sdh1, :u, ip1)
sdh2 = SubDofHandler(dh, getcellset(grid, "right"))
add!(sdh2, :u, ip2)
dh = close!(dh)

a = zeros(ndofs(dh))
apply_analytical!(a, dh, :u, x -> x ⋅ x)

tr = Triangulation(dh, 4)
data = create_data(tr, grid, a, (ip1, ip2))

import GeometryBasics as GB
import CairoMakie as Plt

fig = Plt.Figure()
ax = Plt.Axis(fig[1,1]; aspect=Plt.DataAspect())

nodes = [GB.Point(x.data) for x in tr.nodes]
m = Plt.mesh!(ax, nodes, reshape(tr.triangles, :); color=data)
Plt.Colorbar(fig[1,2], m)
for i in 2:length(tr.tri_edges)
    Plt.lines!(view(nodes, view(tr.edges, tr.tri_edges[i-1]:(tr.tri_edges[i]-1))); color=:black)
end
fig

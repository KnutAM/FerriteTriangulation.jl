
using Ferrite, Tensors
using Gmsh, FerriteGmsh
using FerriteTriangulation: Triangulation

include(joinpath(@__DIR__, "example_utils.jl")) # create_data(tr, grid, a, ips)

function setup_grid(h=0.20)
    # Initialize gmsh
    Gmsh.initialize()
    gmsh.option.set_number("General.Verbosity", 2)

    # Add the points
    o =  gmsh.model.geo.add_point( 0.0,  0.0, 0.0, h)
    p1 = gmsh.model.geo.add_point( 1.0,  0.0, 0.0, h)
    p2 = gmsh.model.geo.add_point( 1.0,  1.0, 0.0, h)
    p3 = gmsh.model.geo.add_point(-1.0,  1.0, 0.0, h)
    p4 = gmsh.model.geo.add_point(-1.0, -1.0, 0.0, h)
    p5 = gmsh.model.geo.add_point( 0.0, -1.0, 0.0, h)

    pts = [o, p1, p2, p3, p4, p5, o]
    # Add the lines
    lines = [gmsh.model.geo.add_line(pts[i-1], pts[i]) for i in 2:length(pts)]
    
    # Create the closed curve loop and the surface
    loop = gmsh.model.geo.add_curve_loop(lines)
    gmsh.model.geo.add_plane_surface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)

    # Save the mesh, and read back in as a Ferrite Grid
    grid = mktempdir() do dir
        path = joinpath(dir, "mesh.msh")
        gmsh.write(path)
        togrid(path)
    end

    # Finalize the Gmsh library
    Gmsh.finalize()

    return grid
end

grid = setup_grid(0.05)
ip = DiscontinuousLagrange{RefTriangle, 1}()^2
dh = close!(add!(DofHandler(grid), :u, ip))

function foo(x::Vec{2}) # Analytical solution
    Δθ = -3π/4 # Rotate to 4th quadrant 
    xp = rotate(x, Δθ)
    r = sqrt(x ⋅ x + eps())
    θ = r ≤ 1e-6 ? zero(eltype(x)) : (atan(xp[2], xp[1]) + Δθ)
    return r^(2//3) * sin(2θ/3)
end

a = zeros(ndofs(dh))

apply_analytical!(a, dh, :u, x -> gradient(foo, x))

tr = Triangulation(dh, 2)
data = create_data(tr, grid, a, (ip,); f = first)

import GeometryBasics as GB
import CairoMakie as Plt

fig = Plt.Figure()
ax = Plt.Axis(fig[1,1]; aspect=Plt.DataAspect())

nodes = [GB.Point(x.data) for x in tr.nodes]
m = Plt.mesh!(ax, nodes, reshape(tr.triangles, :); color=data, 
    colormap=Plt.Makie.wong_colors(), 
    interpolate=false, 
    colorrange = (0, 1.5),
    )
Plt.Colorbar(fig[1,2], m)
#=
for i in 2:length(tr.tri_edges)
    Plt.lines!(view(nodes, view(tr.edges, tr.tri_edges[i-1]:(tr.tri_edges[i]-1))); color=:black)
end # =#
fig
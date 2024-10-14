using Oceananigans, OpenBoundaries, CairoMakie

using Oceananigans.BoundaryConditions: Open

using OpenBoundaries: PertubationAdvectionOpenBoundaryCondition

Re = 900

u = -0.1
l = 1

ν = abs(u) * (l/8) / Re

Lx = 2 * l
Ly = 2 * l

Δx = 1/16

Ny = floor(Int, Ly * 4 / Δx)
Nx = floor(Int, Lx / Ly * Ny)

grid = RectilinearGrid(topology = (Periodic, Bounded, Flat), size = (Nx, Ny), extent = (Lx, Ly))

v_south_boundary  = PertubationAdvectionOpenBoundaryCondition(u; relaxation_timescale = 0.1)#MeanOutflowOpenBoundaryCondition(grid, :west, u; inflow_relaxation_timescale = 10.0, outflow_relaxation_timescale = Inf)
v_north_boundary  = PertubationAdvectionOpenBoundaryCondition(u; relaxation_timescale = 0.1)#MeanOutflowOpenBoundaryCondition(grid, :east, u; inflow_relaxation_timescale = 10.0, outflow_relaxation_timescale = Inf)

u_south_boundary  = GradientBoundaryCondition(0)
u_north_boundary  = GradientBoundaryCondition(0)

u_bcs = FieldBoundaryConditions(south = u_south_boundary, north = u_north_boundary)

v_bcs = FieldBoundaryConditions(south = v_south_boundary, north = v_north_boundary)

closure = ScalarDiffusivity(; ν, κ = ν)

obstruction(x, y) = ifelse((x-l)^2 + (y-Ly/2)^2 < (l/8)^2, 1, 0)

Δt = 0.01#0.2 * min(minimum_xspacing(grid) / 0.1, minimum_xspacing(grid)^2 / ν)

rate = 1/(2*Δt)
u_forcing = Relaxation(; rate, mask = obstruction)
v_forcing = Relaxation(; rate, mask = obstruction)

model = NonhydrostaticModel(; grid, closure, 
                              advection = WENOFifthOrder(grid),
                              boundary_conditions = (; u = u_bcs, v = v_bcs),
                              forcing = (; u = u_forcing, v = v_forcing))

set!(model; v = (x, y) -> (1-obstruction(x, y)) * u, 
            u = (x, y) -> abs(u)/10 * randn() * (1-obstruction(x, y)))

model.velocities.v[:, 0, :] .= u

simulation = Simulation(model; Δt, stop_time = 5 * Lx / abs(u))

@info "Runninng simulation stopping after $(prettytime(simulation.stop_time))"

simulation.output_writers[:velocities] = JLD2OutputWriter(model, model.velocities,
                                                          filename = "2d_cylinder.jld2",
                                                          schedule = TimeInterval(.1),
                                                          overwrite_existing = true)

prog(sim) = @info "Completed $(prettytime(sim)) in $(prettytime(sim.run_wall_time)) with Δt = $(prettytime(simulation.Δt))"
                   #mean left is $(u_west_boundary.classification.matching_scheme.mean_outflow_velocity[]), 
                   #mean right is $(u_east_boundary.classification.matching_scheme.mean_outflow_velocity[]),
                   #should be $(u * cos(time(sim) * 2π / 10))"

add_callback!(simulation, prog, IterationInterval(100))

#conjure_time_step_wizard!(simulation; diffusive_cfl = 0.2)

run!(simulation)

u = FieldTimeSeries("2d_cylinder.jld2", "u")
v = FieldTimeSeries("2d_cylinder.jld2", "v")

n = Observable(1)

u_plt = @lift interior(u[$n], :, :, 1)
v_plt = @lift interior(v[$n], :, :, 1)

title = @lift prettytime(u.times[$n])

fig = Figure(size=(600, 400));

ax1 = Axis(fig[1, 1], title = "u", aspect = DataAspect())
ax2 = Axis(fig[2, 1], title = "v", aspect = DataAspect())

xf, yc, _ = nodes(u)
xc, yf, _ = nodes(v)

hm1 = heatmap!(ax1, xf, yc, u_plt, colorrange = (-0.1, 0.1), colormap = :roma)
hm2 = heatmap!(ax2, xc, yf, v_plt, colorrange = (-0.1, 0.1), colormap = :roma)

Colorbar(fig[1:2, 2], hm1)

#supertitle = Label(fig[0, :], title)

record(fig, "2d_cylinder_both.mp4", 1:length(u.times), framerate=20) do i
    n[] = i
end

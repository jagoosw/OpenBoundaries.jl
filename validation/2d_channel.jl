using Oceananigans, OpenBoundaries, CairoMakie

using OpenBoundaries: PerturbationAdvectionOpenBoundaryCondition

Re = 90

u = 0.1
l = 1

ν = u * l / Re

Lx = 4 * l
Ly = 1 * l

Ny = floor(Int, Ly * 4) * 16
Nx = floor(Int, Lx / Ly * Ny) * 16

grid = RectilinearGrid(topology = (Bounded, Bounded, Flat), size = (4Ny, Ny), extent = (Lx, Ly))

u_west_boundary  = OpenBoundaryCondition(u)#MeanOutflowOpenBoundaryCondition(grid, :west, u; inflow_relaxation_timescale = 10.0, outflow_relaxation_timescale = Inf)
u_east_boundary  = PerturbationAdvectionOpenBoundaryCondition(u; relaxation_timescale = Inf)#MeanOutflowOpenBoundaryCondition(grid, :east, u; inflow_relaxation_timescale = 10.0, outflow_relaxation_timescale = Inf)

u_north_boundary = ValueBoundaryCondition(0)
u_south_boundary = ValueBoundaryCondition(0)

v_west_boundary  = GradientBoundaryCondition(0)
v_east_boundary  = GradientBoundaryCondition(0)

u_bcs = FieldBoundaryConditions(west = u_west_boundary, east = u_east_boundary, 
                                north = u_north_boundary, south = u_south_boundary)

v_bcs = FieldBoundaryConditions(west = v_west_boundary, east = v_east_boundary)

closure = ScalarDiffusivity(; ν, κ = ν)

model = NonhydrostaticModel(; grid, closure, boundary_conditions = (; u = u_bcs, v = v_bcs))

set!(model; u, v = (args...) -> u/10 * randn())

Δt = 0.2 * min(minimum_xspacing(grid) / 0.1, minimum_xspacing(grid)^2 / ν)

simulation = Simulation(model; Δt, stop_time = 3 * Lx / u)

@info "Runninng simulation stopping after $(prettytime(simulation.stop_time))"

simulation.output_writers[:velocities] = JLD2OutputWriter(model, model.velocities,
                                                          filename = "2d_channel.jld2",
                                                          schedule = TimeInterval(1),
                                                          overwrite_existing = true)

prog(sim) = @info "Completed $(prettytime(sim)) in $(prettytime(sim.run_wall_time)) with Δt = $(prettytime(simulation.Δt))"

add_callback!(simulation, prog, IterationInterval(100))

conjure_time_step_wizard!(simulation; diffusive_cfl = 0.2)

run!(simulation)

u = FieldTimeSeries("2d_channel.jld2", "u")
v = FieldTimeSeries("2d_channel.jld2", "v")

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

record(fig, "2d_channel.mp4", 1:length(u.times), framerate=20) do i
    n[] = i
end
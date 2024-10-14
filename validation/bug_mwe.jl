using Oceananigans, Oceananigans.Units, OpenBoundaries

using OpenBoundaries: PerturbationAdvectionOpenBoundaryCondition

using Oceananigans: PartialCellBottom

using Oceananigans.Solvers: ConjugateGradientPoissonSolver, fft_poisson_solver

Lx = Ly = 256
Lz = 16

underlying_grid = RectilinearGrid(size=(64, 64, 16), extent=(Lx, Ly, Lz), topology = (Bounded, Bounded, Bounded))

bathymetry(x, y) = -16 * (1-0.3sin(π * y / Ly)) * x / Lx + 1

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bathymetry))

U = 0
V = 0.1

u_bcs = FieldBoundaryConditions(south = GradientBoundaryCondition(0), 
                                north = GradientBoundaryCondition(0), 
                                #east = OpenBoundaryCondition(U))
                                east = PerturbationAdvectionOpenBoundaryCondition(U, inflow_relaxation_timescale = 100.0,
                                                                                     outflow_relaxation_timescale = 1000.0))

#v_bcs = FieldBoundaryConditions(south = OpenBoundaryCondition(V), north = OpenBoundaryCondition(V), east = GradientBoundaryCondition(0))
v_bcs = FieldBoundaryConditions(south = PerturbationAdvectionOpenBoundaryCondition(V, inflow_relaxation_timescale = 100.0,
                                                                                     outflow_relaxation_timescale = Inf), 
                                north = PerturbationAdvectionOpenBoundaryCondition(V, inflow_relaxation_timescale = 100.0,
                                                                                     outflow_relaxation_timescale = Inf), 
                                east = GradientBoundaryCondition(0))

w_bcs = FieldBoundaryConditions(south = GradientBoundaryCondition(0), north = GradientBoundaryCondition(0), east = GradientBoundaryCondition(0))

pressure_solver = ConjugateGradientPoissonSolver(grid;
                                                 preconditioner = fft_poisson_solver(grid.underlying_grid), maxiter = 1000)

model = NonhydrostaticModel(; grid, 
                              boundary_conditions = (u = u_bcs, v = v_bcs, w = w_bcs), 
                              closure = ScalarDiffusivity(ν = 10^-3),
                              pressure_solver)

# this is a bit of a hard start, we should probably set it so that w = 0 and the flow follows stream lines around the bathymetry, which should produce a steady state with no closures
set!(model, u = U, v = V)

Δt = 0.2 * minimum_yspacing(grid) / abs(V)# / 10

simulation = Simulation(model; Δt, stop_time = Ly / V)

@info "Will stop after $(prettytime(simulation.stop_time))"

simulation.output_writers[:velocities] = JLD2OutputWriter(model, model.velocities; filename = "mwe.jld2", schedule = TimeInterval(10), overwrite_existing = true)

prog(sim) = @info "$(prettytime(sim)) in $(prettytime(sim.run_wall_time)) with Δt = $(prettytime(sim.Δt))"

add_callback!(simulation, prog, IterationInterval(10))

#conjure_time_step_wizard!(simulation, max_change = 2, min_change = 0)

run!(simulation)

results = FieldDataset("mwe.jld2")

times = results["u"].times

using CairoMakie

n = Observable(length(times))

abs_not_nan(x) = ifelse(isfinite(x), abs(x), 0)

u_plt = @lift interior(results["u"][$n], :, :, grid.Nz)
v_plt = @lift interior(results["v"][$n], :, :, grid.Nz)
w_plt = @lift interior(results["w"][$n], :, 16, :)

title = @lift prettytime(times[$n])

fig = Figure()

ax1 = Axis(fig[1:2, 1], aspect = DataAspect(), title = "u")
ax2 = Axis(fig[1:2, 2], aspect = DataAspect(), title = "v")
ax3 = Axis(fig[3, 1:2], aspect = DataAspect(), title = "w")

hm1 = heatmap!(ax1, nodes(results["u"])[1:2]..., u_plt, colorrange = maximum(abs_not_nan, results["u"][:, :, grid.Nz, :]) .* (-1, 1) ./ 2)

hm2 = heatmap!(ax2, nodes(results["v"])[1:2]..., v_plt, colorrange = maximum(abs_not_nan, results["v"][:, :, grid.Nz, :]) .* (-1, 1) ./ 2)

hm3 = heatmap!(ax3, xnodes(results["w"]), znodes(results["w"]), w_plt, colorrange = maximum(abs_not_nan, results["w"][:, 16, :, :]) .* (-1, 1))

supertitle = Label(fig[0, :], title)

fig


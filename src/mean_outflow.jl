using Oceananigans.AbstractOperations: Axᶜᶜᶠ, Ayᶜᶠᶜ, Azᶜᶜᶠ, xspacing, yspacing, zspacing

"""
    MeanOutflow

`MeanOutflow` matching scheme which advects out of the boundary with the mean
boundary normal velocity. This works modaratly well in situations with generally
unidirectional flow and permits eddies to exit the domain fairly unpeturbed.

The flow can also be relaxed to the external value, and if the mean flow is into 
the domain then the flow is simply relaxed with a different time scale to the 
external prescribed values (although this is not very stable).

    ∂u/∂t + U∂u/∂x = {(uₑ - u)/τᵢ for u<0, (uₑ - u)/τ₀ for u>0}

"""
struct MeanOutflow{BF, OV, IT}
                 boundary_flux :: BF
         mean_outflow_velocity :: OV
   inflow_relaxation_timescale :: IT
  outflow_relaxation_timescale :: IT
end

const C = Center()

Adapt.adapt_structure(to, mo::MeanOutflow) = 
    MeanOutflow(nothing, 
                adapt(to, mo.mean_outflow_velocity[]), 
                adapt(to, mo.inflow_relaxation_timescale), 
                adapt(to, mo.outflow_relaxation_timescale))

function MeanOutflowOpenBoundaryCondition(grid, side, val; 
                                          inflow_relaxation_timescale = 1.0, 
                                          outflow_relaxation_timescale = Inf, 
                                          kwargs...)

    boundary_flux = boundary_flux_field(grid, Val(side))

    classifcation = Open(MeanOutflow(boundary_flux, Ref(0.), inflow_relaxation_timescale, outflow_relaxation_timescale))
    
    return BoundaryCondition(classifcation, val; kwargs...)
end

# fields is not yet defiled so have to use an array - we also only need the interior though
boundary_flux_field(grid, ::Union{Val{:west}, Val{:east}}) = on_architecture(architecture(grid), zeros(size(grid, 2), size(grid, 3)))
boundary_flux_field(grid, ::Union{Val{:south}, Val{:north}}) = on_architecture(architecture(grid), zeros(size(grid, 1), size(grid, 3)))
boundary_flux_field(grid, ::Union{Val{:bottom}, Val{:top}}) = on_architecture(architecture(grid), zeros(size(grid, 2), size(grid, 3)))

const MOOBC = BoundaryCondition{<:Open{<:MeanOutflow}}

boundary_normal_velocity(velocities, ::Union{Val{:west}, Val{:east}}) = velocities.u
boundary_normal_velocity(velocities, ::Union{Val{:south}, Val{:north}}) = velocities.v
boundary_normal_velocity(velocities, ::Union{Val{:bottom}, Val{:top}}) = velocities.w

# compute the boundary flux
function update_boundary_condition!(bc::MOOBC, side, field, model)
    ms = bc.classification.matching_scheme

    u = boundary_normal_velocity(model.velocities, side)
    F = ms.boundary_flux

    arch = architecture(model)
    grid = model.grid

    launch!(arch, grid, :yz, _update_boundary_flux!, F, grid, u, side)

    ms.mean_outflow_velocity[] = sum(F) / (grid.Ly * grid.Lz)
end

@kernel function _update_boundary_flux!(F, grid, u, ::Val{:west})
    j, k = @index(Global, NTuple)

    @inbounds F[j, k] = -u[1, j, k] * Axᶜᶜᶠ(1, j, k, grid)
end

@kernel function _update_boundary_flux!(F, grid, u, ::Val{:east})
    j, k = @index(Global, NTuple)

    i = grid.Nx

    @inbounds F[j, k] = u[i, j, k] * Axᶜᶜᶠ(i, j, k, grid)
end

@kernel function _update_boundary_flux!(F, grid, u, ::Val{:south})
    i, k = @index(Global, NTuple)

    @inbounds F[i, k] = -u[i, 1, k] * Ayᶜᶠᶜ(i, 1, k, grid)
end

@kernel function _update_boundary_flux!(F, grid, u, ::Val{:north})
    i, k = @index(Global, NTuple)

    j = grid.Ny

    @inbounds F[i, k] = u[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid)
end

@kernel function _update_boundary_flux!(F, grid, u, ::Val{:bottom})
    i, j = @index(Global, NTuple)

    @inbounds F[i, j] = -u[i, j, 1] * Azᶜᶜᶠ(i, j, 1, grid)
end

@kernel function _update_boundary_flux!(F, grid, u, ::Val{:top})
    i, j = @index(Global, NTuple)

    k = grid.Nz

    @inbounds F[i, j] = u[i, j, k] * Azᶜᶜᶠ(i, j, k, grid)
end

# fill the boundaries

@inline function _fill_west_open_halo!(j, k, grid, ϕ, bc::MOOBC, loc, clock, model_fields)
    ms = bc.classification.matching_scheme

    τᵢ = ms.inflow_relaxation_timescale
    τₒ = ms.outflow_relaxation_timescale

    uⁿ⁻¹ = @inbounds ϕ[1, j, k]
    uⁿ₂  = @inbounds ϕ[2, j, k]
    uₑ   = getbc(bc, j, k, grid, clock, model_fields)

    Δx = xspacing(1, j, k, grid, C, C, C)
    Δt = clock.last_stage_Δt
    
    Ũ = ms.mean_outflow_velocity[] * Δt / Δx

    τ̃ = ifelse(Ũ > 0, τₒ / Δt, τᵢ / Δt)

    # on the first step Δt is Inf so we need to get rid of that

    Ũ = ifelse(isinf(Δt), 0, Ũ)
    τ̃ = ifelse(isinf(Δt), 1, τ̃)

    # also normalise them
    Ũ = min(1, max(0, Ũ))
    τ̃ = max(1, τ̃)

    @inbounds ϕ[1, j, k] = (uⁿ⁻¹ + uⁿ₂ * Ũ + uₑ / τ̃) / (1 + Ũ + 1/τ̃)

    return nothing
end

@inline function _fill_east_open_halo!(j, k, grid, ϕ, bc::MOOBC, loc, clock, model_fields)
    ms = bc.classification.matching_scheme

    τᵢ = ms.inflow_relaxation_timescale
    τₒ = ms.outflow_relaxation_timescale

    i = grid.Nx + 1

    uⁿ⁻¹ = @inbounds ϕ[i, j, k]
    uⁿ₂  = @inbounds ϕ[i-1, j, k]
    uₑ   = getbc(bc, j, k, grid, clock, model_fields)

    Δx = xspacing(i, j, k, grid, C, C, C)
    Δt = clock.last_stage_Δt
    
    Ũ = ms.mean_outflow_velocity[] * Δt / Δx

    τ̃ = ifelse(Ũ > 0, τₒ / Δt, τᵢ / Δt)

    # on the first step Δt is Inf so we need to get rid of that

    Ũ = ifelse(isinf(Δt), 0, Ũ)
    τ̃ = ifelse(isinf(Δt), 1, τ̃)

    # also normalise them
    Ũ = min(1, max(0, Ũ))
    τ̃ = max(1, τ̃)

    @inbounds ϕ[i, j, k] = (uⁿ⁻¹ + uⁿ₂ * Ũ + uₑ / τ̃) / (1 + Ũ + 1/τ̃)

    return nothing
end
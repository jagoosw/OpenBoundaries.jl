using Oceananigans.Grids: xspacing

"""
    PretubationConvection

Assumed that the condition is the mean flow.
"""
struct PertubationConvection{FT}
    relaxation_timescale :: FT
end

function PertubationConvectionOpenBoundaryCondition(val; relaxation_timescale = 10, kwargs...)

    classifcation = Open(PertubationConvection(relaxation_timescale))
    
    return BoundaryCondition(classifcation, val; kwargs...)
end

const PCOBC = BoundaryCondition{<:Open{<:PertubationConvection}}

@inline function _fill_east_open_halo!(j, k, grid, u, bc::PCOBC, loc, clock, model_fields)
    i = grid.Nx + 1

    τ = bc.classification.matching_scheme.relaxation_timescale
    Δt = clock.last_stage_Δt

    Δt = ifelse(isinf(Δt), 0, Δt)

    Δx = xspacing(i, j, k, grid, Face(), Center(), Center())

    ūⁿ⁺¹ = getbc(bc, j, k, grid, clock, model_fields)

    clock.time -= Δt

    ūⁿ   = getbc(bc, j, k, grid, clock, model_fields)

    clock.time += Δt

    u′ᵢⁿ     = @inbounds u[i, j, k] - ūⁿ
    u′ᵢ₋₁ⁿ⁺¹ = @inbounds u[i - 1, j, k] - ūⁿ⁺¹

    u′ᵢⁿ⁺¹ = (u′ᵢⁿ + Δt / Δx * ūⁿ⁺¹ * u′ᵢ₋₁ⁿ⁺¹) / (1 + Δt / τ + Δt / Δx * ūⁿ⁺¹)

    @inbounds u[i, j, k] = ūⁿ⁺¹ + u′ᵢⁿ⁺¹
end

@inline function _fill_west_open_halo!(j, k, grid, u, bc::PCOBC, loc, clock, model_fields)
    τ = bc.classification.matching_scheme.relaxation_timescale
    Δt = clock.last_stage_Δt

    Δt = ifelse(isinf(Δt), 0, Δt)

    Δx = xspacing(0, j, k, grid, Face(), Center(), Center())

    ūⁿ⁺¹ = getbc(bc, j, k, grid, clock, model_fields)

    clock.time -= Δt # this might cause problems from asyncronisity ...

    ūⁿ   = getbc(bc, j, k, grid, clock, model_fields)

    clock.time += Δt

    u′₀ⁿ   = @inbounds u[0, j, k] - ūⁿ
    u′₁ⁿ⁺¹ = @inbounds u[1, j, k] - ūⁿ⁺¹

    u′₀ⁿ⁺¹ = (u′₀ⁿ - Δt / Δx * ūⁿ⁺¹ * u′₁ⁿ⁺¹) / (1 + Δt / τ - Δt / Δx * ūⁿ⁺¹)

    @inbounds u[0, j, k] = ūⁿ⁺¹ + u′₀ⁿ⁺¹
end

@inline function _fill_north_open_halo!(i, k, grid, v, bc::PCOBC, loc, clock, model_fields)
    j = grid.Ny + 1

    τ = bc.classification.matching_scheme.relaxation_timescale
    Δt = clock.last_stage_Δt

    Δt = ifelse(isinf(Δt), 0, Δt)

    Δx = xspacing(i, j, k, grid, Face(), Center(), Center())

    v̄ⁿ⁺¹ = getbc(bc, j, k, grid, clock, model_fields)

    clock.time -= Δt

    v̄ⁿ   = getbc(bc, j, k, grid, clock, model_fields)

    clock.time += Δt

    v′ⱼⁿ     = @inbounds v[i, j, k] - v̄ⁿ
    v′ⱼ₋₁ⁿ⁺¹ = @inbounds v[i, j - 1, k] - v̄ⁿ⁺¹

    v′ⱼⁿ⁺¹ = (v′ⱼⁿ + Δt / Δx * v̄ⁿ⁺¹ * v′ⱼ₋₁ⁿ⁺¹) / (1 + Δt / τ + Δt / Δx * v̄ⁿ⁺¹)

    @inbounds v[i, j, k] = v̄ⁿ⁺¹ + v′ⱼⁿ⁺¹
end

@inline function _fill_south_open_halo!(i, k, grid, v, bc::PCOBC, loc, clock, model_fields)
    τ = bc.classification.matching_scheme.relaxation_timescale
    Δt = clock.last_stage_Δt

    Δt = ifelse(isinf(Δt), 0, Δt)

    Δx = xspacing(i, 0, k, grid, Face(), Center(), Center())

    v̄ⁿ⁺¹ = getbc(bc, i, k, grid, clock, model_fields)

    clock.time -= Δt

    v̄ⁿ   = getbc(bc, i, k, grid, clock, model_fields)

    clock.time += Δt

    v′₀ⁿ   = @inbounds v[i, 0, k] - v̄ⁿ
    v′₁ⁿ⁺¹ = @inbounds v[i, 1, k] - v̄ⁿ⁺¹

    v′₀ⁿ⁺¹ = (v′₀ⁿ - Δt / Δx * v̄ⁿ⁺¹ * v′₁ⁿ⁺¹) / (1 + Δt / τ - Δt / Δx * v̄ⁿ⁺¹)

    @inbounds v[i, 0, k] = v̄ⁿ⁺¹ + v′₀ⁿ⁺¹
end
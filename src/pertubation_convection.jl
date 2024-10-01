using Oceananigans.Grids: xspacing
using Oceananigans.TimeSteppers: Clock

"""
    PretubationConvection

Assumed that the condition is the mean flow.
"""
struct PertubationConvection{FT, C}
    relaxation_timescale :: FT
              last_clock :: C
end

function PertubationConvectionOpenBoundaryCondition(val, FT = Float64; relaxation_timescale = 10, kwargs...)
    last_clock = Clock(; time = zero(FT))

    classifcation = Open(PertubationConvection(relaxation_timescale, last_clock))
    
    return BoundaryCondition(classifcation, val; kwargs...)
end

const PCOBC = BoundaryCondition{<:Open{<:PertubationConvection}}

@inline update_boundary_condition!(bc::PCOBC, side, field, model) = 
    (bc.classification.matching_scheme.last_clock.time = model.clock.time - model.clock.last_stage_Δt)

@inline function _fill_east_open_halo!(j, k, grid, u, bc::PCOBC, loc, clock, model_fields)
    i = grid.Nx + 1

    τ = bc.classification.matching_scheme.relaxation_timescale
    Δt = clock.last_stage_Δt
    tⁿ = bc.classifcation.matching_scheme.last_clock

    Δt = ifelse(isinf(Δt), 0, Δt)

    Δx = xspacing(i, j, k, grid, Face(), Center(), Center())

    ūⁿ⁺¹ = getbc(bc, j, k, grid, clock, model_fields)

    ūⁿ   = getbc(bc, j, k, grid, tⁿ, model_fields)

    u′ᵢⁿ     = @inbounds u[i, j, k] - ūⁿ
    u′ᵢ₋₁ⁿ⁺¹ = @inbounds u[i - 1, j, k] - ūⁿ⁺¹

    U = max(0, min(1, Δt / Δx * ūⁿ⁺¹))

    u′ᵢⁿ⁺¹ = (u′ᵢⁿ + U * u′ᵢ₋₁ⁿ⁺¹) / (1 + Δt / τ + U)

    @inbounds u[i, j, k] = ūⁿ⁺¹ + u′ᵢⁿ⁺¹
end

@inline function _fill_west_open_halo!(j, k, grid, u, bc::PCOBC, loc, clock, model_fields)
    τ = bc.classification.matching_scheme.relaxation_timescale
    Δt = clock.last_stage_Δt
    tⁿ = bc.classifcation.matching_scheme.last_clock

    Δt = ifelse(isinf(Δt), 0, Δt)

    Δx = xspacing(0, j, k, grid, Face(), Center(), Center())

    ūⁿ⁺¹ = getbc(bc, j, k, grid, clock, model_fields)

    ūⁿ   = getbc(bc, j, k, grid, tⁿ, model_fields)

    u′₀ⁿ   = @inbounds u[0, j, k] - ūⁿ
    u′₁ⁿ⁺¹ = @inbounds u[1, j, k] - ūⁿ⁺¹

    U = min(0, max(1, Δt / Δx * ūⁿ⁺¹))

    u′₀ⁿ⁺¹ = (u′₀ⁿ - U * u′₁ⁿ⁺¹) / (1 + Δt / τ - U)

    @inbounds u[0, j, k] = ūⁿ⁺¹ + u′₀ⁿ⁺¹
end

@inline function _fill_north_open_halo!(i, k, grid, v, bc::PCOBC, loc, clock, model_fields)
    j = grid.Ny + 1

    τ = bc.classification.matching_scheme.relaxation_timescale
    Δt = clock.last_stage_Δt
    tⁿ = bc.classifcation.matching_scheme.last_clock

    Δt = ifelse(isinf(Δt), 0, Δt)

    Δx = xspacing(i, j, k, grid, Face(), Center(), Center())

    v̄ⁿ⁺¹ = getbc(bc, j, k, grid, clock, model_fields)

    v̄ⁿ   = getbc(bc, j, k, grid, tⁿ, model_fields)

    v′ⱼⁿ     = @inbounds v[i, j, k] - v̄ⁿ
    v′ⱼ₋₁ⁿ⁺¹ = @inbounds v[i, j - 1, k] - v̄ⁿ⁺¹

    V = max(0, min(1, Δt / Δx * v̄ⁿ⁺¹))

    v′ⱼⁿ⁺¹ = (v′ⱼⁿ + V * v′ⱼ₋₁ⁿ⁺¹) / (1 + Δt / τ + V)

    @inbounds v[i, j, k] = v̄ⁿ⁺¹ + v′ⱼⁿ⁺¹
end

@inline function _fill_south_open_halo!(i, k, grid, v, bc::PCOBC, loc, clock, model_fields)
    τ = bc.classification.matching_scheme.relaxation_timescale
    Δt = clock.last_stage_Δt
    tⁿ = bc.classifcation.matching_scheme.last_clock

    Δt = ifelse(isinf(Δt), 0, Δt)

    Δx = xspacing(i, 0, k, grid, Face(), Center(), Center())

    v̄ⁿ⁺¹ = getbc(bc, i, k, grid, clock, model_fields)

    v̄ⁿ   = getbc(bc, i, k, grid, tⁿ, model_fields)

    v′₀ⁿ   = @inbounds v[i, 0, k] - v̄ⁿ
    v′₁ⁿ⁺¹ = @inbounds v[i, 1, k] - v̄ⁿ⁺¹

    V = min(0, max(1, Δt / Δx * v̄ⁿ⁺¹))

    v′₀ⁿ⁺¹ = (v′₀ⁿ - V * v′₁ⁿ⁺¹) / (1 + Δt / τ - V)

    @inbounds v[i, 0, k] = v̄ⁿ⁺¹ + v′₀ⁿ⁺¹
end
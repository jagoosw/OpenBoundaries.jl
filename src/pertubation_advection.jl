using Oceananigans.Grids: xspacing
using Oceananigans.TimeSteppers: Clock

"""
    PretubationConvection

Assumed that the condition is the mean flow.
"""
struct PertubationAdvection{FT, C}
 outflow_relaxation_timescale :: FT
  inflow_relaxation_timescale :: FT
                   last_clock :: C
end

function PertubationAdvectionOpenBoundaryCondition(val, FT = Float64; 
                                                   outflow_relaxation_timescale = Inf, 
                                                   inflow_relaxation_timescale = 10.0, kwargs...)
    last_clock = Clock(; time = zero(FT))

    classification = Open(PertubationAdvection(outflow_relaxation_timescale, inflow_relaxation_timescale, last_clock))
    
    return BoundaryCondition(classification, val; kwargs...)
end

const PAOBC = BoundaryCondition{<:Open{<:PertubationAdvection}}

@inline function update_boundary_condition!(bc::PAOBC, side, field, model)
    t = model.clock.time
    Δt = model.clock.last_stage_Δt

    Δt = ifelse(isinf(Δt), 0, Δt)
    
    bc.classification.matching_scheme.last_clock.time = t - Δt
    
    return nothing
end

@inline function _fill_east_open_halo!(j, k, grid, u, bc::PAOBC, loc, clock, model_fields)
    i = grid.Nx + 1

    Δt = clock.last_stage_Δt
    tⁿ = bc.classification.matching_scheme.last_clock

    Δt = ifelse(isinf(Δt), 0, Δt)

    Δx = xspacing(i, j, k, grid, Face(), Center(), Center())

    ūⁿ⁺¹ = getbc(bc, j, k, grid, clock, model_fields)

    ūⁿ   = getbc(bc, j, k, grid, tⁿ, model_fields)

    u′ᵢⁿ     = @inbounds u[i, j, k] - ūⁿ
    u′ᵢ₋₁ⁿ⁺¹ = @inbounds u[i - 1, j, k] - ūⁿ⁺¹

    U = max(0, min(1, Δt / Δx * ūⁿ⁺¹))

    τ = ifelse(u′ᵢ₋₁ⁿ⁺¹ + ūⁿ⁺¹ > 0, 
               bc.classification.matching_scheme.outflow_relaxation_timescale, 
               bc.classification.matching_scheme.inflow_relaxation_timescale)

    u′ᵢⁿ⁺¹ = (u′ᵢⁿ + U * u′ᵢ₋₁ⁿ⁺¹) / (1 + Δt / τ + U)

    @inbounds u[i, j, k] = ūⁿ⁺¹ + u′ᵢⁿ⁺¹
end

@inline function _fill_west_open_halo!(j, k, grid, u, bc::PAOBC, loc, clock, model_fields)
    Δt = clock.last_stage_Δt
    tⁿ = bc.classification.matching_scheme.last_clock

    Δt = ifelse(isinf(Δt), 0, Δt)

    Δx = xspacing(1, j, k, grid, Face(), Center(), Center())

    ūⁿ⁺¹ = getbc(bc, j, k, grid, clock, model_fields)

    ūⁿ   = getbc(bc, j, k, grid, tⁿ, model_fields)

    u′₀ⁿ   = @inbounds u[0, j, k] - ūⁿ
    u′₁ⁿ⁺¹ = @inbounds u[2, j, k] - ūⁿ⁺¹

    U = min(0, max(1, Δt / Δx * ūⁿ⁺¹))

    τ = ifelse(u′ᵢ₋₁ⁿ⁺¹ + ūⁿ⁺¹ < 0, 
               bc.classification.matching_scheme.outflow_relaxation_timescale, 
               bc.classification.matching_scheme.inflow_relaxation_timescale)

    u′₀ⁿ⁺¹ = @show (u′₀ⁿ - U * u′₁ⁿ⁺¹) / (1 + Δt / τ - U)

    # this is a temporaty hack because the 1, j, k point is getting stepped during 
    # the timestepping (based on erronious gradients) so we can't just integrate
    # it here
    @inbounds u[1, j, k] = ūⁿ⁺¹ + u′₀ⁿ⁺¹
    @inbounds u[0, j, k] = ūⁿ⁺¹ + u′₀ⁿ⁺¹
end

@inline function _fill_north_open_halo!(i, k, grid, v, bc::PAOBC, loc, clock, model_fields)
    j = grid.Ny + 1

    τ = bc.classification.matching_scheme.relaxation_timescale
    Δt = clock.last_stage_Δt
    tⁿ = bc.classification.matching_scheme.last_clock

    Δt = ifelse(isinf(Δt), 0, Δt)

    Δx = xspacing(i, j, k, grid, Center(), Face(), Center())

    v̄ⁿ⁺¹ = getbc(bc, i, k, grid, clock, model_fields)

    v̄ⁿ   = getbc(bc, i, k, grid, tⁿ, model_fields)

    v′ⱼⁿ     = @inbounds v[i, j, k] - v̄ⁿ
    v′ⱼ₋₁ⁿ⁺¹ = @inbounds v[i, j - 1, k] - v̄ⁿ⁺¹

    V = max(0, min(1, Δt / Δx * v̄ⁿ⁺¹))

    τ = ifelse(v′ⱼ₋₁ⁿ⁺¹ + v̄ⁿ⁺¹ > 0, 
               bc.classification.matching_scheme.outflow_relaxation_timescale, 
               bc.classification.matching_scheme.inflow_relaxation_timescale)
    
    v′ⱼⁿ⁺¹ = (v′ⱼⁿ + V * v′ⱼ₋₁ⁿ⁺¹) / (1 + Δt / τ + V)

    @inbounds v[i, j, k] = v̄ⁿ⁺¹ + v′ⱼⁿ⁺¹
end

@inline function _fill_south_open_halo!(i, k, grid, v, bc::PAOBC, loc, clock, model_fields)
    τ = bc.classification.matching_scheme.relaxation_timescale
    Δt = clock.last_stage_Δt
    tⁿ = bc.classification.matching_scheme.last_clock

    Δt = ifelse(isinf(Δt), 0, Δt)

    Δx = xspacing(i, 1, k, grid, Center(), Face(), Center())

    v̄ⁿ⁺¹ = getbc(bc, i, k, grid, clock, model_fields)

    v̄ⁿ   = getbc(bc, i, k, grid, tⁿ, model_fields)

    v′₀ⁿ   = @inbounds v[i, 0, k] - v̄ⁿ
    v′₁ⁿ⁺¹ = @inbounds v[i, 2, k] - v̄ⁿ⁺¹

    V = min(0, max(1, Δt / Δx * v̄ⁿ⁺¹))

    τ = ifelse(v′ⱼ₋₁ⁿ⁺¹ + v̄ⁿ⁺¹ < 0, 
               bc.classification.matching_scheme.outflow_relaxation_timescale, 
               bc.classification.matching_scheme.inflow_relaxation_timescale)

    v′₀ⁿ⁺¹ = (v′₀ⁿ - V * v′₁ⁿ⁺¹) / (1 + Δt / τ - V)

    # see note above
    @inbounds v[i, 1, k] = v̄ⁿ⁺¹ + v′₀ⁿ⁺¹
    @inbounds v[i, 0, k] = v̄ⁿ⁺¹ + v′₀ⁿ⁺¹
end
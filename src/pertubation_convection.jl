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
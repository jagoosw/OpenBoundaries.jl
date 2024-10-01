struct FixedPhaseSpeed{C, T}
             phase_speed :: C
    relaxation_timescale :: T
end

const FCOBC = BoundaryCondition{<:Open{<:FixedPhaseSpeed}}

@inline function _fill_west_open_halo!(j, k, grid, ϕ, bc::FCOBC, loc, clock, model_fields)
    ms = bc.classification.matching_scheme

    τ = ms.relaxation_timescale

    uⁿ⁻¹ = @inbounds ϕ[1, j, k]
    uⁿ₂  = @inbounds ϕ[2, j, k]
    uₑ   = getbc(bc, j, k, grid, clock, model_fields)

    Δx = xspacing(1, j, k, grid, C, C, C)
    Δt = clock.last_stage_Δt
    
    Ũ = ms.phase_speed * Δt / Δx
    τ̃ = τ / Δt

    Ũ = ifelse(isinf(Δt), 0, Ũ)
    τ̃ = ifelse(isinf(Δt), 1, τ̃)

    Ũ = max(-1, min(0, Ũ))
    τ̃ = max(1, τ̃)

    if Ũ < 0
        @inbounds ϕ[1, j, k] = (uⁿ⁻¹ - uⁿ₂ * Ũ) / (1 - Ũ)
    else
        @inbounds ϕ[1, j, k] = uⁿ⁻¹ + (uₑ - uⁿ⁻¹) / τ̃
    end
    
    return nothing
end

@inline function _fill_east_open_halo!(j, k, grid, ϕ, bc::FCOBC, loc, clock, model_fields)
    ms = bc.classification.matching_scheme

    τ = ms.relaxation_timescale

    i = grid.Nx + 1

    uⁿ⁻¹ = @inbounds ϕ[i, j, k]
    uⁿ₋₁  = @inbounds ϕ[i-1, j, k]
    uₑ   = getbc(bc, j, k, grid, clock, model_fields)

    Δx = xspacing(i, j, k, grid, C, C, C)
    Δt = clock.last_stage_Δt
    
    Ũ = ms.phase_speed * Δt / Δx
    τ̃ = τ / Δt

    Ũ = ifelse(isinf(Δt), 0, Ũ)
    τ̃ = ifelse(isinf(Δt), 1, τ̃)

    Ũ = min(1, max(0, Ũ))
    τ̃ = max(1, τ̃)

    if Ũ > 0
        @inbounds ϕ[i, j, k] = (uⁿ⁻¹ + uⁿ₋₁ * Ũ) / (1 + Ũ)
    else
        @inbounds ϕ[i, j, k] = uⁿ⁻¹ + (uₑ - uⁿ⁻¹) / τ̃
    end

    return nothing
end
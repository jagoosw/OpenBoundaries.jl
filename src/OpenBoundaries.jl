module OpenBoundaries

export MeanOutflowOpenBoundaryCondition, FixedPhaseSpeed

using Adapt, Oceananigans

using KernelAbstractions: @kernel, @index

using Oceananigans.Architectures: on_architecture, architecture
using Oceananigans.BoundaryConditions: Open, getbc
using Oceananigans.Utils: launch!

import Adapt: adapt_structure

import Oceananigans.BoundaryConditions: _fill_west_halo!,
                                        _fill_east_halo!,
                                        _fill_south_halo!,
                                        _fill_north_halo!,
                                        _fill_bottom_halo!,
                                        _fill_top_halo!,
                                        update_boundary_condition!

#include("mean_outflow.jl")
#include("fixed_phase_speed.jl")
include("perturbation_advection.jl")

end # module

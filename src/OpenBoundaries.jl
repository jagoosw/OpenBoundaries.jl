module OpenBoundaries

export MeanOutflowOpenBoundaryCondition, FixedPhaseSpeed

using Adapt, Oceananigans

using KernelAbstractions: @kernel, @index

using Oceananigans.Architectures: on_architecture, architecture
using Oceananigans.BoundaryConditions: Open, getbc
using Oceananigans.Utils: launch!

import Adapt: adapt_structure

import Oceananigans.BoundaryConditions: _fill_west_open_halo!,
                                        _fill_east_open_halo!,
                                        _fill_south_open_halo!,
                                        _fill_north_open_halo!,
                                        _fill_bottom_open_halo!,
                                        _fill_top_open_halo!,
                                        update_boundary_condition!

#include("mean_outflow.jl")
#include("fixed_phase_speed.jl")
include("pertubation_advection.jl")

end # module

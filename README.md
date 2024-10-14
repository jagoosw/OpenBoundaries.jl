# OpenBoundaries.jl

This repository contains the beginnings of a library providing open boundary matching schemes for [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) simulations.

So far, I have only implemented `PerturbationAdvection` where the perturbation component (defined as the boundary adjacent wall normal velocity minus the prescribed boundary velocity) is advected through the domain, and relaxed to the prescribed velocities at different rates for inflow and outflow.
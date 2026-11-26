```@raw html
---
layout: home

hero:
  name: "EKRem.jl"
  text: "Electrokinetic remediation"
  tagline: Multicomponent ion transport in porous media under an electric field, coupled with geochemical reactions
  actions:
    - theme: brand
      text: Get started
      link: /guide/installation
    - theme: alt
      text: Theory
      link: /theory
    - theme: alt
      text: Examples
      link: /examples/1D_trace_electromigration

features:
  - title: Nernst–Planck transport
    details: Diffusion and electromigration of any number of aqueous species, discretised by finite volumes with Scharfetter–Gummel fluxes.
    link: /theory#transport
  - title: Electric potential
    details: Poisson equation, charge conservation in an electroneutral solution, or a fixed field; voltage or current control.
    link: /theory#electric-potential
  - title: Geochemical reactions
    details: Equilibrium and kinetic reactions computed by PhreeqcRM in every node, with PHREEQC databases and input.
    link: /guide/reactions
  - title: Electrode compartments
    details: Well-mixed anolyte and catholyte with water electrolysis and flushing.
    link: /theory#electrode-compartments
---
```

```@raw html
<div class="vp-doc" style="width:80%; margin:auto">
```

## What is EKRem?

EKRem is a Julia package for the simulation and design of electrokinetic remediation
systems. Multicomponent transport is described by the Nernst–Planck equations, the electric
potential by the Poisson or the charge-conservation equation, and geochemical reactions are
computed by [PhreeqcRM](https://github.com/phreeqc-dev/phreeqcrm). The transport equations
are discretised with finite volumes and Scharfetter–Gummel fluxes
([VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl)).

The models are organised by spatial dimension: [`Model1D`](@ref) for columns and laboratory
setups, and `Model2D` and `Model3D` for field applications.

## Quick start

Electromigration of a trace cation through a NaCl solution under an applied voltage of 1 V across
a 10 cm column:

```julia
using EKRem

species = [Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9), Species("T+", 1, 1.0e-9)]
background = Dict("Na+" => 100.0, "Cl-" => 100.0)                 # mol/m³

model = Model1D(;
    length=0.1,                                                     # m
    grid=(h_min=1.0e-5, h_max=1.0e-4, growth=1.1),
    species,
    initial=background,
    left=Reservoir(Dict("Na+" => 100.0, "Cl-" => 100.01, "T+" => 0.01); potential=1.0),
    right=Reservoir(background; potential=0.0),
    porosity=0.4,
    tortuosity=0.5,
)
run!(model, [1.0e4, 1.0e5])                                         # s

x = coordinates(model)
c = concentration(model, "T+", 1.0e5)
φ = electric_potential(model, 1.0e5)
```

## Status

| Feature | Status |
|:--------|:-------|
| `Model1D`: Nernst–Planck transport with reservoirs at both ends | implemented |
| Potential: Poisson, charge conservation, fixed field | implemented |
| Voltage control and current control | implemented |
| Geochemical reactions (PhreeqcRM), sequential non-iterative coupling | implemented |
| Electrode compartments under current or voltage control | implemented |
| Electroosmosis (uniform electroosmotic permeability) | implemented |
| `Model2D`, `Model3D` | not implemented |

The package is at a pre-alpha stage; its interface may change.

```@raw html
</div>
```

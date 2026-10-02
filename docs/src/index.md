```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "HeatExchange.jl"
  text: "Heat and water budgets of organisms"
  tagline: "solve for the body temperature, or for the metabolic rate, of an animal or a leaf exchanging heat with its environment by radiation, convection, conduction and evaporation, with units."
  actions:
    - theme: brand
      text: Get Started
      link: /get_started
    - theme: alt
      text: View on Github
      link: https://github.com/BiophysicalEcology/HeatExchange.jl
    - theme: alt
      text: API Reference
      link: /api

features:
  - title: ⚖️ Heat balance
    details: <a class="highlight-link">Solar and longwave radiation, convection, conduction, evaporation, respiration and metabolism</a>, summed into a budget whose residual is zero at steady state.
    link: /manual/heat_balance
  - title: 🌡️ Temperature or metabolic rate
    details: Find the <a class="highlight-link">body temperature</a> for a given metabolism, or the <a class="highlight-link">metabolic rate</a> that holds a given core temperature, with bare skin or with insulation.
    link: /manual/solvers
  - title: 🧥 Fur, feathers and fat
    details: Heat conducted from the core through <a class="highlight-link">flesh, fat and a porous coat</a> as a stack of radial layers, with radiation and conduction through the fibres.
    link: /manual/radial_layers
  - title: 🐕 Bodies of many parts
    details: A back and a belly, or a trunk, head and limbs, each with <a class="highlight-link">its own heat budget</a>, sharing a core or conducting heat to each other.
    link: /manual/multipart
  - title: 🍃 Leaves
    details: The same heat budget with <a class="highlight-link">stomatal conductance</a> in place of wet skin, for leaf temperature and transpiration.
    link: /tutorials/leaf
  - title: 🦎 For NicheMapR users
    details: How the <a class="highlight-link">ectotherm and endotherm models of NicheMapR</a> map onto this package, and what has changed.
    link: /manual/nichemapr
  - title: 📐 Differentiable
    details: Written so that an optimiser can drive the heat budget with <a class="highlight-link">automatic differentiation</a>, as BiophysicalBehaviour.jl does.
    link: /manual/autodiff
  - title: 📏 Units
    details: Every temperature, heat flow and trait is a <a class="highlight-link">Unitful.jl</a> quantity.
    link: /manual/introduction
---
```

## How to install HeatExchange.jl?

HeatExchange.jl can be installed from the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/BiophysicalEcology/HeatExchange.jl")
```

## Manual

HeatExchange.jl computes the exchange of heat, and of the water that goes with it, between an organism and its
environment. Given the organism's size, shape, surface properties and physiology, and the air temperature, wind,
humidity and radiation around it, it solves the steady-state heat budget for the body temperature, or for the
metabolic rate needed to hold a body temperature.

It is part of the [BiophysicalEcology](https://github.com/BiophysicalEcology) ecosystem for mechanistic niche
modelling. The bodies come from
[BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl), the environments from
[Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl), and the behaviour that changes both from
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl). See the
[Introduction](manual/introduction.md) for the design of the package,
[Environments and the ecosystem](manual/ecosystem.md) for how the packages fit together, and
[For NicheMapR users](manual/nichemapr.md) for its origin in [NicheMapR](https://github.com/mrke/NicheMapR).

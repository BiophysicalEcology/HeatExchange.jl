using Documenter
using DocumenterVitepress
using HeatExchange
using CairoMakie
using Unitful

# Don't output huge svgs for Makie plots
CairoMakie.activate!(type = "png")

# Helpers for the figures, loaded in the examples with `using Main.FigureHelpers`
include("figure_helpers.jl")

# Drawings of bodies, from the documentation of BiophysicalGeometry.jl, loaded with `using Main.GeometryFigures`
include("geometry_figures.jl")

makedocs(
    modules = [HeatExchange],
    sitename = "HeatExchange.jl",
    authors = "Michael Kearney, Urtzi Enriquez Urzelai, Rafael Schouten et al.",
    clean = true,
    doctest = false,
    checkdocs = :exports,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/BiophysicalEcology/HeatExchange.jl", # this must be the full URL!
        devbranch = "main",
        devurl = "dev";
    ),
    source = "src",
    build = "build",
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/BiophysicalEcology/HeatExchange.jl",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)

module FigureHelpers

using CairoMakie
using Markdown
using Unitful
using HeatExchange
using BiophysicalGeometry
using ModelParameters: stripparams

export figure_axis, markdown_table, parameter_table, flow_table, celsius, watts, budget_bars, budget_bars!,
    heat_budget_diagram, layer_budget_diagram, solver_paths_diagram, smoothing_figure, radial_stack_diagram,
    compartment_diagram, radial_network_diagram, FLOW_COLOURS

"""
    figure_axis(xlabel, ylabel; size=(700, 500), kw...)

A `Figure` and `Axis` with minor ticks and grid lines, as used for the figures of the manual.
"""
function figure_axis(xlabel, ylabel; size=(700, 500), kw...)
    fig = Figure(; size)
    ax = Axis(fig[1, 1];
        xlabel, ylabel,
        xminorticksvisible=true, yminorticksvisible=true,
        xminorgridvisible=true, yminorgridvisible=true,
        xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
        kw...,
    )
    return fig, ax
end

"""
    markdown_table(header, rows)

A Markdown table with column names `header` and one row per element of `rows`, with numbers and quantities
rounded to 4 significant digits.
"""
function markdown_table(header, rows)
    lines = ["| " * join(header, " | ") * " |", "|" * repeat(":--|", length(header))]
    append!(lines, ["| " * join(map(_format, row), " | ") * " |" for row in rows])
    return Markdown.parse(join(lines, "\n"))
end

"""
    parameter_table(pars)
    parameter_table(label => pars, ...)

A Markdown table of the fields of one parameter struct, or of several of the same type side by side.
"""
parameter_table(pars) = parameter_table("Value" => pars)
function parameter_table(columns::Pair...)
    objects = map(stripparams ∘ last, columns)
    rows = [("`$name`", (getfield(object, name) for object in objects)...) for name in fieldnames(typeof(first(objects)))]
    return markdown_table(["Parameter", first.(columns)...], rows)
end

"""
    flow_table(flows; header=["Term", "Value"])

A Markdown table of the numbers in a NamedTuple, with powers in W.
"""
function flow_table(flows; header=["Term", "Value"])
    rows = [("`$name`", _canonical(value)) for (name, value) in pairs(flows) if value isa Number]
    return markdown_table(header, rows)
end

"""
    celsius(temperature)

The temperature in °C.
"""
celsius(temperature) = uconvert(u"°C", temperature)

"""
    watts(power)

The power in W, as a number.
"""
watts(power) = ustrip(u"W", power)

_canonical(x) = x
function _canonical(x::Unitful.AbstractQuantity)
    dimension(x) == dimension(1u"W") && return uconvert(u"W", x)
    return x
end

_format(x::Unitful.Quantity) = (y = _canonical(x); string(round(unit(y), y; sigdigits = 4)))
_format(x::Bool) = string(x)
_format(x::Real) = string(round(x; sigdigits = 4))
_format(::Nothing) = "–"
_format(x::Symbol) = "`$x`"
_format(x::AbstractString) = x
_format(x) = string(nameof(typeof(x)))

# ── Colours ──────────────────────────────────────────────────────────────────

const FLOW_COLOURS = (
    solar = RGBf(0.95, 0.70, 0.15), longwave = RGBf(0.80, 0.35, 0.30), metabolism = RGBf(0.55, 0.40, 0.70),
    convection = RGBf(0.35, 0.60, 0.80), conduction = RGBf(0.55, 0.42, 0.35), evaporation = RGBf(0.30, 0.65, 0.60),
    respiration = RGBf(0.60, 0.75, 0.45),
)
const LAYER_COLOURS = (flesh = RGBf(0.88, 0.48, 0.42), fat = RGBf(1.00, 0.93, 0.55), fur = RGBf(0.76, 0.62, 0.42))

# ── Heat budgets as bars ─────────────────────────────────────────────────────

"""
    budget_bars(gains, losses; title="")
    budget_bars!(ax, gains, losses)

Two stacked bars, heat gained and heat lost, from `label => (power, colour)` pairs.
"""
function budget_bars(gains, losses; size=(640, 360), title="")
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; ylabel = "Heat flow (W)", xticks = (1:2, ["Gained", "Lost"]), title)
    budget_bars!(ax, gains, losses)
    Legend(fig[1, 2], ax; framevisible = false)
    return fig
end
function budget_bars!(ax, gains, losses)
    for (x, terms) in enumerate((gains, losses))
        base = 0.0
        for (label, (power, colour)) in terms
            height = watts(power)
            barplot!(ax, [x], [height]; offset = base, color = colour, strokecolor = :black, strokewidth = 0.5,
                     width = 0.6, label)
            base += height
        end
    end
    return ax
end

# ── Schematics ───────────────────────────────────────────────────────────────

_blank_axis(position; kw...) = (ax = Axis(position; aspect = DataAspect(), kw...); hidedecorations!(ax); hidespines!(ax); ax)

function _arrow!(ax, from, to; color = :black, linewidth = 2.5)
    lines!(ax, [from[1], to[1]], [from[2], to[2]]; color, linewidth)
    angle = atan(to[2] - from[2], to[1] - from[1])
    scatter!(ax, [to[1]], [to[2]]; marker = :utriangle, rotation = angle - π / 2, color, markersize = 14)
    return nothing
end

_circle(r; n = 120) = [Point2f(r * cos(θ), r * sin(θ)) for θ in range(0, 2π; length = n)]

"""
    heat_budget_diagram()

The heat flows of an organism with bare skin, after Fig. 2 of Kearney and Porter (2020).
"""
function heat_budget_diagram()
    fig = Figure(size = (700, 400))
    ax = _blank_axis(fig[1, 1])
    poly!(ax, [Point2f(1.6 * p[1], p[2]) for p in _circle(1.0)]; color = (LAYER_COLOURS.flesh, 0.6), strokecolor = :black,
          strokewidth = 1)
    text!(ax, 0, 0.15; text = rich("Q", subscript("met")), align = (:center, :center), fontsize = 18,
          color = FLOW_COLOURS.metabolism)
    text!(ax, 0, -0.35; text = rich("T", subscript("c")), align = (:center, :center), fontsize = 16)
    lines!(ax, [-3.4, 3.4], [-1.05, -1.05]; color = :grey40, linewidth = 2)
    flows = (
        ((-2.6, 2.4), (-1.15, 0.75), rich("Q", subscript("sol")), FLOW_COLOURS.solar, (-2.75, 2.5)),
        ((-0.9, 2.5), (-0.45, 1.05), rich("Q", subscript("IR,in")), FLOW_COLOURS.longwave, (-1.2, 2.75)),
        ((0.45, 1.05), (0.9, 2.5), rich("Q", subscript("IR,out")), FLOW_COLOURS.longwave, (1.25, 2.75)),
        ((1.2, 0.75), (2.6, 1.9), rich("Q", subscript("conv")), FLOW_COLOURS.convection, (2.85, 2.0)),
        ((1.6, 0.1), (3.1, 0.5), rich("Q", subscript("evap")), FLOW_COLOURS.evaporation, (3.35, 0.62)),
        ((-1.6, 0.1), (-3.1, 0.5), rich("Q", subscript("resp")), FLOW_COLOURS.respiration, (-3.45, 0.75)),
        ((0.0, -1.0), (0.0, -1.9), rich("Q", subscript("cond")), FLOW_COLOURS.conduction, (0.55, -1.75)),
    )
    for (from, to, label, colour, at) in flows
        _arrow!(ax, from, to; color = colour)
        text!(ax, at...; text = label, align = (:center, :center), fontsize = 16, color = colour)
    end
    text!(ax, 3.3, -1.3; text = "substrate", align = (:right, :center), fontsize = 12, color = :grey40)
    limits!(ax, -3.8, 3.8, -2.2, 3.0)
    return fig
end

"""
    layer_budget_diagram()

A cross-section of an insulated organism with the three temperatures and the heat flows between them, after the
system diagram of the NicheMapR endotherm model.
"""
function layer_budget_diagram()
    fig = Figure(size = (720, 420))
    ax = _blank_axis(fig[1, 1])
    for (r, colour) in ((2.0, LAYER_COLOURS.fur), (1.5, LAYER_COLOURS.fat), (1.2, LAYER_COLOURS.flesh))
        poly!(ax, _circle(r); color = (colour, 0.75), strokecolor = :black, strokewidth = 1)
    end
    text!(ax, 0, 0.25; text = rich("Q", subscript("gen")), align = (:center, :center), fontsize = 18,
          color = FLOW_COLOURS.metabolism)
    scatter!(ax, [0.0, 1.5, 2.0], [-0.35, 0.0, 0.0]; color = :black, markersize = 9)
    text!(ax, 0, -0.75; text = rich("T", subscript("c")), align = (:center, :center), fontsize = 16)
    text!(ax, 1.3, -0.35; text = rich("T", subscript("s")), align = (:center, :center), fontsize = 16)
    text!(ax, 2.25, -0.3; text = rich("T", subscript("fa")), align = (:center, :center), fontsize = 16)
    for (y, label) in ((1.0, "flesh"), (1.35, "fat"), (1.75, "fur"))
        text!(ax, 0, y; text = label, align = (:center, :center), fontsize = 11, color = :grey20)
    end
    flows = (
        ((-0.9, 0.0), (-3.2, 0.0), rich("Q", subscript("resp")), FLOW_COLOURS.respiration, (-3.7, 0.0)),
        ((-1.1, 1.05), (-2.6, 2.0), rich("Q", subscript("evap")), FLOW_COLOURS.evaporation, (-3.0, 2.2)),
        ((1.45, 1.4), (2.5, 2.3), rich("Q", subscript("conv")), FLOW_COLOURS.convection, (2.9, 2.5)),
        ((0.3, 2.0), (0.6, 3.0), rich("Q", subscript("rad")), FLOW_COLOURS.longwave, (0.75, 3.25)),
        ((-0.9, 3.0), (-0.45, 2.0), rich("Q", subscript("sol")), FLOW_COLOURS.solar, (-1.05, 3.25)),
        ((0.0, -2.0), (0.0, -2.9), rich("Q", subscript("cond")), FLOW_COLOURS.conduction, (0.6, -2.7)),
        ((2.0, 0.6), (3.2, 0.9), rich("Q", subscript("evap,fur")), FLOW_COLOURS.evaporation, (3.9, 1.0)),
    )
    for (from, to, label, colour, at) in flows
        _arrow!(ax, from, to; color = colour)
        text!(ax, at...; text = label, align = (:center, :center), fontsize = 15, color = colour)
    end
    lines!(ax, [-4.2, 4.2], [-2.05, -2.05]; color = :grey40, linewidth = 2)
    limits!(ax, -4.6, 4.6, -3.2, 3.6)
    return fig
end

function _box!(ax, x, y, w, h, label; color = :grey92, fontsize = 12)
    poly!(ax, Rect2f(x - w / 2, y - h / 2, w, h); color, strokecolor = :black, strokewidth = 1)
    text!(ax, x, y; text = label, align = (:center, :center), fontsize)
end

"""
    solver_paths_diagram()

The two ways the residuals of one part are driven to zero.
"""
function solver_paths_diagram()
    fig = Figure(size = (720, 330))
    ax = _blank_axis(fig[1, 1])
    blue, orange = RGBf(0.78, 0.87, 0.95), RGBf(0.98, 0.87, 0.70)
    _box!(ax, 0, 0, 5.2, 0.9, "solve_part_heat_balance\ntemperatures and metabolic rate in, residuals out"; color = :grey88)
    _box!(ax, -3.2, 1.8, 3.6, 0.8, "solve_part_surface\nNewton iteration"; color = blue)
    _box!(ax, 3.2, 1.8, 3.6, 0.8, "part_surface_residuals\nno iteration"; color = orange)
    _box!(ax, -3.2, 3.5, 3.6, 0.8, "solve_metabolic_rate\nsolve_coupled_metabolic_rate"; color = blue)
    _box!(ax, 3.2, 3.5, 3.6, 0.8, "IPOPT, with Enzyme\n(BiophysicalBehaviour.jl)"; color = orange)
    for x in (-3.2, 3.2)
        _arrow!(ax, (x, 3.1), (x, 2.25); linewidth = 1.5)
        _arrow!(ax, (x, 1.4), (sign(x) * 1.6, 0.55); linewidth = 1.5)
    end
    text!(ax, -3.2, 4.25; text = "rules", align = (:center, :center), fontsize = 13, font = :bold)
    text!(ax, 3.2, 4.25; text = "optimisation", align = (:center, :center), fontsize = 13, font = :bold)
    limits!(ax, -5.4, 5.4, -0.7, 4.6)
    return fig
end

"""
    smoothing_figure(; ε=(1e-2, 1e-1, 3e-1))

`safe_abs` and `safe_step` near zero for `HardBound` and for `SmoothBound` of several widths.
"""
function smoothing_figure(; ε = (1e-2, 1e-1, 3e-1))
    fig = Figure(size = (760, 340))
    x = range(-1, 1; length = 801)
    for (i, (f, name)) in enumerate(((safe_abs, "safe_abs(x)"), (safe_step, "safe_step(x)")))
        ax = Axis(fig[1, i]; xlabel = "x", ylabel = name)
        lines!(ax, x, f.(Ref(HardBound()), x); color = :black, linewidth = 2, label = "HardBound()")
        for e in ε
            lines!(ax, x, f.(Ref(SmoothBound(e)), x); linewidth = 2, label = "SmoothBound($e)")
        end
        i == 2 && axislegend(ax; position = :rb, labelsize = 11)
    end
    return fig
end

"""
    radial_stack_diagram(labels, radii, temperatures; colours)

Concentric layers on the left and the temperature at each boundary on the right. `radii` are the outer radii of
the layers, from the centre outwards, and `temperatures` has one more element, starting at the centre.
"""
function radial_stack_diagram(labels, radii, temperatures;
                              colours = (LAYER_COLOURS.flesh, LAYER_COLOURS.fat, LAYER_COLOURS.fur))
    fig = Figure(size = (760, 340))
    r = [ustrip(u"cm", x) for x in radii]
    T = [ustrip(u"°C", x) for x in temperatures]
    ax1 = _blank_axis(fig[1, 1])
    for i in reverse(eachindex(r))
        poly!(ax1, _circle(r[i]); color = (colours[i], 0.8), strokecolor = :black, strokewidth = 1)
    end
    ax2 = Axis(fig[1, 2]; xlabel = "Distance from the centre (cm)", ylabel = "Temperature (°C)")
    edges = [0.0; r]
    for i in eachindex(r)
        vspan!(ax2, edges[i], edges[i + 1]; color = (colours[i], 0.5))
    end
    scatterlines!(ax2, edges, T; color = :black, linewidth = 2, markersize = 9)
    span = maximum(T) - minimum(T)
    ylims!(ax2, minimum(T) - 0.3 * span, maximum(T) + 0.05 * span)
    for i in eachindex(r)
        narrow = edges[i + 1] - edges[i] < 0.12 * r[end]
        text!(ax2, (edges[i] + edges[i + 1]) / 2, minimum(T) - 0.27 * span; text = labels[i], fontsize = 11,
              align = narrow ? (:left, :center) : (:center, :bottom), rotation = narrow ? π / 2 : 0.0)
    end
    colsize!(fig.layout, 1, Relative(0.38))
    return fig
end

"""
    compartment_diagram(graph; edges=(), positions=nothing)

The parts of a `CompartmentGraph`, coloured by compartment. `edges` are `(part, part, style)` with style
`:shared` or `:conductive`.
"""
function compartment_diagram(graph; edges = (), positions = nothing, size = (640, 300))
    names = compartment_part_names(graph)
    positions = positions === nothing ? Dict(name => Point2f(1.6 * (i - 1), isodd(i) ? 0.0 : 0.6) for (i, name) in enumerate(names)) :
                                    Dict(name => Point2f(xy...) for (name, xy) in positions)
    palette = [RGBf(0.30, 0.55, 0.75), RGBf(0.90, 0.62, 0.20), RGBf(0.35, 0.68, 0.45), RGBf(0.80, 0.40, 0.45),
               RGBf(0.58, 0.45, 0.72), RGBf(0.55, 0.42, 0.35), RGBf(0.85, 0.55, 0.75), RGBf(0.50, 0.50, 0.50)]
    fig = Figure(; size)
    ax = _blank_axis(fig[1, 1])
    for (a, b, style) in edges
        pa, pb = positions[a], positions[b]
        lines!(ax, [pa, pb]; color = :black, linewidth = style == :shared ? 5 : 2,
               linestyle = style == :shared ? :solid : :dash)
    end
    for name in names
        colour = palette[mod1(compartment_of(graph, name), length(palette))]
        scatter!(ax, [positions[name]]; color = colour, markersize = 46, strokecolor = :black, strokewidth = 1)
        left = positions[name][1] <= sum(first, values(positions)) / length(positions)
        text!(ax, positions[name] + Point2f(left ? -0.3 : 0.3, 0); text = string(name),
              align = (left ? :right : :left, :center), fontsize = 12)
        text!(ax, positions[name]; text = string(compartment_of(graph, name)), align = (:center, :center), fontsize = 14,
              color = :white, font = :bold)
    end
    xs, ys = first.(values(positions)), last.(values(positions))
    limits!(ax, minimum(xs) - 1.4, maximum(xs) + 1.4, minimum(ys) - 0.6, maximum(ys) + 0.6)
    return fig
end


"""
    radial_network_diagram(; ground=true, sources=true)

The heat balance of one part as a network of temperature nodes joined by thermal resistances, from the core to
the environment, with the branch through compressed fur to the substrate and the heat sources at each node.
"""
function radial_network_diagram(; ground = true, sources = true)
    fig = Figure(size = (760, ground ? 400 : 300))
    ax = _blank_axis(fig[1, 1])
    nodes = (core = (0.0, 0.0), boundary = (2.2, 0.0), skin = (4.4, 0.0), surface = (6.6, 0.0), environment = (8.8, 0.0))
    resistor(a, b, label; dy = 0.0) = begin
        lines!(ax, [a[1], b[1]], [a[2], b[2]]; color = :black, linewidth = 2)
        mid = ((a[1] + b[1]) / 2, (a[2] + b[2]) / 2)
        poly!(ax, Rect2f(mid[1] - 0.45, mid[2] - 0.17, 0.9, 0.34); color = :white, strokecolor = :black, strokewidth = 1.5)
        text!(ax, mid[1], mid[2] + dy; text = label, align = (:center, :center), fontsize = 12)
    end
    resistor(nodes.core, nodes.boundary, rich("R", subscript("flesh")))
    resistor(nodes.boundary, nodes.skin, rich("R", subscript("fat")))
    resistor(nodes.skin, nodes.surface, rich("R", subscript("fur")))
    resistor(nodes.surface, nodes.environment, "conv, rad")
    if ground
        compressed, substrate = (6.6, -1.6), (8.8, -1.6)
        lines!(ax, [nodes.skin[1], nodes.skin[1], compressed[1] - 1.1], [nodes.skin[2], compressed[2], compressed[2]]; color = :black, linewidth = 2)
        resistor((4.4 + 0.6, -1.6), compressed, rich("R", subscript("fur,comp")))
        resistor(compressed, substrate, "contact")
        scatter!(ax, [Point2f(compressed)]; color = LAYER_COLOURS.fur, markersize = 30, strokecolor = :black, strokewidth = 1)
        scatter!(ax, [Point2f(substrate)]; color = :grey80, markersize = 30, strokecolor = :black, strokewidth = 1)
        text!(ax, compressed[1], compressed[2] - 0.45; text = "compressed fur", align = (:center, :top), fontsize = 11)
        text!(ax, substrate[1], substrate[2] - 0.45; text = "substrate", align = (:center, :top), fontsize = 11)
    end
    colours = (core = LAYER_COLOURS.flesh, boundary = :white, skin = RGBf(0.95, 0.80, 0.70), surface = LAYER_COLOURS.fur, environment = :grey80)
    labels = (core = rich("T", subscript("c")), boundary = "", skin = rich("T", subscript("s")), surface = rich("T", subscript("fa")), environment = "")
    names = (core = "core", boundary = "flesh–fat
(eliminated)", skin = "skin", surface = "fur surface", environment = "air, sky,
ground")
    for key in keys(nodes)
        scatter!(ax, [Point2f(nodes[key])]; color = colours[key], markersize = key == :boundary ? 18 : 34, strokecolor = key == :boundary ? :grey50 : :black,
                 strokewidth = 1)
        text!(ax, nodes[key]...; text = labels[key], align = (:center, :center), fontsize = 13)
        text!(ax, nodes[key][1], nodes[key][2] + 0.5; text = names[key], align = (:center, :bottom), fontsize = 11)
    end
    if sources
        for (x, label, colour, up) in ((0.0, rich("+Q", subscript("gen"), "  −Q", subscript("resp")), FLOW_COLOURS.metabolism, false),
                                       (4.4, rich("−Q", subscript("evap")), FLOW_COLOURS.evaporation, true),
                                       (6.6, rich("+Q", subscript("sol"), "  −Q", subscript("evap,fur")), FLOW_COLOURS.solar, true))
            y = up ? 1.55 : -0.75
            text!(ax, x, y; text = label, align = (:center, :center), fontsize = 12, color = colour)
        end
    end
    limits!(ax, -1.2, 10.0, ground ? -2.6 : -1.3, 2.0)
    return fig
end

end

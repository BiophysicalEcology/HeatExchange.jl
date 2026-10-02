module GeometryFigures

using CairoMakie
using Markdown
using Unitful
using BiophysicalGeometry
import BiophysicalGeometry: Sphere, Cylinder, Ellipsoid, Cone, Plate, Half, Top, Bottom
import BiophysicalGeometry: AbstractCylindrical, AbstractEllipsoidal, AbstractSpherical

const BG = BiophysicalGeometry

export temperature_views, body_axis, single, draw_parts!, draw_layers!, shape_gallery, surface_diagram,
    surface_diagram!, layer_diagram, layer_diagram!, layer_sections, long_section!, body_graph, body_graph!, composite_views,
    silhouette_panel!, PART_COLOURS, LAYER_COLOURS

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

_format(x::Unitful.Quantity) = string(round(unit(x), x; sigdigits = 4))
_format(x::Real) = string(round(x; sigdigits = 4))
_format(::Nothing) = "–"
_format(x) = string(x)

# ── Colours ──────────────────────────────────────────────────────────────────

const LAYER_COLOURS = (flesh = RGBf(0.88, 0.48, 0.42), fat = RGBf(1.00, 0.93, 0.55), fur = RGBf(0.76, 0.62, 0.42))
const PART_COLOURS = [RGBf(0.30, 0.55, 0.75), RGBf(0.90, 0.62, 0.20), RGBf(0.35, 0.68, 0.45), RGBf(0.80, 0.40, 0.45),
                      RGBf(0.58, 0.45, 0.72), RGBf(0.55, 0.42, 0.35), RGBf(0.85, 0.55, 0.75), RGBf(0.50, 0.50, 0.50)]
const SURFACE_COLOURS = [RGBf(0.30, 0.55, 0.75), RGBf(0.90, 0.62, 0.20), RGBf(0.35, 0.68, 0.45), RGBf(0.80, 0.40, 0.45),
                         RGBf(0.58, 0.45, 0.72), RGBf(0.55, 0.42, 0.35)]

_cm(x) = ustrip(u"cm", x)

# ── 3-D axes and composite drawing ───────────────────────────────────────────

"""
    body_axis(position; decorations=true, azimuth=1.25π, elevation=π/7, kw...)

An `Axis3` with equal scaling on all axes, in cm, for drawing bodies.
"""
function body_axis(position; decorations=true, azimuth=1.25π, elevation=π/7, kw...)
    ax = Axis3(position; aspect=:data, perspectiveness=0.0, viewmode=:fit, azimuth, elevation,
        xlabel="x (cm)", ylabel="y (cm)", zlabel="z (cm)",
        xlabelsize=11, ylabelsize=11, zlabelsize=11,
        xticklabelsize=9, yticklabelsize=9, zticklabelsize=9, protrusions=(30, 10, 10, 10), kw...)
    if !decorations
        hidedecorations!(ax)
        hidespines!(ax)
    end
    return ax
end

"""
    single(body)

A one-part `CompositeBody` holding `body`, for the functions that need a composite.
"""
single(body::Body) = CompositeBody(; parts = (; body), joins = ())
single(body::CompositeBody) = body

# ── Drawing surfaces ─────────────────────────────────────────────────────────
#
# Drawing is done by the Makie extension of the package, which puts everything for one axis into a single mesh so
# that its faces are sorted by depth. `tiles` is a vector of `(X, Y, Z, colour)` grids.

const EXT = Base.get_extension(BiophysicalGeometry, :BiophysicalGeometryMakieExt)

_view_direction(ax) = EXT._view_direction(ax.azimuth[], ax.elevation[])
_draw_tiles!(ax, tiles) = EXT._mesh_tiles!(ax, tiles, ax.azimuth[], ax.elevation[])

function _part_tiles(body, colors, sc)
    composite = single(body)
    names = propertynames(composite.parts)
    used = map(enumerate(names)) do (i, name)
        colors isa NamedTuple ? getfield(colors, name) : colors[mod1(i, length(colors))]
    end
    tiles = []
    for (name, col) in zip(names, used)
        part = getfield(composite.parts, name)
        pose = getfield(composite.poses, name)
        for grid in BG._part_outer_meshes(part.shape, part, sc)
            push!(tiles, (BG._transform_mesh(grid..., pose, sc)..., col))
        end
    end
    return tiles, NamedTuple{names}(Tuple(used))
end

"""
    draw_parts!(ax, body; colors=PART_COLOURS, sc=100.0)

Draw the outer surface of every part of `body` in its pose, one colour per part. `colors` is a vector, cycled
over the parts, or a `NamedTuple` keyed by part name. Returns the colours used, by part name.
"""
function draw_parts!(ax, body; colors=PART_COLOURS, sc=100.0)
    tiles, used = _part_tiles(body, colors, sc)
    _draw_tiles!(ax, tiles)
    return used
end

# ── Cutaways ─────────────────────────────────────────────────────────────────

_layer_colour(body::Body) = BG.outer_insulation(body.insulation) isa FibrousLayer ? LAYER_COLOURS.fur :
                            LAYER_COLOURS.flesh

"""
    draw_layers!(ax, body; cut=π/2, sc=100.0)

Draw `body` with its layers of flesh, fat and fibres, and the part of the outer layers facing the viewer cut away
over the angle `cut`. Layers too thin to see are not cut away.
"""
function draw_layers!(ax, body::Body; cut=π / 2, sc=100.0)
    if insulation_radius(body) - flesh_radius(body) < 0.03 * insulation_radius(body)
        _draw_tiles!(ax, first(_part_tiles(body, [_layer_colour(body)], sc)))
    else
        tiles = EXT._cutaway_tiles(body.shape, body, sc, LAYER_COLOURS,
            range(ax.azimuth[] + cut / 2, ax.azimuth[] + 2π - cut / 2; length=73))
        _draw_tiles!(ax, tiles)
    end
    return ax
end

"""
    shape_gallery(label => body, ...; ncols=3, size)

A grid of drawings, one per `label => body`: bodies with their layers cut away, composites by part.
"""
function shape_gallery(items::Pair...; ncols=3, size=(260 * min(ncols, length(items)), 250 * cld(length(items), ncols)),
                       decorations=false, kw...)
    fig = Figure(; size)
    for (i, (label, body)) in enumerate(items)
        row, col = fldmod1(i, ncols)
        ax = body_axis(fig[row, col]; title=label, titlesize=13, decorations, kw...)
        body isa CompositeBody ? draw_parts!(ax, body) : draw_layers!(ax, body)
    end
    return fig
end

"""
    composite_views(body; views=(:oblique, :side, :front, :top), titles, colors, legend=true, size)

The parts of `body` from several directions, one colour per part.
"""
function composite_views(body; views=(:oblique, :side, :front, :top), titles=string.(views), colors=PART_COLOURS,
                         legend=true, size=(230 * length(views), legend ? 300 : 250))
    angles = (oblique = (1.25π, π / 7), three_quarter = (-0.25π, π / 7), side = (-π / 2, 0.0), front = (0.0, 0.0), top = (-π / 2, π / 2),
              back = (π, 0.0))
    fig = Figure(; size)
    used = nothing
    for (i, view) in enumerate(views)
        azimuth, elevation = getfield(angles, view)
        ax = body_axis(fig[1, i]; decorations=false, azimuth, elevation, title=titles[i], titlesize=12)
        used = draw_parts!(ax, body; colors)
    end
    if legend
        Legend(fig[2, 1:length(views)], [PolyElement(color=c) for c in values(used)], collect(string.(keys(used)));
            orientation=:horizontal, framevisible=false, labelsize=11, nbanks=cld(length(used), 7))
    end
    return fig
end

# ── Named surfaces ───────────────────────────────────────────────────────────
#
# For each shape, the named attachment surfaces with the mesh tiles that show them, at skin level in cm.

function _surface_tiles(::Cylinder, body)
    r = _cm(skin_radius(body)); L = _cm(body.geometry.length.length_skin)
    [EndA() => [BG._cylinder_cap(r, 0.0)], EndB() => [BG._cylinder_cap(r, L)], Lateral() => [BG._cylinder_tube(r, L)]]
end
function _surface_tiles(sh::Cone, body)
    r = _cm(skin_radius(body)); L = _cm(body.geometry.length.length_skin); t = sh.top_ratio
    tiles = Pair[EndA() => [BG._cylinder_cap(r, 0.0)], Lateral() => [BG._cone_tube(r, t * r, L)]]
    t > 0 && push!(tiles, EndB() => [BG._cylinder_cap(t * r, L)])
    tiles
end
function _surface_tiles(::Sphere, body)
    r = _cm(skin_radius(body))
    [Radial() => [BG._ellipsoid_mesh(r, r)]]
end
function _surface_tiles(sh::Plate, body)
    l = body.geometry.length
    hl, hw, hh = _cm(l.length_skin) / 2, _cm(l.width_skin) / 2, _cm(l.height_skin) / 2
    [Top() => [BG._box_face_z(-hl, hl, -hw, hw, hh)], Bottom() => [BG._box_face_z(-hl, hl, -hw, hw, -hh)],
     SideA() => [BG._box_face_x(hl, -hw, hw, -hh, hh)], SideB() => [BG._box_face_x(-hl, -hw, hw, -hh, hh)],
     SideC() => [BG._box_face_y(-hl, hl, hw, -hh, hh)], SideD() => [BG._box_face_y(-hl, hl, -hw, -hh, hh)]]
end
function _surface_tiles(::Half{<:AbstractCylindrical}, body)
    r = _cm(skin_radius(body)); L = _cm(body.geometry.length.length_skin)
    [EndA() => [BG._cylinder_cap(r, 0.0; θ_end=π)], EndB() => [BG._cylinder_cap(r, L; θ_end=π)],
     Lateral() => [BG._cylinder_tube(r, L; θ_end=π)], Flat() => [BG._half_cylinder_flat(r, L)]]
end
function _surface_tiles(sh::Half{<:Union{AbstractEllipsoidal,AbstractSpherical}}, body)
    a, b, _ = _cm.(BG._domed_semiaxes(sh, body))
    [Dome() => [BG._ellipsoid_mesh(a, b; φ_end=π / 2)], Flat() => [BG._half_ellipsoid_flat_mesh(a, b)]]
end

_label(loc) = string(nameof(typeof(loc)))

function _extent(tiles)
    lo = fill(Inf, 3); hi = fill(-Inf, 3)
    for (_, meshes) in tiles, mesh in meshes, k in 1:3
        lo[k] = min(lo[k], minimum(mesh[k])); hi[k] = max(hi[k], maximum(mesh[k]))
    end
    return lo, hi
end

function _local_axes!(ax, lo, hi)
    reach = 0.22 * maximum(hi .- lo)
    for (k, name) in enumerate(("x", "y", "z"))
        tip = zeros(3); tip[k] = max(hi[k], 0.0) + reach
        lines!(ax, [Point3f(0, 0, 0), Point3f(tip...)]; color=:grey30, linewidth=1.5)
        text!(ax, Point3f(tip...); text=name, fontsize=13, font=:italic, color=:grey20, align=(:center, :bottom),
            overdraw=true)
    end
end

function _surface_label!(ax, sh, body, loc, col, offset)
    c = _cm.(BG.surface_centroid(sh, body, loc)); n = BG.surface_centroid_normal(sh, body, loc)
    p = Point3f((c .+ offset .* n)...)
    lines!(ax, [Point3f(c...), p]; color=:black, linewidth=1, overdraw=true)
    scatter!(ax, [Point3f(c...)]; color=:black, markersize=5, overdraw=true)
    d = _view_direction(ax)
    right = n[1] * -d[2] + n[2] * d[1]      # is the label to the right of the shape as seen?
    align = abs(right) < 0.3 ? (:center, n[3] < 0 ? :top : :bottom) : (right > 0 ? :left : :right, :center)
    text!(ax, p; text=_label(loc), fontsize=13, font=:bold, color=col, align, overdraw=true)
end

"""
    surface_diagram!(ax, body)

Draw `body` with each named attachment surface in its own colour, labelled, with the local axes.
"""
function surface_diagram!(ax, body::Body)
    sh = body.shape
    surfaces = _surface_tiles(sh, body)
    lo, hi = _extent(surfaces)
    offset = 0.3 * maximum(hi .- lo)
    colour(i) = SURFACE_COLOURS[mod1(i, length(SURFACE_COLOURS))]
    _draw_tiles!(ax, [(grid..., colour(i)) for (i, (_, grids)) in enumerate(surfaces) for grid in grids])
    _local_axes!(ax, lo, hi)
    for (i, (loc, _)) in enumerate(surfaces)
        _surface_label!(ax, sh, body, loc, colour(i), offset)
    end
    return ax
end

# An ellipsoid's surfaces are two poles and a ring, not faces.
function surface_diagram!(ax, body::Body{<:Ellipsoid})
    sh = body.shape
    l = body.geometry.length
    a, b = _cm(l.a_semi_major_skin), _cm(l.b_semi_minor_skin)
    grey = RGBf(0.85, 0.87, 0.89)
    x_ratio = 1 - sh.pole_a_truncation
    if sh.pole_a_truncation == 0
        _draw_tiles!(ax, [(BG._ellipsoid_mesh(a, b)..., grey)])
    else
        _draw_tiles!(ax, [(BG._ellipsoid_mesh_truncated(a, b, b, x_ratio)..., grey),
                          (BG._ellipsoid_pole_a_cap(a, b, b, x_ratio)..., SURFACE_COLOURS[1])])
    end
    _local_axes!(ax, [-a, -b, -b], [a, b, b])
    ts = range(0, 2π; length=100)
    lines!(ax, [Point3f(0, b * cos(t), b * sin(t)) for t in ts]; color=SURFACE_COLOURS[3], linewidth=3, overdraw=true)
    scatter!(ax, [Point3f(a * x_ratio, 0, 0)]; color=SURFACE_COLOURS[1], markersize=13, overdraw=true)
    scatter!(ax, [Point3f(-a, 0, 0)]; color=SURFACE_COLOURS[2], markersize=13, overdraw=true)
    offset = 0.3 * 2a
    for (i, loc) in enumerate((PoleA(), PoleB(), Equator()))
        _surface_label!(ax, sh, body, loc, SURFACE_COLOURS[i], offset)
    end
    return ax
end

"""
    surface_diagram(label => body, ...; ncols, size)

A grid of `surface_diagram!` drawings.
"""
function surface_diagram(items::Pair...; ncols=length(items), size=(330 * min(ncols, length(items)),
                         300 * cld(length(items), ncols)), kw...)
    fig = Figure(; size)
    for (i, (label, body)) in enumerate(items)
        row, col = fldmod1(i, ncols)
        ax = body_axis(fig[row, col]; decorations=false, title=label, titlesize=13, kw...)
        surface_diagram!(ax, body)
    end
    return fig
end
surface_diagram(body::Body; kw...) = surface_diagram("" => body; kw...)

# ── Layers ───────────────────────────────────────────────────────────────────

_ring(r; n=200) = [Point2f(r * cos(t), r * sin(t)) for t in range(0, 2π; length=n)]

"""
    layer_diagram!(ax, body; labels=true)

The transverse section of `body`, with its layers and the flesh, skin and insulation radii marked.
"""
function layer_diagram!(ax, body::Body; labels=true)
    flesh, skin, ins = _cm(flesh_radius(body)), _cm(skin_radius(body)), _cm(insulation_radius(body))
    has_fat = skin > flesh * (1 + 1e-9); has_fur = ins > skin * (1 + 1e-9)
    has_fur && poly!(ax, _ring(ins); color=LAYER_COLOURS.fur, strokecolor=:grey30, strokewidth=1)
    has_fat && poly!(ax, _ring(skin); color=LAYER_COLOURS.fat)
    poly!(ax, _ring(flesh); color=LAYER_COLOURS.flesh)
    lines!(ax, _ring(skin); color=:black, linewidth=2)
    if labels
        marks = [(flesh, 200.0, "flesh_radius")]
        has_fat && push!(marks, (skin, 160.0, "skin_radius"))
        !has_fat && !has_fur && (marks[1] = (flesh, 200.0, "all three radii"))
        !has_fat && has_fur && (marks[1] = (flesh, 200.0, "flesh and skin radius"))
        has_fur && push!(marks, (ins, 120.0, "insulation_radius"))
        for (r, angle, name) in marks
            tip = Point2f(r * cosd(angle), r * sind(angle))
            lines!(ax, [Point2f(0, 0), tip]; color=:black, linewidth=1.5)
            scatter!(ax, [tip]; color=:black, markersize=7)
            out = Point2f(1.25 * ins * cosd(angle), 1.25 * ins * sind(angle))
            lines!(ax, [tip, out]; color=:black, linewidth=0.7, linestyle=:dot)
            text!(ax, out; text=name, fontsize=11, font=:bold, align=(:right, :center), offset=(-3, 0))
        end
        text!(ax, Point2f(0.45 * flesh, -0.3 * flesh); text="flesh", fontsize=12, align=(:center, :center))
        has_fat && text!(ax, Point2f(1.25 * ins * cosd(20), 1.25 * ins * sind(20)); text="fat", fontsize=12,
            align=(:left, :center), offset=(3, 0))
        has_fat && lines!(ax, [Point2f((flesh + skin) / 2 * cosd(20), (flesh + skin) / 2 * sind(20)),
            Point2f(1.25 * ins * cosd(20), 1.25 * ins * sind(20))]; color=:black, linewidth=0.7)
        text!(ax, Point2f(1.25 * ins * cosd(-10), 1.25 * ins * sind(-10)); text="skin", fontsize=12,
            align=(:left, :center), offset=(3, 0))
        lines!(ax, [Point2f(skin * cosd(-10), skin * sind(-10)), Point2f(1.25 * ins * cosd(-10), 1.25 * ins * sind(-10))];
            color=:black, linewidth=0.7)
        has_fur && text!(ax, Point2f(1.25 * ins * cosd(-40), 1.25 * ins * sind(-40)); text="fibres", fontsize=12,
            align=(:left, :center), offset=(3, 0))
        has_fur && lines!(ax, [Point2f((skin + ins) / 2 * cosd(-40), (skin + ins) / 2 * sind(-40)),
            Point2f(1.25 * ins * cosd(-40), 1.25 * ins * sind(-40))]; color=:black, linewidth=0.7)
    end
    limits!(ax, -2.6 * ins, 1.9 * ins, -1.3 * ins, 1.5 * ins)
    hidedecorations!(ax); hidespines!(ax)
    return ax
end

"""
    layer_diagram(label => body, ...; size)

A row of `layer_diagram!` drawings.
"""
function layer_diagram(items::Pair...; size=(360 * length(items), 280), kw...)
    fig = Figure(; size)
    for (i, (label, body)) in enumerate(items)
        ax = Axis(fig[1, i]; aspect=DataAspect(), title=label, titlesize=13)
        layer_diagram!(ax, body; kw...)
    end
    return fig
end

_ellipse(a, b; n=200) = [Point2f(a * cos(t), b * sin(t)) for t in range(0, 2π; length=n)]

# Outline of each layer in a section along the long axis, which is drawn left to right.
function _long_outlines(body::Body{<:Union{Cylinder,Cone}})
    (; L, pad, taper, rf, rs, ri) = EXT._axial_layers(body, 100.0)
    box(r, z0, z1) = [Point2f(z0, -r * taper(z0)), Point2f(z1, -r * taper(z1)), Point2f(z1, r * taper(z1)),
                      Point2f(z0, r * taper(z0)), Point2f(z0, -r * taper(z0))]
    (flesh = box(rf, 0.0, L), skin = box(rs, 0.0, L), fur = box(ri, -pad, L + pad))
end
function _long_outlines(body::Body{<:Union{Sphere,Ellipsoid}})
    ratio = body.shape isa Ellipsoid ? Float64(body.shape.axis_ratio_b) : 1.0
    bf, bs, bi = _cm(flesh_radius(body)), _cm(skin_radius(body)), _cm(insulation_radius(body))
    af = ratio * bf; as = af + (bs - bf); ai = as + (bi - bs)
    (flesh = _ellipse(af, bf), skin = _ellipse(as, bs), fur = _ellipse(ai, bi))
end

"""
    long_section!(ax, body)

The section of `body` along its long axis, with its layers.
"""
function long_section!(ax, body::Body)
    outlines = _long_outlines(body)
    has_fat = skin_radius(body) > flesh_radius(body) * (1 + 1e-9)
    has_fur = insulation_radius(body) > skin_radius(body) * (1 + 1e-9)
    has_fur && poly!(ax, outlines.fur; color=LAYER_COLOURS.fur, strokecolor=:grey30, strokewidth=1)
    has_fat && poly!(ax, outlines.skin; color=LAYER_COLOURS.fat)
    poly!(ax, outlines.flesh; color=LAYER_COLOURS.flesh)
    lines!(ax, outlines.skin; color=:black, linewidth=2)
    hidedecorations!(ax); hidespines!(ax)
    return ax
end

"""
    layer_sections(label => body, ...; size)

For each body, its section across the long axis above its section along it.
"""
function layer_sections(items::Pair...; size=(250 * length(items) + 60, 420))
    fig = Figure(; size)
    Label(fig[1, 0], "Across"; rotation=π / 2, fontsize=13, font=:bold, tellheight=false)
    Label(fig[2, 0], "Along"; rotation=π / 2, fontsize=13, font=:bold, tellheight=false)
    for (i, (label, body)) in enumerate(items)
        across = Axis(fig[1, i]; aspect=DataAspect(), title=label, titlesize=13)
        layer_diagram!(across, body; labels=false)
        r = 1.15 * _cm(insulation_radius(body))
        limits!(across, -r, r, -r, r)
        along = Axis(fig[2, i]; aspect=DataAspect())
        long_section!(along, body)
    end
    Legend(fig[3, 1:length(items)], [PolyElement(color=LAYER_COLOURS.flesh), PolyElement(color=LAYER_COLOURS.fat),
        LineElement(color=:black, linewidth=2), PolyElement(color=LAYER_COLOURS.fur)],
        ["flesh", "fat", "skin", "fibres"]; orientation=:horizontal, framevisible=false, labelsize=12)
    return fig
end

# ── Graphs ───────────────────────────────────────────────────────────────────

# Tree layout: depth to the right, leaves stacked downwards, each parent centred on its children.
function _tree_layout(names, edges)
    root = first(names)
    neighbours = Dict(n => Symbol[] for n in names)
    for (p, c) in edges
        push!(neighbours[p], c); push!(neighbours[c], p)
    end
    position = Dict{Symbol,Point2f}(); parent = Dict{Symbol,Symbol}()
    next_leaf = Ref(0.0)
    function place(node, depth)
        children = [n for n in neighbours[node] if n != get(parent, node, :_none_)]
        if isempty(children)
            position[node] = Point2f(depth, -next_leaf[]); next_leaf[] += 1
        else
            for child in children
                parent[child] = node
                place(child, depth + 1)
            end
            position[node] = Point2f(depth, sum(position[c][2] for c in children) / length(children))
        end
    end
    place(root, 0)
    return position
end

"""
    body_graph!(ax, body; colors=PART_COLOURS, edge_labels=true)

Draw the parts of `body` as nodes and its joins as edges, root on the left. Each edge is labelled with the
surfaces it joins, parent side first.
"""
function body_graph!(ax, body::CompositeBody; colors=PART_COLOURS, edge_labels=true, markersize=44)
    names = collect(propertynames(body.parts))
    edges = [join_partners(j) for j in body.joins]
    position = _tree_layout(names, edges)
    for (j, (p, c)) in zip(body.joins, edges)
        a, b = position[p], position[c]
        lines!(ax, [a, b]; color=:grey40, linewidth=2)
        if edge_labels
            label = _label(j.parent_attachment.location) * " – " * _label(j.child_attachment.location)
            for halo in (true, false)   # white behind the text, so the edge does not run through it
                text!(ax, a + 0.5f0 * (b - a); text=label, fontsize=10, align=(:center, :center),
                    color=halo ? :white : :grey20, strokecolor=:white, strokewidth=halo ? 5 : 0)
            end
        end
    end
    for (i, name) in enumerate(names)
        col = colors isa NamedTuple ? getfield(colors, name) : colors[mod1(i, length(colors))]
        scatter!(ax, [position[name]]; color=col, markersize, strokecolor=i == 1 ? :black : :white,
            strokewidth=i == 1 ? 3 : 1)
        text!(ax, position[name]; text=string(name), fontsize=11, font=i == 1 ? :bold : :regular,
            align=(:center, :top), offset=(0, -markersize / 2 - 2))
    end
    xs = [p[1] for p in values(position)]; ys = [p[2] for p in values(position)]
    limits!(ax, minimum(xs) - 0.5, maximum(xs) + 0.5, minimum(ys) - 0.8, maximum(ys) + 0.6)
    hidedecorations!(ax); hidespines!(ax)
    return ax
end

"""
    body_graph(body; size, kw...)

A `Figure` with the body on the left and its graph on the right, in matching colours.
"""
function body_graph(body::CompositeBody; size=(760, 340), azimuth=1.25π, elevation=π / 7, kw...)
    fig = Figure(; size)
    ax3 = body_axis(fig[1, 1]; decorations=false, azimuth, elevation)
    draw_parts!(ax3, body)
    ax = Axis(fig[1, 2])
    body_graph!(ax, body; kw...)
    colsize!(fig.layout, 1, Relative(0.45))
    return fig
end

# ── Silhouettes ──────────────────────────────────────────────────────────────

"""
    silhouette_panel!(ax, body, direction; resolution=300)

Draw the silhouette of `body` seen from `direction` and return its area.
"""
function silhouette_panel!(ax, body, direction; resolution=300, color=RGBf(0.2, 0.2, 0.25))
    result = silhouette_rasterized(single(body), direction; resolution, return_image=true)
    xs = range(result.x_range[1] * 100, result.x_range[2] * 100; length=resolution)
    ys = range(result.y_range[1] * 100, result.y_range[2] * 100; length=resolution)
    heatmap!(ax, xs, ys, Float32.(result.bitmap); colormap=[RGBAf(1, 1, 1, 0), color], colorrange=(0, 1))
    ax.aspect = DataAspect()
    return result.area
end


# ── Parts coloured by a quantity ─────────────────────────────────────────────

"""
    temperature_views(body, temperatures; views=(:oblique, :front), label, colormap, colorrange, size)

The parts of `body` from several directions, each coloured by its value in the `NamedTuple` `temperatures`
(in °C), with a colour bar.
"""
function temperature_views(body, temperatures::NamedTuple; views=(:oblique, :front), titles=string.(views),
                           label="Skin temperature (°C)", colormap=:plasma,
                           colorrange=extrema(values(temperatures)), size=(230 * length(views) + 110, 320))
    angles = (oblique = (1.25π, π / 7), three_quarter = (-0.25π, π / 7), side = (-π / 2, 0.0), front = (0.0, 0.0),
              top = (-π / 2, π / 2), back = (π, 0.0))
    scheme = Makie.to_colormap(colormap)
    lo, hi = colorrange
    colour(t) = scheme[clamp(round(Int, 1 + (t - lo) / max(hi - lo, eps()) * (length(scheme) - 1)), 1, length(scheme))]
    colors = map(colour, temperatures)
    fig = Figure(; size)
    for (i, view) in enumerate(views)
        azimuth, elevation = getfield(angles, view)
        ax = body_axis(fig[1, i]; decorations=false, azimuth, elevation, title=titles[i], titlesize=12)
        draw_parts!(ax, body; colors)
    end
    Colorbar(fig[1, length(views) + 1]; colormap, limits=(lo, hi), label)
    return fig
end

end

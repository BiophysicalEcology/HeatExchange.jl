const _W = 26

# Strip Param wrapper if present, reconstructing units from metadata when needed.
_val(x) = x
function _val(p::ModelParameters.Param)
    v = p.val
    NT = typeof(p).parameters[2]  # the NamedTuple type holding metadata
    :units ∉ fieldnames(NT) && return v
    u = p.units
    isnothing(u) ? v : v * u
end

# Generic field-by-field printer for any parameter struct.
function _show_fields(io, obj)
    fnames = fieldnames(typeof(obj))
    n = length(fnames)
    for (i, fname) in enumerate(fnames)
        v = _val(getfield(obj, fname))
        label = "  $(fname):"
        i < n ? println(io, rpad(label, _W), v) : print(io, rpad(label, _W), v)
    end
end

# ── Param-based parameter structs ─────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", p::AbstractMorphologyParameters)
    println(io, nameof(typeof(p)))
    _show_fields(io, p)
end

function Base.show(io::IO, ::MIME"text/plain", p::AbstractPhysiologyParameters)
    println(io, nameof(typeof(p)))
    _show_fields(io, p)
end

function Base.show(io::IO, ::MIME"text/plain", p::AbstractModelParameters)
    println(io, nameof(typeof(p)))
    _show_fields(io, p)
end

function Base.show(io::IO, ::MIME"text/plain", p::AbstractEnvironmentalPars)
    println(io, nameof(typeof(p)))
    _show_fields(io, p)
end

# ── FibreProperties ───────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", f::FibreProperties)
    println(io, "FibreProperties")
    _show_fields(io, f)
end

# ── InsulationParameters ──────────────────────────────────────────────────────
# InsulationParameters <: AbstractMorphologyParameters but needs nested display.

function Base.show(io::IO, ::MIME"text/plain", p::InsulationParameters)
    println(io, "InsulationParameters")
    print(io, "  Dorsal:  ")
    show(io, MIME"text/plain"(), p.dorsal)
    println(io)
    print(io, "  Ventral: ")
    show(io, MIME"text/plain"(), p.ventral)
    println(io)
    println(io, rpad("  depth_compressed:", _W), _val(p.depth_compressed))
    print(io,   rpad("  longwave_depth_fraction:", _W), _val(p.longwave_depth_fraction))
end

# ── EnvironmentalVars ─────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", e::EnvironmentalVars)
    println(io, "EnvironmentalVars")
    _show_fields(io, e)
end

# ── HeatExchangeTraits ────────────────────────────────────────────────────────

function Base.show(io::IO, mime::MIME"text/plain", t::HeatExchangeTraits)
    header = "HeatExchangeTraits"
    println(io, header)
    println(io, "─" ^ length(header))
    print(io, "Shape:      ")
    show(io, mime, t.shape_pars)
    println(io)
    println(io)
    print(io, "Insulation: ")
    show(io, mime, t.insulation_pars)
    println(io)
    println(io)
    for fname in (
        :conduction_pars_external, :conduction_pars_internal, :radiation_pars,
        :convection_pars, :evaporation_pars, :hydraulic_pars,
        :respiration_pars, :metabolism_pars, :options,
    )
        component = getfield(t, fname)
        println(io, fname, ":")
        _show_fields(io, component)
        println(io)
    end
end

"""
    SmoothingStrategy

Policy for handling kinks (`abs`, `max(x, 0)`, branchy `> 0` masks, `clamp`) so that
forward simulation gets exact piecewise behaviour while automatic differentiation gets
well-defined gradients.

Concrete subtypes:
- [`HardBound`](@ref): exact `abs`, `max`, `clamp`, ternary masks.
- [`SmoothBound`](@ref): replace each kink with a `sqrt(x² + ε²)`-based smoothing.

Generic helpers dispatching on the strategy:
- [`safe_abs`](@ref) — smooth `|x|`
- [`safe_relu`](@ref) — smooth `max(x, 0)`
- [`safe_step`](@ref) — smooth Heaviside step at 0 (returns ≈1 for x>0, ≈0 for x<0)
- [`safe_max`](@ref) / [`safe_min`](@ref) / [`safe_clamp`](@ref)
"""
abstract type SmoothingStrategy end

"""
    HardBound() <: SmoothingStrategy

Use exact `abs`, `max`, `clamp` and ternary masks. Correct values, fast, but introduces
gradient kinks at singular points that break Enzyme's reverse-mode pass.
"""
struct HardBound <: SmoothingStrategy end

"""
    SmoothBound(ε) <: SmoothingStrategy

`ε` is a **dimensionless relative tolerance** in roughly `[0, 1]`. Each helper accepts
an optional `scale` keyword giving the natural magnitude of the input at that call
site; the actual smoothing window has size `ε · scale`. So `ε` is a single knob the
user tunes (default `1e-3`) while every call site bakes in its own physical scale.

Sweet spot: `ε = 1e-3` to `1e-2`. At `ε ≪ 1e-6` the smoothing collapses below Float64
precision near the kink; at `ε ≳ 1e-1` the smoothing window starts to bias
`safe_step` / `safe_relu` outputs noticeably away from `{0, 1}`.
"""
struct SmoothBound{T<:Real} <: SmoothingStrategy
    ε::T
end

SmoothBound() = SmoothBound(1.0e-3)

# --- helpers ---------------------------------------------------------------
#
# Each helper takes `scale` as a kwarg with default `oneunit(x)` (or `oneunit(a-b)`
# for two-arg helpers). Pass an explicit `scale` when the input's natural magnitude
# differs from `oneunit` (e.g. an unboxed Float64 wetness in [0, 0.05]; a Quantity
# in metres whose typical nonzero value is 1e-4 m). Smoothing window = ε · scale.

"""
    safe_abs(s::SmoothingStrategy, x; scale=oneunit(x))

`abs(x)` (HardBound) or `sqrt(x² + (ε·scale)²)` (SmoothBound).
The smoothing only deviates from `|x|` when `|x| ≲ ε·scale`, so `safe_abs` is
robust to a small mismatch between `scale` and the actual magnitude of `x`.
"""
safe_abs(::HardBound, x; scale=oneunit(x))    = abs(x)
safe_abs(s::SmoothBound, x; scale=oneunit(x)) = sqrt(x*x + (s.ε * scale)^2)

"""
    safe_relu(s::SmoothingStrategy, x; scale=oneunit(x))

`max(x, 0)` (HardBound) or `(x + safe_abs(s, x; scale)) / 2` (SmoothBound).
"""
safe_relu(::HardBound, x; scale=oneunit(x))    = max(x, zero(x))
safe_relu(s::SmoothBound, x; scale=oneunit(x)) = (x + safe_abs(s, x; scale)) / 2

"""
    safe_step(s::SmoothingStrategy, x; scale=oneunit(x))

Heaviside step at 0 — `1.0` if `x > 0` else `0.0` (HardBound), or a smooth
transition centred on `x = 0` with window width `ε·scale` (SmoothBound).
Pick `scale` so that `ε·scale` is well below the smallest "meaningfully positive"
value of `x` you expect.
"""
safe_step(::HardBound, x; scale=oneunit(x))    = x > zero(x) ? 1.0 : 0.0
safe_step(s::SmoothBound, x; scale=oneunit(x)) = (1.0 + x / safe_abs(s, x; scale)) / 2

"""
    safe_max(s::SmoothingStrategy, a, b; scale=oneunit(a-b))

`max(a, b)` (HardBound) or `(a + b + safe_abs(s, a - b; scale)) / 2` (SmoothBound).
"""
safe_max(::HardBound, a, b; scale=oneunit(a-b))    = max(a, b)
safe_max(s::SmoothBound, a, b; scale=oneunit(a-b)) = (a + b + safe_abs(s, a - b; scale)) / 2

"""
    safe_min(s::SmoothingStrategy, a, b; scale=oneunit(a-b))

`min(a, b)` (HardBound) or `(a + b - safe_abs(s, a - b; scale)) / 2` (SmoothBound).
"""
safe_min(::HardBound, a, b; scale=oneunit(a-b))    = min(a, b)
safe_min(s::SmoothBound, a, b; scale=oneunit(a-b)) = (a + b - safe_abs(s, a - b; scale)) / 2

"""
    safe_clamp(s::SmoothingStrategy, x, lo, hi; scale=oneunit(x))

`clamp(x, lo, hi)` (HardBound) or smoothed via `safe_max(s, lo, safe_min(s, hi, x))`.
"""
safe_clamp(s::SmoothingStrategy, x, lo, hi; scale=oneunit(x)) =
    safe_max(s, lo, safe_min(s, hi, x; scale); scale)

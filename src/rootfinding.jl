# Re-export Roots.jl solvers so users of HeatExchange can choose any solver
# without an additional Roots dependency.
using Roots: find_zero, Bisection, A42, AlefeldPotraShi, FalsePosition, Order0

export find_zero, Bisection, A42, AlefeldPotraShi, FalsePosition, Order0

"""
    zbrac(f, a, b; factor=1.6, maxiter=50) → (a, b)

Expand the interval `[a, b]` until `f(a)` and `f(b)` have opposite signs (i.e.
bracket a root), returning the bracketing pair.

Ported from NicheMapR `ZBRAC.f` (Brent 2002, Dover).  At each iteration the
end-point with the *smaller* function magnitude is moved outward by `factor`,
which is more efficient than symmetric expansion.

Works with plain `Float64` values.  Throws an error if no bracket is found
within `maxiter` iterations.
"""
function zbrac(f, a::Real, b::Real; factor::Real=1.6, maxiter::Int=50)
    fa = f(a)
    fb = f(b)
    for _ in 1:maxiter
        fa * fb < 0 && return (a, b)
        if abs(fa) < abs(fb)
            a  += factor * (a - b)
            fa  = f(a)
        else
            b  += factor * (b - a)
            fb  = f(b)
        end
    end
    error("zbrac: failed to bracket root in $maxiter iterations")
end

"""
    zbrent(f, a, b, tol; maxiter=300, eps=3e-8) → root

Brent's method root-finder for `f(x) = 0` on `[a, b]`, where `f(a)` and
`f(b)` must have opposite signs.

Ported from NicheMapR `ZBRENT.f` (Brent 2002, Dover), with the Porter (2003)
early-exit when `|f(b)| ≤ tol` after the first iteration.

Works with plain `Float64` values.  For Unitful quantities, strip units before
calling and restore them after:

```julia
T_K = zbrent(T -> ustrip(u"W", heat_balance(T * u"K")), 273.15, 343.15, 1e-3) * u"K"
```

Returns `b` if convergence is not achieved within `maxiter` iterations.
"""
function zbrent(f, a::Real, b::Real, tol::Real; maxiter::Int=300, eps::Float64=3e-8)
    qa = f(a)
    qb = f(b)

    c  = zero(b)
    d  = zero(b)
    e  = zero(b)
    qc = qb

    @inbounds for i in 1:maxiter
        if qb * qc > 0
            c  = a;   qc = qa
            d  = b - a
            e  = d
        end
        if abs(qc) < abs(qb)
            a = b;  b = c;  c = a
            qa = qb;  qb = qc;  qc = qa
        end

        tol1 = 2eps * abs(b) + tol / 2
        xm   = (c - b) / 2

        (abs(xm) ≤ tol1 || qb == 0) && return b
        abs(qb) ≤ tol1 && i > 1     && return b

        if abs(e) ≥ tol1 && abs(qa) > abs(qb)
            s = qb / qa
            if a == c
                p = 2xm * s
                q = 1 - s
            else
                q = qa / qc
                r = qb / qc
                p = s * (2xm * q * (q - r) - (b - a) * (r - 1))
                q = (q - 1) * (r - 1) * (s - 1)
            end
            p > 0 && (q = -q)
            p = abs(p)
            if 2p < min(3xm * q - abs(tol1 * q), abs(e * q))
                e = d;  d = p / q
            else
                d = xm;  e = d
            end
        else
            d = xm;  e = d
        end

        a = b;  qa = qb
        b += abs(d) > tol1 ? d : sign(xm) * tol1
        qb = f(b)
    end

    return b
end

# =============================================================================
#  2.  Full thermodynamics
# =============================================================================

"""
    thermodynamics!(m; N = 60, nphi = 400, spread = (-2, 5), lambda = 0.4,
                    tol = 1e-6, maxit = 20_000, progress = true, findTc = true)

Sweep phiH from ~0 up to `phiMax` by continuation and store the thermodynamic
curve in `m.thermo` as a NamedTuple of vectors

    phiH, T, S, Lambda, t = T/Lambda, s = S/Lambda^3, p = (T,s) free energy,
    converged

Points where Newton fails are kept as `NaN` (never silently dropped, so the
index still lines up with `phiH`).  `p` is obtained from dp = s dT integrated
along the branch, which is what makes the swallowtail visible.  With
`findTc = true` the first-order critical temperature is located from the
self-intersection of the (T, p) curve and stored in `m.Tc`.

The phiH grid is `phiMax*(tanh(range(spread...)) + 1)/2`: clustered near
`phiMax` where the branch turns, without collapsing points onto it.
"""
function thermodynamics!(m::Bea_Mateos; N::Int = 60, nphi::Int = 400,
                         spread = (-2, 5), lambda::Real = 0.4, tol::Real = 1e-6,
                         maxit::Int = 20_000, progress::Bool = true,
                         findTc::Bool = true)
    x, D, D2 = cheby(N)
    phiH_list = m.phiMax .* (tanh.(range(spread[1], spread[2], length = nphi)) .+ 1) ./ 2

    T = fill(NaN, nphi); S = fill(NaN, nphi); Lam = fill(NaN, nphi)
    conv = falses(nphi)
    F0 = vcat(1 .- x.^4, ones(N), phiH_list[1] .* x)

    iter = progress ? Progress(nphi; desc = "thermodynamics: ") : nothing
    for (i, phiH) in enumerate(phiH_list)
        f, g, ph, res = _newton(m, x, D, D2, phiH, F0;
                                lambda = lambda, tol = tol, maxit = maxit)
        cand = vcat(f, g, ph)
        ok = isfinite(res) && res < tol && !any(isnan, cand) && g[end] > 0
        if ok
            T[i]   = -1/(4π) * (D[end:end, :] * f)[1] * sqrt(g[end])
            S[i]   = π * g[end]^(3/2)
            Lam[i] = ph[1]
            conv[i] = true
            F0 = cand                      # advance the guess only on a good solve
        end
        progress && next!(iter)
    end
    nbad = count(!, conv)
    nbad > 0 && @info "thermodynamics!: $nbad of $nphi points did not converge (stored as NaN)"

    t = T ./ Lam
    s = S ./ Lam.^3
    p = _free_energy(t, s)                  # p/Lambda^4
    F = -p                                  # free energy density, F = -p
    e = s .* t .- p                         # energy density, eps + p = sT

    m.thermo = (phiH = phiH_list, T = T, S = S, Lambda = Lam, t = t, s = s,
                p = p, F = F, e = e, converged = conv, N = N)

    if findTc
        try
            m.Tc = critical_temperature(s, t, phiH_list)
            @info "critical temperature" Tc = m.Tc
        catch err
            m.Tc = nothing
            @warn "Tc not found (the sampled branch may not close a swallowtail)" err
        end
    end
    return m.thermo
end

"""
    _free_energy(t, s) -> p

Pressure along the branch from `dp = s dT`, integrated along the sweep with the
trapezoid rule; the free energy is `F = -p` and the energy density
`eps = sT - p`.

`p` is fixed to zero at the *coldest* point of the sweep.  Two remarks on that:

  * Only an additive constant is undetermined -- pinning it properly needs the
    renormalized on-shell action (holographic counterterms), which `(T, s)`
    alone cannot give.  So `F` here is a free energy up to a constant, which is
    all the swallowtail and the phase transition depend on.
  * The anchor is deliberately at the cold end rather than the hot end.
    Anchoring at high T with the conformal value `p = sT/4` looks more physical
    but is useless numerically: `p ~ T^4`, so the integration path down to Tc
    spans many orders of magnitude and the accumulated quadrature error swamps
    the O(0.1) structure near the transition.  Measured: `F` at fixed
    `T/Lambda = 0.3923` came out as -23.6, +2.7e5, +1.0e9, +3.8e12 as the hot
    end was pushed from `T/Lambda = 8` to 3400, while the spread of `F` across
    the transition window stayed at 0.16-0.19 throughout.
"""
function _free_energy(t, s)
    p = fill(NaN, length(t))
    idx = [i for i in eachindex(t) if isfinite(t[i]) && isfinite(s[i])]
    isempty(idx) && return p
    p[idx[end]] = 0.0                          # coldest sampled point
    for k in (length(idx)-1):-1:1
        i, j = idx[k], idx[k+1]
        p[i] = p[j] + 0.5 * (s[i] + s[j]) * (t[i] - t[j])
    end
    return p
end

# --- swallowtail Tc finder -----------------------------------------------------

function _seg(x1, y1, x2, y2, x3, y3, x4, y4)
    d = (x2-x1)*(y4-y3) - (y2-y1)*(x4-x3)
    iszero(d) && return nothing
    t = ((x3-x1)*(y4-y3) - (y3-y1)*(x4-x3)) / d
    u = ((x3-x1)*(y2-y1) - (y3-y1)*(x2-x1)) / d
    (0 <= t <= 1 && 0 <= u <= 1) ? (t, u) : nothing
end

"""
    critical_temperature(s, T, phiH; iters = 12, tol = 1e-14, mingap_rel = 1e-4,
                         refine = true)

First-order critical temperature: the self-intersection of the (T, p) curve,
with p from s dT.  Points are thinned to a physical spacing first, because a
tanh-clustered phiH grid piles thousands of near-degenerate knots at phiMax that
both ill-condition the spline and flood the seed search with fake crossings.

With `refine = false` the crossing of the piecewise-linear (T, p) curve is
returned as is, with p from the trapezoid rule -- no splines, no Newton step.
"""
function critical_temperature(s::AbstractVector, T::AbstractVector, phiH::AbstractVector;
                              iters::Int = 12, tol::Real = 1e-14, mingap_rel::Real = 1e-4,
                              refine::Bool = true)
    @assert length(T) == length(s) == length(phiH) "s, T, phiH must have equal length"

    o = sortperm(phiH)
    φ = collect(float.(phiH[o])); Tv = collect(float.(T[o])); sv = collect(float.(s[o]))
    fin = isfinite.(φ) .& isfinite.(Tv) .& isfinite.(sv)
    φ, Tv, sv = φ[fin], Tv[fin], sv[fin]

    span = maximum(φ) - minimum(φ)
    keepidx = [1]
    for i in 2:length(φ)
        (φ[i] - φ[keepidx[end]]) >= mingap_rel * span && push!(keepidx, i)
    end
    φ, Tv, sv = φ[keepidx], Tv[keepidx], sv[keepidx]

    n = length(φ)
    n < 4 && error("Too few usable points after cleaning ($n).")

    # p at the nodes.  With `refine` this is the exact per-segment quadrature of
    # s dT/dphi on cubic splines; otherwise the trapezoid rule on the same nodes.
    p = zeros(eltype(Tv), n)
    if refine
        Ts, ss = CubicSpline(Tv, φ), CubicSpline(sv, φ)
        Tp(z) = DataInterpolations.derivative(Ts, z)
        for i in 2:n
            p[i] = p[i-1] + first(quadgk(z -> ss(z) * Tp(z), φ[i-1], φ[i]))
        end
    else
        for i in 2:n
            p[i] = p[i-1] + 0.5 * (sv[i] + sv[i-1]) * (Tv[i] - Tv[i-1])
        end
    end

    # Coarse seed: piecewise-linear self-intersection of (T, p).  Take the FIRST
    # (lowest-phi) crossing -- that is the physical swallowtail; residual noise
    # sits at higher phi.
    seed = nothing; Tseed = NaN
    @inbounds for i in 1:n-1, j in i+2:n-1
        r = _seg(Tv[i], p[i], Tv[i+1], p[i+1], Tv[j], p[j], Tv[j+1], p[j+1])
        if r !== nothing
            t, u = r
            seed = (φ[i] + t*(φ[i+1]-φ[i]), φ[j] + u*(φ[j+1]-φ[j]))
            Tseed = Tv[i] + t*(Tv[i+1]-Tv[i])
            break
        end
    end
    seed === nothing && error("No self-intersection: the branch does not close a swallowtail.")
    refine || return Tseed

    # Refine: solve T(a) = T(b), p(a) = p(b) by 2-D Newton, with
    # J = [T'(a) -T'(b); s(a)T'(a) -s(b)T'(b)].
    ps = CubicSpline(p, φ)
    a, b = seed
    for _ in 1:iters
        g1, g2 = Ts(a) - Ts(b), ps(a) - ps(b)
        (abs(g1) < tol && abs(g2) < tol) && break
        Ta, Tb = Tp(a), Tp(b)
        det = Ta*(-ss(b)*Tb) - (-Tb)*(ss(a)*Ta)
        iszero(det) && break
        a -= ((-ss(b)*Tb)*g1 - (-Tb)*g2) / det
        b -= (-(ss(a)*Ta)*g1 + Ta*g2) / det
    end
    return Ts(a)
end

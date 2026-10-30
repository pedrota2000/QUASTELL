# =============================================================================
#  4.  Quasinormal modes
# =============================================================================

"""
    singular_band(m, b, q) -> (lo, hi)

Range of `|omega|/2piT` over which the *principal symbol* of the master system
vanishes somewhere in `0 <= u <= 1`.

The coefficient of the second derivative is `A0(u) + omega^2 A2(u)`, and it has
a zero at `omega^2 = a(u)` with

    a(u) = k^2 [ f (3(2g - u g')^2 + 4 chi u^2 g^2 phi'^2) + 24 g^2 - 8 chi g^2 V ]
           / ( 9 (2g - u g')^2 ),        k = 2 pi q T .

For `omega` in this band the ODE system acquires an *interior* singular point,
and collocation responds by producing one essentially-real eigenvalue pair per
grid point there.  Those are artifacts of the master variables, not modes: they
do not move when N changes, so no convergence test can remove them, but the band
is known analytically and can simply be excluded.

At pure AdS (`f = 1-u^4`, `g = 1`, `phi = 0`) this reduces to
`a/k^2 = 1 - u^4/3`, so the band is `|lambda|/q` in `[sqrt(2/3), 1]`.
"""
function principal_zeros(m::Bea_Mateos, b::Background, q::Real, T::Real)
    u = b.x; f = b.f0; g = b.g0; dg = b.D * b.g0
    dphi = m.nu .* u.^(m.nu - 1) .* b.phi0 .+ u.^m.nu .* (b.D * b.phi0)
    VV = m.V.(u.^m.nu .* b.phi0)
    G = (2 .* g .- u .* dg).^2
    k2 = (2π * q * T)^2
    return k2 .* (f .* (3 .* G .+ 4 .* m.chi .* u.^2 .* g.^2 .* dphi.^2) .+
                  24 .* g.^2 .- 8 .* m.chi .* g.^2 .* VV) ./ (9 .* G)   # = a(u_i), in omega^2
end

function singular_band(m::Bea_Mateos, b::Background, q::Real, T::Real)
    a = principal_zeros(m, b, q, T)
    s = sqrt.(max.(a, 0.0)) ./ (2π * T)           # in units of omega/2piT
    return minimum(s), maximum(s)
end

"""
    _spectrum(m, phiH, q, N; cutoff) -> omegas, T

Raw spectrum on one grid: every finite eigenvalue of the companion pencil with
|omega| < cutoff.  Most of these are numerical artifacts.
"""
function _spectrum(m::Bea_Mateos, phiH::Real, q::Real, N::Int; cutoff::Real = Inf,
                   balance::Bool = true, rescale = (2, 0), divide_f::Bool = true,
                   bcexpand::Bool = false)
    _, _, Mt, b = mtildes(m, phiH; q = q, N = N, rescale = rescale, divide_f = divide_f,
                          bcexpand = bcexpand)
    sc = 2π * b.T

    if balance
        # Two exact equivalences that make the pencil far better conditioned.
        # (i) solve for lambda = omega/2piT instead of omega: Mt[n] carries wildly
        #     different powers of T (Mt[1] ~ q^2 T^2, Mt[4] ~ T^0), which spreads
        #     the pencil over many orders of magnitude;
        # (ii) divide each row by its largest entry across all Mt[n] -- a left
        #     diagonal scaling, so the spectrum is untouched.
        A = [ComplexF64.(Mt[n]) .* sc^(n-1) for n in eachindex(Mt)]
        for i in axes(A[1], 1)
            s = maximum(maximum(abs, view(X, i, :)) for X in A)
            s > 0 && for X in A
                @views X[i, :] ./= s
            end
        end
        M0, M1 = companion(A)
        w = eigvals(M0, -M1) .* sc               # back to omega
    else
        M0, M1 = companion(Mt)
        w = eigvals(M0, -M1)
    end

    keep = [z for z in w if isfinite(z) && abs(z) < cutoff]
    return ComplexF64.(keep), b.T
end

"""
    qnms(m, phiH; q = 0.0, Ns = (30, 40, 50), rtol = 1e-4, cutoff = Inf,
         verbose = true) -> QNMResult

Quasinormal frequencies at horizon value `phiH` and dimensionless momentum `q`.

The spectrum is computed independently on **every** grid size in `Ns` — the
background is re-solved on each one, never interpolated, so the only thing that
changes between runs is the resolution.  A frequency is accepted only if it
appears, to within `rtol` relative distance, in all of them.  Physical modes sit
still under a change of N; artifacts move or vanish.  This is the convergence
test of Jansen, arXiv:1709.09178, eq. (3.2), and with three grids instead of two
it also protects against an artifact that happens to survive one comparison.

Results are sorted by decreasing `imag(lambda)` (least damped first): any mode
with `imag > 0` is an instability.  `r.lambda = r.omega / 2 pi T` is the
dimensionless frequency; `r.err` is the residual drift between the two finest
grids.  The run is appended to `m.qnm_runs`.
"""
function qnms(m::Bea_Mateos, phiH::Real; q::Real = 0.0, Ns = (30, 40, 50),
              rtol::Real = 1e-4, cutoff::Real = Inf, verbose::Bool = true,
              drop_band::Bool = true, bandtol::Real = 1e-2, origintol::Real = 1e-3,
              rescale = (2, 0), divide_f::Bool = true, progress::Bool = false,
              bcexpand::Bool = false)
    Nlist = sort(collect(Int, Ns))
    length(Nlist) >= 2 || error("give at least two grid sizes to separate physical from spurious modes")

    raw = Vector{Vector{ComplexF64}}(undef, length(Nlist))
    Tref = NaN
    bar = progress ? Progress(length(Nlist); desc = "qnms: ") : nothing
    for (i, N) in enumerate(Nlist)
        raw[i], Tref = _spectrum(m, phiH, q, N; cutoff = cutoff, rescale = rescale,
                                 bcexpand = bcexpand,
                                 divide_f = divide_f)
        verbose && @printf("  N = %3d : %4d finite eigenvalues\n", N, length(raw[i]))
        progress && _tick!(bar, [(:N, N), (:eigenvalues, length(raw[i]))])
    end

    # Artifacts of the master variables: essentially-real eigenvalues sitting in
    # the band where the principal symbol degenerates (see `singular_band`).
    # They survive the resolution test, so they have to go by their location.
    if drop_band
        ndrop = 0
        if q == 0
            # At k = 0 the band a(u) ~ k^2 collapses onto omega = 0, so k^2 is no
            # yardstick and every artifact piles up at the origin (362 of them at
            # N up to 90) and passes the resolution test.  A fixed small cut does
            # the job: measured at N = 90, the artifacts all sit at
            # |lambda| <= 2.5e-7 while the next eigenvalue is at 0.945, so
            # origintol = 1e-3 sits in a six-order-wide gap.
            cut = origintol * 2π * Tref
            for i in eachindex(Nlist)
                ndrop += count(w -> abs(w) < cut, raw[i])
                raw[i] = filter(w -> abs(w) >= cut, raw[i])
            end
            @warn "q = 0: discarding the degenerate band at omega = 0 (|omega/2piT| < " *
                  "$origintol).  The gapped tower is unaffected, but the hydrodynamic " *
                  "mode also sits exactly at omega = 0 here and cannot be separated " *
                  "from it -- use a small nonzero q for that."
            verbose && @printf("  dropped %d eigenvalues at the degenerate band omega = 0\n", ndrop)
        else
            k2 = (2π * q * Tref)^2
            for i in eachindex(Nlist)
                a = principal_zeros(m, background(m, phiH; N = Nlist[i]), q, Tref)
                inband(w) = minimum(abs.(w^2 .- a)) < bandtol * k2
                ndrop += count(inband, raw[i])
                raw[i] = filter(!inband, raw[i])
            end
            if verbose
                lo, hi = singular_band(m, background(m, phiH; N = Nlist[end]), q, Tref)
                @printf("  dropped %d eigenvalues at omega^2 = a(u_i), band |lambda|/q in [%.4f, %.4f]\n",
                        ndrop, lo / q, hi / q)
            end
        end
    end

    # Candidates come from the finest grid; each must be matched on all others.
    fine = raw[end]
    others = raw[1:end-1]
    acc = ComplexF64[]; accerr = Float64[]
    for w in fine
        scale = rtol * max(1.0, abs(w))
        dists = [isempty(o) ? Inf : minimum(abs.(o .- w)) for o in others]
        all(d -> d <= scale, dists) || continue
        push!(acc, w); push!(accerr, dists[end])   # drift vs the second-finest grid
    end

    ord = sortperm(imag.(acc), rev = true)
    acc, accerr = acc[ord], accerr[ord]
    lam = acc ./ (2π * Tref)

    if verbose
        @printf("  converged on all %d grids: %d modes  (T = %.6g)\n",
                length(Nlist), length(acc), Tref)
        nunst = count(z -> imag(z) > 0, lam)
        nunst > 0 && @printf("  *** %d mode(s) with Im(omega) > 0: UNSTABLE ***\n", nunst)
    end

    r = QNMResult(float(phiH), float(q), Tref, Nlist, raw, acc, lam, accerr)
    push!(m.qnm_runs, r)
    return r
end

"""
    show_modes(r; nmax = 10)

Table of the converged dimensionless frequencies lambda = omega/2piT.
"""
function show_modes(r::QNMResult; nmax::Int = 10)
    @printf("phiH = %g, q = %g, T = %.8g, grids %s\n",
            r.phiH, r.q, r.T, string(r.Ns))
    @printf("%4s  %-22s %-22s %-10s\n", "n", "Re(w/2piT)", "Im(w/2piT)", "err(w)")
    for n in 1:min(nmax, length(r.lambda))
        @printf("%4d  %-22.14g %-22.14g %-10.3g\n",
                n, real(r.lambda[n]), imag(r.lambda[n]), r.err[n])
    end
    length(r.lambda) > nmax && @printf("...  (%d modes total)\n", length(r.lambda))
    return nothing
end

"""
    dispersion(m, phiH; qs = 0:0.25:2, nmodes = 3, kwargs...)
        -> (q = qs, lambda = Matrix, results = [QNMResult])

The `nmodes` least-damped modes as a function of momentum.  `lambda[i, n]` is
mode `n` at `qs[i]`, `NaN + NaN im` where fewer than `n` modes converged.
Extra keyword arguments are forwarded to `qnms`.
"""
function dispersion(m::Bea_Mateos, phiH::Real; qs = 0.1:0.1:1.0, nmodes::Int = 3,
                    verbose::Bool = false, progress::Bool = true, foldtol::Real = 1e-5,
                    desc::AbstractString = "dispersion: ", kwargs...)
    qv = sort(collect(float.(qs)))
    if any(iszero, qv)
        @warn "dropping q = 0: at k = 0 the omega^0 block of the master system degenerates"
        qv = filter(!iszero, qv)
    end

    lam  = fill(ComplexF64(NaN, NaN), length(qv), nmodes)
    runs = QNMResult[]
    bar  = progress ? Progress(length(qv); desc = desc) : nothing

    for (i, q) in enumerate(qv)
        r = qnms(m, phiH; q = q, verbose = verbose, kwargs...)
        push!(runs, r)
        progress && _tick!(bar, [(:q, q), (:modes, length(r.lambda))])

        # Fold the +-Re(lambda) symmetry, keeping one representative per pair, so a
        # branch is a single-valued function of q.  Purely imaginary modes have
        # Re = 0 and are all kept: they are distinct modes, not a mirror pair.
        #
        # The tolerance has to scale with the numerical noise on Re, not be a hard
        # 1e-9.  A purely imaginary mode only has Re = 0 up to the eigenvalue error,
        # and with `divide_f = true` the pencil is conditioned ~1e12, so Re lands
        # around 1e-7 rather than 1e-13.  A fixed -1e-9 threshold discards such a
        # mode -- and if that happens at the first q, the branch never initialises
        # and the continuation below tracks something else for the entire sweep,
        # which is how an inverted parabola turns into a flat damped line.  Any
        # genuine mirror pair has |Re| far larger than this, so only one of its two
        # members still survives.
        cand = sort(filter(z -> real(z) >= -foldtol * max(1.0, abs(z)), r.lambda),
                    by = z -> -imag(z))

        if i == 1
            for n in 1:min(nmodes, length(cand))
                lam[i, n] = cand[n]
            end
        else
            # Continue each branch by nearest neighbour to a linear extrapolation
            # in q, and consume the match so two branches cannot claim one mode.
            for n in 1:nmodes
                prev = lam[i-1, n]
                (isnan(real(prev)) || isempty(cand)) && continue
                pred = prev
                if i > 2 && !isnan(real(lam[i-2, n]))
                    pred = prev + (prev - lam[i-2, n]) *
                                  ((qv[i] - qv[i-1]) / (qv[i-1] - qv[i-2]))
                end
                k = argmin(abs.(cand .- pred))
                lam[i, n] = cand[k]
                deleteat!(cand, k)
            end
        end
    end
    return (q = qv, lambda = lam, results = runs)
end

"""
    sound_speed(d; nmax = 3) -> (cs, residual)

Speed of sound from the least-damped branch of a `dispersion` result, fitted over
the `nmax` smallest momenta.

The sound dispersion is `lambda = cs q - i Gamma q^2 / 2 + O(q^3)`, so the real
part is odd: `Re lambda = cs q + b q^3`.  Both terms are fitted (falling back to
the pure slope if only one momentum is available) -- fitting a straight line
instead biases `cs` by the cubic term, which at q = 0.15 is already a 0.5%
error.

Because `lambda = omega/2piT` and `q = k/2πT`, the slope *is* `c_s` in units of
the speed of light -- `1/sqrt(3) = 0.57735` for a conformal plasma.  On a branch
with negative specific heat `c_s^2 < 0`, the sound mode turns purely imaginary:
`Re lambda = 0`, so this returns ~0 and the information moves to `Im lambda`,
which then grows linearly in q.
"""
function sound_speed(d::NamedTuple; nmax::Int = 3)
    q, l = d.q, d.lambda[:, 1]
    k = [i for i in 1:min(nmax, length(q)) if !isnan(real(l[i]))]
    isempty(k) && return (NaN, NaN)
    x, y = q[k], real.(l[k])
    if length(k) == 1
        cs = y[1] / x[1]
        return (cs, 0.0)
    end
    A = hcat(x, x.^3)                    # Re lambda = cs*q + b*q^3
    c = A \ y
    return (c[1], sqrt(sum((y .- A*c).^2) / max(length(k) - 2, 1)))
end

"""
    growth_rate(m, phiH; qs = 0.05:0.1:1.2, Nbg = 60, kwargs...) -> NamedTuple

Growth rate of the unstable sound mode on a spinodal background, in the units of
the phase-transition literature: `gamma/Lambda` against `k/Lambda`, where
`gamma = Im omega`.

Inside the spinodal region `c_s^2 < 0`, so `c_s = i|c_s|` and the sound mode is
purely imaginary.  For small momenta

    gamma(k) = |c_s| k - Gamma k^2 / 2 + O(k^3) ,

an inverted parabola whose slope at the origin is `|c_s|` and which crosses zero
at `k = 2|c_s|/Gamma`: modes with `k` below that grow, above it they decay.
Returned fields: `kL`, `gL`, `TL = T/Lambda`, the fitted `cs` and `Gamma`, the
zero crossing `kmax`, and the underlying `dispersion` result `d`.

The conversion is `k/Lambda = 2 pi (T/Lambda) q` and likewise for gamma, since
`q = k/2πT`.
"""
function growth_rate(m::Bea_Mateos, phiH::Real; qs = 0.05:0.1:1.2, Nbg::Int = 60,
                     nfit::Int = 3, kwargs...)
    d = dispersion(m, phiH; qs = qs, nmodes = 1, kwargs...)
    b = background(m, phiH; N = Nbg)
    TL = b.T / b.Lambda
    kL = 2π * TL .* d.q
    gL = 2π * TL .* imag.(d.lambda[:, 1])

    # Fit only the smallest momenta: (2.16) is the *hydrodynamic* limit, and far
    # out the branch leaves the unstable regime and eventually collides with
    # another mode, which a global fit would happily absorb.
    ok = [i for i in 1:min(nfit, length(kL)) if isfinite(gL[i])]
    A = hcat(kL[ok], -0.5 .* kL[ok].^2)              # gamma = |cs| k - Gamma k^2/2
    c = length(ok) >= 2 ? A \ gL[ok] : [gL[ok[1]]/kL[ok[1]], NaN]
    return (kL = kL, gL = gL, TL = TL, cs = c[1], Gamma = c[2],
            kmax = 2c[1]/c[2], d = d)
end

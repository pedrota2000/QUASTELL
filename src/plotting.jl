# =============================================================================
#  Plotting
# =============================================================================

"""Potential V(phi) on [0, phiMax], with the minimum marked."""
function plot_potential(m::Bea_Mateos; npts::Int = 400)
    φ = range(0, 1.1 * m.phiMax, length = npts)
    p = plot(φ, m.V.(φ), lw = 2, label = "V(phi)",
             xlabel = "phi", ylabel = "V",
             title = "Bea-Mateos potential, phiM = $(m.phiM), phiQ = $(m.phiQ)")
    scatter!(p, [m.phiMax], [m.V(m.phiMax)], ms = 5, label = "phiMax")
    return p
end

"""The three background functions f, g and phi = u^nu phi0 against u."""
function plot_solution(m::Bea_Mateos, b::Background)
    p = plot(xlabel = "u", title = @sprintf("phiH = %g, N = %d, T/Lambda = %.6g",
                                            b.phiH, b.N, b.T / b.Lambda))
    plot!(p, b.x, b.f0, lw = 2, label = "f(u)")
    plot!(p, b.x, b.g0, lw = 2, label = "g(u)")
    plot!(p, b.x, b.x.^m.nu .* b.phi0, lw = 2, label = "phi(u)")
    return p
end
plot_solution(m::Bea_Mateos, phiH::Real; kwargs...) =
    plot_solution(m, background(m, phiH; kwargs...))

"""
    plot_thermodynamics(m)

Four panels: T/Lambda and s/Lambda^3 against phiH, the equation of state
s/Lambda^3 vs T/Lambda, and the swallowtail p vs T that fixes Tc (marked).
"""
function plot_thermodynamics(m::Bea_Mateos; size = (1100, 750), zoom = nothing)
    m.thermo === nothing && error("call thermodynamics!(m) first")
    th = m.thermo
    ok = th.converged .& isfinite.(th.t) .& isfinite.(th.s)

    # window to display: everything, or a relative half-width `zoom` around Tc
    sel = ok
    tlo = thi = 0.0
    if zoom !== nothing
        m.Tc === nothing && error("no Tc available to zoom around; run thermodynamics! with findTc = true")
        w = zoom === true ? 0.15 : float(zoom)
        tlo, thi = m.Tc*(1 - w), m.Tc*(1 + w)
        sel = ok .& (th.t .>= tlo) .& (th.t .<= thi)
        count(sel) < 3 && error("only $(count(sel)) points inside the window; widen `zoom`")
    end

    Fok = sel .& isfinite.(th.F)

    p1 = plot(th.phiH[ok], th.t[ok], lw = 2, legend = false,
              xlabel = "phiH", ylabel = "T/Lambda", title = "temperature")
    p2 = plot(th.phiH[ok], th.s[ok], lw = 2, legend = false,
              xlabel = "phiH", ylabel = "s/Lambda^3", title = "entropy density")
    p3 = plot(th.t[ok], th.s[ok], lw = 2, legend = false,
              xlabel = "T/Lambda", ylabel = "s/Lambda^3", title = "equation of state")
    p4 = plot(th.t[ok], th.F[ok], lw = 2, legend = false,
              xlabel = "T/Lambda", ylabel = "F/Lambda^4 (up to a constant)",
              title = "free energy (swallowtail)")

    pl = _pad(th.F[Fok])
    if zoom !== nothing
        plot!(p3, xlims = (tlo, thi), ylims = _pad(th.s[sel]))
        plot!(p4, xlims = (tlo, thi), ylims = pl)
        plot!(p1, xlims = _pad(th.phiH[sel], 0.02), ylims = (tlo, thi))
        plot!(p2, xlims = _pad(th.phiH[sel], 0.02), ylims = _pad(th.s[sel]))
    end

    if m.Tc !== nothing
        vline!(p3, [m.Tc], ls = :dash, c = :black, label = "")
        vline!(p4, [m.Tc], ls = :dash, c = :black, label = "")
        hline!(p1, [m.Tc], ls = :dash, c = :black, label = "")
        annotate!(p4, m.Tc, pl[1] + 0.08*(pl[2] - pl[1]),
                  text(@sprintf(" Tc = %.5g", m.Tc), 8, :left))
    end
    return plot(p1, p2, p3, p4, layout = (2, 2), size = size,
                left_margin = 6Plots.mm, bottom_margin = 4Plots.mm)
end

"padded (min, max) of `v`, with `frac` relative margin"
function _pad(v, frac = 0.05)
    isempty(v) && return (0.0, 1.0)
    lo, hi = minimum(v), maximum(v)
    d = hi - lo
    pad = d > 0 ? frac*d : max(abs(hi), 1.0)*frac
    return (lo - pad, hi + pad)
end

"""
    plot_qnms(r; showraw = true)

Converged modes in the complex lambda = omega/2piT plane.  With `showraw` the
unfiltered eigenvalues of the coarsest grid are drawn behind them in grey: those
grey points with no black point on top are exactly the artifacts the
multi-resolution test removed.  The horizontal line is Im lambda = 0, above
which a mode is unstable.
"""
function plot_qnms(r::QNMResult; showraw::Bool = true, size = (700, 550))
    p = plot(xlabel = "Re ω/2πT", ylabel = "Im ω/2πT", size = size,
             title = @sprintf("phiH = %g, q = %g  (grids %s)", r.phiH, r.q, string(r.Ns)))
    if showraw
        scale = 2π * r.T
        for (i, N) in enumerate(r.Ns)
            scatter!(p, real.(r.raw[i]) ./ scale, imag.(r.raw[i]) ./ scale,
                     ms = 3, msw = 0, alpha = 0.35, label = "raw, N = $N")
        end
    end
    st = filter(z -> imag(z) <= 0, r.lambda)
    un = filter(z -> imag(z) >  0, r.lambda)
    scatter!(p, real.(st), imag.(st), ms = 6, c = :black, msw = 0,
             label = "converged ($(length(st)))")
    isempty(un) || scatter!(p, real.(un), imag.(un), ms = 10, c = :red, m = :star5,
                            label = "UNSTABLE ($(length(un)))")
    hline!(p, [0.0], ls = :dash, c = :gray, label = "")
    return p
end

"""
    plot_convergence(r)

|Im lambda| of the raw spectrum at each resolution, side by side.  Physical
modes line up in a column across the panels; artifacts do not.
"""
function plot_convergence(r::QNMResult; size = (800, 550))
    scale = 2π * r.T
    p = plot(xlabel = "N (collocation points)", ylabel = "Im ω/2πT", size = size,
             title = @sprintf("resolution test, phiH = %g, q = %g", r.phiH, r.q))
    for (i, N) in enumerate(r.Ns)
        y = imag.(r.raw[i]) ./ scale
        scatter!(p, fill(N, length(y)), y, ms = 3, msw = 0, alpha = 0.4,
                 label = i == 1 ? "raw" : "")
    end
    for w in r.lambda
        hline!(p, [imag(w)], ls = :dot, c = :black, alpha = 0.6, label = "")
    end
    return p
end

"""
    plot_growth_rates(m, phiHs; qs = 0.05:0.1:1.2, hydro = true, kwargs...)

`gamma/Lambda` against `k/Lambda` for several spinodal states, with the fitted
hydrodynamic parabola of `growth_rate` overlaid as a dashed line.
"""
function plot_growth_rates(m::Bea_Mateos, phiHs; qs = 0.05:0.1:1.2,
                           hydro::Bool = true, size = (740, 540), kwargs...)
    p = plot(xlabel = "k/Lambda", ylabel = "gamma/Lambda", size = size,
             legend = :bottomleft, title = "growth rate of the unstable mode",
             left_margin = 5Plots.mm, bottom_margin = 4Plots.mm)
    gmax = 0.0
    for (i, phiH) in enumerate(phiHs)
        g = growth_rate(m, phiH; qs = qs,
                        desc = @sprintf("phiH = %.4f: ", phiH), kwargs...)
        pk = maximum(g.gL[isfinite.(g.gL)]; init = 0.0)
        gmax = max(gmax, pk)
        # Truncate where the tracked branch leaves the unstable regime: past the
        # zero crossing it collides with other modes and drops far below -pk,
        # which would otherwise squash the parabolas flat.
        last = findfirst(x -> x < -pk, g.gL)
        keep = last === nothing ? eachindex(g.kL) : (1:max(last-1, 1))
        plot!(p, g.kL[keep], g.gL[keep], lw = 2, m = :circle, ms = 3, c = i,
              label = @sprintf("T/Lambda = %.4f  (|cs| = %.3f)", g.TL, g.cs))
        if hydro && isfinite(g.Gamma)
            kk = range(0, maximum(g.kL[keep]), length = 200)
            plot!(p, kk, g.cs .* kk .- 0.5 * g.Gamma .* kk.^2, ls = :dash, c = i,
                  alpha = 0.6, label = "")
        end
    end
    hline!(p, [0.0], c = :black, label = "")
    gmax > 0 && plot!(p, ylims = (-0.6gmax, 1.25gmax))
    return p
end

"""
    plot_dispersion(d)  /  plot_dispersion(m, phiH; kwargs...)

Re and Im of omega/2piT against the dimensionless momentum q.
"""
function plot_dispersion(d::NamedTuple; size = (1000, 420), sound = true)
    nmodes = Base.size(d.lambda, 2)      # `size` is the figure-size keyword here
    p1 = plot(xlabel = "q = k/2πT", ylabel = "Re ω/2πT", title = "real part",
              legend = :topleft)
    p2 = plot(xlabel = "q = k/2πT", ylabel = "Im ω/2πT", title = "imaginary part",
              legend = :bottomleft)
    for n in 1:nmodes
        plot!(p1, d.q, real.(d.lambda[:, n]), lw = 2, m = :circle, ms = 3, label = "n = $n")
        plot!(p2, d.q, imag.(d.lambda[:, n]), lw = 2, m = :circle, ms = 3, label = "n = $n")
    end
    if sound
        cs, res = sound_speed(d)
        isfinite(cs) && plot!(p1, [0.0; d.q], cs .* [0.0; d.q], ls = :dash, c = :black,
                              label = @sprintf("cs = %.5f  (1/sqrt3 = %.5f)", cs, 1/sqrt(3)))
    end
    hline!(p2, [0.0], ls = :dash, c = :red, label = "")
    return plot(p1, p2, layout = (1, 2), size = size,
                left_margin = 5Plots.mm, bottom_margin = 4Plots.mm)
end
plot_dispersion(m::Bea_Mateos, phiH::Real; kwargs...) =
    plot_dispersion(dispersion(m, phiH; kwargs...))

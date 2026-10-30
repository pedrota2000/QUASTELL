# Regenerate the README animation:
#
#     julia --project assets/make_gif.jl
#
# A black brane is moved along its equation of state at small backreaction while the
# gapped QNMs at zero momentum are followed in the complex frequency plane.  The
# lightest gapped pair falls onto the imaginary axis, collides, and one of the two
# modes crosses into Im > 0 -- and stays there across a long stretch where ds/dT > 0,
# i.e. where the equation of state alone reports a stable plasma.
ENV["GKSwstype"] = "100"                     # draw off-screen
using QUASTELL, Plots, Printf

const m    = Bea_Mateos(0.8, 10; chi = 0.3)
const NS   = (40, 50, 60)
# Hot branch, middle branch, and -- past the second fold of T(phiH) at 1.760 -- the
# cold stable branch, which is squeezed against phiMax = 1.7708: uniform steps up to
# 1.70, then logarithmic in phiMax - phiH down to 5e-5.
const PHIS = vcat(collect(range(0.05, 1.70, length = 67)),
                  m.phiMax .- 10 .^ range(log10(m.phiMax - 1.71), -4.3, length = 30))
const OUT  = joinpath(@__DIR__, "gapmode.gif")
const LMAX = 1.45                            # the lightest gapped pair and its partner

th = thermodynamics!(m; N = 40, nphi = 400, progress = false, findTc = false)

# Cold starts fail above phiH ~ 1.6, so walk up the branch on every grid, seeding each
# background with the previous one.  `qnms` then finds them in the cache.
for N in NS
    prev = nothing
    for phiH in PHIS
        b = background(m, phiH; N = N, warm = prev)
        b.converged || (b = background(m, phiH; N = N, force = true))
        prev = b
    end
end

# Gapped modes at q = 0.  origintol = 1e-6 keeps light modes near the crossing.
usable(phiH) = all(N -> haskey(m.bg, (N, phiH)) && m.bg[(N, phiH)].converged, NS)
frames = map(filter(usable, PHIS)) do phiH
    r = qnms(m, phiH; q = 0.0, Ns = NS, rtol = 1e-4, cutoff = 30.0, origintol = 1e-6,
             bcexpand = true, verbose = false)
    light = [r.lambda[i] for i in eachindex(r.lambda)
             if r.err[i] < 1e-2 * abs(r.lambda[i]) && abs(r.lambda[i]) < LMAX]
    b = background(m, phiH; N = NS[end])
    @printf("phiH = %.3f  T/Lambda = %.5f  %d modes%s\n", phiH, b.T / b.Lambda,
            length(light), any(z -> imag(z) > 0, light) ? "  (unstable)" : "")
    (phiH = phiH, t = b.T / b.Lambda, s = b.S / b.Lambda^3, modes = light)
end

# Equation of state over the same stretch, with the (short) ds/dT < 0 part in red.
ok = findall(i -> th.converged[i] && PHIS[1] - 0.02 <= th.phiH[i] <= PHIS[end] + 0.02,
             eachindex(th.phiH))
t, s, ph = th.t[ok], th.s[ok], th.phiH[ok]
dsdT = [(s[min(i+1, end)] - s[max(i-1, 1)]) / (t[min(i+1, end)] - t[max(i-1, 1)])
        for i in eachindex(t)]
s_neg = [d < 0 ? s[i] : NaN for (i, d) in enumerate(dsdT)]
dsdT_at(phiH) = dsdT[argmin(abs.(ph .- phiH))]

function draw(k)
    f = frames[k]
    unstable = filter(z -> imag(z) > 0, f.modes)
    stable   = filter(z -> imag(z) <= 0, f.modes)

    p1 = plot(t, s; c = :gray30, lw = 2, label = "", xlabel = "T/Λ", ylabel = "s/Λ³",
              title = "equation of state", xscale = :log10, yscale = :log10,
              legend = :topleft)
    plot!(p1, t, s_neg; c = :red, lw = 3, label = "ds/dT < 0")
    scatter!(p1, [f.t], [f.s]; ms = 9, c = isempty(unstable) ? :dodgerblue : :red,
             msw = 2, msc = :black, label = "")

    p2 = plot(xlabel = "Re ω/2πT", ylabel = "Im ω/2πT", title = "gapped QNMs at k = 0",
              xlims = (-1.5, 1.5), ylims = (-1.5, 0.8), legend = :bottomright)
    hline!(p2, [0.0]; c = :black, lw = 1, label = "")
    for j in 1:k-1                                              # trails
        scatter!(p2, real.(frames[j].modes), imag.(frames[j].modes); ms = 2, msw = 0,
                 c = :gray70, alpha = 0.6, label = "")
    end
    scatter!(p2, real.(stable), imag.(stable); ms = 7, c = :dodgerblue, msw = 1,
             label = "damped")
    isempty(unstable) || scatter!(p2, real.(unstable), imag.(unstable); ms = 12,
                                  m = :star5, c = :red, msw = 1, label = "unstable")
    annotate!(p2, -1.45, 0.68, text(@sprintf("φH = %.4f   T/Λ = %.3f", f.phiH, f.t), 9, :left))
    if !isempty(unstable) && dsdT_at(f.phiH) > 0
        annotate!(p2, -1.45, 0.5, text("unstable although ds/dT > 0", 9, :left, :red))
    end

    plot(p1, p2; layout = (1, 2), size = (960, 400), left_margin = 5Plots.mm,
         bottom_margin = 5Plots.mm)
end

anim = Animation()
for k in vcat(1:length(frames), length(frames)-1:-1:2)         # there and back
    frame(anim, draw(k))
end
gif(anim, OUT; fps = 12)

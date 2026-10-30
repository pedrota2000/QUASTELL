# Regenerate the README animation:
#
#     julia --project assets/make_gif.jl
#
# A black brane is carried across the first-order transition -- from the hot branch,
# over the folded middle branch, onto the cold stable branch -- while the gapped QNMs
# at zero momentum are followed in the complex frequency plane.  At k = 0 energy conservation
# freezes the sound mode, so what moves is the order parameter: its gapped mode turns
# unstable exactly between the two zeros of ds/dT on the middle branch, where c_v > 0
# and the equation of state alone would call the plasma stable.
ENV["GKSwstype"] = "100"                     # draw off-screen
using QUASTELL, Plots, Printf

const m    = Bea_Mateos(0.8, 10; chi = 1.0)
const NS   = (40, 50, 60)
const OUT  = joinpath(@__DIR__, "gapmode.gif")
const WMAX = 2.0                             # window in omega/2piTc: the lightest gapped
                                             # pair and its partner
const IMTOL = 1e-5                           # Im above this is unstable; the nearly-real
                                             # cold-branch pair has |Im| ~ 1e-8 roundoff

th = thermodynamics!(m; N = 40, nphi = 1500, progress = false)
const Tc = m.Tc                              # in units of Lambda

# Around the transition: from T = 1.15 Tc on the hot branch, over the folded middle
# branch, down to T = 0.7 Tc on the cold one.  Uniform steps to 1.70, then logarithmic
# in phiMax - phiH, since the cold branch is squeezed against phiMax.
const TWIN = (0.65, 1.35)                    # T/Tc shown in the left panel
i2 = findfirst(i -> th.converged[i] && th.t[i] <= 1.15Tc, eachindex(th.phiH))
const PHIS = vcat(collect(range(th.phiH[i2], 1.70, length = 70)),
                  m.phiMax .- 10 .^ range(log10(m.phiMax - 1.71), -4.3, length = 26))

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

# Gapped modes at k = 0, as omega/2piTc = lambda T/Tc: in units of T a fixed gap
# runs off to infinity down the cold branch, in units of Tc it stays put.
# origintol = 1e-6 keeps light modes near the crossing.
usable(phiH) = all(N -> haskey(m.bg, (N, phiH)) && m.bg[(N, phiH)].converged, NS)
Tof(phiH) = (b = m.bg[(NS[end], phiH)]; b.T / b.Lambda / Tc)
inwindow(phiH) = phiH < 1.66 || Tof(phiH) >= 0.7         # stop at 0.7 Tc on the cold branch
frames = map(filter(p -> usable(p) && inwindow(p), PHIS)) do phiH
    r = qnms(m, phiH; q = 0.0, Ns = NS, rtol = 1e-4, cutoff = 30.0, origintol = 1e-6,
             bcexpand = true, verbose = false)
    b = background(m, phiH; N = NS[end])
    tt = b.T / b.Lambda / Tc
    light = [r.lambda[i] * tt for i in eachindex(r.lambda)
             if r.err[i] < 1e-2 * abs(r.lambda[i]) && abs(r.lambda[i]) * tt < WMAX]
    @printf("phiH = %.5f  T/Tc = %.4f  %d modes%s\n", phiH, tt,
            length(light), any(z -> imag(z) > IMTOL, light) ? "  (unstable)" : "")
    (phiH = phiH, t = tt, s = b.S / b.Lambda^3 / Tc^3, modes = light)
end

# Equation of state in units of Tc over the same stretch, ds/dT < 0 in red.
ok = findall(i -> th.converged[i] && PHIS[1] - 0.01 <= th.phiH[i], eachindex(th.phiH))
t, s, ph = th.t[ok] ./ Tc, th.s[ok] ./ Tc^3, th.phiH[ok]
dsdT = [(s[min(i+1, end)] - s[max(i-1, 1)]) / (t[min(i+1, end)] - t[max(i-1, 1)])
        for i in eachindex(t)]
s_neg = [d < 0 ? s[i] : NaN for (i, d) in enumerate(dsdT)]
dsdT_at(phiH) = dsdT[argmin(abs.(ph .- phiH))]
smax = maximum(f.s for f in frames)

function draw(k)
    f = frames[k]
    unstable = filter(z -> imag(z) > IMTOL, f.modes)
    stable   = filter(z -> imag(z) <= IMTOL, f.modes)

    p1 = plot(t, s; c = :gray30, lw = 2, label = "", xlabel = "T/Tc", ylabel = "s/Tc³",
              title = "equation of state", xlims = TWIN, ylims = (0, 1.1smax),
              legend = :topleft)
    plot!(p1, t, s_neg; c = :red, lw = 3, label = "ds/dT < 0")
    vline!(p1, [1.0]; ls = :dash, c = :black, label = "Tc")
    scatter!(p1, [f.t], [f.s]; ms = 9, c = isempty(unstable) ? :dodgerblue : :red,
             msw = 2, msc = :black, label = "")

    p2 = plot(xlabel = "Re ω/2πTc", ylabel = "Im ω/2πTc", title = "gapped QNMs at k = 0",
              xlims = (-2.0, 2.0), ylims = (-1.6, 0.8), legend = :bottomright)
    hline!(p2, [0.0]; c = :black, lw = 1, label = "")
    for j in 1:k-1                                              # trails
        scatter!(p2, real.(frames[j].modes), imag.(frames[j].modes); ms = 2, msw = 0,
                 c = :gray70, alpha = 0.6, label = "")
    end
    scatter!(p2, real.(stable), imag.(stable); ms = 7, c = :dodgerblue, msw = 1,
             label = "damped")
    isempty(unstable) || scatter!(p2, real.(unstable), imag.(unstable); ms = 12,
                                  m = :star5, c = :red, msw = 1, label = "unstable")
    annotate!(p2, -1.93, 0.68, text(@sprintf("φH = %.4f   T/Tc = %.3f", f.phiH, f.t), 9, :left))
    if !isempty(unstable) && dsdT_at(f.phiH) > 0
        annotate!(p2, -1.93, 0.5, text("unstable although ds/dT > 0", 9, :left, :red))
    end

    plot(p1, p2; layout = (1, 2), size = (960, 400), left_margin = 5Plots.mm,
         bottom_margin = 5Plots.mm)
end

anim = Animation()
for k in vcat(1:length(frames), length(frames)-1:-1:2)         # there and back
    frame(anim, draw(k))
end
gif(anim, OUT; fps = 12)

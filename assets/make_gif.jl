# Regenerate the README animation:
#
#     julia --project assets/make_gif.jl
#
# A black brane is moved along its equation of state, across the spinodal branch,
# while the light quasinormal modes at small momentum are followed in the complex
# frequency plane.  On the spinodal branch the sound pair collides on the
# imaginary axis and one of the two modes grows.
ENV["GKSwstype"] = "100"                     # draw off-screen
using QUASTELL, Plots, Printf

const m    = Bea_Mateos(1.0, 10; chi = 1.2)
const q    = 0.05
const NS   = (40, 50, 60)
const PHIS = range(0.70, 1.50, length = 61)
const OUT  = joinpath(@__DIR__, "spinodal.gif")

th = thermodynamics!(m; N = 40, nphi = 400, progress = false)

# Light modes along the sweep.  `err < 1e-2 |lambda|` also drops the exact zeros
# that `bcexpand = true` leaves in the spectrum at q != 0.
frames = map(PHIS) do phiH
    r = qnms(m, phiH; q = q, Ns = NS, rtol = 1e-4, bcexpand = true, verbose = false)
    light = [r.lambda[i] for i in eachindex(r.lambda)
             if r.err[i] < 1e-2 * abs(r.lambda[i]) && abs(r.lambda[i]) < 0.1]
    b = background(m, phiH; N = NS[end])
    @printf("phiH = %.3f  T/Lambda = %.5f  %d light modes%s\n", phiH, b.T / b.Lambda,
            length(light), any(z -> imag(z) > 0, light) ? "  (unstable)" : "")
    (phiH = phiH, t = b.T / b.Lambda, s = b.S / b.Lambda^3, modes = light)
end

# Equation of state over the same stretch of the branch, spinodal part in red.
ok = findall(i -> th.converged[i] && PHIS[1] - 0.05 <= th.phiH[i] <= PHIS[end] + 0.02,
             eachindex(th.phiH))
t, s = th.t[ok], th.s[ok]
dsdT = [(s[min(i+1, end)] - s[max(i-1, 1)]) / (t[min(i+1, end)] - t[max(i-1, 1)])
        for i in eachindex(t)]
s_spinodal = [d < 0 ? s[i] : NaN for (i, d) in enumerate(dsdT)]

allmodes = reduce(vcat, [f.modes for f in frames])
xmax = 1.15 * maximum(abs.(real.(allmodes)))
ylo, yhi = extrema(imag.(allmodes))
ypad = 0.12 * (yhi - ylo)

function draw(k)
    f = frames[k]
    unstable = filter(z -> imag(z) > 0, f.modes)
    stable   = filter(z -> imag(z) <= 0, f.modes)

    p1 = plot(t, s; c = :gray30, lw = 2, label = "", xlabel = "T/Λ", ylabel = "s/Λ³",
              title = "equation of state")
    plot!(p1, t, s_spinodal; c = :red, lw = 3, label = "spinodal branch")
    vline!(p1, [m.Tc]; ls = :dash, c = :black, label = "Tc")
    scatter!(p1, [f.t], [f.s]; ms = 9, c = isempty(unstable) ? :dodgerblue : :red,
             msw = 2, msc = :black, label = "")

    p2 = plot(xlabel = "Re ω/2πT", ylabel = "Im ω/2πT",
              title = @sprintf("light QNMs at k/2πT = %.2f", q),
              xlims = (-xmax, xmax), ylims = (ylo - ypad, yhi + ypad), legend = :bottomright)
    hline!(p2, [0.0]; c = :black, lw = 1, label = "")
    for j in 1:k-1                                              # trails
        scatter!(p2, real.(frames[j].modes), imag.(frames[j].modes); ms = 2, msw = 0,
                 c = :gray70, alpha = 0.6, label = "")
    end
    scatter!(p2, real.(stable), imag.(stable); ms = 7, c = :dodgerblue, msw = 1,
             label = "damped")
    isempty(unstable) || scatter!(p2, real.(unstable), imag.(unstable); ms = 12,
                                  m = :star5, c = :red, msw = 1, label = "unstable")
    annotate!(p2, -0.95xmax, yhi + 0.5ypad,
              text(@sprintf("φH = %.2f   T/Λ = %.4f", f.phiH, f.t), 9, :left))

    plot(p1, p2; layout = (1, 2), size = (960, 400), left_margin = 5Plots.mm,
         bottom_margin = 5Plots.mm)
end

anim = Animation()
for k in vcat(1:length(frames), length(frames)-1:-1:2)         # there and back
    frame(anim, draw(k))
end
println("frames in ", anim.dir)
gif(anim, OUT; fps = 12)

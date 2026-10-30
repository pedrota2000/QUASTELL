using QUASTELL
using Test

# All checks use phiM = 1, phiQ = 10, chi = 1.2, whose spinodal branch (ds/dT < 0)
# is phiH in [1.033, 1.271].  Grids stay at N <= 60: beyond that the background
# Newton residual floors above its tolerance and every solve runs to `maxit`.
const m = Bea_Mateos(1.0, 10; chi = 1.2)
const NS = (40, 50, 60)

# At q = 0 `qnms` drops the degenerate band at omega = 0 and warns that the
# hydrodynamic mode goes with it.  The warning is expected, so check for it
# instead of printing it.
qnms_q0(phiH; kw...) =
    @test_logs (:warn, r"q = 0") match_mode = :any qnms(m, phiH; q = 0.0, verbose = false, kw...)

@testset "QUASTELL" begin

    @testset "Chebyshev grid" begin
        x, D, D2 = cheby(24)
        @test x[1] == 0 && x[end] == 1
        @test issorted(x)
        @test D * x.^5 ≈ 5 .* x.^4
        @test D2 * x.^5 ≈ 20 .* x.^3
    end

    @testset "model" begin
        @test m.phiMax > 0
        @test abs(m.Vp(m.phiMax)) < 1e-5            # phiMax is the minimum of V
        @test m.V(0.0) == 0                          # m.V is the potential + 3
        @test m.Vpp(0.0) == -3                       # m^2 = -3, Delta = 3
    end

    @testset "background" begin
        bg = background(m, 0.3; N = 40)
        @test bg.converged
        @test bg.res < 1e-10
        @test bg.T      ≈ 0.31942632 rtol = 1e-7
        @test bg.S      ≈ 2.964797   rtol = 1e-6
        @test bg.Lambda ≈ 0.35152238 rtol = 1e-7
        @test bg.f0[1] ≈ 1 && bg.g0[1] ≈ 1 && abs(bg.f0[end]) < 1e-12
        @test background(m, 0.3; N = 40) === bg      # cached
    end

    @testset "fluctuation matrices" begin
        M0, M1, Mt, b = mtildes(m, 0.3; q = 0.3, N = 40)
        @test length(Mt) == 4                        # cubic in omega
        @test all(size(X) == (80, 80) for X in Mt)   # 2 fields x N
        @test size(M0) == size(M1) == (240, 240)     # companion form
    end

    @testset "AdS5-Schwarzschild limit" begin
        r = qnms_q0(1e-3; Ns = NS, rtol = 1e-4, cutoff = 40.0)
        refs = [1.099407283 - 0.879767314im,         # Zphi, Delta = 3, n = 1
                1.559725796 - 1.373337872im,         # Z2,   Delta = 4, n = 1
                2.105948627 - 1.887444155im,         # Zphi, n = 2
                2.584760    - 2.381785im]            # Z2,   n = 2
        for z in refs
            @test minimum(abs.(r.lambda .- z)) < 1e-5
        end
        @test all(imag.(r.lambda) .< 0)              # nothing unstable
    end

    @testset "conformal sound mode" begin
        q = 0.05
        r = qnms(m, 1e-3; q = q, Ns = NS, rtol = 1e-4, cutoff = 40.0, verbose = false)
        z = r.lambda[argmin(abs.(r.lambda .- (q/sqrt(3) - im*q^2/3)))]
        # lambda = q/sqrt(3) - i q^2/3 + (3 - 2 log 2)/(6 sqrt(3)) q^3 + O(q^4)
        @test real(z) ≈ q/sqrt(3) + (3 - 2log(2))/(6sqrt(3)) * q^3 atol = 1e-6
        @test imag(z) ≈ -q^2/3 atol = 1e-6
    end

    @testset "thermodynamics and Tc" begin
        mt = Bea_Mateos(1.0, 10; chi = 1.2)
        th = thermodynamics!(mt; N = 40, nphi = 400, progress = false)
        @test all(th.converged)
        i0 = findfirst(th.converged)                 # hottest point: nearly conformal
        @test th.e[i0] / (3 * th.p[i0]) ≈ 1 atol = 5e-3
        @test mt.Tc ≈ 0.37356750 rtol = 1e-6        # N = 50, nphi = 2000 gives 0.3735674956

        dsdT(phiH) = (i = argmin(abs.(th.phiH .- phiH));
                      (th.s[i+1] - th.s[i-1]) / (th.t[i+1] - th.t[i-1]))
        @test dsdT(1.16) < 0                         # spinodal branch
        @test dsdT(1.40) > 0                         # stable branch
    end

    @testset "spinodal instability" begin
        # Inside the spinodal region c_s^2 < 0: the sound mode is purely imaginary and grows.
        r = qnms(m, 1.16; q = 0.05, Ns = NS, rtol = 1e-4, bcexpand = true, verbose = false)
        lead = r.lambda[1]
        @test imag(lead) ≈ 0.0110679 atol = 1e-6
        @test abs(real(lead)) < 1e-8
        @test count(z -> imag(z) > 0, r.lambda) == 1

        # Outside it the same mode is ordinary, damped sound.
        r = qnms(m, 1.40; q = 0.05, Ns = NS, rtol = 1e-4, bcexpand = true, verbose = false)
        @test count(z -> imag(z) > 0, r.lambda) == 0
        sound = r.lambda[findfirst(z -> abs(z) > 1e-8, r.lambda)]
        @test abs(real(sound)) ≈ 0.0193923 atol = 1e-6
        @test imag(sound) < 0
    end

    @testset "plots build" begin
        r = qnms(m, 0.3; q = 0.3, Ns = (40, 50), rtol = 1e-4, verbose = false)
        @test plot_qnms(r) isa QUASTELL.Plots.Plot
        @test plot_potential(m) isa QUASTELL.Plots.Plot
    end
end

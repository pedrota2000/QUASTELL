# =============================================================================
#  1.  Background solution for a given phiH
# =============================================================================

"""
    _newton_matrix(m, x, D, D2, F0, phiH) -> MMop, v

One Newton step of the background problem: the linearized operator acting on
(delta f, delta g, delta phi0) and the residual source.  Rows 1, N+1, 2N+1 carry
the boundary conditions f(0) = g(0) = 1 and phi0(1) = phiH; rows N, 2N, 3N the
horizon conditions f(1) = 0 and the two regularity constraints.
"""
function _newton_matrix(m::Bea_Mateos, x, D, D2, F0, phiH)
    N = length(x)
    f0 = F0[1:N]; g0 = F0[N+1:2N]; phi0 = F0[2N+1:end]
    Id = I(N)
    nu = m.nu; chi = m.chi

    VV   = m.V.(x.^nu .* phi0)
    VVp  = m.Vp.(x.^nu .* phi0)
    VVpp = m.Vpp.(x.^nu .* phi0)

    A11 = x .* g0
    B11 = -3 .* g0 .+ 2 .* x .* (D * g0)
    C11 = zeros(N); A12 = zeros(N)
    B12 = 2 .* x .* (D * f0)
    C12 = -3 .* (D * f0) .+ x .* (D2 * f0)
    A13 = zeros(N); B13 = zeros(N); C13 = zeros(N)
    S1 = 2 .* x .* (D * f0) .* (D * g0) .+ g0 .* (-3 .* (D * f0) .+ x .* (D2 * f0))

    A21 = zeros(N); B21 = zeros(N); C21 = zeros(N)
    A22 = -3 .* x .* g0; A23 = zeros(N)
    B22 = -3 .* g0 .+ 6 .* x .* (D * g0)
    B23 = -8 .* x.^(2nu) .* chi .* g0.^2 .* (nu .* phi0 .+ x .* (D * phi0))
    C22 = -8 .* x.^(-1 + 2nu) .* chi .* g0 .* (nu .* phi0 .+ x .* (D * phi0)).^2 .-
          3 .* ((D * g0) .+ x .* (D2 * g0))
    C23 = -8 .* chi .* nu .* x.^(-1 + 2nu) .* g0.^2 .* (nu .* phi0 .+ x .* (D * phi0))
    S2 = 3 .* x .* (D * g0).^2 .-
         4 .* x.^(-1 + 2nu) .* chi .* g0.^2 .* (nu .* phi0 .+ x .* (D * phi0)).^2 .-
         3 .* g0 .* ((D * g0) .+ x .* (D2 * g0))

    A31 = zeros(N); A32 = zeros(N)
    A33 = x.^(2 + nu) .* f0
    B31 = x.^(1 + nu) .* (nu .* phi0 .+ x .* (D * phi0))
    B32 = (2 .* x.^(1 + nu) .* f0 .* (nu .* phi0 .+ x .* (D * phi0))) ./ g0
    B33 = x.^(1 + nu) .* (x .* (D * f0) .+ f0 .* ((-3 + 2nu) .+ (2 .* x .* (D * g0)) ./ g0))
    C31 = x.^nu .* (2 .* x .* (D * g0) .* (nu .* phi0 .+ x .* (D * phi0)) .+
          g0 .* ((-4 + nu) .* nu .* phi0 .+
                 x .* ((-3 + 2nu) .* (D * phi0) .+ x .* (D2 * phi0)))) ./ g0
    C32 = -(2 .* x.^(1 + nu) .* f0 .* (D * g0) .* (nu .* phi0 .+ x .* (D * phi0))) ./ g0.^2
    C33 = x.^nu ./ g0 .* (nu .* f0 .* ((-4 + nu) .* g0 .+ 2 .* x .* (D * g0)) .+
                          g0 .* (x .* nu .* (D * f0) .- VVpp))
    S3 = x.^(1 + nu) .* nu .* phi0 .* (D * f0) .- VVp .+
         x.^(2 + nu) .* (D * f0) .* (D * phi0) .+
         (x.^nu .* f0) ./ g0 .* (2 .* x .* (D * g0) .* (nu .* phi0 .+ x .* (D * phi0)) .+
             g0 .* ((-4 + nu) .* nu .* phi0 .+
                    x .* ((-3 + 2nu) .* (D * phi0) .+ x .* (D2 * phi0))))

    BH11 = nu .* phi0[end] .+ D[end:end, :] * phi0
    BH12 = 0.0; BH13 = D[end:end, :] * f0
    CH11 = 0.0; CH12 = 0.0
    CH13 = nu .* (D[end:end, :] * f0) .- VVpp[end]
    SH1 = -VVp[end] .+ (D[end:end, :] * f0) .* (nu .* phi0[end] .+ D[end:end, :] * phi0)

    BH21 = 3/4 * (-2 .+ (D[end:end, :] * g0) ./ g0[end])
    BH22 = (3 .* D[end:end, :] * f0) ./ (4 .* g0[end])
    BH23 = 0.0; CH21 = 0.0
    CH22 = .-(3 .* (D[end:end, :] * f0) .* (D[end:end, :] * g0)) ./ (4 .* g0[end].^2)
    CH23 = 2 .* chi .* VVp[end]
    SH2 = -6 .+ 2 * chi * VV[end] .+
          (3 .* (D[end:end, :] * f0) .* (-2 .* g0[end] .+ D[end:end, :] * g0)) ./ (4 .* g0[end])

    MMf = hcat(Diagonal(A11)*D2 .+ Diagonal(B11)*D .+ Diagonal(C11),
               Diagonal(A12)*D2 .+ Diagonal(B12)*D .+ Diagonal(C12),
               Diagonal(A13)*D2 .+ Diagonal(B13)*D .+ Diagonal(C13))
    MMg = hcat(Diagonal(A21)*D2 .+ Diagonal(B21)*D .+ Diagonal(C21),
               Diagonal(A22)*D2 .+ Diagonal(B22)*D .+ Diagonal(C22),
               Diagonal(A23)*D2 .+ Diagonal(B23)*D .+ Diagonal(C23))
    MMphi = hcat(Diagonal(A31)*D2 .+ Diagonal(B31)*D .+ Diagonal(C31),
                 Diagonal(A32)*D2 .+ Diagonal(B32)*D .+ Diagonal(C32),
                 Diagonal(A33)*D2 .+ Diagonal(B33)*D .+ Diagonal(C33))
    MMop = vcat(MMf, MMg, MMphi)

    MMop[1, :] .= 0.0;   MMop[1, 1] = 1.0
    MMop[N+1, :] .= 0.0; MMop[N+1, N+1] = 1.0
    MMop[2N+1, :] .= 0.0; MMop[2N+1, 3N] = 1.0
    MMop[N, :] .= 0.0;   MMop[N, N] = 1.0
    MMop[2N:2N, :] .= hcat(BH11 .* D[end:end, :] .+ CH11 .* Id[end:end, :],
                           BH12 .* D[end:end, :] .+ CH12 .* Id[end:end, :],
                           BH13 .* D[end:end, :] .+ CH13 .* Id[end:end, :])
    MMop[3N:3N, :] .= hcat(BH21 .* D[end:end, :] .+ CH21 .* Id[end:end, :],
                           BH22 .* D[end:end, :] .+ CH22 .* Id[end:end, :],
                           BH23 .* D[end:end, :] .+ CH23 .* Id[end:end, :])

    v = vcat(vcat(1 - f0[1],      -S1[2:end-1], -f0[end]),
             vcat(1 - g0[1],      -S2[2:end-1], -SH1),
             vcat(phiH - phi0[end], -S3[2:end-1], -SH2))
    return MMop, v
end

"""
    _newton(m, x, D, D2, phiH, F0; lambda, tol, maxit)

Damped Newton iteration at fixed `phiH`, warm-started from `F0`.
"""
function _newton(m::Bea_Mateos, x, D, D2, phiH, F0; lambda = 1.0, tol = 1e-10, maxit = 2000)
    N = length(x)
    f = copy(F0[1:N]); g = copy(F0[N+1:2N]); ph = copy(F0[2N+1:end])
    _, v = _newton_matrix(m, x, D, D2, F0, phiH)
    res = maximum(abs.(v)); it = 1
    while res >= tol && it <= maxit
        MMop, v = _newton_matrix(m, x, D, D2, vcat(f, g, ph), phiH)
        step = MMop \ v
        f  .+= lambda .* step[1:N]
        g  .+= lambda .* step[N+1:2N]
        ph .+= lambda .* step[2N+1:end]
        res = maximum(abs.(v)); it += 1
        isfinite(res) || break
    end
    return f, g, ph, res
end

"""
    background(m, phiH; N = 40, nsteps = 20, lambda = 1.0, tol = 1e-10,
               maxit = 2000, warm = nothing, force = false) -> Background

Solve the background for horizon scalar `phiH` on an `N`-point grid, by
continuation in phiH from pure AdS (`nsteps` stages) unless a warm start `warm`
(a `Background` or a stacked `[f; g; phi0]`) is given.  Results are cached, so
calling it twice with the same `(N, phiH)` is free; pass `force = true` to redo.

`T`, `S`, `Lambda` follow the background solver: T = -f'(1) sqrt(g(1)) / 4pi,
S = pi g(1)^(3/2), Lambda = phi0(0).
"""
function background(m::Bea_Mateos, phiH::Real; N::Int = 40, nsteps::Int = 20,
                    lambda::Real = 1.0, tol::Real = 1e-10, maxit::Int = 2000,
                    warm = nothing, force::Bool = false)
    key = (N, float(phiH))
    (!force && haskey(m.bg, key)) && return m.bg[key]

    x, D, D2 = cheby(N)

    F0 = if warm isa Background
        vcat(warm.f0, warm.g0, warm.phi0)
    elseif warm isa AbstractVector
        collect(float.(warm))
    else
        nothing
    end

    ramp = F0 === nothing ? collect(range(phiH / nsteps, phiH, length = nsteps)) : [float(phiH)]
    F0 === nothing && (F0 = vcat(1 .- x.^4, ones(N), (phiH / nsteps) .* x))

    f = g = ph = nothing; res = Inf
    for target in ramp
        f, g, ph, res = _newton(m, x, D, D2, target, F0;
                                lambda = lambda, tol = tol, maxit = maxit)
        isfinite(res) || break
        F0 = vcat(f, g, ph)
    end

    T = -1/(4π) * (D[end:end, :] * f)[1] * sqrt(g[end])
    S = π * g[end]^(3/2)
    Lam = ph[1]
    conv = isfinite(res) && res < max(tol, 1e-8) && !any(isnan, F0)
    conv || @warn "background did not converge" phiH N res

    bgsol = Background(N, x, D, D2, f, g, ph, float(phiH), T, S, Lam, res, conv)
    m.bg[key] = bgsol
    return bgsol
end

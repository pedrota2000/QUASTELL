# =============================================================================
#  3.  The Mtilde matrices
# =============================================================================

# Generated from Mathematica: the coefficients of the two coupled fluctuation
# equations, P = 3 (highest power of omega) and Nf = 2 master fields.  Each `cc`
# is one coefficient evaluated on the grid; it multiplies Dm[1] = I, Dm[2] = D or
# Dm[3] = D2 according to which derivative of the field it belongs to.  Regenerate
# this block with the exporter in mathematica/Derive_the_equations.nb if the
# equations change.
function build_Mtildes(x, D, D2, f0, g0, phi0, Vfuncs, q, T, chi; nu = 1.0)
    N = length(x); Nf = 2; P = 3
    phi  = x.^nu .* phi0
    dphi = nu .* x.^(nu-1) .* phi0 .+ x.^nu .* (D * phi0)
    df = D*f0; dg = D*g0
    VV = Vfuncs[1].(phi); VVp = Vfuncs[2].(phi); VVpp = Vfuncs[3].(phi)
    Dm = (I(N), D, D2)
    rows(a) = (a-1)*N+1 : a*N
    Mt = [zeros(ComplexF64, Nf*N, Nf*N) for _ in 0:P]

    cc = @. (f0*g0^2*pi^2*q^2*T^2*(-36*dg^3*f0*pi^2*q^2*T^2*x^5 + 9*dg^2*f0*g0*x^3*(-(dg*VVpp) + 24*pi^2*q^2*T^2*x) - 6*dg*g0^2*x^2*(-16*pi^2*q^2*T^2*(-3 + chi*VV)*x + f0*(dg*(-9*VVpp + 16*chi*dphi*VVp*x) + 8*pi^2*q^2*T^2*x*(9 + chi*dphi^2*x^2))) - 8*g0^4*(6*(-3 + chi*VV)*VVpp - 3*f0*(3*VVpp + chi*dphi^2*(-48 + 16*chi*VV + VVpp)*x^2) + 4*chi*dphi*VVp*x*(6 - 2*chi*VV + f0*(12 + chi*dphi^2*x^2))) + 12*g0^3*x*(-2*(-3 + chi*VV)*(-(dg*VVpp) + 8*pi^2*q^2*T^2*x) + f0*(8*pi^2*q^2*T^2*x*(3 + chi*dphi^2*x^2) - dg*(9*VVpp - 32*chi*dphi*VVp*x + chi*dphi^2*(-48 + 16*chi*VV + VVpp)*x^2)))))/2
    Mt[1][rows(1), rows(1)] .+= Diagonal(cc) * Dm[1]

    cc = @. (f0*g0^3*pi^2*q^2*T^2*x*(9*dg^3*f0^2*x^3 - 36*dg*f0*g0^2*x*(6 - 2*chi*VV + f0*(-3 + chi*dphi^2*x^2)) + 6*dg^2*f0*g0*x^2*(12 - 4*chi*VV + f0*(-9 + 2*chi*dphi^2*x^2)) + 8*g0^3*(8*(-3 + chi*VV)^2 - 2*f0*(-3 + chi*VV)*(3 + 4*chi*dphi^2*x^2) + f0^2*(-9 + 3*chi*dphi^2*x^2 + 2*chi^2*dphi^4*x^4))))/2
    Mt[1][rows(1), rows(1)] .+= Diagonal(cc) * Dm[2]

    cc = @. (3*f0^2*g0^3*pi^2*q^2*T^2*x^2*(-2*g0 + dg*x)*(-12*dg*f0*g0*x + 3*dg^2*f0*x^2 + 4*g0^2*(6 - 2*chi*VV + f0*(3 + chi*dphi^2*x^2))))/2
    Mt[1][rows(1), rows(1)] .+= Diagonal(cc) * Dm[3]

    cc = @. (g0^4*x^2*(6 - 2*chi*VV + chi*dphi^2*f0*x^2)*(-3*dg*VVp*x + g0*(6*VVp - 8*dphi*(-3 + chi*VV)*x)))/2
    Mt[1][rows(1), rows(2)] .+= Diagonal(cc) * Dm[1]

    cc = @. (3*f0*g0^3*x^3*(-2*g0 + dg*x)*(3*dg*VVp*x + g0*(-6*VVp + 8*dphi*(-3 + chi*VV)*x)))/4
    Mt[1][rows(1), rows(2)] .+= Diagonal(cc) * Dm[2]

    cc = @. 16*chi*g0^(5/2)*pi^4*q^4*T^4*(12*dg*f0*g0*x - 3*dg^2*f0*x^2 + 2*g0^2*(6 - 2*chi*VV + f0*(-6 + chi*dphi^2*x^2)))*(3*dg^2*f0*VVp*x^2 + 12*dg*f0*g0*x*(-VVp + 2*dphi*(-3 + chi*VV)*x) + 4*g0^2*(-12*dphi*f0*(-3 + chi*VV)*x + VVp*(6 - 2*chi*VV + f0*(3 + chi*dphi^2*x^2))))
    Mt[1][rows(2), rows(1)] .+= Diagonal(cc) * Dm[1]

    cc = @. (pi^2*q^2*T^2*x^2*(27*dg^6*f0^2*x^6 - 108*dg^4*f0*g0*x^5*(3*dg*f0 + pi^2*q^2*T^2*x) + 64*g0^6*(4*chi*dphi^2*(-3 + chi*VV)^2*x^2 + f0^2*(-3 + chi*dphi^2*x^2)^2*(3 + chi*dphi^2*x^2) - 2*f0*(-3 + chi*VV)*(-18 - 3*chi*dphi^2*x^2 + 2*chi^2*dphi^4*x^4)) + 36*dg^3*f0*g0^2*x^4*(24*pi^2*q^2*T^2*x + dg*(4*(-3 + chi*VV) + f0*(45 - chi*dphi^2*x^2))) - 192*g0^5*x*(-6*pi^2*q^2*T^2*(-3 + chi*VV)*x - dg*f0^2*(-27 + 6*chi*dphi^2*x^2 + chi^2*dphi^4*x^4) + f0*(3*pi^2*q^2*T^2*x*(3 + chi*dphi^2*x^2) + 2*dg*(-3 + chi*VV)*(12 + chi*dphi^2*x^2))) - 144*dg^2*g0^3*x^3*(-2*pi^2*q^2*T^2*(-3 + chi*VV)*x - 2*dg*f0^2*(-15 + chi*dphi^2*x^2) + f0*(8*dg*(-3 + chi*VV) + pi^2*q^2*T^2*x*(18 + chi*dphi^2*x^2))) + 48*dg*g0^4*x^2*(-24*pi^2*q^2*T^2*(-3 + chi*VV)*x - dg*f0^2*(-135 + 18*chi*dphi^2*x^2 + chi^2*dphi^4*x^4) + 2*f0*(6*pi^2*q^2*T^2*x*(6 + chi*dphi^2*x^2) + dg*(-3 + chi*VV)*(36 + chi*dphi^2*x^2)))))/(4*sqrt(g0))
    Mt[1][rows(2), rows(2)] .+= Diagonal(cc) * Dm[1]

    cc = @. (3*sqrt(g0)*pi^2*q^2*T^2*x^3*(2*g0 - dg*x)*(135*dg^3*f0^2*g0*x^3 - 18*dg^4*f0^2*x^4 + 12*dg*f0*g0^3*x*(26*(-3 + chi*VV) + 3*f0*(13 + chi*dphi^2*x^2)) - 6*dg^2*f0*g0^2*x^2*(12*(-3 + chi*VV) + f0*(63 + 2*chi*dphi^2*x^2)) + 8*g0^4*(-42*f0*(-3 + chi*VV) - 8*(-3 + chi*VV)^2 + f0^2*(-27 - 3*chi*dphi^2*x^2 + 2*chi^2*dphi^4*x^4))))/4
    Mt[1][rows(2), rows(2)] .+= Diagonal(cc) * Dm[2]

    cc = @. (9*f0*g0^(3/2)*pi^2*q^2*T^2*x^4*(-2*g0 + dg*x)^2*(-12*dg*f0*g0*x + 3*dg^2*f0*x^2 + 4*g0^2*(6 - 2*chi*VV + f0*(3 + chi*dphi^2*x^2))))/4
    Mt[1][rows(2), rows(2)] .+= Diagonal(cc) * Dm[3]

    cc = @. (9*f0*g0^(3/2)*im*pi^2*q^2*T^2*x*(-2*g0 + dg*x)^2*(-12*dg*f0*g0*x + 3*dg^2*f0*x^2 + 4*g0^2*(6 - 2*chi*VV + f0*(3 + chi*dphi^2*x^2))))/4
    Mt[2][rows(1), rows(1)] .+= Diagonal(cc) * Dm[1]

    cc = @. -3*f0*g0^(5/2)*im*pi^2*q^2*T^2*x^2*(2*g0 - dg*x)*(-12*dg*f0*g0*x + 3*dg^2*f0*x^2 + 4*g0^2*(6 - 2*chi*VV + f0*(3 + chi*dphi^2*x^2)))
    Mt[2][rows(1), rows(1)] .+= Diagonal(cc) * Dm[2]

    cc = @. (3*g0^(5/2)*im*x^3*(-2*g0 + dg*x)*(3*dg*VVp*x + g0*(-6*VVp + 8*dphi*(-3 + chi*VV)*x)))/4
    Mt[2][rows(1), rows(2)] .+= Diagonal(cc) * Dm[1]

    cc = @. (3*im*pi^2*q^2*T^2*x^3*(2*g0 - dg*x)*(504*dg^3*f0*g0*x^3 - 63*dg^4*f0*x^4 + 48*dg*g0^3*x*(10*(-3 + chi*VV) + 3*f0*(14 + chi*dphi^2*x^2)) - 12*dg^2*g0^2*x^2*(10*(-3 + chi*VV) + 3*f0*(42 + chi*dphi^2*x^2)) + 16*g0^4*(-2*(-3 + chi*VV)*(15 + 4*chi*dphi^2*x^2) + f0*(-63 - 9*chi*dphi^2*x^2 + 4*chi^2*dphi^4*x^4))))/8
    Mt[2][rows(2), rows(2)] .+= Diagonal(cc) * Dm[1]

    cc = @. (9*g0*im*pi^2*q^2*T^2*x^4*(-2*g0 + dg*x)^2*(-12*dg*f0*g0*x + 3*dg^2*f0*x^2 + 4*g0^2*(6 - 2*chi*VV + f0*(3 + chi*dphi^2*x^2))))/2
    Mt[2][rows(2), rows(2)] .+= Diagonal(cc) * Dm[2]

    cc = @. (-3*f0*g0^2*(2*g0 - dg*x)*(36*dg^2*pi^2*q^2*T^2*x^4 + 9*dg*g0*x^2*(dg*VVpp - 16*pi^2*q^2*T^2*x) + 4*g0^3*(9*VVpp - 24*chi*dphi*VVp*x + 16*chi*dphi^2*(-3 + chi*VV)*x^2) + 12*g0^2*x*(12*pi^2*q^2*T^2*x + dg*(-3*VVpp + 4*chi*dphi*VVp*x))))/8
    Mt[3][rows(1), rows(1)] .+= Diagonal(cc) * Dm[1]

    cc = @. (-9*f0*g0^3*x*(-2*g0 + dg*x)^2*(3*dg*f0*x + g0*(24 - 8*chi*VV + f0*(-6 + 4*chi*dphi^2*x^2))))/8
    Mt[3][rows(1), rows(1)] .+= Diagonal(cc) * Dm[2]

    cc = @. (-27*f0^2*g0^3*x^2*(-2*g0 + dg*x)^3)/8
    Mt[3][rows(1), rows(1)] .+= Diagonal(cc) * Dm[3]

    cc = @. 12*chi*g0^(3/2)*pi^2*q^2*T^2*(2*g0 - dg*x)*(12*dg*f0*g0*x - 3*dg^2*f0*x^2 + 2*g0^2*(6 - 2*chi*VV + f0*(-6 + chi*dphi^2*x^2)))*(-6*dg^2*dphi*f0*x^3 + 3*dg*g0*x*(VVp + 8*dphi*f0*x) + g0^2*(-6*VVp + 4*dphi*f0*x*(-6 + chi*dphi^2*x^2)))
    Mt[3][rows(2), rows(1)] .+= Diagonal(cc) * Dm[1]

    cc = @. (27*x^2*(-2*g0 + dg*x)^4*(-4*g0^2*(-6 + 3*f0 + 2*chi*VV) - 3*dg^2*f0*x^2 + 12*g0*x*(dg*f0 + pi^2*q^2*T^2*x)))/(16*sqrt(g0))
    Mt[3][rows(2), rows(2)] .+= Diagonal(cc) * Dm[1]

    cc = @. (27*sqrt(g0)*x^3*(2*g0 - dg*x)^3*(27*dg*f0*g0*x - 6*dg^2*f0*x^2 + g0^2*(24 - 8*chi*VV + f0*(-30 + 4*chi*dphi^2*x^2))))/16
    Mt[3][rows(2), rows(2)] .+= Diagonal(cc) * Dm[2]

    cc = @. (-81*f0*g0^(3/2)*x^4*(-2*g0 + dg*x)^4)/16
    Mt[3][rows(2), rows(2)] .+= Diagonal(cc) * Dm[3]

    cc = @. (-81*f0*g0^(3/2)*im*x*(-2*g0 + dg*x)^4)/16
    Mt[4][rows(1), rows(1)] .+= Diagonal(cc) * Dm[1]

    cc = @. (-27*f0*g0^(5/2)*im*x^2*(-2*g0 + dg*x)^3)/4
    Mt[4][rows(1), rows(1)] .+= Diagonal(cc) * Dm[2]

    cc = @. (81*im*x^3*(-2*g0 + dg*x)^5)/32
    Mt[4][rows(2), rows(2)] .+= Diagonal(cc) * Dm[1]

    cc = @. (-81*g0*im*x^4*(-2*g0 + dg*x)^4)/8
    Mt[4][rows(2), rows(2)] .+= Diagonal(cc) * Dm[2]

    d = Nf*N
    Idb = Matrix{ComplexF64}(I, d, d)
    M0 = zeros(ComplexF64, P*d, P*d)
    M1 = zeros(ComplexF64, P*d, P*d)
    blk(i,j) = ((i-1)*d+1:i*d, (j-1)*d+1:j*d)
    for j in 1:P
        M0[blk(1,j)...] .= Mt[j]
    end
    for i in 2:P
        M0[blk(i,i)...]   .=  Idb
        M1[blk(i,i-1)...] .= -Idb
    end
    M1[blk(1,P)...] .= Mt[P+1]
    return M0, M1, Mt
end


"""
    mtildes(m, phiH; q = 0.0, N = 40, bg = nothing) -> M0, M1, Mt, bg

Build the fluctuation matrices (see `build_Mtildes` above) on the background at
`(N, phiH)`.  `Mt[n+1]` is the coefficient of omega^n in

    M(omega) = Mt[1] + omega Mt[2] + ... + omega^P Mt[P+1],   P = 3,

each a `Nf*N x Nf*N` block matrix over the two master fields, assembled the same
way as the background operator: `Diagonal(coef) * D2 + Diagonal(coef) * D +
Diagonal(coef)`.  `M0`, `M1` are the companion linearization of that polynomial
pencil, so that the quasinormal frequencies are the eigenvalues of

    M0 v = omega (-M1) v .

`q` is the dimensionless momentum k/2piT.
"""
function mtildes(m::Bea_Mateos, phiH::Real; q::Real = 0.0, N::Int = 40, bg = nothing,
                bc::Bool = true, rescale = (2, 0), divide_f::Bool = true,
                bcexpand::Bool = true)
    b = bg === nothing ? background(m, phiH; N = N) : bg
    # A diverged Newton solve returns a Background full of NaN, which otherwise
    # only surfaces much later as LAPACKException(360) out of `eigvals`.  The
    # cold-start ramp in `background` gives up somewhere above phiH ~ 1.6 for
    # typical potentials; pass `warm = <nearby Background>` to walk past that.
    b.converged || error("background at phiH = $phiH, N = $N did not converge " *
                         "(res = $(b.res)); no fluctuation matrices can be built. " *
                         "Try a warm start: background(m, $phiH; N = $N, warm = b_nearby).")
    _, _, Mt = build_Mtildes(b.x, b.D, b.D2, b.f0, b.g0, b.phi0, m.Vfuncs,
                             float(q), b.T, m.chi; nu = m.nu)
    divide_f && divide_eq1_by_f!(Mt, b)      # before the rescaling -- see its docstring
    rescale_fields!(Mt, b, rescale)          # before the boundary rows
    bc && boundary_rows!(Mt, b, m, float(q); horizon = !divide_f)
    # no-source on the EXPANSION rather than a Dirichlet value; the row it uses depends on
    # the rescaling exponent of Zphi -- see no_source_row!.  rescale = (0, 0) is the cheap case.
    bcexpand && no_source_row!(Mt, b, rescale[1])
    M0, M1 = companion(Mt)
    return M0, M1, Mt, b
end

"""
    divide_eq1_by_f!(Mt, b) -> maximum extraction residual

Divide the first equation by `f(u)`, removing the spurious overall factor the
export carries, so its collocation row at the horizon comes back to life and no
ingoing condition has to be imposed by hand.

Why it is needed: as exported, every coefficient of equation 1 vanishes one power
too fast at the horizon -- `(1-u)^2, (1-u), (1-u)` where a regular singular point
wants `(1-u), 1, 1` -- so that row reads `0 = 0` and has to be replaced.  The
replacement supplied by the derivation degenerates: its derivative part factors as
`(f'(1) sqrt(g(1)) + 4 i pi T lambda) * (...)`, and `f'(1) sqrt(g(1)) = -4 pi T`
identically, so it vanishes at `lambda = -i` for every background and momentum and
plants a spurious mode there.  Dividing the equation out instead removes both
problems at once.

Two subtleties:

  * `f` is *not* an explicit factor of every coefficient -- for the omega^0 Z2
    coefficient the vanishing follows from the background horizon conditions
    instead.  So the quotient is only finite on-shell, and the u = 1 value has to
    come from a limit: `c = f*chat` gives `chat(1) = c'(1)/f'(1)`.
  * This must run *before* `rescale_fields!`.  Afterwards the Zphi block is
    `D2*diag(u^2)`, which is not of the form `a D2 + b D + c I`, and the
    coefficient extraction used here would silently mangle it (~4% error).

ON BY DEFAULT.  Note that a symbolic division in the exporter is *not* available
as an alternative: with `f, g, phi, V` free symbols the omega^0 Z2 coefficient
genuinely does not vanish at u = 1, so `Cancel[eq1/f[u]]` cannot succeed without
first substituting the background equations of motion -- which the master
equations were already derived with.  Nor can a row operation help: the u = 1 row
is identically zero, so any combination of rows lies in the existing row space and
rank is invariant.  The missing condition `E1(1) = 0` is genuinely new
information, and the limit below is the way to get it.

Verified at N = 60, 90, 120 (the multi-resolution filter in `qnms` is part of
what makes it safe -- see the third bullet):

  * full rank *after* row balancing: 120 of 120 at N = 60, 180 of 180 at N = 90.
    Before balancing it reads 119 of 120, but the near-null direction is row
    N+2 -- equation 2 at the second gridpoint, whose row norm is ~1e-9 of the
    largest because of *its* own overall factor of `u` -- and has nothing to do
    with this division.  `_spectrum` balances, so the pencil actually handed to
    `eigvals` is nonsingular.
  * `lambda = -i` is gone and the Zphi tower comes back.  At phiH = 1e-3,
    q = 0.05 it resolves `Zphi n = 1, 2, 3` (1.100993 - 0.879452i,
    2.106942 - 1.887211i, 3.108251 - 2.888019i against the AdS5 Delta = 3 values
    1.099407 - 0.879767i, 2.105949 - 1.887444i, 3.107 - 2.891i), where the
    hand-imposed horizon row resolves only n = 1.  At u = 1 the divided
    coefficients come out as `ahat = 0`, `bhat = +0.7200`, `chat = -0.5382`:
    the second-derivative term drops, i.e. the equation degenerates to the
    first-order ingoing relation by itself, exactly as it should in EF
    coordinates.
  * it finds the spinodal instability under a *stricter* filter than the hand row
    can survive.  At phiH = 1.3126, q = 0.1 both give the physical mode at
    +0.10437i at rtol = 1e-3, but at rtol = 1e-4 only this one still has it.

REQUIRES `rtol <= 1e-4` FOR ANY STABILITY CLAIM.  This is the one genuine cost.
Dividing by a quantity that vanishes at the horizon takes cond(M) from ~2e5 to
~1.5e12 at N = 90, so eigenvalues carry an error of order 1e-4 and a nearly-real
mode can have its imaginary part smeared onto the wrong side of zero.  At
phiH = 0.3, q = 0.3 -- T/Lambda = 0.905, far above Tc, so certainly stable -- a
spurious pair at +-0.3389 + 9e-5i survives the three-grid filter at rtol = 1e-3
and looser, and is rejected at the default rtol = 1e-4.  (It is not a
singular-band artifact: it sits 0.28 k^2 from the band, well outside `bandtol`.)
An earlier revision of this docstring recorded that pair as a reason to leave the
factor in; it is instead a reason not to loosen `rtol`.  Measured at
Ns = (50,70,90), q = 0.3, phiH = 0.3:

    rtol    1e-2   1e-3   1e-4   1e-5
    fake       2      2      0      0

THE FAKE PAIR HAS TWO FACES, and only the near-real one is described above.  Its
`omega^2` is real and changes sign, so the same branch appears as:

  * `omega^2 > 0` : a near-real pair `+-R + i*eps` with `|eps| ~ 1e-4` of random
    sign -- the case above, a fake weakly damped or weakly growing oscillation;
  * `omega^2 < 0` : a pair `+-i y` on the imaginary axis, whose upper member is a
    LARGE FAKE INSTABILITY.  At phiM = 0.8, phiQ = 10, phiH = 1.0, chi = 0.6,
    q = 0 it sits at +0.52384i -- a bigger growth rate than the physical mode
    (-0.01126i) at the same point, so it captures any "least damped mode" seed.

Measured signatures of the branch (Ns = (60,70,80) unless stated), at
phiM = 0.8, phiQ = 10, q = 0:

  * the two members are related by `z -> -z`, NOT by the physical `z -> -conj(z)`.
    On the imaginary axis that is decisive: a physical purely imaginary mode is
    its OWN partner under `-conj`, so a partner at `-i y` means spurious.  A
    physical mirror pair (e.g. +-0.50056 - 0.02822i at phiH = 1.65, chi = 1) has
    both members at the same Im and is correctly untouched by that test.
  * `lambda^2` is real: -0.2744 (chi 0.6), -0.1502 (0.8), +0.0420 (1.05),
    +0.2824 (1.30) at phiH = 1.0, crossing zero at chi ~ 0.995 -- and the sign
    change is a square-root branch point, `+-0.92 sqrt(chi - chi_c)`, with the
    same coefficient on the imaginary side.
  * on the real branch `Im` does not converge, it goes to 0 with N:
    -6.31e-5, -4.95e-5, -4.02e-5, -3.37e-5 for finest N = 60, 70, 80, 90.
  * it is absent with `divide_f = false`, while the physical mode there agrees to
    1e-4 (-0.125884i vs -0.125985i at phiH = 1.0, chi = 0.8).

WHAT IT ACTUALLY IS (settled by looking at the eigenFUNCTION, not the eigenvalue):
a BOUNDARY-LAYER artifact pinned to the second gridpoint.  At phiM = 0.8, phiQ = 10,
phiH = 0.3, chi = 0.6, q = 0 the +0.2194i member has

    N     eigenvalue      peak at        |Z2| at gridpoints 2..6
   60    +0.21937i     i=2, u=7.09e-4    1.00 0.250 0.112 0.063 0.040
   70    +0.21936i     i=2, u=5.18e-4    1.00 0.250 0.111 0.063 0.040
   80    +0.21935i     i=2, u=3.95e-4    1.00 0.250 0.111 0.063 0.040

against a physical mode (+1.0888-0.8631i) which peaks at u ~ 0.57 for every N with all
neighbouring values ~0.52, i.e. smooth and resolved.

Two things follow.  The profile 1, 1/4, 1/9, 1/16, 1/25 combined with the Chebyshev-Lobatto
clustering u_i ~ (i-1)^2 near the edge (u3/u2 = 3.996) means the eigenfunction is exactly
`Z2 ~ 1/u`: the NON-NORMALISABLE branch at the AdS boundary, the source solution that the
no-source condition is meant to exclude.  Collocation cannot fully exclude it because 1/u is
only ever sampled at u2 > 0, so it survives as a spike of amplitude 1/u2.  And its support
shrinks as 1/N^2 while its eigenvalue is N-independent to 5 digits -- so it has no continuum
limit even though it "converges" beautifully.

THEREFORE NO EIGENVALUE-BASED TEST CAN REMOVE IT.  It is the same class as the singular_band
artifacts ("they do not move when N changes"), with a boundary rather than interior origin.
The discriminant has to be the eigenFUNCTION: reject modes whose amplitude peaks in the first
few gridpoints, or whose support shrinks with N.  That needs `eigen` rather than `eigvals` in
`_spectrum`.  Interim workarounds, all with known failure modes:

  * the `z -> -z` pairing test -- valid ONLY for on-axis modes at q = 0.  At q != 0 it
    misidentifies the SOUND MODE: at phiH = 1e-3, q = 0.05 the pair +-0.0289 - 0.0008i is
    physical (0.0289/0.05 = 0.578 vs the conformal 1/sqrt(3) = 0.5774) and gets flagged,
    because for small |Im| the maps -z and -conj(z) coincide.
  * intersecting the `divide_f = true` and `false` spectra -- 0 leakage measured at
    phiH = 1.0 and 1.35, but it discards the Zphi tower (present only in `true`), leaves
    holes wherever `false` drops a light mode, and degrades where the spurious branch
    crosses a physical one (phiH = 1.35, chi ~ 0.5: the two formulations disagree even on
    the sign, -0.43064i vs +0.43196i).

`rtol <= 1e-4` removes it, as the table says, but is NOT a blanket fix: the
acceptance test is `err < rtol*max(1,|lambda|)`, i.e. ABSOLUTE for |lambda| < 1,
and the fake pair's absolute drift lies in the same 1e-4..1e-3 band as genuine
light modes.  At phiH = 1.25, chi = 0.9 rtol = 1e-4 deletes the fake pair
(+-0.244) AND the physical mode (+0.16239).  A relative-drift cut does not
separate them either (fake 2e-4..8e-4, physical light modes 2.6e-3..5.7e-3 -- the
physical ones are WORSE).  Use the `z -> -z` pairing test above, or require the
mode to survive `divide_f = false`.

Float64 also runs out of digits near N ~ 120, so prefer the largest `Ns` that
still converges over pushing the grid further.  At rtol = 1e-4 the deepest Z2
overtone drops out (Z2 2 of 3, against 3 of 3 at rtol = 1e-2) -- that trade is
worth it, and the dispersion and growth-rate helpers keep a loose `rtol` on
purpose, since they track one branch by position rather than by sign of Im.
"""
function divide_eq1_by_f!(Mt::Vector{<:AbstractMatrix}, b::Background)
    N = b.N; D = b.D; D2 = b.D2; f = b.f0
    fp1 = (D * f)[end]
    resid = 0.0
    for n in eachindex(Mt), j in 1:2
        cols = (j-1)*N+1 : j*N
        a = zeros(ComplexF64, N); bb = zeros(ComplexF64, N); c = zeros(ComplexF64, N)
        for i in 1:N
            ei = zeros(N); ei[i] = 1.0
            v = [D2[i, :] D[i, :] ei] \ Vector(Mt[n][i, cols])
            a[i], bb[i], c[i] = v
            resid = max(resid, norm(v[1].*D2[i, :] .+ v[2].*D[i, :] .+ v[3].*ei .-
                                    Mt[n][i, cols]))
        end
        ah = a ./ f; bh = bb ./ f; ch = c ./ f          # f[N] = 0 -> fixed next
        ah[N] = (D*a)[N] / fp1
        bh[N] = (D*bb)[N] / fp1
        ch[N] = (D*c)[N] / fp1
        for i in 1:N
            ei = zeros(N); ei[i] = 1.0
            Mt[n][i, cols] .= ah[i].*D2[i, :] .+ bh[i].*D[i, :] .+ ch[i].*ei
        end
    end
    return resid
end

"""
    no_source_row!(Mt, b, s; Deltap = 3) -> Mt

Impose the *no-source* condition on Zphi as a condition on its boundary EXPANSION,
in place of a Dirichlet value at u = 0.

Why the Dirichlet row is not enough.  The two near-boundary branches are
`Zphi ~ A u^Delta_- + B u^Delta_+` with `Delta_- = d - Delta_+ = 1`, `A` the source
(non-normalizable) and `B` the vev (normalizable).  After the rescaling
`Zphi = u^s Ztilde` these become `A u^(1-s)` and `B u^(3-s)`.  For `s = 2`,
`rescale_fields!` argues the sourced one "diverges and is unrepresentable" -- true
in the continuum, false on a grid: collocation carries N *values*, and a vector
that is zero at u_1 = 0 and goes as 1/u from u_2 onwards satisfies
`Ztilde(u_1) = 0` while being pure sourced branch.  It survives as a spike at the
second gridpoint of amplitude 1/u_2 ~ N^2 whose eigenvalue is N-INDEPENDENT, so no
convergence test can see it (see `divide_eq1_by_f!`).

The condition that does exclude it is `A = 0`.  Using the two leading branches at
the first two interior points, `p = Delta_+ - s`,

    u_3^p Ztilde(u_2) - u_2^p Ztilde(u_3) = 0        <=>   Ztilde ~ u^p   <=>   A = 0

Row 1 -- equation 1 at u = 0 -- is identically zero (measured 2.6e-23 relative to the
largest row), so the condition goes into a FREE row and no collocation point is lost.
One row per field, exactly as at the horizon.

For `s >= 1` the column scaling `diag(u^s)` also sends column 1 to zero, so
`Ztilde(u_1)` appears in no equation at all.  That does NOT require a second row:
the pencil is then structurally singular and simply returns a few extra non-finite
eigenvalues, which `_spectrum` already discards via `isfinite`.  Measured at
phiH = 1.3, chi = 1, N = 60: 346 finite of 360 for `s = 2` against 349 of 360 for
`s = 0`, with the physical spectrum identical to twelve digits.  An earlier revision
of this docstring claimed the orphan had to be pinned and the condition pushed into
row 2, sacrificing equation 1 at u_2; that was wrong.

Measured (phiM = 0.8, phiQ = 10, q = 0), spread of a mode over N = 50..90:

    mode                          Dirichlet   s=0        s=2
    switch-off, chi .6 phiH 1.0    9.49e-06    2.09e-10   1.85e-10
    on-axis,    chi 1  psiH 1.10   5.14e-05    1.12e-10   1.11e-10
    unstable,   chi 1  psiH 1.30   1.38e-05    3.59e-10   3.58e-10
    AdS tower n=1                  1.99e-08    1.87e-08   1.87e-08

so the light modes gain four to five orders, and `s = 0` and `s = 2` are equivalent
(`s = 2` marginally better on the lightest mode).  Keep the default `s = 2`.  `max|Re lambda|` of an on-axis mode drops from
3.8e-4 to 7e-13 (it must be 0, so this is the error monitor).  Unchanged: the AdS5
Delta = 3 tower (error 3.5e-7 and 4.7e-7 against 1.099407-0.879767i and
2.105949-1.887444i), the hydrodynamic sound mode (c_s = 0.57774 vs the conformal
0.57735 at phiH = 1e-3, q = 0.05), and five psiH reference modes at chi = 1, all to
six digits.  The spurious `+-i y` / `+-R` pair documented in `divide_eq1_by_f!`
disappears -- it *is* the sourced branch: at phiH = 1.0 it is present for
chi = 0.4, 0.6, 0.8 with Dirichlet and absent for every chi with this row.

`cond(M0)` cannot be used to compare the variants: the companion linearisation is
singular by construction (infinite eigenvalues), so it is `Inf` for all of them.
The convergence spread above is the usable proxy.

Caveats.  (i) Two branches only, so the row is accurate to `O(u_2^(2))`; a
three-point version against `A u^(1-s) + B u^(3-s) + C u^(5-s)` would be sharper.
(ii) `Z2` is not rescaled and its sourced branch is a constant, which plain
Dirichlet does exclude -- consistent with the artifact living entirely in Zphi
(`|Z2| = 5e-17` in its eigenvector), though Z2's branch powers were not
re-derived here.  (iii) It changes every light-mode value in the fourth digit.
ON by default since 2026-11-02: without it, at phiM = 0.7, chi = 1, phiH = 1.30,
q = 0 the N = 30 grid puts the unstable gap mode at 0.3183i instead of 0.3119i,
so `qnms` on `Ns = (30, 40, 50)` dropped it for every phiH in [1.20, 1.70].  With
it, N = 30 already gives seven digits.
"""
function no_source_row!(Mt::Vector{<:AbstractMatrix}, b::Background, s::Integer;
                        Deltap::Integer = 3)
    u  = b.x
    pw = Deltap - s               # Ztilde ~ A u^(1-s) + B u^(pw);  A = 0 is imposed
    for X in Mt                   # row 1 (eq 1 at u = 0) is identically zero: use it
        X[1, :] .= 0
    end
    Mt[1][1, 2] =  u[3]^pw
    Mt[1][1, 3] = -u[2]^pw
    return Mt
end

function rescale_fields!(Mt::Vector{<:AbstractMatrix}, b::Background, s)
    N = b.N
    for (j, sj) in enumerate(s)
        sj == 0 && continue
        cols = (j-1)*N+1 : j*N
        Sj = Diagonal(b.x .^ sj)
        for X in Mt
            X[:, cols] .= X[:, cols] * Sj
        end
    end
    return Mt
end

"""
    boundary_rows!(Mt, b, m, q) -> Mt

Overwrite the three collocation rows that the bulk equations leave empty.

At u = 0 the second equation degenerates to `0 = 0` (all of its coefficients
vanish like u^4, u^3, u^2), and at u = 1 the first equation does the same (its
coefficients vanish like (1-u)^2, (1-u), (1-u)).  Without this surgery
`det M(omega) = 0` identically and the pencil has no spectrum at all.

The replacements are the ones the derivation supplies:

  * row 1     -> `Zphi(0) = 0`, row N+1 -> `Z2(0) = 0`, both omega-independent;
  * row N     -> the horizon regularity condition, which is *cubic* in omega and
    so is spread across all four `Mt[n]`.  It is first order (no second
    derivative), and it does not involve `Z2'`.

Convention conversions from the derivation to the exported coefficients: its
frequency variable is the dimensionless `lambda = omega/2piT`, hence the
`(2 pi T)^(n-1)`; and its potential is `V -> chi*V - 3`, so `V' -> sqrt(chi) V'`
and `V'' -> V''`.
"""
function boundary_rows!(Mt::Vector{<:AbstractMatrix}, b::Background,
                        m::Bea_Mateos, q::Real; horizon::Bool = true)
    N = b.N; T = b.T; chi = m.chi; D = b.D

    # --- u = 0 : Dirichlet on both fields -----------------------------------
    for X in Mt
        X[1, :]   .= 0
        X[N+1, :] .= 0
    end
    Mt[1][1, 1]     = 1
    Mt[1][N+1, N+1] = 1

    # --- u = 1 : horizon regularity, one condition, cubic in omega ----------
    # Only needed when equation 1 still carries its spurious factor of f, which
    # leaves this row dead.  With `divide_eq1_by_f!` applied the row is alive and
    # supplies the ingoing condition itself, so imposing this one would only put
    # back the spurious mode at lambda = -i.
    horizon || return Mt
    phiH = b.x[end]^m.nu * b.phi0[end]
    VH   = chi * m.V(phiH) - 3
    VpH  = sqrt(chi) * m.Vp(phiH)
    VppH = m.Vpp(phiH)
    df   = (D * b.f0)[end]                 # f'(1)
    gH   = b.g0[end]

    BH = ComplexF64[ 8π^2*q^2*T^2 * df^4 * gH^(3/2) * VH,
                    32im*π^3*q^2*T^3 * df^3 * gH     * VH,
                    64π^2*T^2       * df^2 * gH^(3/2) * VH^2,
                   256im*π^3*T^3    * df   * gH      * VH^2]

    W = -gH*VpH^2 + VH*(4π^2*q^2*T^2 - df*gH + gH*VppH)
    C1 = ComplexF64[ -8π^2*q^2*T^2 * df^3 * sqrt(gH) * W,
                     32im*π^3*q^2*T^3 * df^2 * gH * (df - 2VH) * VH,
                    -64π^2*T^2 * df * sqrt(gH) * VH * W,
                   -256im*π^3*T^3 * gH * VH^2 * (-df + 2VH)]

    U = VpH*(VpH^2 - VH*VppH)
    C2 = ComplexF64[df*sqrt(gH)*U, -4im*π*T*U, 0, 0]

    eN = zeros(N); eN[end] = 1.0
    for n in 1:length(Mt)
        s = (2π*T)^(n-1)
        Mt[n][N, :] .= 0
        Mt[n][N, 1:N]        .= (BH[n] .* D[end, :] .+ C1[n] .* eN) ./ s
        Mt[n][N, N+1:2N]     .= (C2[n] .* eN) ./ s        # no Z2' in the condition
    end
    return Mt
end

"""
    companion(Mt) -> M0, M1

Companion linearization of `Mt[1] + omega Mt[2] + ... + omega^P Mt[P+1]`, so
that the frequencies are the eigenvalues of `M0 v = omega (-M1) v`.
"""
function companion(Mt::Vector{<:AbstractMatrix})
    P = length(Mt) - 1
    d = size(Mt[1], 1)
    P == 0 && error("the equations do not depend on omega")
    P == 1 && return Matrix{ComplexF64}(Mt[1]), Matrix{ComplexF64}(Mt[2])
    Idb = Matrix{ComplexF64}(I, d, d)
    M0 = zeros(ComplexF64, P*d, P*d)
    M1 = zeros(ComplexF64, P*d, P*d)
    blk(i, j) = ((i-1)*d+1:i*d, (j-1)*d+1:j*d)
    for j in 1:P
        M0[blk(1, j)...] .= Mt[j]
    end
    for i in 2:P
        M0[blk(i, i)...]   .=  Idb
        M1[blk(i, i-1)...] .= -Idb
    end
    M1[blk(1, P)...] .= Mt[P+1]
    return M0, M1
end

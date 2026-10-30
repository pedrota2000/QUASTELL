"""
    QUASTELL -- QUAsinormal Spectra of ThErmaL hoLography

One object, `Bea_Mateos(phiM, phiQ)`, that owns the whole pipeline:

    m = Bea_Mateos(0.8, 10; chi = 1.2)

    bg  = background(m, 0.3; N = 40)        # 1. one background solution
    th  = thermodynamics!(m; N = 60)        # 2. the full thermodynamic curve
    M0, M1, Mt = mtildes(m, 0.3; q = 0.5)   # 3. the Mtilde / companion matrices
    r   = qnms(m, 0.3; q = 0.5)             # 4. the quasinormal frequencies

    plot_potential(m); plot_solution(m, bg); plot_thermodynamics(m)
    plot_qnms(r); plot_convergence(r); plot_dispersion(m, 0.3)

Conventions
-----------
* Radial coordinate u in [0,1], x[1] = 0 (boundary), x[end] = 1 (horizon),
  on the Chebyshev-Gauss-Lobatto grid of `cheby`.
* The scalar is stored rescaled: phi(u) = u^nu * phi0(u), nu = 1.
* Background units: f(0) = g(0) = 1, so the source Lambda = phi0(0) sets the
  scale and T, s are measured in it (T/Lambda, s/Lambda^3).
* The QNM eigenvalue returned by the Mtilde pencil is omega itself, in those
  same units; `q` is the *dimensionless* momentum k/2piT (that is how it
  enters the generated coefficients, as k = 2 pi q T).  The physical output
  is the dimensionless frequency lambda = omega / 2 pi T.

Physical modes vs numerical artifacts: `qnms` never trusts a single grid.  It
recomputes background + spectrum at every N in `Ns` and keeps only the
eigenvalues that reappear, unmoved, in all of them (Jansen, arXiv:1709.09178,
eq. 3.2).  Everything else is discarded as an artifact.
"""
module QUASTELL

using LinearAlgebra, Printf
using Optim, Plots, ProgressMeter
using DataInterpolations, QuadGK

export Bea_Mateos, Background, QNMResult
export cheby, background, thermodynamics!, critical_temperature
export mtildes, build_Mtildes, companion,
       divide_eq1_by_f!, no_source_row!, rescale_fields!, boundary_rows!
export principal_zeros, singular_band, qnms, show_modes,
       dispersion, sound_speed, growth_rate
export plot_potential, plot_solution, plot_thermodynamics, plot_qnms,
       plot_convergence, plot_dispersion, plot_growth_rates

include("model.jl")            # containers and the Bea_Mateos model
include("grid.jl")             # Chebyshev grid, progress-bar helper
include("background.jl")       # 1. background at one phiH
include("thermodynamics.jl")   # 2. s(T), free energy, Tc
include("fluctuations.jl")     # 3. the Mtilde matrices and boundary surgery
include("qnms.jl")             # 4. spectrum, multi-grid filter, dispersion
include("plotting.jl")

end # module

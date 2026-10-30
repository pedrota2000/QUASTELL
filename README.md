# QUASTELL

```text
    \o/                                                                               \o/
     |                                                                                 |
   o/ \o        ___    _   _      _      ____    _____   _____   _       _           o/ \o
  /|\ /|\      / _ \  | | | |    / \    / ___|  |_   _| | ____| | |     | |         /|\ /|\
 o/ \o/ \o    | | | | | | | |   / _ \   \___ \    | |   |  _|   | |     | |        o/ \o/ \o
/|\ /|\ /|\   | |_| | | |_| |  / ___ \   ___) |   | |   | |___  | |___  | |___    /|\ /|\ /|\
/ \ / \ / \    \__\_\  \___/  /_/   \_\ |____/    |_|   |_____| |_____| |_____|   / \ / \ / \
                          QUAsinormal Spectra of ThErmaL hoLography
```

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**QUA**sinormal **S**pectra of **T**h**E**rma**L** ho**L**ography.

The name also reads as *qua + castell*: QNMs come in towers, like the Catalan human towers.

<p align="center">
  <img src="assets/gapmode.gif" alt="The gapped order-parameter mode of a black brane turning unstable across a first-order transition" width="900">
</p>

*A black brane is carried across the first-order transition (left, φ<sub>M</sub> = 0.8, χ = 1),
and its gapped modes are followed at zero momentum (right). Between the two zeros of ds/dT the
order-parameter mode crosses into Im ω > 0. It is unstable there although c<sub>v</sub> > 0,
and becomes stable again on the cold branch.*

QUASTELL is a Julia pipeline for the quasinormal modes (QNMs) of five-dimensional holographic
black branes with a self-interacting scalar field. It goes from the bulk potential to the
background geometry, then to the thermodynamics and its first-order phase transition, and
finally to the QNM spectrum and its instabilities:

```
V(φ)  ──►  background  ──►  thermodynamics, Tc  ──►  fluctuation matrices  ──►  QNMs
            (Newton)          (swallowtail)            (cubic pencil)          (multi-grid filter)
```

Everything is pseudospectral (Chebyshev collocation). Spurious eigenvalues are separated from
physical modes by re-solving the problem on several grids and by removing the known numerical
artifacts of the master equations (see [Physical modes and artifacts](#physical-modes-and-artifacts)).

## The model

QUASTELL implements the model of Bea and Mateos
([arXiv:1805.01806](https://arxiv.org/abs/1805.01806)): Einstein gravity in AdS<sub>5</sub>
coupled to a scalar, with potential derived from the superpotential (their eqs. 2.2–2.3, with L = 1)

$$
W(\phi) = -\frac{3}{2} - \frac{\phi^2}{2} - \frac{\phi^4}{4\phi_M^2} + \frac{\phi^6}{\phi_Q},
\qquad
V(\phi) = -\frac{4}{3}W(\phi)^2 + \frac{1}{2}W'(\phi)^2 .
$$

The parameters `phiM` and `phiQ` set the potential. Expanding around φ = 0 gives
V = −3 − (3/2)φ² + …, so the scalar has m² = −3 (L = 1) and is dual to an operator of
dimension Δ = 3. Its source Λ is a relevant deformation that breaks conformal invariance.
For suitable `phiM`, the equation of state s(T) is multivalued.
That gives a first-order transition at Tc, with a spinodal branch where the plasma is
dynamically unstable.

In the code, `m.V` is V(φ) + 3, and the potential enters the equations as `chi*V(phi) - 3`.
The parameter `chi` controls the backreaction of the scalar on the geometry; `chi = 1` is the
model above.

**Conventions**

| | |
|---|---|
| radial coordinate | u ∈ [0, 1], u = 0 the AdS boundary, u = 1 the horizon (ingoing Eddington–Finkelstein) |
| scalar | φ(u) = u<sup>ν</sup> φ₀(u), ν = 1; branches are labelled by the horizon value `phiH` |
| units | f(0) = g(0) = 1, source Λ = φ₀(0); thermodynamics is reported as T/Λ and s/Λ³ |
| frequencies | λ = ω/2πT |
| momentum | q = k/2πT |
| fields | two gauge-invariant master fields, Z<sub>φ</sub> (scalar channel) and Z₂ |

## Installation

QUASTELL needs Julia ≥ 1.7. Clone the repository and install its dependencies:

```bash
git clone https://github.com/<your-user>/QUASTELL.git
cd QUASTELL
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Then, from Julia started with `julia --project` in that folder:

```julia
using QUASTELL
```

The first `using` precompiles Plots and takes a few minutes. Later sessions start much faster.

In the REPL, in Jupyter and in VS Code notebooks, loading the package greets you with a castell. To skip it, set
`ENV["QUASTELL_BANNER"] = "false"` before `using QUASTELL`.

## Quick start

```julia
using QUASTELL

m  = Bea_Mateos(1.0, 10; chi = 1.2)            # phiM, phiQ, backreaction chi

bg = background(m, 0.3; N = 40)                # one black brane, labelled by phiH
bg.T / bg.Lambda                               # 0.9087

th = thermodynamics!(m; N = 50, nphi = 2000)   # the s(T) curve; m.Tc = 0.3735675
plot_thermodynamics(m)

# phiH = 1.16 lies on the spinodal branch (ds/dT < 0 for phiH in [1.033, 1.271])
r  = qnms(m, 1.16; q = 0.05, Ns = (40, 50, 60), rtol = 1e-4, bcexpand = true)
r.lambda[1]                                    # 6.8e-12 + 0.011068im
plot_qnms(r)
```

Modes are sorted least damped first, so `r.lambda[1]` is the most unstable one. Here it is
purely imaginary with Im λ > 0: in the spinodal region c<sub>s</sub>² < 0, and the sound mode
grows instead of propagating. At `phiH = 1.40`, on the stable branch, the same call gives
ordinary damped sound, λ = ±0.01939 − 0.00096 i.

[`notebooks/QUASTELL_demo.ipynb`](notebooks/QUASTELL_demo.ipynb) walks through every function,
one per cell. It also includes the validation against AdS<sub>5</sub>-Schwarzschild, the
dispersion relations and the growth rate of the instability.

## What is in the package

Every function has a docstring with the details (`?qnms` at the REPL).

| Step | Function | Result |
|---|---|---|
| model | `Bea_Mateos(phiM, phiQ; chi = 1, nu = 1)` | V, V′, V″, the bound `phiMax`; caches backgrounds, thermodynamics and QNM runs |
| 1. background | `background(m, phiH; N, warm)` | `Background`: f, g, φ₀ on the grid, plus T, S, Λ |
| 2. thermodynamics | `thermodynamics!(m; N, nphi)` | phiH, T/Λ, s/Λ³, pressure, free energy, energy density; `m.Tc` |
| | `critical_temperature(s, T, phiH)` | Tc from the self-intersection of the free energy (swallowtail) |
| 3. fluctuations | `mtildes(m, phiH; q, N)` | `Mt` (coefficients of M(ω) = Σ ωⁿ Mₙ, n ≤ 3) and the companion pair `M0`, `M1` |
| 4. QNMs | `qnms(m, phiH; q, Ns, rtol, cutoff)` | `QNMResult`: converged λ sorted by Im λ, their drift `err`, and the raw spectra |
| | `dispersion(m, phiH; qs, nmodes)` | the `nmodes` least-damped branches λₙ(q), followed by continuation |
| | `sound_speed(d)` | c<sub>s</sub> fitted from a dispersion result |
| | `growth_rate(m, phiH; qs)` | γ(k)/Λ of the unstable mode, with the hydrodynamic fit \|c<sub>s</sub>\|k − Γk²/2 |
| | `singular_band(m, b, q, T)` | the band of ω where the master equations acquire an interior singular point |
| plots | `plot_potential`, `plot_solution`, `plot_thermodynamics`, `plot_qnms`, `plot_convergence`, `plot_dispersion`, `plot_growth_rates` | Plots.jl figures |

The lower-level pieces of the fluctuation problem are also exported: `build_Mtildes`,
`divide_eq1_by_f!`, `rescale_fields!`, `boundary_rows!`, `no_source_row!`, `companion`.

## Deriving the equations

[`mathematica/Derive_the_equations.nb`](mathematica/Derive_the_equations.nb) derives the
equations for the two master fields at finite ω and k, simplified on-shell with the background
equations of motion. It works first in Fefferman–Graham coordinates and then in general
Eddington–Finkelstein coordinates, with a varying backreaction parameter. Its last part exports
the coefficients of the master equations as Julia code (`Mtildes.jl`, written next to the
notebook).

The coefficient block of `build_Mtildes` in [`src/fluctuations.jl`](src/fluctuations.jl) is
that export, pasted in. If the equations change, re-run the notebook and replace the block.

## Physical modes and artifacts

The raw spectrum of a collocation pencil is mostly numerical artifacts. `qnms` removes them in
several layers.

1. **Multi-resolution test.** The background and the spectrum are recomputed from scratch on
   every grid in `Ns` (never interpolated). A frequency is kept only if it reappears, within
   `rtol * max(1, |λ|)`, on all of them. This is the test of Jansen's
   [QNMspectral](https://github.com/APJansen/QNMspectral)
   ([arXiv:1709.09178](https://arxiv.org/abs/1709.09178), eq. (3.2)).
2. **Singular band.** At q ≠ 0, the coefficient of the second derivative vanishes at
   ω² = a(u) for some u. In that band, collocation produces one nearly real pair per grid
   point that does *not* move with N. These pairs are removed by location (`drop_band`,
   `bandtol`). At q = 0 the band collapses onto ω = 0, so `|λ| < origintol` is dropped. The
   hydrodynamic mode also sits at ω = 0 there, so use a small q > 0 to see it.
3. **Spurious factor of f in equation 1** (`divide_f = true`, the default). As exported,
   equation 1 carries an extra overall factor of f(u), which kills its row at the horizon. It
   is divided out, with the horizon value taken as a limit. This removes a fake mode at λ = −i
   and restores the Z<sub>φ</sub> tower. The cost is conditioning, so **only trust the sign of
   Im λ with `rtol <= 1e-4`**.
4. **AdS boundary.** Z<sub>φ</sub> is rescaled as u² Z̃ (`rescale = (2, 0)`). With
   `bcexpand = true`, the no-source condition is imposed on the boundary *expansion* instead of
   as a Dirichlet value. That removes an N-independent artifact (the non-normalisable 1/u
   branch, living on the second gridpoint) and improves the convergence of the light modes by
   four to five orders of magnitude. It is off by default, but recommended for light modes.
   The rescaling zeroes the u = 0 column of the pencil, which plants an exact eigenvalue
   λ = 0 on every grid. The raw spectrum drops it (|λ| < 10⁻¹²), so it can't pass the test
   or capture a branch in `dispersion`.
5. **Balancing.** The pencil is solved for λ rather than ω, with its rows equilibrated. Both
   are exact transformations that leave the spectrum unchanged.

Practical notes:

- Float64 runs out of digits around N ≈ 120.
- From N ≈ 70 the background Newton residual floors above its default `tol = 1e-10` (about
  2 × 10⁻¹⁰ at N = 70 and 10⁻⁹ at N = 90). Each continuation step then runs to `maxit`. The
  solution is still accepted, since `converged` only asks for 10⁻⁸, but it takes minutes
  instead of seconds. Grids of N ≤ 60 avoid this.
- Cold-start backgrounds fail above `phiH ≈ 1.6`. For larger values, pass a nearby solution
  as `warm = ...` and move in small steps.
- `plot_convergence(r)` shows the raw spectra of all grids side by side. Physical modes line
  up across the panels; artifacts do not.

## Validation

In the limit `phiH → 0` the scalar decouples, and each channel has to reproduce its
AdS<sub>5</sub>-Schwarzschild tower: Z₂ the Δ = 4 one and Z<sub>φ</sub> the Δ = 3 one
(`phiH = 1e-3`, q = 0, `Ns = (60, 70, 90)`):

| mode | QUASTELL | reference | \|difference\| |
|---|---|---|---|
| Z<sub>φ</sub>, n = 1 | 1.099407474 − 0.879767076 i | 1.099407283 − 0.879767314 i | 3.1 × 10⁻⁷ |
| Z₂, n = 1 | 1.559725860 − 1.373337493 i | 1.559725796 − 1.373337872 i | 3.8 × 10⁻⁷ |
| Z<sub>φ</sub>, n = 2 | 2.105948634 − 1.887443607 i | 2.105948627 − 1.887444155 i | 5.5 × 10⁻⁷ |
| Z₂, n = 2 | 2.584760685 − 2.381784179 i | 2.584760 − 2.381785 i | 1 × 10⁻⁶ (reference to 6 digits) |

The hydrodynamic sound mode must approach λ = q/√3 − i q²/3. At q = 0.05 the package gives
λ = 0.028886922 − 0.000833145 i, against 0.028867513 − 0.000833333 i. The residual difference
is the O(q³) correction.

## Tests

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

The test suite checks the following:

- the Chebyshev grid;
- a background solution and its thermodynamics;
- the structure of the fluctuation matrices;
- the AdS<sub>5</sub> towers above;
- the conformal sound mode;
- the critical temperature, and the extent of the spinodal branch;
- the growing sound mode inside the spinodal branch, and ordinary damped sound outside it.

It takes about six minutes, half of it precompilation.

## Repository layout

```
QUASTELL/
├── LICENSE
├── Project.toml
├── src/
│   ├── QUASTELL.jl           module definition and exports
│   ├── model.jl            Bea_Mateos, Background, QNMResult
│   ├── grid.jl             Chebyshev–Gauss–Lobatto grid, progress-bar helper
│   ├── background.jl       damped Newton solver for the black brane
│   ├── thermodynamics.jl   s(T), free energy, Tc
│   ├── fluctuations.jl     fluctuation matrices, boundary conditions, companion form
│   ├── qnms.jl             spectrum, multi-grid filter, dispersion, growth rate
│   └── plotting.jl         figures (Plots.jl)
├── test/
│   └── runtests.jl
├── notebooks/
│   └── QUASTELL_demo.ipynb   tour of the whole pipeline
├── assets/
│   ├── gapmode.gif           the animation at the top of this README
│   └── make_gif.jl           regenerates it
└── mathematica/
    └── Derive_the_equations.nb   derivation of the equations, and the Julia exporter
```

## Acknowledgements and references

QUASTELL builds on the following work, and if you use it, please cite these as well.

- **The model.** Y. Bea and D. Mateos, *Heating up Exotic RG Flows with Holography*,
  JHEP **08** (2018) 034, [arXiv:1805.01806](https://arxiv.org/abs/1805.01806).

- **QNMspectral**, the Mathematica package by Aron Jansen for computing quasinormal modes with
  spectral methods ([github.com/APJansen/QNMspectral](https://github.com/APJansen/QNMspectral)).
  The multi-resolution convergence test in `qnms` is the one introduced there.

  A. Jansen, *Overdamped modes in Schwarzschild-de Sitter and a Mathematica package for the
  numerical computation of quasinormal modes*, Eur. Phys. J. Plus **132**, 546 (2017),
  [arXiv:1709.09178](https://arxiv.org/abs/1709.09178).

- The work of Mikel Sánchez-Garitaonandia and Carlos Hoyos on the dynamics near the spinodal
  point of holographic first-order phase transitions.

  A. Caddeo, O. Henriksson, C. Hoyos and M. Sánchez-Garitaonandia, *Spinodal slowing down and
  scaling in a holographic model*, [arXiv:2406.15297](https://arxiv.org/abs/2406.15297).

## License

QUASTELL is released under the [MIT License](LICENSE).

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

Julia code for the quasinormal modes (QNMs) of five-dimensional holographic black branes with a
self-interacting scalar. The name also reads *qua + castell*: QNMs come in towers, like the
Catalan human towers.

<p align="center">
  <img src="assets/gapmode.gif" alt="The gapped order-parameter mode of a black brane turning unstable across a first-order transition" width="900">
</p>

*φ<sub>M</sub> = 0.8, χ = 1. Left: the equation of state across the first-order transition.
Right: the gapped QNMs at k = 0. Between the two zeros of ds/dT the order-parameter mode is
unstable, even though c<sub>v</sub> > 0.*

```
V(φ)  ──►  background  ──►  thermodynamics, Tc  ──►  fluctuation matrices  ──►  QNMs
            (Newton)          (swallowtail)            (cubic pencil)          (multi-grid filter)
```

Everything is Chebyshev collocation. Spurious eigenvalues are removed by solving on several
grids and by filtering the known artifacts of the master equations
([below](#spurious-modes)).

## The model

The model of Bea and Mateos ([arXiv:1805.01806](https://arxiv.org/abs/1805.01806)):
Einstein gravity in AdS<sub>5</sub> coupled to a scalar with potential (L = 1)

$$
W(\phi) = -\frac{3}{2} - \frac{\phi^2}{2} - \frac{\phi^4}{4\phi_M^2} + \frac{\phi^6}{\phi_Q},
\qquad
V(\phi) = -\frac{4}{3}W(\phi)^2 + \frac{1}{2}W'(\phi)^2 .
$$

V = −3 − (3/2)φ² + …, so φ is dual to a Δ = 3 operator, whose source Λ breaks conformal
invariance. For suitable `phiM` the equation of state s(T) is multivalued: there is a
first-order transition at Tc and a spinodal branch where the plasma is unstable.

The metric is written in ingoing Eddington–Finkelstein coordinates, as in the Mathematica
notebook:

$$
ds^2 = \frac{1}{u^2}\left[-f(u) \, g(u) \, dt^2 - 2\sqrt{g(u)} \, dt \, du + g(u)\left(dx^2 + dy^2 + dz^2\right)\right],
\qquad
\phi = u \, \phi_0(u).
$$

The boundary is at u = 0 (f = g = 1) and the horizon at u = 1 (f = 0), so
T = −f′(1)√g(1)/4π. Branches are labelled by `phiH` = φ₀(1), and Λ = φ₀(0).

`chi` multiplies the scalar stress tensor, G<sub>μν</sub> − 6g<sub>μν</sub> = χT<sub>μν</sub>;
`chi = 1` is the model above. In the code `m.V` is V + 3.

The QNMs are computed for two gauge-invariant master fields, Z<sub>φ</sub> and Z₂, in units
λ = ω/2πT and q = k/2πT. Thermodynamics is in units of Λ.

## Installation

```bash
git clone https://github.com/<your-user>/QUASTELL.git
cd QUASTELL
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Needs Julia ≥ 1.7. The first `using QUASTELL` precompiles Plots and takes a few minutes. Set
`ENV["QUASTELL_BANNER"] = "false"` before loading to hide the banner.

## Quick start

```julia
using QUASTELL

m  = Bea_Mateos(1.0, 10; chi = 1.2)            # phiM, phiQ, chi
bg = background(m, 0.3; N = 40)                # black brane with phiH = 0.3
bg.T / bg.Lambda                               # 0.9087

thermodynamics!(m; N = 50, nphi = 2000)        # s(T); sets m.Tc = 0.3735675
plot_thermodynamics(m)

# phiH = 1.16 is on the spinodal branch, phiH in [1.033, 1.271]
r = qnms(m, 1.16; q = 0.05, Ns = (40, 50, 60), rtol = 1e-4, bcexpand = true)
r.lambda[1]                                    # 6.8e-12 + 0.011068im
plot_qnms(r)
```

Modes are sorted by Im λ, so `r.lambda[1]` is the most unstable one: sound with
c<sub>s</sub>² < 0, which grows instead of propagating. At `phiH = 1.40` the same call gives
damped sound, λ = ±0.01939 − 0.00096i.

For the eigenfunctions, add `eigenfunctions = true`:

```julia
r = qnms(m, 1.16; q = 0.05, Ns = (40, 50, 60), rtol = 1e-4, bcexpand = true,
         eigenfunctions = true)
plot_eigenfunctions(r; modes = 1:3)
```

The fields are stored in `r.Zphi` and `r.Z2` on the grid `r.u`.

[`notebooks/QUASTELL_demo.ipynb`](notebooks/QUASTELL_demo.ipynb) goes through every function.

## Functions

See the docstrings for details (`?qnms`).

| | Function | Returns |
|---|---|---|
| model | `Bea_Mateos(phiM, phiQ; chi = 1, nu = 1)` | V, V′, V″; caches backgrounds, thermodynamics and QNM runs |
| background | `background(m, phiH; N, warm)` | f, g, φ₀ on the grid; T, S, Λ |
| thermodynamics | `thermodynamics!(m; N, nphi)` | s(T), pressure, free energy, energy density; sets `m.Tc` |
| | `critical_temperature(s, T, phiH)` | Tc from the swallowtail |
| fluctuations | `mtildes(m, phiH; q, N)` | the matrices of M(ω) = Σ ω<sup>n</sup> M<sub>n</sub>, n ≤ 3, and its companion form |
| QNMs | `qnms(m, phiH; q, Ns, rtol, cutoff, eigenfunctions)` | converged λ sorted by Im λ, their error, the raw spectra; optionally Z<sub>φ</sub>(u), Z₂(u) |
| | `dispersion(m, phiH; qs, nmodes)` | λ(q) of the least damped modes |
| | `sound_speed(d)` | c<sub>s</sub> from a dispersion |
| | `growth_rate(m, phiH; qs)` | γ(k) of the unstable sound mode, with the fit \|c<sub>s</sub>\|k − Γk²/2 |
| | `singular_band(m, b, q, T)` | the band of ω where the master equations have an interior singular point |
| plots | `plot_potential`, `plot_solution`, `plot_thermodynamics`, `plot_qnms`, `plot_convergence`, `plot_dispersion`, `plot_growth_rates`, `plot_eigenfunctions` | |

The lower-level steps are exported too: `build_Mtildes`, `divide_eq1_by_f!`,
`rescale_fields!`, `boundary_rows!`, `no_source_row!`, `companion`.

## Deriving the equations

[`mathematica/Derive_the_equations.nb`](mathematica/Derive_the_equations.nb) derives the master
equations at finite ω and k and exports their coefficients as Julia (`Mtildes.jl`). That export
is the coefficient block of `build_Mtildes` in [`src/fluctuations.jl`](src/fluctuations.jl). If
the equations change, re-run the notebook and paste it in again.

## Spurious modes

Most eigenvalues of the discretised problem are artifacts. `qnms` removes them as follows.

1. **Several grids.** A frequency is kept only if it appears on every grid in `Ns`, within
   `rtol * max(1, |λ|)`. This is the test of Jansen's
   [QNMspectral](https://github.com/APJansen/QNMspectral).
2. **Singular band.** For q ≠ 0, collocation puts N-independent fake modes where the master
   equations have an interior singular point. They are removed by location (`bandtol`). At
   q = 0 the band sits at ω = 0, so |λ| < `origintol` is dropped. Use a small q > 0 to see the
   hydrodynamic mode.
3. **Factor of f.** Equation 1 carries a spurious factor f(u), which gives a fake mode at
   λ = −i. `divide_f = true` (the default) removes it. Trust the sign of Im λ only with
   `rtol <= 1e-4`.
4. **Boundary condition.** `bcexpand = true` imposes the no-source condition on the boundary
   expansion. This removes a boundary artifact and improves the light modes by four to five
   orders of magnitude. It is off by default but recommended.

Notes:

- Keep N ≤ 60. Above N ≈ 70 the background solver stalls (minutes per step), and Float64 runs
  out of digits around N ≈ 120.
- Cold starts fail above `phiH ≈ 1.6`. Pass a nearby solution as `warm`.
- `plot_convergence(r)` shows the raw spectra of all grids. Physical modes line up; artifacts
  don't.
- A QNM has Z<sub>φ</sub> ~ u³ at the boundary. The artifact of step 4 has Z<sub>φ</sub> ~ u,
  which its `eferr` doesn't show.

## Validation

As `phiH → 0` the scalar decouples, and the two channels must reproduce the
AdS<sub>5</sub>-Schwarzschild towers: Δ = 3 for Z<sub>φ</sub> and Δ = 4 for Z₂. With
`phiH = 1e-3`, q = 0 and `Ns = (60, 70, 90)`:

| mode | QUASTELL | reference | difference |
|---|---|---|---|
| Z<sub>φ</sub>, n = 1 | 1.099407474 − 0.879767076i | 1.099407283 − 0.879767314i | 3.1 × 10⁻⁷ |
| Z₂, n = 1 | 1.559725860 − 1.373337493i | 1.559725796 − 1.373337872i | 3.8 × 10⁻⁷ |
| Z<sub>φ</sub>, n = 2 | 2.105948634 − 1.887443607i | 2.105948627 − 1.887444155i | 5.5 × 10⁻⁷ |
| Z₂, n = 2 | 2.584760685 − 2.381784179i | 2.584760 − 2.381785i | ~10⁻⁶ (6-digit reference) |

Sound at q = 0.05 gives λ = 0.028886922 − 0.000833145i. Hydrodynamics, q/√3 − iq²/3, gives
0.028867513 − 0.000833333i. The difference is O(q³).

## Tests

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

About six minutes, half of it precompilation.

## Repository layout

```
QUASTELL/
├── src/
│   ├── QUASTELL.jl         module and exports
│   ├── model.jl            Bea_Mateos, Background, QNMResult
│   ├── grid.jl             Chebyshev grid
│   ├── background.jl       Newton solver for the black brane
│   ├── thermodynamics.jl   s(T), free energy, Tc
│   ├── fluctuations.jl     fluctuation matrices and boundary conditions
│   ├── qnms.jl             spectrum, multi-grid filter, dispersion, growth rate
│   └── plotting.jl
├── test/runtests.jl
├── notebooks/QUASTELL_demo.ipynb
├── assets/                 gapmode.gif and the script that makes it
└── mathematica/Derive_the_equations.nb
```

## References

If you use QUASTELL, please also cite:

- Y. Bea and D. Mateos, *Heating up Exotic RG Flows with Holography*, JHEP **08** (2018) 034,
  [arXiv:1805.01806](https://arxiv.org/abs/1805.01806). This is the model.
- A. Jansen, *Overdamped modes in Schwarzschild-de Sitter and a Mathematica package for the
  numerical computation of quasinormal modes*, Eur. Phys. J. Plus **132**, 546 (2017),
  [arXiv:1709.09178](https://arxiv.org/abs/1709.09178). This is
  [QNMspectral](https://github.com/APJansen/QNMspectral), whose convergence test `qnms` uses.
- A. Caddeo, O. Henriksson, C. Hoyos and M. Sánchez-Garitaonandia, *Spinodal slowing down and
  scaling in a holographic model*, [arXiv:2406.15297](https://arxiv.org/abs/2406.15297). This is
  Mikel Sánchez-Garitaonandia and Carlos Hoyos's work on the dynamics near the spinodal point.

## License

[MIT](LICENSE).

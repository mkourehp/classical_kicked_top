# classical_kicked_top

Trajectories and phase space of the classical kicked top. The computation is in Fortran; Python only plots the data it writes.

## Model

Hamiltonian $H(t) = \frac{k}{2j}J_z^2 + p\,J_y\sum_n\delta(t-n)$: a continuous torsion about $z$ and, at integer times, a kick that rotates about $y$. A unit vector $(X, Y, Z)$ on the sphere, observed just before each kick. One period:

1. kick, rotation about $y$ by $p$:  $X' = X\cos p + Z\sin p,\ Y' = Y,\ Z' = -X\sin p + Z\cos p$
2. free torsion about $z$ by $kZ'$: $X'' = X'\cos kZ' - Y'\sin kZ',\ Y'' = X'\sin kZ' + Y'\cos kZ',\ Z'' = Z'$

Phase space is shown in the canonical (area-preserving) coordinates $(\phi, Z = \cos\theta)$.

## Layout

```
config.nml                      classical parameters (namelist)
config_quantum.nml              quantum parameters (namelist)
Makefile                        gfortran + OpenMP, build artifacts in build/
src/common/kinds_mod.f90        precision and constants
src/common/sphere_mod.f90       unit vector from (theta, phi), equal-area grid in (phi, Z)
src/common/io_utils_mod.f90     namelist-read check, output directory
src/classical/config_mod.f90    namelist reading and validation
src/classical/map_mod.f90       one period (kick + torsion) and its tangent map
src/classical/energy_mod.f90    energy terms pY, (k/2)Z^2 and their sum
src/classical/init_mod.f90      initial conditions: single / grid / random
src/classical/evolve_mod.f90    iteration of the map over all initial conditions
src/classical/lyapunov_mod.f90  largest Lyapunov exponent (Benettin algorithm)
src/classical/io_mod.f90        ASCII output
src/classical/main.f90          driver -> build/kicked_top
src/quantum/quantum_config_mod.f90  namelist reading and validation
src/quantum/spin_mod.f90        J_y, <J>, spin coherent states
src/quantum/floquet_mod.f90     U = exp(-i k/(2j) Jz^2) exp(-i p Jy) (LAPACK zheev)
src/quantum/husimi_mod.f90      Husimi distribution Q(theta, phi)
src/quantum/spectrum_mod.f90    quasienergies per parity sector (zgeev), spacing statistics
src/quantum/quantum_io_mod.f90  ASCII output
src/quantum/quantum_main.f90    driver -> build/kicked_top_quantum
scripts/post_processor.py          classical plots
scripts/quantum_post_processor.py  quantum plots
```

New features go in their own module with their own namelist group, called from `main.f90`.

The physics and the code are explained in the lecture notes `docs/lecture_note/kicked_top_notes.pdf` (`make doc lecture_note` rebuilds them, `make doc lecture_note compute` first recomputes their figures). Code comments cite the equations they implement by number and LaTeX label, e.g. `eq. (20) [eq:map]`.

## Usage

```bash
make                                    # builds all programs: kicked_top, kicked_top_quantum, kicked_top_spin_chain, kicked_top_combined
make classical                          # builds and runs ./build/kicked_top config_classical.nml
make spin_chain                         # qubit chain alone: ./build/kicked_top_spin_chain config_spin_chain.nml
make combined                           # classical sections + qubit chain: ./build/kicked_top_combined config_combined.nml
python scripts/post_processor.py        # output/phase_space.png
python scripts/post_processor.py --sphere --show   # also 3D trajectory, open windows
python scripts/post_processor.py --lyapunov        # also output/lyapunov.png (needs &lyapunov enabled)
python scripts/spin_chain_post_processor.py        # output_spin_chain/*.png (Husimi, fidelity, total and single-qubit spins)
python scripts/combined_post_processor.py          # output_combined/*.png (the same plus the classical sections)
```

Python needs `numpy` and `matplotlib`. The quantum code needs LAPACK and BLAS (Ubuntu: `sudo apt install liblapack-dev libblas-dev`).

The notebook `notebooks/kam_theory.ipynb` explains KAM theory, the Poincaré–Birkhoff theorem, the Diophantine condition and golden-mean tori, using the standard map as the simplest example. It is standalone (no Fortran build or output needed) and needs `numpy`, `matplotlib` and `jupyter`: `jupyter notebook notebooks/kam_theory.ipynb`.

## Parameters (`config.nml`)

| group | name | meaning |
|---|---|---|
| `&system` | `p`, `k` | kick rotation angle and torsion strength |
| `&simulation` | `n_kicks` | kicks recorded after the transient |
| | `n_transient` | kicks discarded before recording |
| | `save_every` | store every `save_every`-th kick |
| `&initial` | `mode` | `'single'`, `'grid'` or `'random'` |
| | `theta0`, `phi0` | single initial point (radians) |
| | `n_theta`, `n_phi` | grid size, uniform in $Z=\cos\theta$ and $\phi$ |
| | `n_random`, `seed` | random points, uniform on the sphere |
| `&output` | `out_dir` | output directory |
| `&lyapunov` | `enabled` | compute the largest Lyapunov exponent of every initial condition |
| | `n_renorm` | rescale the tangent vector every `n_renorm` kicks |

## Output

ASCII, one block per initial condition separated by a blank line; first line is a `# p = ..., k = ...` header.

- `trajectory.dat`: `ic n X Y Z`
- `phase_space.dat`: `ic n phi Z`, with $\phi \in (-\pi, \pi]$
- `energy.dat` (modes `single` and `grid`): `pY  k/2*Z^2  total`, row-aligned with the two files above
- `lyapunov.dat` (if `&lyapunov` is enabled): `ic phi0 Z0 Lambda`, one line per initial condition

## Quantum kicked top

$U = e^{-i\frac{k}{2j}J_z^2}e^{-ipJ_y}$ (kick, then torsion) on the $2j+1$ states $|j,m\rangle$. A spin coherent state centred at $(\theta_0, \phi_0)$ is evolved kick by kick; optionally the quasienergies of $U$ are computed separately in the two parity sectors of $e^{-i\pi J_y}$.

```bash
make quantum                                       # builds and runs ./build/kicked_top_quantum config_quantum.nml
python scripts/quantum_post_processor.py           # output_quantum/*.png
python scripts/quantum_post_processor.py --classical output/trajectory.dat --sphere
```

`--classical` overlays a classical run started from the same point (`mode = 'single'`, same `theta0`, `phi0`, `p`, `k`).

| group | name | meaning |
|---|---|---|
| `&system` | `p`, `k`, `j` | kick rotation angle, torsion strength, spin (integer or half-integer) |
| `&simulation` | `n_kicks` | number of kicks |
| `&initial` | `theta0`, `phi0` | centre of the initial coherent state |
| `&husimi` | `enabled`, `n_theta`, `n_phi`, `kicks` | Husimi snapshots on the equal-area grid at the listed kicks |
| | `n_husimi` | if > 0, also `n_husimi` snapshots evenly spaced from kick 0 to `n_kicks` (merged with `kicks`, sorted, duplicates removed) |
| `&spectrum` | `enabled`, `n_bins`, `s_max` | quasienergies and spacing histogram on $[0, s_{max}]$ |
| `&output` | `out_dir` | output directory (default `output_quantum`) |

Output (ASCII, `# p = ..., k = ..., j = ...` header):

- `expectation.dat`: `n Jx/j Jy/j Jz/j norm spread`, spread $=(\langle J^2\rangle-|\langle J\rangle|^2)/j^2$
- `husimi.dat`: `n phi Z X Y Q`, one block per snapshot
- `quasienergies.dat`: `sector eps`
- `spacings.dat`: `s P(s) Poisson COE CUE`
- `summary.dat`: correctness checks (unitarity, parity commutator, norm, Husimi normalisation, mean gap ratio)

For level statistics use a generic $p$ (e.g. 1.7): at $p=\pi/2$ one parity sector still mixes two independent spectra.

# classical_kicked_top

Trajectories and phase space of the classical kicked top. The computation is in Fortran; Python only plots the data it writes.

## Model

A unit vector $(X, Y, Z)$ on the sphere. One kick (Haake convention):

1. rotation about $y$ by $p$:  $X' = X\cos p + Z\sin p,\ Y' = Y,\ Z' = -X\sin p + Z\cos p$
2. torsion about $z$ by $kZ'$: $X'' = X'\cos kZ' - Y'\sin kZ',\ Y'' = X'\sin kZ' + Y'\cos kZ',\ Z'' = Z'$

Phase space is shown in the canonical (area-preserving) coordinates $(\phi, Z = \cos\theta)$.

## Layout

```
config.nml               parameters (namelist)
Makefile                 gfortran + OpenMP, build artifacts in build/
src/kinds_mod.f90        precision and constants
src/config_mod.f90       namelist reading and validation
src/map_mod.f90          one kick (rotation + torsion) and its tangent map
src/energy_mod.f90       energy terms pY, (k/2)Z^2 and their sum
src/init_mod.f90         initial conditions: single / grid / random
src/evolve_mod.f90       iteration of the map over all initial conditions
src/lyapunov_mod.f90     largest Lyapunov exponent (Benettin algorithm)
src/io_mod.f90           ASCII output
src/main.f90             driver
scripts/post_processor.py  plots
notebooks/kam_theory.ipynb KAM theory explainer (standalone)
```

New features go in their own module with their own namelist group, called from `main.f90`.

The physics and the code are explained in the lecture notes `docs/kicked_top_notes.pdf` (`make docs` rebuilds them, `make figures` recomputes their figures). Code comments cite the equations they implement by number and LaTeX label, e.g. `eq. (20) [eq:map]`.

## Usage

```bash
make                                    # builds build/kicked_top
./build/kicked_top config.nml           # or: make run
python scripts/post_processor.py        # output/phase_space.png
python scripts/post_processor.py --sphere --show   # also 3D trajectory, open windows
python scripts/post_processor.py --lyapunov        # also output/lyapunov.png (needs &lyapunov enabled)
```

Python needs `numpy` and `matplotlib`.

The notebook `notebooks/kam_theory.ipynb` explains KAM theory, the Poincaré–Birkhoff theorem, the Diophantine condition and golden-mean tori, using the standard map as the simplest example. It is standalone (no Fortran build or output needed) and needs `numpy`, `matplotlib` and `jupyter`: `jupyter notebook notebooks/kam_theory.ipynb`.

## Parameters (`config.nml`)

| group | name | meaning |
|---|---|---|
| `&system` | `p`, `k` | rotation angle and torsion strength |
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
- `energy.dat`: `pY  k/2*Z^2  total`, row-aligned with the two files above
- `lyapunov.dat` (if `&lyapunov` is enabled): `ic phi0 Z0 Lambda`, one line per initial condition

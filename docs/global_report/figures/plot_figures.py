"""Plots the figures of the lecture notes from the data written by `make figures`."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = Path(__file__).resolve().parent
DATA_DIR = FIG_DIR / "data"
PORTRAIT_K = ["0.5", "1.5", "2.5", "3.0", "4.5", "6.0"]

plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "savefig.bbox": "tight"})


def load(run, name):
    """Numeric data of output file name of run (a subdirectory of data/)."""
    return np.loadtxt(DATA_DIR / run / name)


def phase_axes(ax):
    """Axes limits, ticks and labels of a (phi, Z) phase portrait."""
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-1, 1)
    ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi], [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.set_xlabel(r"$\phi$")
    ax.set_ylabel(r"$Z=\cos\theta$")


def portraits():
    """Fig. 5 [fig:portraits]: phase portraits for p = pi/2 and six values of k."""
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.5), sharex=True, sharey=True)
    for ax, k in zip(axes.flat, PORTRAIT_K):
        d = load(f"portrait_k{k}", "phase_space.dat")
        ax.scatter(d[:, 2], d[:, 3], c=d[:, 0], s=0.05, cmap="turbo", rasterized=True)
        phase_axes(ax)
        ax.set_title(rf"$k = {k}$")
        ax.label_outer()
    fig.savefig(FIG_DIR / "portraits.pdf", dpi=150)


def near_flow():
    """Fig. 3 [fig:near-flow]: small p, k; colour = E of eq. (23)."""
    d = load("near_flow", "phase_space.dat")
    e = load("near_flow", "energy.dat")
    fig, ax = plt.subplots(figsize=(7, 4))
    sc = ax.scatter(d[:, 2], d[:, 3], c=e[:, 2], s=0.1, cmap="viridis", rasterized=True)
    phase_axes(ax)
    fig.colorbar(sc, ax=ax, label=r"$E = pY + \frac{k}{2}Z^2$")
    fig.savefig(FIG_DIR / "near_flow.pdf", dpi=150)


def sphere():
    """Fig. 7 [fig:sphere]: the k = 3 portrait on the unit sphere."""
    d = load("portrait_k3.0", "trajectory.dat")
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(projection="3d")
    ax.scatter(d[:, 2], d[:, 3], d[:, 4], c=d[:, 0], s=0.05, cmap="turbo", rasterized=True)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.view_init(elev=20, azim=35)
    fig.savefig(FIG_DIR / "sphere.pdf", dpi=150)


def energy():
    """Fig. 4 [fig:energy]: energy terms of eq. (23) along a regular and a chaotic orbit."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.5), sharey=True)
    for ax, run, title in zip(axes, ["single_regular", "single_chaotic"], ["regular orbit", "chaotic orbit"]):
        e = load(run, "energy.dat")
        n = np.arange(len(e))
        ax.plot(n, e[:, 0], lw=0.8, label=r"$pY$")
        ax.plot(n, e[:, 1], lw=0.8, label=r"$\frac{k}{2}Z^2$")
        ax.plot(n, e[:, 2], "k", lw=1.0, label="total")
        ax.set_xlim(0, 100)
        ax.set_xlabel("kick $n$")
        ax.set_title(title)
    axes[0].set_ylabel("energy")
    axes[0].legend(loc="lower left", ncol=3, fontsize=8)
    fig.savefig(FIG_DIR / "energy.pdf")


def divergence():
    """Fig. 8 [fig:divergence]: distance between two nearby orbits."""
    fig, ax = plt.subplots(figsize=(6, 4))
    for run, label in [("single_regular", "regular"), ("single_chaotic", "chaotic")]:
        a = load(run, "trajectory.dat")
        b = load(f"{run}_shifted", "trajectory.dat")
        ax.semilogy(a[:, 1], np.linalg.norm(a[:, 2:] - b[:, 2:], axis=1), label=label)
    ax.set_xlim(0, 60)
    ax.set_xlabel("kick $n$")
    ax.set_ylabel(r"distance $|\mathbf{s}_n - \mathbf{s}'_n|$")
    ax.legend()
    fig.savefig(FIG_DIR / "divergence.pdf")


def single_orbits():
    """Fig. 6 [fig:single-orbits]: one regular and one chaotic orbit for k = 3."""
    d = load("portrait_k3.0", "phase_space.dat")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.scatter(d[:, 2], d[:, 3], c="0.8", s=0.05, rasterized=True)
    for run, color in [("single_regular", "C0"), ("single_chaotic", "C3")]:
        o = load(run, "phase_space.dat")
        ax.scatter(o[:, 2], o[:, 3], c=color, s=2, rasterized=True)
        ax.plot(o[0, 2], o[0, 3], "k*", ms=10)
    phase_axes(ax)
    fig.savefig(FIG_DIR / "single_orbits.pdf", dpi=150)


def lyapunov_map():
    """Fig. 9 [fig:lyapunov-map]: Lyapunov exponent of each initial condition, k = 3."""
    d = load("lyapunov_k3.0", "lyapunov.dat")
    fig, ax = plt.subplots(figsize=(7, 4))
    sc = ax.scatter(d[:, 1], d[:, 2], c=d[:, 3], s=3, marker="s", cmap="magma", rasterized=True)
    phase_axes(ax)
    ax.set_xlabel(r"$\phi_0$")
    ax.set_ylabel(r"$Z_0=\cos\theta_0$")
    fig.colorbar(sc, ax=ax, label=r"$\Lambda$")
    fig.savefig(FIG_DIR / "lyapunov_map.pdf", dpi=150)


def quantum_husimi():
    """Fig. [fig:quantum-husimi]: Husimi snapshots of a coherent state in an island and in the chaotic sea."""
    fig, axes = plt.subplots(2, 4, figsize=(13, 5.5), sharex=True, sharey=True)
    for row, run in zip(axes, ["quantum_island", "quantum_chaotic"]):
        d = load(run, "husimi.dat")
        kicks = np.unique(d[:, 0]).astype(int)
        n_phi = len(np.unique(d[d[:, 0] == kicks[0], 1]))
        for ax, n in zip(row, kicks):
            b = d[d[:, 0] == n]
            ax.pcolormesh(b[:, 1].reshape(-1, n_phi), b[:, 2].reshape(-1, n_phi), b[:, 5].reshape(-1, n_phi),
                          shading="nearest", cmap="viridis", rasterized=True)
            phase_axes(ax)
            ax.set_title(f"$n = {n}$")
            ax.label_outer()
    fig.savefig(FIG_DIR / "quantum_husimi.pdf", dpi=150)


def quantum_expectation():
    """Fig. [fig:quantum-expectation]: <J>/j of the coherent state against the classical orbit."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.5), sharey=True)
    for ax, run, title in zip(axes, ["quantum_island", "quantum_chaotic"], ["island", "chaotic sea"]):
        q = load(run, "expectation.dat")
        c = load(f"{run}_classical", "trajectory.dat")
        for col, name, color in [(1, "x", "C0"), (2, "y", "C1"), (3, "z", "C2")]:
            ax.plot(q[:, 0], q[:, col], color=color, label=rf"$\langle J_{name}\rangle/j$")
            ax.plot(c[:, 1], c[:, col + 1], "--", color=color, lw=0.8)
        ax.set_xlabel("kick $n$")
        ax.set_title(title)
    axes[0].legend(fontsize=8, ncol=3, loc="lower left")
    fig.savefig(FIG_DIR / "quantum_expectation.pdf")


def quantum_spacings():
    """Fig. [fig:quantum-spacings]: quasienergy spacing histogram against Poisson, COE and CUE."""
    d = load("quantum_spectrum", "spacings.dat")
    fig, ax = plt.subplots(figsize=(6, 3.8))
    ax.bar(d[:, 0], d[:, 1], width=d[1, 0] - d[0, 0], color="0.75", label="quasienergies")
    for col, label in [(2, "Poisson"), (3, "COE"), (4, "CUE")]:
        ax.plot(d[:, 0], d[:, col], label=label)
    ax.set_xlabel("$s$")
    ax.set_ylabel("$P(s)$")
    ax.legend()
    fig.savefig(FIG_DIR / "quantum_spacings.pdf")


if __name__ == "__main__":
    portraits()
    near_flow()
    sphere()
    energy()
    divergence()
    single_orbits()
    lyapunov_map()
    quantum_husimi()
    quantum_expectation()
    quantum_spacings()

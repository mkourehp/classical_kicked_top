"""Plots the figures of docs/three_regimes.tex from the data written by `make figures`
(kicked_top_regimes with docs/figures/configs_regimes/three_regimes.nml)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = Path(__file__).resolve().parent
DATA_DIR = FIG_DIR / "data" / "three_regimes"
K_CHIRIKOV = 0.9716

plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "savefig.bbox": "tight"})


def load(i, name):
    """Numeric data of output file name of run i (1-based)."""
    return np.loadtxt(DATA_DIR / f"h_{i}" / name, ndmin=2)


def runs():
    """(i, p, k, n_q) of every run."""
    return [(int(r[0]), r[1], r[2], int(r[4])) for r in np.loadtxt(DATA_DIR / "regimes.dat", ndmin=2)]


def phase_axes(ax):
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-1, 1)
    ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi], [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.set_xlabel(r"$\phi$")


def sections():
    """Fig. [fig:sections]: Poincare sections; dashed: |Z| = Z_c where k p sqrt(1-Z^2) = K_c."""
    rs = runs()
    fig, axes = plt.subplots(1, len(rs), figsize=(4.2 * len(rs), 3.8), sharey=True)
    for ax, (i, p, k, _) in zip(axes, rs):
        s = load(i, "classical/section/phase_space.dat")
        o = load(i, "classical/orbit/phase_space.dat")
        ax.plot(s[:, 2], s[:, 3], ",", color="k", alpha=0.6, rasterized=True)
        ax.plot(o[:, 2], o[:, 3], ".", color="C3", ms=1.5, rasterized=True)
        if K_CHIRIKOV < k * p:
            zc = np.sqrt(1 - (K_CHIRIKOV / (k * p)) ** 2)
            for z in (zc, -zc):
                ax.axhline(z, color="C0", ls="--", lw=1)
        phase_axes(ax)
        ax.set_title(rf"$p = {p:g}$, $kp = {k * p:g}$")
    axes[0].set_ylabel(r"$Z=\cos\theta$")
    fig.savefig(FIG_DIR / "three_regimes_sections.pdf", dpi=200)


def husimi():
    """Fig. [fig:husimi]: Husimi distribution of the chain at t = 0, t_r, 2 t_r."""
    rs = runs()
    d0 = load(1, "chain/husimi.dat")
    steps = np.unique(d0[:, 0]).astype(int)
    n_phi = len(np.unique(d0[d0[:, 0] == steps[0], 1]))
    fig, axes = plt.subplots(len(rs), len(steps), figsize=(3.6 * len(steps), 3.0 * len(rs)),
                             sharex=True, sharey=True)
    for row, (i, p, _, _) in zip(axes, rs):
        d = load(i, "chain/husimi.dat")
        for ax, t in zip(row, steps):
            b = d[d[:, 0] == t]
            shape = (-1, n_phi)
            ax.pcolormesh(b[:, 1].reshape(shape), b[:, 2].reshape(shape), b[:, 5].reshape(shape),
                          shading="nearest", cmap="viridis", rasterized=True)
            phase_axes(ax)
            ax.set_title(rf"$p = {p:g}$, $t = {t}$")
            ax.label_outer()
        row[0].set_ylabel(r"$Z=\cos\theta$")
    fig.savefig(FIG_DIR / "three_regimes_husimi.pdf", dpi=200)


def echo():
    """Fig. [fig:echo]: fidelity F(t) and single-qubit Bloch-vector length r(t) = |<J>|/j."""
    fig, (ax_f, ax_r) = plt.subplots(2, 1, figsize=(7, 5.5), sharex=True)
    for i, p, _, n_q in runs():
        f = load(i, "chain/fidelity.dat")
        s = load(i, "chain/total_spin.dat")
        ax_f.plot(f[:, 0], f[:, 1], "o-", ms=2.5, lw=1, label=rf"$p = {p:g}$")
        ax_r.plot(s[:, 0], np.linalg.norm(s[:, 1:4], axis=1) / (n_q / 2), "o-", ms=2.5, lw=1)
    t_r = (len(f) - 1) // 2
    for ax in (ax_f, ax_r):
        ax.axvline(t_r, color="0.5", ls="--", lw=0.8)
        ax.set_ylim(-0.02, 1.02)
    ax_f.axhline(1 / (n_q + 1), color="0.5", ls=":", lw=0.8)
    ax_f.set_ylabel(r"$F(t)$")
    ax_r.set_ylabel(r"$r(t)$")
    ax_r.set_xlabel(r"$t$")
    ax_f.legend(loc="center right", fontsize=8)
    fig.savefig(FIG_DIR / "three_regimes_echo.pdf")


if __name__ == "__main__":
    sections()
    husimi()
    echo()

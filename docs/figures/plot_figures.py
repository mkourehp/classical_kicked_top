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
    return np.loadtxt(DATA_DIR / run / name)


def phase_axes(ax):
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-1, 1)
    ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi], [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.set_xlabel(r"$\phi$")
    ax.set_ylabel(r"$Z=\cos\theta$")


def portraits():
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.5), sharex=True, sharey=True)
    for ax, k in zip(axes.flat, PORTRAIT_K):
        d = load(f"portrait_k{k}", "phase_space.dat")
        ax.scatter(d[:, 2], d[:, 3], c=d[:, 0], s=0.05, cmap="turbo", rasterized=True)
        phase_axes(ax)
        ax.set_title(rf"$k = {k}$")
        ax.label_outer()
    fig.savefig(FIG_DIR / "portraits.pdf", dpi=150)


def near_flow():
    d = load("near_flow", "phase_space.dat")
    e = load("near_flow", "energy.dat")
    fig, ax = plt.subplots(figsize=(7, 4))
    sc = ax.scatter(d[:, 2], d[:, 3], c=e[:, 2], s=0.1, cmap="viridis", rasterized=True)
    phase_axes(ax)
    fig.colorbar(sc, ax=ax, label=r"$E = pY + \frac{k}{2}Z^2$")
    fig.savefig(FIG_DIR / "near_flow.pdf", dpi=150)


def sphere():
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
    d = load("portrait_k3.0", "phase_space.dat")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.scatter(d[:, 2], d[:, 3], c="0.8", s=0.05, rasterized=True)
    for run, color in [("single_regular", "C0"), ("single_chaotic", "C3")]:
        o = load(run, "phase_space.dat")
        ax.scatter(o[:, 2], o[:, 3], c=color, s=2, rasterized=True)
        ax.plot(o[0, 2], o[0, 3], "k*", ms=10)
    phase_axes(ax)
    fig.savefig(FIG_DIR / "single_orbits.pdf", dpi=150)


if __name__ == "__main__":
    portraits()
    near_flow()
    sphere()
    energy()
    divergence()
    single_orbits()

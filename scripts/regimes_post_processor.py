"""Plots the plot-ready data written by kicked_top_regimes.

All computation is done in Fortran (build/kicked_top_regimes); this script only reads
the ASCII files in the output directory and plots them, one column per value of h:
the classical Poincare sections, the Husimi snapshots of the qubit chain at t = 0,
t_r and 2 t_r, the fidelity, the total spin and the single-qubit spins.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


class RegimesPostProcessor:
    """Reads the output in out_dir (layout described in src/regimes/regimes_io_mod.f90)
    and writes PNG figures into it."""

    def __init__(self, out_dir="output_regimes"):
        self.out_dir = Path(out_dir)
        self.regimes = np.loadtxt(self.out_dir / "regimes.dat", ndmin=2)

    def _load(self, i, name):
        """Numeric data of output file name of run i (1-based)."""
        return np.loadtxt(self.out_dir / f"h_{i}" / name, ndmin=2)

    def _runs(self):
        """(i, h, k) of every run."""
        return [(int(r[0]), r[1], r[2]) for r in self.regimes]

    def _t_r(self):
        """Reversal step t_r, read from the length of the chain time series (2 t_r + 1 steps)."""
        return (len(self._load(1, "chain/fidelity.dat")) - 1) // 2

    def plot_sections(self):
        """Classical Poincare sections in (phi, Z), with the orbit of the initial state in red."""
        runs = self._runs()
        fig, axes = plt.subplots(1, len(runs), figsize=(4.5 * len(runs), 4), squeeze=False)
        for ax, (i, h, k) in zip(axes[0], runs):
            s = self._load(i, "classical/section/phase_space.dat")
            o = self._load(i, "classical/orbit/phase_space.dat")
            ax.plot(s[:, 2], s[:, 3], ",", color="k", alpha=0.5, rasterized=True)
            ax.plot(o[:, 2], o[:, 3], ".", color="C3", ms=2, rasterized=True)
            ax.set_xlim(-np.pi, np.pi)
            ax.set_ylim(-1, 1)
            ax.set_xlabel(r"$\phi$")
            ax.set_title(f"$h = p = {h:g}$, $k = {k:g}$")
        axes[0, 0].set_ylabel(r"$Z = \cos\theta$")
        fig.suptitle("Classical kicked top, Poincare sections")
        fig.tight_layout()
        fig.savefig(self.out_dir / "sections.png", dpi=200)
        return fig

    def plot_husimi(self):
        """Husimi distribution of the chain state in (phi, Z): one row per h, one column per snapshot."""
        runs = self._runs()
        d0 = self._load(1, "chain/husimi.dat")
        steps = np.unique(d0[:, 0]).astype(int)
        n_phi = len(np.unique(d0[d0[:, 0] == steps[0], 1]))
        fig, axes = plt.subplots(len(runs), len(steps), figsize=(4 * len(steps), 3.4 * len(runs)),
                                 squeeze=False)
        for row, (i, h, k) in zip(axes, runs):
            d = self._load(i, "chain/husimi.dat")
            for ax, t in zip(row, steps):
                b = d[d[:, 0] == t]
                shape = (-1, n_phi)
                ax.pcolormesh(b[:, 1].reshape(shape), b[:, 2].reshape(shape), b[:, 5].reshape(shape),
                              shading="nearest", cmap="viridis", rasterized=True)
                ax.set_xlim(-np.pi, np.pi)
                ax.set_ylim(-1, 1)
                ax.set_title(f"$h = {h:g}$, $t = {t}$")
            row[0].set_ylabel(r"$Z = \cos\theta$")
        for ax in axes[-1]:
            ax.set_xlabel(r"$\phi$")
        fig.suptitle("Qubit chain, Husimi distribution")
        fig.tight_layout()
        fig.savefig(self.out_dir / "husimi.png", dpi=200)
        return fig

    def plot_fidelity(self):
        """Fidelity F(t) = |<psi0|psi(t)>|^2 through the time reversal at t_r."""
        fig, ax = plt.subplots(figsize=(7, 4))
        for i, h, k in self._runs():
            d = self._load(i, "chain/fidelity.dat")
            ax.plot(d[:, 0], d[:, 1], "o-", ms=3, label=f"$h = {h:g}$")
        ax.axvline(self._t_r(), color="0.5", ls="--", lw=0.8)
        ax.set_xlabel("step $t$")
        ax.set_ylabel(r"$F(t)$")
        ax.set_ylim(-0.02, 1.02)
        ax.legend()
        ax.set_title(f"Qubit chain, time reversal at $t_r = {self._t_r()}$")
        fig.tight_layout()
        fig.savefig(self.out_dir / "fidelity.png", dpi=200)
        return fig

    def plot_total_spin(self):
        """<S_a> (top) and <S_a^2>, <S^2> (bottom) of the total spin, one column per h."""
        runs = self._runs()
        fig, axes = plt.subplots(2, len(runs), figsize=(4.5 * len(runs), 6), sharex=True, squeeze=False)
        for col, (i, h, k) in enumerate(runs):
            d = self._load(i, "chain/total_spin.dat")
            for c, a in [(1, "x"), (2, "y"), (3, "z")]:
                axes[0, col].plot(d[:, 0], d[:, c], label=rf"$\langle S_{a}\rangle$")
                axes[1, col].plot(d[:, 0], d[:, c + 3], label=rf"$\langle S_{a}^2\rangle$")
            axes[1, col].plot(d[:, 0], d[:, 7], "k", label=r"$\langle S^2\rangle$")
            for ax in axes[:, col]:
                ax.axvline(self._t_r(), color="0.5", ls="--", lw=0.8)
            axes[0, col].set_title(f"$h = {h:g}$")
            axes[1, col].set_xlabel("step $t$")
        axes[0, 0].legend(fontsize=8)
        axes[1, 0].legend(fontsize=8)
        fig.suptitle("Qubit chain, total spin")
        fig.tight_layout()
        fig.savefig(self.out_dir / "total_spin.png", dpi=200)
        return fig

    def plot_site_spins(self):
        """<sigma_a^i>(t) of every qubit i as a colour map: rows a = x, y, z, one column per h."""
        runs = self._runs()
        fig, axes = plt.subplots(3, len(runs), figsize=(4.5 * len(runs), 7), sharex=True, squeeze=False)
        for col, (i, h, k) in enumerate(runs):
            d = self._load(i, "chain/site_spins.dat")
            n_q = (d.shape[1] - 1) // 3
            for row, a in enumerate("xyz"):
                ax = axes[row, col]
                im = ax.pcolormesh(d[:, 0], np.arange(1, n_q + 1), d[:, 1 + row * n_q:1 + (row + 1) * n_q].T,
                                   shading="nearest", cmap="RdBu_r", vmin=-1, vmax=1)
                ax.axvline(self._t_r(), color="k", ls="--", lw=0.8)
                ax.set_ylabel(rf"qubit $i$, $\langle\sigma_{a}^i\rangle$")
            axes[0, col].set_title(f"$h = {h:g}$")
            axes[-1, col].set_xlabel("step $t$")
        fig.colorbar(im, ax=axes, shrink=0.6)
        fig.suptitle("Qubit chain, single-qubit spins")
        fig.savefig(self.out_dir / "site_spins.png", dpi=200)
        return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", default="output_regimes", help="directory with the kicked_top_regimes output")
    parser.add_argument("--show", action="store_true", help="open the figures in a window")
    args = parser.parse_args()

    rp = RegimesPostProcessor(args.out_dir)
    rp.plot_sections()
    rp.plot_husimi()
    rp.plot_fidelity()
    rp.plot_total_spin()
    rp.plot_site_spins()
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

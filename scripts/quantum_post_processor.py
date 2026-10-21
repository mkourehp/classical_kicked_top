"""Plots the plot-ready data written by the Fortran quantum kicked top code.

All computation is done in Fortran (build/kicked_top_quantum); this script only reads
the ASCII files in the output directory and plots them. Section numbers refer to the
lecture notes, docs/kicked_top_notes.pdf.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


class QuantumPostProcessor:
    """Reads the quantum output in out_dir and writes PNG figures next to it."""

    def __init__(self, out_dir="output_quantum"):
        self.out_dir = Path(out_dir)

    def _load(self, name):
        """Returns the header title (first line) and the numeric data of an output file."""
        path = self.out_dir / name
        with path.open() as f:
            title = f.readline().lstrip("#").strip()
        return title, np.loadtxt(path, ndmin=2)

    def plot_expectation(self, classical=None):
        """<J>/j and the spread of the state per kick; optionally the classical orbit
        (trajectory.dat of the classical code, started at the same point) for comparison."""
        title, d = self._load("expectation.dat")
        fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, height_ratios=[2, 1])
        for col, label, color in [(1, "X", "C0"), (2, "Y", "C1"), (3, "Z", "C2")]:
            ax.plot(d[:, 0], d[:, col], color=color, label=rf"$\langle J_{label.lower()}\rangle/j$")
        if classical is not None:
            c = np.loadtxt(classical, ndmin=2)
            c = c[c[:, 0] == 1]
            for col, color in [(2, "C0"), (3, "C1"), (4, "C2")]:
                ax.plot(c[:, 1], c[:, col], "--", color=color, lw=0.8)
            ax.plot([], [], "k--", lw=0.8, label="classical")
        ax.set_ylabel("expectation value")
        ax.legend(ncol=4, fontsize=8)
        ax.set_title(f"Quantum kicked top, {title}")
        ax2.plot(d[:, 0], d[:, 5], "k")
        ax2.set_xlabel("kick $n$")
        ax2.set_ylabel(r"$(\langle J^2\rangle - |\langle J\rangle|^2)/j^2$")
        fig.tight_layout()
        fig.savefig(self.out_dir / "expectation.png", dpi=200)
        return fig

    def plot_husimi(self):
        """Husimi distribution Q in the canonical coordinates (phi, Z), one panel per snapshot."""
        title, d = self._load("husimi.dat")
        kicks = np.unique(d[:, 0]).astype(int)
        n_phi = len(np.unique(d[d[:, 0] == kicks[0], 1]))
        fig, axes = plt.subplots(1, len(kicks), figsize=(4 * len(kicks), 3.6), squeeze=False)
        for ax, n in zip(axes[0], kicks):
            b = d[d[:, 0] == n]
            shape = (-1, n_phi)
            ax.pcolormesh(b[:, 1].reshape(shape), b[:, 2].reshape(shape), b[:, 5].reshape(shape),
                          shading="nearest", cmap="viridis", rasterized=True)
            ax.set_xlim(-np.pi, np.pi)
            ax.set_ylim(-1, 1)
            ax.set_xlabel(r"$\phi$")
            ax.set_title(f"$n = {n}$")
        axes[0, 0].set_ylabel(r"$Z = \cos\theta$")
        fig.suptitle(f"Husimi distribution, {title}")
        fig.tight_layout()
        fig.savefig(self.out_dir / "husimi.png", dpi=200)
        return fig

    def plot_husimi_sphere(self):
        """Last Husimi snapshot drawn on the unit sphere (X, Y, Z)."""
        title, d = self._load("husimi.dat")
        b = d[d[:, 0] == d[-1, 0]]
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(projection="3d")
        ax.scatter(b[:, 3], b[:, 4], b[:, 2], c=b[:, 5], s=2, cmap="viridis")
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title(f"Husimi distribution at n = {int(b[0, 0])}, {title}")
        fig.tight_layout()
        fig.savefig(self.out_dir / "husimi_sphere.png", dpi=200)
        return fig

    def plot_spacings(self):
        """Histogram of unfolded quasienergy spacings against the Poisson, COE and CUE curves."""
        title, d = self._load("spacings.dat")
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(d[:, 0], d[:, 1], width=d[1, 0] - d[0, 0], color="0.75", label="quasienergies")
        for col, label in [(2, "Poisson"), (3, "COE"), (4, "CUE")]:
            ax.plot(d[:, 0], d[:, col], label=label)
        ax.set_xlabel("spacing $s$")
        ax.set_ylabel("$P(s)$")
        ax.set_title(title)
        ax.legend()
        fig.tight_layout()
        fig.savefig(self.out_dir / "spacings.png", dpi=200)
        return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", default="output_quantum", help="directory with the quantum Fortran output")
    parser.add_argument("--classical", help="trajectory.dat of a classical run from the same initial point, to overlay")
    parser.add_argument("--sphere", action="store_true", help="also draw the last Husimi snapshot on the sphere")
    parser.add_argument("--show", action="store_true", help="open the figures in a window")
    args = parser.parse_args()

    qp = QuantumPostProcessor(args.out_dir)
    qp.plot_expectation(args.classical)
    if (qp.out_dir / "husimi.dat").exists():
        qp.plot_husimi()
        if args.sphere:
            qp.plot_husimi_sphere()
    if (qp.out_dir / "spacings.dat").exists():
        qp.plot_spacings()
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

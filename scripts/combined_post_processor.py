"""Plots the plot-ready data written by kicked_top_combined.

All computation is done in Fortran (build/kicked_top_combined); this script only reads
the ASCII files in the output directory and plots them, one column per value of p:
the classical Poincare sections, the Husimi snapshots of the qubit chain at t = 0,
t_r/2, t_r, 3 t_r/2 and 2 t_r, the fidelity, the total spin and the single-qubit spins.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from spin_chain_post_processor import SpinChainPostProcessor


class CombinedPostProcessor(SpinChainPostProcessor):
    """Reads the output in out_dir (layout described in src/combined/combined_io_mod.f90)
    and writes PNG figures into it: the chain plots of SpinChainPostProcessor plus the
    classical Poincare sections."""

    def __init__(self, out_dir="output_combined"):
        super().__init__(out_dir)

    def plot_sections(self):
        """Classical Poincare sections in (phi, Z), with the orbit of the initial state in red."""
        runs = self._runs()
        fig, axes = plt.subplots(1, len(runs), figsize=(4.5 * len(runs), 4), squeeze=False)
        for ax, (i, p, k) in zip(axes[0], runs):
            s = self._load(i, "classical/section/phase_space.dat")
            o = self._load(i, "classical/orbit/phase_space.dat")
            ax.plot(s[:, 2], s[:, 3], ",", color="k", alpha=0.5, rasterized=True)
            ax.plot(o[:, 2], o[:, 3], ".", color="C3", ms=2, rasterized=True)
            ax.set_xlim(-np.pi, np.pi)
            ax.set_ylim(-1, 1)
            ax.set_xlabel(r"$\phi$")
            ax.set_title(f"$p = {p:g}$, $k = {k:g}$")
        axes[0, 0].set_ylabel(r"$Z = \cos\theta$")
        fig.suptitle("Classical kicked top, Poincare sections")
        fig.tight_layout()
        fig.savefig(self.out_dir / "sections.png", dpi=200)
        return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", default="output_combined", help="directory with the kicked_top_combined output")
    parser.add_argument("--show", action="store_true", help="open the figures in a window")
    args = parser.parse_args()

    rp = CombinedPostProcessor(args.out_dir)
    rp.plot_sections()
    rp.plot_husimi()
    rp.plot_fidelity()
    rp.plot_total_spin()
    rp.plot_site_spins()
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

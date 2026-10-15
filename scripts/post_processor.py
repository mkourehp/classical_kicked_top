"""Plots the plot-ready data written by the Fortran kicked top code."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


class PostProcessor:
    def __init__(self, out_dir="output"):
        self.out_dir = Path(out_dir)

    def _load(self, name):
        path = self.out_dir / name
        with path.open() as f:
            title = f.readline().lstrip("#").strip()
        return title, np.loadtxt(path)

    def plot_phase_space(self):
        title, data = self._load("phase_space.dat")
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.scatter(data[:, 2], data[:, 3], c=data[:, 0], s=0.2, cmap="turbo", rasterized=True)
        ax.set_xlim(-np.pi, np.pi)
        ax.set_ylim(-1, 1)
        ax.set_xlabel(r"$\phi$")
        ax.set_ylabel(r"$Z = \cos\theta$")
        ax.set_title(f"Phase space, {title}")
        fig.tight_layout()
        fig.savefig(self.out_dir / "phase_space.png", dpi=200)
        return fig

    def plot_trajectory(self):
        title, data = self._load("trajectory.dat")
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(projection="3d")
        ax.scatter(data[:, 2], data[:, 3], data[:, 4], c=data[:, 0], s=0.2, cmap="turbo")
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title(f"Trajectory on the sphere, {title}")
        fig.tight_layout()
        fig.savefig(self.out_dir / "trajectory.png", dpi=200)
        return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="output", help="directory with the Fortran output")
    parser.add_argument("--sphere", action="store_true", help="also plot the 3D trajectory on the sphere")
    parser.add_argument("--show", action="store_true", help="open the figures in a window")
    args = parser.parse_args()

    pp = PostProcessor(args.out_dir)
    pp.plot_phase_space()
    if args.sphere:
        pp.plot_trajectory()
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

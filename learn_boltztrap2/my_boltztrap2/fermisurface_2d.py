# -*- coding: utf-8 -*
#    my_boltztrap2, a program for interpolating band structures and calculating
#    semi-classical transport coefficients.
#    Copyright (C) 2017-2025 Georg K. H. Madsen <georg.madsen@tuwien.ac.at>
#    Copyright (C) 2017-2025 Jesús Carrete <jesus.carrete.montana@tuwien.ac.at>
#    Copyright (C) 2017-2025 Matthieu J. Verstraete <matthieu.verstraete@ulg.ac.be>
#    Copyright (C) 2018-2019 Genadi Naydenov <gan503@york.ac.uk>
#    Copyright (C) 2020 Gavin Woolman <gwoolma2@staffmail.ed.ac.uk>
#    Copyright (C) 2020 Roman Kempt <roman.kempt@tu-dresden.de>
#    Copyright (C) 2022 Robert Stanton <stantor@clarkson.edu>
#    Copyright (C) 2024 Haoyu (Daniel) Yang <yanghaoyu97@outlook.com>
#
#    This file is part of my_boltztrap2.
#
#    my_boltztrap2 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    my_boltztrap2 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with my_boltztrap2. If not, see <http://www.gnu.org/licenses/>.


import argparse
import colorsys
import itertools
import sys

import matplotlib
import matplotlib.colors as mc
import matplotlib.gridspec as gridspec
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.linalg as la
import scipy.spatial
from matplotlib.widgets import Slider

import my_boltztrap2 as bt2
import my_boltztrap2.fite
import my_boltztrap2.serialization
from my_boltztrap2.units import eV


def plot_fermisurface_2d(
    data, equivalences, coeffs, mu, unit="eV", broadening=0.05
):
    assert unit in [
        "eV",
        "Hartree",
        "a.u.",
        "hartree",
        "au",
    ], "unit not recognized"
    # Check that this is indeed a 2D system oriented perpendicular to OZ.
    if not np.allclose(0.0, data.kpoints[:, 2]):
        sys.exit("Error: This file does not seem to describe a 2D system.")
    lattvec = data.get_lattvec()
    if not (
        np.allclose(0.0, lattvec[2, :2]) and np.allclose(0.0, lattvec[:2, 2])
    ):
        sys.exit("Error: The system is not perpendicular to the OZ axis.")

    # Rebuild the bands
    ebands_in = bt2.fite.getBTPbands(equivalences, coeffs, lattvec)[0]

    # Get the extent of the regular grid in reciprocal space
    hdims = np.max(np.abs(np.vstack(equivalences)), axis=0)
    dims = 2 * hdims + 1

    # Compute the coordinates of the points of the grid.
    rlattvec = 2.0 * np.pi * la.inv(lattvec).T
    lattvec = lattvec[:2, :2]
    rlattvec = rlattvec[:2, :2]
    a0, b0 = np.meshgrid(
        np.arange(dims[0], dtype=np.float64) / dims[0],
        np.arange(dims[1], dtype=np.float64) / dims[1],
        indexing="ij",
    )
    xy0 = (
        a0[..., np.newaxis] * rlattvec[:, 0]
        + b0[..., np.newaxis] * rlattvec[:, 1]
    )

    # Note that the bands are still expressed on a 3D grid, and in general that
    # grid contains more than one point along the third reciprocal axis. We
    # therefore take a slice to sample the energies on our 2D grid.
    ebands_in = ebands_in.reshape([-1] + dims.tolist())[..., hdims[-1]]

    # Extend the array with the positions of three more images of the grid
    # and create an array of band energies corresponding to all the points
    # in the extended set.
    ebands = []
    xy = []
    for i_a, i_b in itertools.product(range(-1, 1), repeat=2):
        xy.append(xy0 + i_a * rlattvec[:2, 0] + i_b * rlattvec[:2, 1])
        ebands.append(ebands_in)
    # Play with the shapes a little bit to make the indexing simpler and to try
    # to avoid weird teleconnections in the contour plot.
    xy = np.stack(xy).reshape((2, 2, dims[0], dims[1], 2))
    xy = xy.transpose((0, 2, 1, 3, 4)).reshape(
        (2 * xy0.shape[0], 2 * xy0.shape[1], 2)
    )
    ebands = np.stack(ebands).reshape((2, 2, -1, dims[0], dims[1]))
    ebands = ebands.transpose((2, 0, 3, 1, 4)).reshape(
        (-1, 2 * xy0.shape[0], 2 * xy0.shape[1])
    )

    # Center the energies at the Fermi level.
    ebands -= data.fermi

    # I'm a chemist, I prefer having stuff in eV
    if unit == "eV":
        ebands /= eV

    # Obtain the first Brillouin zone as the Voronoi polygon of the
    # Gamma point.
    points = []
    for ij0 in itertools.product(range(5), repeat=2):
        ij = [i if i <= 2 else i - 5 for i in ij0]
        points.append(rlattvec[:2, :2] @ np.array(ij))
    voronoi = scipy.spatial.Voronoi(points)
    region_index = voronoi.point_region[0]
    vertex_indices = voronoi.regions[region_index]
    vertices = voronoi.vertices[vertex_indices, :]

    # Zooming in the axes limits on the polygon
    xmin = np.min(np.vstack(vertices)[:, 0]) * 1.05
    xmax = np.max(np.vstack(vertices)[:, 0]) * 1.05
    ymin = np.min(np.vstack(vertices)[:, 1]) * 1.05
    ymax = np.max(np.vstack(vertices)[:, 1]) * 1.05

    # All the elements are ready. Now, create the figure.

    fig = plt.figure(constrained_layout=True, figsize=(5, 5))
    spec = gridspec.GridSpec(
        ncols=1, nrows=2, figure=fig, height_ratios=[15, 1]
    )
    ax = fig.add_subplot(spec[0])
    axbar = fig.add_subplot(spec[1])

    def adjust_lightness(color, amount=0.5):
        """Helper function to adjust brightness of a color."""

        try:
            c = mc.cnames[color]
        except KeyError:
            c = color
        c = colorsys.rgb_to_hls(*mc.to_rgb(c))
        return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

    def plot_fermi(ax, vertices, ebands, mu):
        # Add the first BZ as a polygon.
        polygon = matplotlib.patches.Polygon(
            vertices, facecolor="none", edgecolor="black", linewidth=2
        )
        ax.add_patch(polygon)
        # Plot some contours and use the same polygon as a clipping path.

        # Generating a linear color map for the bands. The bands at the fermi-level are going to fall into the middle of the color map.
        cmaps = matplotlib.colors.LinearSegmentedColormap.from_list(
            "",
            [
                "black",
                "purple",
                "green",
                "royalblue",
                "crimson",
                "orange",
                "black",
            ],
        )

        # Adding a broadening around the chemical potential of 50 meV.
        levels = sorted([mu - broadening, mu + broadening])

        for i in range(ebands.shape[0]):
            # Each band gets a different color which is shaded according to the broadening, giving a sense of depth.
            j = i / ebands.shape[0]
            c1, c2 = adjust_lightness(cmaps(j), 0.75), adjust_lightness(
                cmaps(j), 1.25
            )
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "", [c1, c2]
            )
            contours = ax.contour(
                xy[..., 0],
                xy[..., 1],
                ebands[i, ...],
                levels=np.arange(levels[0], levels[1], 0.01),
                cmap=cmap,
            )
            for c in contours.collections:
                c.set_linestyle("solid")
                c.set_clip_path(polygon)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_frame_on(False)
        ax.set_aspect("equal")
        ax.set_xlim((xmin, xmax))
        ax.set_ylim((ymin, ymax))

    plot_fermi(ax, vertices, ebands, mu)

    spot = Slider(
        axbar, r"$\mu$ [eV]", -3, +3, valinit=0, valstep=0.01, valfmt="% 2.2f"
    )

    def update(ax):
        nu = spot.val
        ax.clear()
        plot_fermi(ax, vertices, ebands, nu)
        ax.set_frame_on(False)
        fig.canvas.draw_idle()

    spot.on_changed(lambda x: update(ax))
    plt.show()

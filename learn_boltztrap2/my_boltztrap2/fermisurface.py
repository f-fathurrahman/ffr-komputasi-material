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

import itertools

import matplotlib
import matplotlib.colors as colors
import numpy as np

# Since we only use the module for shifts in this this part of the code,
# we do not bother finding the highest-performance implementation here.
import scipy as sp
import scipy.linalg as la
import scipy.spatial

from my_boltztrap2.misc import TimerContext, info, warning

# If VTK cannot be imported, the module-level variable "available" will be set
# to False and the features of this submodule will fail to work.
try:
    import vtk

    available = True
except ImportError:
    warning("vtk not found. The 'fermisurface' command will be disabled")
    available = False

# Colors for all elements other than the surface itself.
OTHER_COLORS = dict(
    slider_tube="#2e3436",
    slider_slider="#a40000",
    slider_cap="#babdb6",
    slider_selected="#a40000",
    slider_title="#2e3436",
    brillouin_zone="#204a87",
    sphere="#f57900",
)


def plot_fermisurface(
    data, equivalences, ebands_in, mu, color_cycle=None, edge_thickness=1.0
):
    """Launch an interactive VTK representation of the Fermi surface.

    Make sure to check the module-level variable "available" before calling
    this function.

    Args:
        data: a DFTData object
        equivalences: list of k-point equivalence classes
        ebands_in: (nbands, nkpoints) array with the band energies
        mu: initial value of the energy at which the surface will be plotted,
            with respect to data.fermi. This can later be changed by the user
            using an interactive slider.
        color_cycle: sequence of colors to loop over, as hex strings. If it is
            None, use matplotlib's default.
        edge_thickness: thickness of the edges of the first BZ in the plot. A
            non-positive value causes the edges to be hidden.

    Returns:
        None.
    """
    lattvec = data.get_lattvec()
    rlattvec = 2.0 * np.pi * la.inv(lattvec).T

    if color_cycle is None:
        color_cycle = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]

    # Find the ratio between the circumscribed and inscribed radii of the BZ.
    bz_volume = abs(la.det(rlattvec))
    bz_max_area = max(
        la.norm(np.cross(rlattvec[:, i], rlattvec[:, j]))
        for (i, j) in itertools.product(range(3), repeat=2)
    )
    bz_r_in = 0.5 * bz_volume / bz_max_area
    bz_corners = np.asarray(
        [
            np.zeros(3),
            rlattvec[:, 0],
            rlattvec[:, 1],
            rlattvec[:, 2],
            rlattvec[:, 0] + rlattvec[:, 1],
            rlattvec[:, 0] + rlattvec[:, 2],
            rlattvec[:, 1] + rlattvec[:, 2],
            rlattvec[:, 0] + rlattvec[:, 1] + rlattvec[:, 2],
        ]
    )
    distances = sp.spatial.distance.pdist(bz_corners)
    bz_r_out = 0.5 * distances.max()
    MAX_IMAGE = 2 * (int(np.ceil(bz_r_out / bz_r_in)) // 2) + 1

    # Obtain the first Brillouin zone as the Voronoi polyhedron of Gamma
    points = []
    RANGE_LIMIT = 2 * MAX_IMAGE + 1
    for ijk0 in itertools.product(range(RANGE_LIMIT), repeat=3):
        ijk = [i if i <= MAX_IMAGE else i - RANGE_LIMIT for i in ijk0]
        points.append(rlattvec @ np.array(ijk))
    voronoi = scipy.spatial.Voronoi(points)
    region_index = voronoi.point_region[0]
    vertex_indices = voronoi.regions[region_index]
    vertices = voronoi.vertices[vertex_indices, :]

    # Compute a center and an outward-pointing normal for each of the facets
    # of the BZ
    facets = []
    for ridge in voronoi.ridge_vertices:
        if all(i in vertex_indices for i in ridge):
            facets.append(ridge)
    centers = []
    normals = []
    for f in facets:
        corners = np.array([voronoi.vertices[i, :] for i in f])
        center = corners.mean(axis=0)
        v1 = corners[0, :]
        for i in range(1, corners.shape[0]):
            v2 = corners[i, :]
            prod = np.cross(v1 - center, v2 - center)
            if not np.allclose(prod, 0.0):
                break
        if np.dot(center, prod) < 0.0:
            prod = -prod
        centers.append(center)
        normals.append(prod)

    # Get the extent of the regular grid in reciprocal space
    hdims = np.max(np.abs(np.vstack(equivalences)), axis=0)
    dims = 2 * hdims + 1

    class FermiInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
        """Custom interaction style enabling the user to pick points on
        the screen and overriding some default VTK keyboard bindings.
        """

        def __init__(self, parent=None):
            """Simple constructor that adds an observer to the middle mouse
            button press and another one to "char events" (originating from key
            presses).
            """
            self.AddObserver("MiddleButtonPressEvent", self.pick_point)
            # Remove the original char event handler before installing our own.
            self.RemoveObservers("CharEvent")
            self.AddObserver("CharEvent", self.process_keypress)

        def pick_point(self, obj, event):
            """Get the coordinates of the point selected with the middle mouse
            button, find the nearest data point, and print its direct
            coordinates.
            """
            interactor = self.GetInteractor()
            picker = interactor.GetPicker()
            pos = interactor.GetEventPosition()
            picker.Pick(
                pos[0],
                pos[1],
                0,
                interactor.GetRenderWindow().GetRenderers().GetFirstRenderer(),
            )
            picked = np.array(picker.GetPickPosition())
            # Move the sphere to the new coordinates and make it visible
            sphere.SetCenter(*picked.tolist())
            sphere_actor.VisibilityOn()
            picked = la.solve(rlattvec, picked)
            print("Point picked:", picked)
            self.OnMiddleButtonDown()

        def process_keypress(self, obj, event):
            """Capture 's', which should not do anything, and 'w', which should
            only toggle between surface and wireframe modes for the Fermi
            surface and not for the BZ polyhedron. Any other event is passed
            on to the handler in the base class: for instance, pressing 'e'
            will still close the window and exit the program, as usual.
            """
            interactor = self.GetInteractor()
            key = interactor.GetKeySym()
            if key == "s":
                return
            elif key == "w":
                if len(fermiactors) == 0:
                    return
                representation = (
                    fermiactors[0].GetProperty().GetRepresentationAsString()
                )
                if representation == "Wireframe":
                    for f in fermiactors:
                        f.GetProperty().SetRepresentationToSurface()
                else:
                    for f in fermiactors:
                        f.GetProperty().SetRepresentationToWireframe()
                window.Render()
            else:
                super().OnChar()

    # Create the VTK representation of the grid.
    # Create an array with the direct coordinates of the points in the
    # grid, from 0 to 1.
    a0, b0, c0 = np.meshgrid(
        np.arange(dims[0], dtype=np.float64) / dims[0],
        np.arange(dims[1], dtype=np.float64) / dims[1],
        np.arange(dims[2], dtype=np.float64) / dims[2],
        indexing="ij",
    )
    xyz0 = (
        a0[..., np.newaxis] * rlattvec[:, 0]
        + b0[..., np.newaxis] * rlattvec[:, 1]
        + c0[..., np.newaxis] * rlattvec[:, 2]
    )
    # Extend the array with the positions of seven more images of the grid
    # and create an array of band energies corresponding to all the points
    # in the extended set.
    ebands = []
    xyz = []
    for i_a, i_b, i_c in itertools.product(range(-1, 1), repeat=3):
        xyz.append(
            xyz0
            + i_a * rlattvec[:, 0]
            + i_b * rlattvec[:, 1]
            + i_c * rlattvec[:, 2]
        )
        ebands.append(ebands_in)
    # Play with the axes of these arrays to massage them into a format VTK
    # likes.
    xyz = np.stack(xyz).reshape((2, 2, 2, dims[0], dims[1], dims[2], 3))
    # Note that for VTK the "a" coordinate changes faster than the "b"
    # coordinate, which in turns changes faster than the "c" coordinate.
    # The convention used in my_boltztrap2 is the inverse.
    xyz = xyz.transpose((2, 5, 1, 4, 0, 3, 6))
    xyz = xyz.reshape((-1, 3))
    ebands = np.stack(ebands).reshape((2, 2, 2, -1, dims[0], dims[1], dims[2]))
    ebands = ebands.transpose((3, 2, 6, 1, 5, 0, 4))
    ebands = ebands.reshape((ebands.shape[0], -1))

    sgrid = vtk.vtkStructuredGrid()
    sgrid.SetDimensions(2 * dims[0], 2 * dims[1], 2 * dims[2])
    spoints = vtk.vtkPoints()
    for p in xyz:
        spoints.InsertNextPoint(*p.tolist())
    sgrid.SetPoints(spoints)
    ebands -= data.fermi
    emax = ebands.max(axis=1)
    emin = ebands.min(axis=1)

    # Find the shortest distance between points to compute a good
    # radius for the selector sphere later.
    dmin = np.inf
    for i in range(3):
        abc = np.zeros(3)
        abc[i] = 1.0 / dims[i]
        cart = rlattvec @ abc
        if la.norm(cart) > 0:
            dmin = min(dmin, la.norm(cart))

    # Create a 2D chemical potential slider
    slider = vtk.vtkSliderRepresentation2D()
    slider.SetMinimumValue(emin.min())
    slider.SetMaximumValue(emax.max())
    slider.SetValue(mu)
    slider.SetTitleText("Chemical potential")
    slider.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider.GetPoint1Coordinate().SetValue(0.1, 0.9)
    slider.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider.GetPoint2Coordinate().SetValue(0.9, 0.9)
    slider.GetTubeProperty().SetColor(
        *colors.hex2color(OTHER_COLORS["slider_tube"])
    )
    slider.GetSliderProperty().SetColor(
        *colors.hex2color(OTHER_COLORS["slider_slider"])
    )
    slider.GetCapProperty().SetColor(
        *colors.hex2color(OTHER_COLORS["slider_cap"])
    )
    slider.GetSelectedProperty().SetColor(
        *colors.hex2color(OTHER_COLORS["slider_selected"])
    )
    slider.GetTitleProperty().SetColor(
        *colors.hex2color(OTHER_COLORS["slider_title"])
    )

    # Find all the isosurfaces with energy equal to the threshold
    allcontours = []
    with TimerContext() as timer:
        fermiactors = []
        for iband, band in enumerate(ebands):
            band_color = color_cycle[iband % len(color_cycle)]
            sgridp = vtk.vtkStructuredGrid()
            sgridp.DeepCopy(sgrid)
            # Feed the energies to VTK
            scalar = vtk.vtkFloatArray()
            for i in band:
                scalar.InsertNextValue(i)
            sgridp.GetPointData().SetScalars(scalar)
            # Estimate the isosurfaces
            contours = vtk.vtkMarchingContourFilter()
            contours.SetInputData(sgridp)
            contours.UseScalarTreeOn()
            contours.SetValue(0, mu)
            contours.ComputeNormalsOff()
            contours.ComputeGradientsOff()
            allcontours.append(contours)

            # Create the set of clipping planes representing the boundaries
            # of the first Brillouin zone and apply the clipping.
            clipping_centers = vtk.vtkPoints()
            clipping_normals = vtk.vtkDoubleArray()
            clipping_normals.SetNumberOfComponents(3)
            for c, n in zip(centers, normals):
                clipping_centers.InsertNextPoint(*c)
                clipping_normals.InsertNextTuple3(*n)
            clipping_planes = vtk.vtkPlanes()
            clipping_planes.SetPoints(clipping_centers)
            clipping_planes.SetNormals(clipping_normals)
            clipper = vtk.vtkClipPolyData()
            clipper.SetClipFunction(clipping_planes)
            clipper.InsideOutOn()
            clipper.SetInputConnection(contours.GetOutputPort())

            # Compute the normals to the surfaces to obtain better lighting
            normals_calculator = vtk.vtkPolyDataNormals()
            normals_calculator.SetInputConnection(clipper.GetOutputPort())
            normals_calculator.ComputeCellNormalsOn()
            normals_calculator.ComputePointNormalsOn()

            # Create a mapper and an actor for the surfaces
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(normals_calculator.GetOutputPort())
            mapper.ScalarVisibilityOff()
            fermiactors.append(vtk.vtkActor())
            fermiactors[-1].SetMapper(mapper)
            fermiactors[-1].GetProperty().SetColor(
                *colors.hex2color(band_color)
            )
        deltat = timer.get_deltat()
        info("building the VTK surfaces took {:.3g} s".format(deltat))

    # Represent the BZ as a polyhedron in VTK
    points = vtk.vtkPoints()
    for v in voronoi.vertices:
        points.InsertNextPoint(*v)
    fids = vtk.vtkIdList()
    fids.InsertNextId(len(facets))
    for f in facets:
        fids.InsertNextId(len(f))
        for i in f:
            fids.InsertNextId(i)
    fgrid = vtk.vtkUnstructuredGrid()
    fgrid.SetPoints(points)
    fgrid.InsertNextCell(vtk.VTK_POLYHEDRON, fids)

    # Create an actor and a mapper for the BZ
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(fgrid)
    bzactor = vtk.vtkActor()
    bzactor.SetMapper(mapper)
    bzactor.GetProperty().SetRepresentationToWireframe()
    bzactor.GetProperty().SetColor(
        *colors.hex2color(OTHER_COLORS["brillouin_zone"])
    )
    if edge_thickness > 0.0:
        bzactor.GetProperty().SetLineWidth(edge_thickness)

    # Create a visual representation of the selected point, and hide
    # it for the time being.
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(dmin / 2.0)
    sphere_mapper = vtk.vtkPolyDataMapper()
    sphere_mapper.SetInputConnection(sphere.GetOutputPort())
    sphere_mapper.ScalarVisibilityOff()
    sphere_actor = vtk.vtkActor()
    sphere_actor.SetMapper(sphere_mapper)
    sphere_actor.GetProperty().SetColor(
        *colors.hex2color(OTHER_COLORS["sphere"])
    )
    sphere_actor.VisibilityOff()

    # Create a VTK window and other elements of an interactive scene
    renderer = vtk.vtkRenderer()
    renderer.AddActor(bzactor)
    renderer.AddActor(sphere_actor)
    for f in fermiactors:
        renderer.AddActor(f)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Zoom(5.0)
    renderer.SetBackground(1.0, 1.0, 1.0)
    if edge_thickness <= 0.0:
        # If we need to hide the edges of the BZ, do it here after if has been
        # taken into account when zooming.
        bzactor.VisibilityOff()

    window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetInteractorStyle(FermiInteractorStyle())
    interactor.SetRenderWindow(window)

    # Add a set of axes
    axes = vtk.vtkAxesActor()
    assembly = vtk.vtkPropAssembly()
    assembly.AddPart(axes)
    marker = vtk.vtkOrientationMarkerWidget()
    marker.SetOrientationMarker(assembly)
    marker.SetInteractor(interactor)
    marker.SetEnabled(1)
    marker.InteractiveOff()

    def callback(obj, ev):
        """Update the isosurface with a new value."""
        mu = obj.GetRepresentation().GetValue()
        for e, E, c, a in zip(emin, emax, allcontours, fermiactors):
            visible = e <= mu and E >= mu
            a.SetVisibility(visible)
            if visible:
                c.SetValue(0, mu)

    # Add the slider widget
    widget = vtk.vtkSliderWidget()
    widget.SetInteractor(interactor)
    widget.SetRepresentation(slider)
    widget.SetAnimationModeToJump()
    widget.EnabledOn()
    widget.AddObserver(vtk.vtkCommand.InteractionEvent, callback)

    # Launch the visualization
    interactor.Initialize()
    window.Render()
    interactor.Start()

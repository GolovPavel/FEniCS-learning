import numpy as np
import pygmsh


class Cylinder:
    def __init__(self, radius, height, mesh_size):
        self.radius = radius
        self.height = height
        self.mesh_size = mesh_size


class GmshElementary:
    def __init__(self, points, circles, circle_ll, lines, surface_ll, surfaces, plane_surfaces, surface_loop, volume):
        self.points = points
        self.circles = circles
        self.circles_ll = circle_ll
        self.lines = lines
        self.surface_ll = surface_ll
        self.surfaces = surfaces
        self.plane_surfaces = plane_surfaces
        self.surface_loop = surface_loop
        self.volume = volume


def create_geometry(cylinders):
    geometry = pygmsh.built_in.Geometry()
    prev_elementary = None
    for c in cylinders:
        prev_elementary = __add_cylinder(geometry, c.radius, c.mesh_size, c.height, prev_elementary)
    return geometry


def __add_cylinder(geom, radius, lcar, H, prev_elementary):
    base_coords = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.0, 0.0, H]),
    ]

    coords = [
        np.array([radius, 0.0, 0.0]),
        np.array([0.0, radius, 0.0]),
        np.array([-radius, 0.0, 0.0]),
        np.array([0.0, -radius, 0.0]),
        np.array([radius, 0.0, H]),
        np.array([0.0, radius, H]),
        np.array([-radius, 0.0, H]),
        np.array([0.0, -radius, H])
    ]

    circle_points = [
        np.array([2, 0, 3]),
        np.array([3, 0, 4]),
        np.array([4, 0, 5]),
        np.array([5, 0, 2]),

        np.array([6, 1, 7]),
        np.array([7, 1, 8]),
        np.array([8, 1, 9]),
        np.array([9, 1, 6])
    ]

    ll_circles = [
        np.array([0, 1, 2, 3]),
        np.array([4, 5, 6, 7]),
    ]

    line_points = [
        np.array([6, 2]),
        np.array([7, 3]),
        np.array([8, 4]),
        np.array([9, 5]),
    ]

    ll_surfaces = [
        np.array([0, 0, 1, 4]),
        np.array([1, 1, 2, 5]),
        np.array([2, 2, 3, 6]),
        np.array([3, 3, 0, 7])
    ]

    points = [geom.add_point(point, lcar=0.0) for point in base_coords] + \
             [geom.add_point(point, lcar=lcar) for point in coords]
    circles = [geom.add_circle_arc(points[idx[0]], points[idx[1]], points[idx[2]]) for idx in circle_points]
    circle_ll = [geom.add_line_loop([circles[idx[0]], circles[idx[1]], circles[idx[2]], circles[idx[3]]])
                 for idx in ll_circles]
    lines = [geom.add_line(points[idx[0]], points[idx[1]]) for idx in line_points]
    surface_ll = [geom.add_line_loop([lines[idx[0]], circles[idx[1]], -lines[idx[2]], -circles[idx[3]]])
                  for idx in ll_surfaces]
    surfaces = [geom.add_surface(surface) for surface in surface_ll]
    if prev_elementary is None:
        plane_surfaces = [geom.add_plane_surface(circle_loop) for circle_loop in circle_ll]
        surface_loop = geom.add_surface_loop([plane_surfaces[0], *surfaces, plane_surfaces[1]])
    else:
        plane_surfaces = [geom.add_plane_surface(circle_loop, [prev_elementary.circles_ll[i]])
                          for i, circle_loop in enumerate(circle_ll)]
        surface_loop = geom.add_surface_loop(
            [plane_surfaces[0], *surfaces, plane_surfaces[1], *prev_elementary.surfaces])
    volume = geom.add_volume(surface_loop)

    return GmshElementary(points, circles, circle_ll, lines, surface_ll, surfaces, plane_surfaces, surface_loop, volume)

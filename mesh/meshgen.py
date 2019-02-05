import meshio
import pygmsh


class Cylinder:
    def __init__(self, radius, height, mesh_size):
        self.radius = radius
        self.height = height
        self.mesh_size = mesh_size


def create_2d_geometry(cylinders):
    geometry = pygmsh.opencascade.Geometry()
    for cylinder in cylinders:
        geometry.add_cylinder(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, cylinder.height],
            cylinder.radius,
            char_length=cylinder.mesh_size
        )
    return geometry


cylinders = [
    Cylinder(100, 2000, 15),
    Cylinder(1000, 2000, 150),
    Cylinder(5000, 2000, 700),
    Cylinder(30000, 2000, 3000)
]

if __name__ == '__main__':
    geom = create_2d_geometry(cylinders)
    points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom)
    meshio.write_points_cells("cylinder_2d.vtk", points, cells)

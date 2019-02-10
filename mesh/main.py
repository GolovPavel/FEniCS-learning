import meshio

from meshgen import *

cylinders = [
    Cylinder(100, 2000, 15),
    Cylinder(1000, 2000, 150),
    Cylinder(5000, 2000, 700),
    Cylinder(30000, 2000, 3000)
]

GEO_FILE_NAME = "30000_2000h.geo"
MESH_2D_FILE_NAME = "30000_2000h_2d.vtk"
MESH_3D_FILE_NAME = "30000_2000h_3d.vtk"

if __name__ == '__main__':
    geom = create_geometry(cylinders)

    points2d, cells2d, _, _, _ = pygmsh.generate_mesh(geom, geo_filename=GEO_FILE_NAME, dim=2)
    points3d, cells3d, _, _, _ = pygmsh.generate_mesh(geom, dim=3)

    meshio.write_points_cells(MESH_2D_FILE_NAME, points2d, cells2d)
    meshio.write_points_cells(MESH_3D_FILE_NAME, points3d, cells3d)

# About physical points
# http://onelab.info/pipermail/gmsh/2010/005511.html

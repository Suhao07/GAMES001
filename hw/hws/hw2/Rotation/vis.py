import pyvista as pv
import numpy as np
import os

import rot.euler
import rot.rotvec
from rot.quat import Quaternion

pv.global_theme.lighting_params.ambient = 0.15
plotter = pv.Plotter()
plotter.enable_anti_aliasing()
plotter.camera_position = [(-4, 2, -4), (0, 0, 0), (0, 1, 0)]
plotter.show_axes()

cube = pv.Cube(center=(0, 0, 0), x_length=2, y_length=2, z_length=2)
plotter.add_mesh(cube, color="black", style="wireframe")

spot = pv.read("./assets/mesh/spot/spot_subdivision.obj")
spot_tex = pv.read_texture("./assets/mesh/spot/spot_texture.png")
spot.compute_normals(point_normals=True, inplace=True)
spot.points -= spot.center
rest_pos = spot.points.copy()
plotter.add_mesh(spot, smooth_shading=True, texture=spot_tex)

n_frame = 60
os.makedirs("output", exist_ok=True)

rot0 = -np.pi / 2 * 3 * np.array([1, 0, 0])
rot1 = -np.pi / 2 * 3 * np.array([0, 1, 0])

plotter.open_gif("output/01-direct_lerp.gif", fps=30)
for t in np.linspace(0, 1, n_frame):
    spot.points = rest_pos @ rot.rotvec.direct_lerp(rot0, rot1, t).T
    plotter.write_frame()

plotter.open_gif("output/02-rel_lerp.gif", fps=30)
for t in np.linspace(0, 1, n_frame):
    spot.points = rest_pos @ rot.rotvec.rel_lerp(rot0, rot1, t).T
    plotter.write_frame()

qrot0 = Quaternion([0, 1, 0, 0])
qrot1 = Quaternion([0, -1, 1, 0])

plotter.open_gif("output/03-nlerp.gif", fps=30)
for t in np.linspace(0, 1, n_frame):
    spot.points = rest_pos @ Quaternion.Nlerp(qrot0, qrot1, t).to_matrix().T
    plotter.write_frame()

plotter.open_gif("output/04-slerpv1.gif", fps=30)
for t in np.linspace(0, 1, n_frame):
    spot.points = rest_pos @ Quaternion.Slerp(qrot0, qrot1, t, False).to_matrix().T
    plotter.write_frame()

plotter.open_gif("output/05-slerpv2.gif", fps=30)
for t in np.linspace(0, 1, n_frame):
    spot.points = rest_pos @ Quaternion.Slerp(qrot0, qrot1, t).to_matrix().T
    plotter.write_frame()

plotter.close()

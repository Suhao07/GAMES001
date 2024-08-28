import numpy as np
import matplotlib.pyplot as plt
from utils import *
import os

os.makedirs('output', exist_ok=True)

resolution = 64
grid_size = (resolution, resolution)

u = np.zeros(grid_size)
boundary_value = [-1, 1, -1, 1]
u[0, 1:-1] = boundary_value[0]
u[1:-1, 0] = boundary_value[1]
u[-1, 1:-1] = boundary_value[2]
u[1:-1, -1] = boundary_value[3]

iters = 300
result = []
residual = 0
residual_init = compute_laplace_residual(u)

for iter in range(iters + 1):
    iter_jacobi(u)
    if iter % 60 == 0:
        residual = compute_laplace_residual(u)
        print("{:.3e}".format(residual / residual_init))
        result.append(u.copy())


fig, axes = plt.subplots(2, 3)
for i, ax in enumerate(axes.flatten()):
    ax.imshow(result[i], cmap='bwr', extent=[0, 1, 0, 1], origin='lower')
    ax.set_title('iter={}'.format(i * 60))
plt.tight_layout()
plt.savefig('output/laplace.png')

def Ai(n):
    if n % 2 == 0:
        return -4 / ((n+1) * np.pi * np.sinh((n+1)*np.pi))
    
    return 0

def phi1(k = 3):
    nx = ny = resolution
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)

    xv, yv = np.meshgrid(x, y)

    phi = np.zeros((nx, ny))
    for i in range(k):
        n = 2*i
        sinh = np.sinh((n+1) * np.pi * (xv-1))
        sin  = np.sin((n+1) * np.pi * yv)

        phi += Ai(n) * sinh * sin

    return phi

def phi2(k = 3):
    nx = ny = resolution
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)

    xv, yv = np.meshgrid(x, y)

    phi = np.zeros((nx, ny))
    for i in range(k):
        n = 2*i
        sinh = np.sinh((n+1) * np.pi * (yv-1))
        sin  = np.sin((n+1) * np.pi * xv)

        phi += -Ai(n) * sinh * sin

    return phi

def phi3(k = 3):
    nx = ny = resolution
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)

    xv, yv = np.meshgrid(x, y)

    phi = np.zeros((nx, ny))
    for i in range(k):
        n = 2*i
        sinh = np.sinh((n+1) * np.pi * xv)
        sin  = np.sin((n+1) * np.pi * yv)

        phi += -Ai(n) * sinh * sin

    return phi

def phi4(k = 3):
    nx = ny = resolution
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)

    xv, yv = np.meshgrid(x, y)

    phi = np.zeros((nx, ny))
    for i in range(k):
        n = 2*i
        sinh = np.sinh((n+1) * np.pi * yv)
        sin  = np.sin((n+1) * np.pi * xv)

        phi += Ai(n) * sinh * sin

    return phi

k = 100
f1 = phi1(k)
f2 = phi2(k)
f3 = phi3(k)
f4 = phi4(k)
f  = f1+f2+f3+f4

res = [f1, f2, f3, f4]
titles = [r"$\phi_1 (x,y)$", r"$\phi_2 (x,y)$", r"$\phi_3 (x,y)$", r"$\phi_4 (x,y)$"]

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, axes = plt.subplots(2, 2)
for i, ax in enumerate(axes.flatten()):
    im = ax.imshow(res[i], cmap='bwr', vmax=1, vmin=-1, extent=[0, 1, 0, 1], origin='lower')
    ax.set_title(titles[i])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

plt.tight_layout()
plt.savefig('output/phi_grid.png')


fig, axes = plt.subplots(1, 1)
plt.gca().set_aspect('equal', adjustable='box')
plt.imshow(f, vmax=1, vmin=-1, cmap='bwr', extent=[0, 1, 0, 1])
plt.colorbar()
plt.title(r"$\phi (x, y)$")
plt.tight_layout()
plt.savefig('output/phi.png')

def relative_err(pred, gt):
    return np.abs(gt - pred) / (np.abs(gt) + 1e-7)

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, axes = plt.subplots(2, 3, figsize=(16,9))
for i, ax in enumerate(axes.flatten()):
    im = ax.imshow(relative_err(result[i], f), vmax=1.0, vmin=0.0, cmap='bwr', extent=[0, 1, 0, 1], origin='lower')
    ax.set_title('iter={}'.format(i * 60))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

plt.tight_layout()
plt.savefig('output/relative_err.png')
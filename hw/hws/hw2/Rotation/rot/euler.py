import numpy as np
import re


def euler2matrix(seq, angles, degrees=True) -> np.ndarray:
    num_axes = len(seq)
    if num_axes != 3:
        raise ValueError("Wrong euler sequence")
    if re.match(r"^[xyzXYZ]{3}$", seq) is None:
        raise ValueError("Wrong euler sequence")
    seq = seq.lower()
    if degrees:
        angles = np.deg2rad(angles)

    mat = np.eye(3)
    # TODO: your code here
    # 请将 mat 补全，使得函数返回正确的旋转矩阵
    # raise NotImplementedError

    for axis, theta in zip(seq, angles):
        if axis == "x":
            mat = RotX(theta) @ mat
        elif axis == "y":
            mat = RotY(theta) @ mat
        elif axis == "z":
            mat = RotZ(theta) @ mat

    return mat

def RotX(theta):
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([
                    [1, 0,  0],
                    [0, c, -s],
                    [0, s,  c]
    ])

def RotY(theta):
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([
                    [ c, 0, s],
                    [ 0, 1, 0],
                    [-s, 0, c]
    ])

def RotZ(theta):
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([
                    [c, -s, 0],
                    [s,  c, 0],
                    [0,  0, 1]
    ])
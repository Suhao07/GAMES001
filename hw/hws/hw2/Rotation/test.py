import rot.euler
import rot.rotvec
from rot.quat import Quaternion
from scipy.spatial.transform import Rotation as R
import numpy as np


def gen_rand_euler():
    seq = list("xyz")
    np.random.shuffle(seq)
    return ("".join(seq), np.random.rand(3) * 360)


if __name__ == "__main__":
    N = 10
    print("Testing euler2matrix... ", end=" ")
    for i in range(N):
        seq, angles = gen_rand_euler()
        a = R.from_matrix(rot.euler.euler2matrix(seq, angles)).as_quat()
        b = R.from_euler(seq, angles, degrees=True).as_quat()
        print(np.abs(1 - np.abs(np.dot(a, b))) < 1e-5, end=" ")
    print()
    print("Testing rotvec2matrix...", end=" ")
    for i in range(N):
        b = R.random()
        a = R.from_matrix(rot.rotvec.rotvec2matrix(b.as_rotvec())).as_quat()
        print(np.abs(1 - np.abs(np.dot(a, b.as_quat()))) < 1e-5, end=" ")
    print()
    print("Testing quat_mul...     ", end=" ")
    for i in range(N):
        b1 = R.random()
        b2 = R.random()
        q1 = b1.as_quat()
        q2 = b2.as_quat()
        q3 = (b1 * b2).as_quat()
        a = (Quaternion([q1[3], *q1[:3]]) * Quaternion([q2[3], *q2[:3]])).quat
        print(np.abs(1 - np.abs(np.dot(a, np.array([q3[3], *q3[:3]])))) < 1e-5, end=" ")
    print()
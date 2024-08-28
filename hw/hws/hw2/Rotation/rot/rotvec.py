import numpy as np


def rotvec2matrix(vec, degrees=False) -> np.ndarray:
    # TODO: your code here
    # 利用罗德里格旋转公式返回轴角表示对应的旋转矩阵
    # raise NotImplementedError

    norm    = np.sum(vec*vec)
    theta   = np.sqrt(norm)
    
    if np.abs(theta) < 1e-5:
        return np.eye(3)
    
    u       = vec / theta
    cross_u = np.array([
                        [    0, -u[2],  u[1]],
                        [ u[2],     0, -u[0]],
                        [-u[1],  u[0],     0]
                    ])

    R = np.eye(3) + (1 - np.cos(theta)) * (cross_u @ cross_u) + np.sin(theta) * cross_u
    
    return R


def matrix2rotvec(matrix) -> np.ndarray:
    angle = np.arccos((np.trace(matrix) - 1) / 2)
    u = np.array(
        [
            matrix[2, 1] - matrix[1, 2],
            matrix[0, 2] - matrix[2, 0],
            matrix[1, 0] - matrix[0, 1],
        ]
    ) / (2 * np.sin(angle))
    return u * angle


def direct_lerp(a, b, t) -> np.ndarray:
    # TODO: your code here
    # 返回使用轴角表示直接插值的结果
    # raise NotImplementedError
    # return rotvec2matrix(b, False)

    ## linear interpolation of rotation vectors
    vt = (1. - t) * a + t * b
    Rt = rotvec2matrix(vt, False)
    
    return Rt


def rel_lerp(a, b, t) -> np.ndarray:
    # TODO: your code here
    # 返回使用轴角表示相对插值的结果
    # raise NotImplementedError

    # return rotvec2matrix(b, False)

    ## relative rotation
    R0 = rotvec2matrix(a, False)
    R1 = rotvec2matrix(b, False)
    dR = R1 @ R0.T

    ## rotation matrix to axis-angle
    vec = matrix2rotvec(dR)

    norm  = np.sum(vec*vec)
    theta = np.sqrt(norm)

    if np.abs(theta) < 1e-5:
        return np.eye(3)
    
    u     = vec / theta

    theta_t = t * theta

    ## axis-angle to rotation matrix
    dRt = rotvec2matrix(theta_t * u)

    Rt = dRt @ R0
    return Rt
import numpy as np


class Quaternion:
    # (w, ix, jy, kz)
    def __init__(self, quat=[0, 1, 0, 0]) -> None:
        self.quat = self.normalize(quat)

    def normalize(self, quat) -> np.ndarray:
        return quat / np.linalg.norm(quat)

    def to_matrix(self) -> np.ndarray:
        w, x, y, z = self.quat
        return np.array(
            [
                [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
            ]
        )

    def from_rotvec(self, vec):
        vec = np.array(vec)
        angle = np.linalg.norm(vec)
        self.quat = np.array([np.cos(angle / 2), *np.sin(angle / 2) * vec / angle])
        return self

    def __mul__(self, other):
        # TODO: your code here
        # 返回两个四元数的乘积
        # raise NotImplementedError
        a, b, c, d = self.quat
        e, f, g, h = other.quat

        return Quaternion([
                            a*e - b*f - c*g - d*h,
                            b*e + a*f - d*g + c*h,
                            c*e + d*f + a*g - b*h,
                            d*e - c*f + b*g + a*h
                        ])

    @staticmethod
    def Nlerp(a, b, t) -> "Quaternion":
        # TODO: your code here
        # 返回 Nlerp 插值结果
        # raise NotImplementedError
        # return b
    
        q0 = a.quat
        q1 = b.quat

        qt = np.array([(1-t)*x + t*y for x,y in zip(q0, q1)])
        qt = qt / np.linalg.norm(qt)

        return Quaternion(qt)

    @staticmethod
    def Slerp(a, b, t, short_path=True) -> "Quaternion":
        # TODO: your code here
        # short_path 为 True 时，返回最小弧的 Slerp 插值结果
        # 返回 Slerp 插值结果
        # raise NotImplementedError
        # return b

        q0 = a.quat
        q1 = b.quat

        theta = np.arccos(np.sum(q0 * q1))
        cos = np.cos(theta)

        if cos < 0 and short_path:
            bb = Quaternion([-x for x in q1])
            return Quaternion.Slerp(a, bb, t, False)
        
        sin   = np.sin(theta)
        alpha = np.sin((1-t) * theta) / sin
        beta  = np.sin(t * theta) / sin

        qt = np.array([alpha*x + beta*y for x,y in zip(q0, q1)])

        return Quaternion(qt)

import numpy as np


class Quaternion(np.ndarray):
    def __new__(cls, q1=0, q2=0, q3=0, q4=0, dtype=np.float64, *args, **kwargs):
        return np.asarray((q1, q2, q3,q4), dtype=dtype).view(cls)

    def __init__(self, q1=0, q2=0, q3=0, q4=0, dtype=np.float64):
        pass

    def __mul__(self, other):
        m = np.array([
            [ self[3], -self[2],  self[1],  self[0]],
            [ self[2],  self[3], -self[0],  self[1]],
            [-self[1],  self[0],  self[3],  self[2]],
            [-self[0], -self[1], -self[2],  self[3]]
        ])
        r = np.dot(m, other)
        return Quaternion(r[0], r[1], r[2], r[3])

    @property
    def conjugate(self):
        return Quaternion(-self[0], -self[1], -self[2], self[3])

    def apply(self, point):
        q = Quaternion(point[0], point[1], point[2], 0)
        temp1 = self * q
        temp2 = temp1 * self.conjugate
        return temp2[:3]

    def from_axis_angle(self, axis, angle):
        self[0] = axis[0] * np.sin(angle/2.0)
        self[1] = axis[1] * np.sin(angle / 2.0)
        self[2] = axis[2] * np.sin(angle / 2.0)
        self[3] = np.cos(angle / 2.0)
        return self


if __name__ == "__main__":
    q1 = Quaternion(1, 0, 0, 0)
    q2 = Quaternion(0, 1, 0, 0)
    q3 = Quaternion(0, 0, 1, 0)
    q4 = Quaternion(0, 0, 0, 1)

    print
    print np.all(q1 * q1 == -q4)
    print np.all(q1 * q2 == q3)
    print np.all(q1 * q3 == -q2)
    print np.all(q1 * q4 == q1)
    print
    print np.all(q2 * q1 == -q3)
    print np.all(q2 * q2 == -q4)
    print np.all(q2 * q3 == q1)
    print np.all(q2 * q4 == q2)
    print
    print np.all(q3 * q1 == q2)
    print np.all(q3 * q2 == -q1)
    print np.all(q3 * q3 == -q4)
    print np.all(q3 * q4 == q3)
    print
    print np.all(q4 * q1 == q1)
    print np.all(q4 * q2 == q2)
    print np.all(q4 * q3 == q3)
    print np.all(q4 * q4 == q4)

    r = Quaternion().from_axis_angle((0, 1, 0), np.pi/2)
    p = np.array([1, 0, 0])
    print np.allclose(r.apply(p), np.array([0,0,-1]))

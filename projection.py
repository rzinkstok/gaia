import numpy as np
from scipy.optimize import newton


def cartesian2spherical(x, y, z):
    # Theta is right ascension, phi is 90 - declination
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arccos(z/r)
    theta = np.arctan2(y, x)
    return r, theta, phi


def mollweide(longitude, latitude, central_meridian=0, radius=1):
    def f(l):
        return 2 * l + np.sin(2 * l) - np.pi * np.sin(latitude)

    def df(l):
        return 2 + 2 * np.cos(2 * l)

    theta = newton(f, latitude, df)

    x = radius * (2 * np.sqrt(2) / np.pi) * (longitude - central_meridian) * np.cos(theta)
    y = radius * np.sqrt(2) * np.sin(theta)
    return x, y


def aitoff(longitude, latitude, central_meridian=0):
    l = (longitude - central_meridian)/2.0
    p = latitude
    aa = np.arccos(np.cos(p) * np.cos(l))/np.pi
    a = np.sinc(aa)
    x = 2 * np.cos(p) * np.sin(l) / a
    y = np.sin(p) / a
    return x, y

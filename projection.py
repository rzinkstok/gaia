import numpy as np
from scipy.optimize import newton


def cartesian2spherical(x, y, z):
    # Theta is right ascension, phi is 90 - declination
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arccos(z/r)
    theta = np.arctan2(y, x)
    return r, theta, phi


def mollweide(longitude, latitude, central_meridian=0):
    R = 1

    def f(x):
        return 2 * x + np.sin(2*x) - np.pi * np.sin(latitude)

    def df(x):
        return 2 + 2 * np.cos(2*x)

    theta = newton(f, latitude, df)
    x = R * (2 * np.sqrt(2) / np.pi) * (longitude - central_meridian) * np.cos(theta)
    y = R * np.sqrt(2) * np.sin(theta)
    return (x, y)


def aitoff(longitude, latitude, central_meridian=0):
    l = (longitude - central_meridian)/2.0
    p = latitude
    aa = np.arccos(np.cos(p) * np.cos(l))/np.pi
    a = np.sinc(aa)
    x = 2 * np.cos(p) * np.sin(l) / a
    y = np.sin(p) / a
    return (x, y)
import numpy as np
import quaternion
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import datetime
import time

from projection import aitoff, cartesian2spherical
from sun import solar_apparent_longitude, days_since_epoch



PRECESSION_SPEED_CONSTANT = 4.223
SOLAR_ASPECT_ANGLE = np.deg2rad(45.0)
INERTIAL_SPIN_RATE = np.deg2rad(60.0/3600) * 86400 # rad per day


def scan_directions2(solar_longitude, nu, omega):
    N = solar_longitude.shape[0]

    # Build arrays with the quaternion components
    a1 = np.tile(np.array([0, 0, 1]), (N, 1)) * solar_longitude.reshape(N, 1)
    a2 = np.tile(np.array([1, 0, 0]), (N, 1)) * (np.pi/2.0 - nu).reshape(N, 1)
    a3 = np.tile(np.array([0, 1, 0]), (N, 1)) * (np.pi/2.0 - SOLAR_ASPECT_ANGLE)
    a4 = np.tile(np.array([0, 0, 1]), (N, 1)) * omega.reshape(N, 1)

    # Convert to quaternions
    q1 = quaternion.from_rotation_vector(a1)
    q2 = quaternion.from_rotation_vector(a2)
    q3 = quaternion.from_rotation_vector(a3)
    q4 = quaternion.from_rotation_vector(a4)

    # Compute composite rotation quaternions
    q12 = q1 * q2
    q123 = q12 * q3
    q = q123 * q4

    # Apply to initial vector (x axis)
    v = np.array([1, 0, 0])
    rv1 = quaternion.rotate_vectors(q123, v)
    rv2 = quaternion.rotate_vectors(q, v)

    # Convert to spherical coordinates
    r1, theta1, phi1 = cartesian2spherical(rv1)
    r2, theta2, phi2 = cartesian2spherical(rv2)

    return theta1, phi1, theta2, phi2


def calculate_scan_directions(ts, nus, omegas, filename):
    t1 = time.clock()
    solar_longitudes = solar_apparent_longitude(ts)[1]
    t2 = time.clock()
    print("Solar longitude calculation: {:.2f}".format(t2-t1))

    theta1, phi1, theta2, phi2 = scan_directions2(solar_longitudes, nus, omegas)
    t3 = time.clock()
    print("Scan direction calculation: {:.2f}s".format(t3 - t2))

    n = ts.shape[0]
    ts = ts.reshape(n, 1)
    nus = nus.reshape(n, 1)
    omegas = omegas.reshape(n, 1)
    solar_longitudes = solar_longitudes.reshape(n, 1)
    theta1 = theta1.reshape(n, 1)
    phi1 = phi1.reshape(n, 1)
    theta2 = theta2.reshape(n, 1)
    phi2 = phi2.reshape(n, 1)
    t4 = time.clock()
    print("Array juggling time: {:.2f}s".format(t4 - t3))

    if filename is not None:
        res = np.hstack((ts, solar_longitudes, nus, omegas, theta1, phi1, theta2, phi2))
        np.savetxt(filename, res, delimiter=",")
        t5 = time.clock()
        print("Write time: {:.2f}s".format(t5 - t4))
    return ts.reshape(n), theta2.reshape(n), phi2.reshape(n) - 0.5 * np.pi


def nsl_derivative(y, t):
    _, sl, dsl = solar_apparent_longitude(t)

    d1 = dsl * (np.sqrt(PRECESSION_SPEED_CONSTANT ** 2 - np.cos(y[0])**2) + np.cos(SOLAR_ASPECT_ANGLE) * np.sin(y[0])) / np.sin(SOLAR_ASPECT_ANGLE)
    d2 = INERTIAL_SPIN_RATE - np.cos(SOLAR_ASPECT_ANGLE) * d1 - np.sin(SOLAR_ASPECT_ANGLE) * np.sin(y[0]) * sl
    return np.array([d1, d2])


def calculate(filename, startdate):
    initial_precession_phase = 0.0
    initial_spin_phase = 0.0
    y = np.array([initial_precession_phase, initial_spin_phase])

    t0 = days_since_epoch(startdate)
    dt = 0.0001
    tmax = 5*365.25
    nt = int(round(tmax/dt))+1
    t = np.linspace(t0, tmax, nt)

    print("Solving ODEs")
    sol = odeint(nsl_derivative, y, t)

    print("Calculating scan directions")
    return calculate_scan_directions(t, sol[:, 0], sol[:, 1], filename)


def plot_from_file(filename):
    print("Reading file")
    t1 = time.time()
    with open(filename, "r") as fp:
        thetas = np.array([float(x.split(',')[6]) for x in fp])
    with open(filename, "r") as fp:
        phis = np.array([float(x.split(',')[7]) for x in fp]) - 0.5 * np.pi
    t2 = time.time()
    print("Reading time: {:.2f}s".format(t2-t1))
    plot(thetas, phis)


def plot(thetas, phis):
    plot_density_x, plot_density_y = aitoff(thetas, phis)

    print("Plotting histogram")
    plt.figure(figsize=(7.8, 4))
    plt.hist2d(plot_density_x, plot_density_y, (1200, 940), cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlim(-np.pi, np.pi)
    plt.ylim(-np.pi/2, np.pi/2)
    #ax = plt.axes().set_aspect(4/np.pi)
    plt.savefig("results.png")
    plt.show()



if __name__ == "__main__":
    filename = "results.csv"
    #ts, theta2, phi2 = calculate(filename, datetime.datetime(2000, 1, 1, 0, 0, 0))
    if filename:
        plot_from_file(filename)
    else:
        plot(theta2, phi2)

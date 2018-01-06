import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import datetime


from projection import aitoff, cartesian2spherical
from sun import solar_apparent_longitude
from quaternion import Quaternion


PRECESSION_SPEED_CONSTANT = 4.223
SOLAR_ASPECT_ANGLE = np.deg2rad(45.0)
INERTIAL_SPIN_RATE = np.deg2rad(60.0/3600) * 86400 # rad per day


def scan_directions(solar_longitude, nu, omega):
    # rotate around axis 3, angle solar_longitude
    q1 = Quaternion().from_axis_angle((0, 0, 1), solar_longitude)

    # rotate around axis 1, angle 90 - nu
    q2 = Quaternion().from_axis_angle((1, 0, 0), np.pi/2.0 - nu)

    # rotate around axis 2, angle 90 - SOLAR_ASPECT_ANGLE
    q3 = Quaternion().from_axis_angle((0, 1, 0), np.pi/2.0 - SOLAR_ASPECT_ANGLE)
    #q3 = Quaternion().from_axis_angle((0, -1, 0), SOLAR_ASPECT_ANGLE)

    # rotate around axis 3, angle omega
    q4 = Quaternion().from_axis_angle((0, 0, 1), omega)

    q12 = q1 * q2
    q123 = q12 * q3
    q = q123 * q4

    p1 = q123.apply(np.array([1, 0, 0]))
    p2 = q.apply(np.array([1, 0, 0]))

    # Calculate the precession direction
    r1, theta1, phi1 = cartesian2spherical(*p1)

    # Calculate the scan direction
    r2, theta2, phi2 = cartesian2spherical(*p2)

    return theta1, phi1, theta2, phi2


def calculate_scan_directions(ts, nus, omegas, start_datetime, filename):
    with open(filename, "w") as fp:
        for t, nu, omega in zip(ts, nus, omegas):
            tt = start_datetime + datetime.timedelta(days=t)
            solar_longitude = solar_apparent_longitude(tt)[1]
            theta1, phi1, theta2, phi2 = scan_directions(solar_longitude, nu, omega)
            fp.write("{},{},{},{},{},{},{},{}\n".format(t, solar_longitude, nu, omega, theta1, phi1, theta2, phi2))


def nsl_derivative(y, t, start_datetime):
    dt = start_datetime + datetime.timedelta(days=t)
    _, sl, dsl = solar_apparent_longitude(dt)

    d1 = dsl * (np.sqrt(PRECESSION_SPEED_CONSTANT ** 2 - np.cos(y[0])**2) + np.cos(SOLAR_ASPECT_ANGLE) * np.sin(y[0])) / np.sin(SOLAR_ASPECT_ANGLE)
    d2 = INERTIAL_SPIN_RATE - np.cos(SOLAR_ASPECT_ANGLE) * d1 - np.sin(SOLAR_ASPECT_ANGLE) * np.sin(y[0]) * sl
    return np.array([d1, d2])


def calculate(filename):
    initial_precession_phase = 0.0
    initial_spin_phase = 0.0
    y = np.array([initial_precession_phase, initial_spin_phase])

    t0 = 0.0
    dt = 0.0001
    tmax = 5*365.25
    nt = int(round(tmax/dt))+1
    t = np.linspace(t0, tmax, nt)
    start_dt = datetime.datetime(2000, 1, 1, 0, 0, 0)

    print "Solving ODEs"
    sol = odeint(nsl_derivative, y, t, args=(start_dt,))

    print "Calculating scan directions"
    calculate_scan_directions(t, sol[:, 0], sol[:, 1], start_dt, filename)


def plot(filename):
    print "Plotting scan directions"
    plot_density_x = []
    plot_density_y = []

    print "Reading file"
    with open(filename, "r") as fp:
        lines = fp.readlines()
        nlines = len(lines)
        for i, l in enumerate(lines):
            if i%100000 == 0:
                print "Progress: {:.2f}%".format(100.0*float(i)/nlines)
            parts = [float(x) for x in l.split(",")]
            theta, phi = parts[6:8]
            phi -= 0.5 * np.pi
            mx, my = aitoff(theta, phi)
            plot_density_x.append(mx)
            plot_density_y.append(my)

    print "Plotting histogram"
    plt.figure(figsize=(7.8, 4))
    plt.hist2d(plot_density_x, plot_density_y, (1200, 940), cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlim(-np.pi, np.pi)
    plt.ylim(-np.pi/2, np.pi/2)
    plt.axes().set_aspect(4/np.pi)
    plt.savefig("results.png")
    plt.show()



if __name__ == "__main__":
    filename = "results.csv"
    #calculate(filename)
    plot(filename)
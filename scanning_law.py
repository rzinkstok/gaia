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


def scan_direction(solar_longitude, nu, omega):
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

    p = q.apply(np.array([1, 0, 0]))
    r, theta, phi = cartesian2spherical(*p)
    return theta, phi


def calculate_scan_directions(ts, nus, omegas, start_datetime):
    solar_longitudes = [solar_apparent_longitude(start_datetime + datetime.timedelta(days=t))[1] for t in ts]
    return [scan_direction(solar_longitude, nu, omega) for solar_longitude, nu, omega in zip(solar_longitudes, nus, omegas)]


def density_plot_aitoff(spherical_points):
    plot_density_x = []
    plot_density_y = []
    for theta, phi in spherical_points:
        phi -= 0.5 * np.pi
        mx, my = aitoff(theta, phi)
        plot_density_x.append(mx)
        plot_density_y.append(my)
    plt.hist2d(plot_density_x, plot_density_y, (1200, 600), cmap=plt.cm.jet)
    plt.colorbar()
    #plt.axes().set_aspect('equal', 'datalim')

    plt.show()


def nsl_derivative(y, t, start_datetime):
    dt = start_datetime + datetime.timedelta(days=t)
    _, sl, dsl = solar_apparent_longitude(dt)

    d1 = dsl * (np.sqrt(PRECESSION_SPEED_CONSTANT ** 2 - np.cos(y[0])**2) + np.cos(SOLAR_ASPECT_ANGLE) * np.sin(y[0])) / np.sin(SOLAR_ASPECT_ANGLE)
    d2 = INERTIAL_SPIN_RATE - np.cos(SOLAR_ASPECT_ANGLE) * d1 - np.sin(SOLAR_ASPECT_ANGLE) * np.sin(y[0]) * sl
    return np.array([d1, d2])


def run():
    initial_precession_phase = 0.0
    initial_spin_phase = 0.0
    y = np.array([initial_precession_phase, initial_spin_phase])

    t0 = 0.0
    dt = 0.0001
    tmax = 5*365.25
    nt = int(round(tmax/dt))+1
    print nt
    t = np.linspace(t0, tmax, nt)
    start_dt = datetime.datetime(2000, 1, 1, 0, 0, 0)

    print "Solving ODEs"
    sol = odeint(nsl_derivative, y, t, args=(start_dt,))

    print "Calculating scan directions"
    scan_directions = calculate_scan_directions(t, sol[:, 0], sol[:, 1], start_dt)

    print "Plotting scan directions"
    density_plot_aitoff(scan_directions)


if __name__ == "__main__":
    run()
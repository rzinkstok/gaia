import numpy as np
import datetime
import time
import sunpy.sun as sun


EPOCH_J2000 = datetime.datetime(2000, 1, 1, 12, 0, 0)


def days_since_epoch(d, epoch=EPOCH_J2000):
    delta = d-epoch
    return delta.days + delta.seconds/86400.0


def solar_apparent_longitude(t):
    # From Lindgren SAG-LL-35, p. 10
    # t is days since epoch J2000.0

    e = 0.016709
    a0 = np.deg2rad(280.458)
    a1 = np.deg2rad(0.98560911)
    g0 = np.deg2rad(357.528)
    g1 = np.deg2rad(0.98560028)

    a = a0 + a1 * t
    g = g0 + g1 * t

    solar_longitude = a + 2 * e * np.sin(g) + 1.25 * e**2 * np.sin(2 * g)
    solar_longitude = np.divmod(solar_longitude, 2 * np.pi)[1]
    derivative_solar_longitude = a1 + (2 * e * np.cos(g) + 2 * 1.25 * e**2 * np.cos(2 * g)) * g1
    solar_distance = 1 - e * np.cos(g + e * np.sin(g) + 0.5 * e**2 * np.sin(2 * g))

    return solar_distance, solar_longitude, derivative_solar_longitude


def dms2rad(dms):
    return  np.deg2rad(dms[0] + dms[1] / 60.0 + dms[2] / (60 * 60.0))


def solar_apparent_longitude2(dt):
    """Use sunpy for the apparent longitude; needs a datetime date."""

    deltat = datetime.timedelta(days=1)
    l = dms2rad(sun.apparent_longitude(dt).dms)
    l1 = dms2rad(sun.apparent_longitude(dt - deltat).dms)
    l2 = dms2rad(sun.apparent_longitude(dt + deltat).dms)
    dl = 0.5 * (l2 - l1)
    return None, l, dl


if __name__ == "__main__":
    dt = 0.0001
    dt = 0.0001
    tmax = 5 * 365.25
    nt = int(round(tmax / dt)) + 1
    t = np.linspace(0, tmax, nt)
    t1 = time.clock()
    solar_apparent_longitude(t)
    print time.clock() - t1


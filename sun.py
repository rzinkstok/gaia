import numpy as np
import datetime


EPOCH_J2000 = datetime.datetime(2000, 1, 1, 12, 0, 0)


def days_since_epoch(d, epoch):
    delta = d-epoch
    return delta.days + delta.seconds/86400.0


def solar_apparent_longitude(dt):
    # From Lindgren SAG-LL-35, p. 10
    t = days_since_epoch(dt, EPOCH_J2000)

    e = 0.016709
    a0 = np.deg2rad(280.458)
    a1 = np.deg2rad(0.98560911)
    g0 = np.deg2rad(357.528)
    g1 = np.deg2rad(0.98560028)

    a = a0 + a1 * t
    g = g0 + g1 * t

    solar_longitude = a + 2 * e * np.sin(g) + 1.25 * e**2 * np.sin(2 * g)
    solar_longitude = divmod(solar_longitude, 2 * np.pi)[1]
    derivative_solar_longitude = a1 + (2 * e * np.cos(g) + 2 * 1.25 * e**2 * np.cos(2 * g)) * g1
    solar_distance = 1 - e * np.cos(g + e * np.sin(g) + 0.5 * e**2 * np.sin(2 * g))

    return solar_distance, solar_longitude, derivative_solar_longitude


if __name__ == "__main__":
    import sunpy.sun as sun
    d = datetime.datetime.now() + datetime.timedelta(days=0)
    print solar_apparent_longitude(d)[1]
    l = sun.apparent_longitude(d).dms
    print np.deg2rad(l[0] + l[1]/60.0 + l[2]/(60*60.0))

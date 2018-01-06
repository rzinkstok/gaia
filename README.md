# gaia

Simple simulations of ESA's Gaia star mapping satellite.

## Nominal Scanning Law 

The Gaia satellite orbits the L<sub>2</sub> Lagrange point of Earth. It scans the sky using a clever combination of three rotations:
- a rotation Ω (spin) around an axis perpendicular to the lines of sight of its two telescopes, with a period of 6 hours;
- a precession ν of this rotation axis around the line between the spacecraft and the Sun, with a period of 63 days;
- the yearly motion of the Sun as seen from Earth (ecliptic latitude β<sub>s</sub> = 0, ecliptic latitude λ<sub>s</sub>).

The combination of these rotations yields a full coverage of the sky: this is called the GAIA nominal scanning law.

### Calculation
Lennart Lindegren's technical note 
[Calculating the GAIA Nominal Scanning law](http://www.astro.lu.se/~lennart/Astrometry/TN/Gaia-LL-035-20010219-Calculating-the-GAIA-Nominal-Scanning-Law.pdf)
describes a method to crunch the numbers. It involves the integration of two coupled differential equations that relate 
the spin phase Ω(t) and the precession phase ν(t) with their derivatives. The nominal longitude of the Sun λ<sub>s</sub>(t) 
is considered a given, and can be approximated using the first terms of a series expansion of Kepler's equation. 

So given the apparent solar longitude λ<sub>s</sub>, the angle between the direction to the Sun and the precession axis ξ,
the precession rate S, the spin rate ω, and the initial values Ω<sub>0</sub> and ν<sub>0</sub>, the scanning law is completely
determined.

The values of Ω(t) and ν(t) can be used to compute the direction of the line of sight; this is done quite easily using quaternions.


### Results
When the instantaneous scan directions are projected using an Aitoff projection, and binned in 1200 by 940 bins, this is 
the resulting density map (using ecliptic coordinates):
![GAIA Nominal Scanning Law: sky coverage of Gaia](results.png?raw=true "The GAIA Nominal Scanning Law: sky coverage of Gaia")
Figure 1: The GAIA Nominal Scanning Law: sky coverage of Gaia.

### Limitations
Formally, the data should be binned before projecting, so this is only an approximation. Furthermore, the field of view of the
imaging telescopes is not taken into account. And of course the solar longitude is only an approximation.
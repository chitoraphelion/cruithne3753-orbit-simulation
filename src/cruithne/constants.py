"""Physical constants, masses (SI), and conversion factors."""

from __future__ import annotations

G = 6.67430e-11  # m^3 kg^-1 s^-2

M_SUN = 1.98847e30
M_MERCURY = 3.3011e23
M_VENUS = 4.8675e24
M_EARTH = 5.97219e24
M_MOON = 7.342e22
M_MARS = 6.4171e23
M_JUPITER = 1.89813e27
M_SATURN = 5.6834e26
M_URANUS = 8.6810e25
M_NEPTUNE = 1.02413e26

GM_SUN = G * M_SUN
GM_MERCURY = G * M_MERCURY
GM_VENUS = G * M_VENUS
GM_EMB = G * (M_EARTH + M_MOON)
GM_MARS = G * M_MARS
GM_JUPITER = G * M_JUPITER
GM_SATURN = G * M_SATURN
GM_URANUS = G * M_URANUS
GM_NEPTUNE = G * M_NEPTUNE

AU = 1.495978707e11  # metres
DAY = 86400.0  # seconds
YEAR = 365.25 * DAY

PLANETS = (
    ("Sun", "10", "majorbody", M_SUN),
    ("Mercury", "1", "majorbody", M_MERCURY),
    ("Venus", "2", "majorbody", M_VENUS),
    ("EMB", "3", "majorbody", M_EARTH + M_MOON),
    ("Mars", "4", "majorbody", M_MARS),
    ("Jupiter", "5", "majorbody", M_JUPITER),
    ("Saturn", "6", "majorbody", M_SATURN),
    ("Uranus", "7", "majorbody", M_URANUS),
    ("Neptune", "8", "majorbody", M_NEPTUNE),
)
CRUITHNE_ID = ("3753", "smallbody")

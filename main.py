import numpy as np
import Measurements
from pathlib import Path


def s_ac(u_a, u_s, f):
    w = 2 * np.pi * f
    # w = f
    v = 3e11 / 1.4

    temp = (3*u_s/2) * (-u_a + u_a*np.sqrt(1+(w**2/(v**2*u_a**2))))
    # temp *=

    return -temp**0.5


def s_ph(u_a, u_s, f):
    w = 2 * np.pi * f
    # w = f
    v = 3e11 / 1.4

    temp = (3 * u_s / 2) * (u_a + u_a*np.sqrt(1+(w**2/(v**2*u_a**2))))
    # temp = temp

    return temp**0.5


measurementPath = Path("")
measurement1 = Measurements.DualSlopePhantomMeasurement(location=measurementPath, common='detector')
measurement2 = Measurements.DualSlopePhantomMeasurement(location=measurementPath, common='source')




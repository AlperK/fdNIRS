import numpy as np
import Measurements
import fdNIRS
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


measurementPath = Path("2023-11-15", "DUAL-SLOPE-BOTH-3", "2")
measurement1 = Measurements.DualSlopePhantomMeasurement(location=measurementPath,
                                                        common='detector')
measurement2 = Measurements.DualSlopePhantomMeasurement(location=measurementPath,
                                                        common='source')

# print(measurement1.phase_slopes[:, 1].mean(axis=1)[23])
print('------------------------------')
print(fdNIRS.compute_optical_parameters(measurement1.dual_amplitude_slopes_color1,
                                        measurement1.dual_phase_slopes_color1,
                                        measurement1.modulation_frequency)[0])
print('------------------------------')
print(fdNIRS.compute_optical_parameters(measurement1.dual_amplitude_slopes_color2,
                                        measurement1.dual_phase_slopes_color2,
                                        measurement1.modulation_frequency)[0])
absorption_coefficient_830 = fdNIRS.compute_optical_parameters(
    measurement1.dual_amplitude_slopes_color1,
    measurement1.dual_phase_slopes_color1,
    measurement1.modulation_frequency)[0].mean(axis=1)
absorption_coefficient_685 = fdNIRS.compute_optical_parameters(
    measurement1.dual_amplitude_slopes_color2,
    measurement1.dual_phase_slopes_color2,
    measurement1.modulation_frequency)[0].mean(axis=1)
print('------------------------------')
# print(fdNIRS.compute_hemoglobins(absorption_coefficient_830,
#                                  absorption_coefficient_685))
# print('------------------------------')
print(fdNIRS.compute_hemoglobin_concentrations(absorption_coefficient_830,
                                absorption_coefficient_685).mean(axis=1))
print('------------------------------')

import numpy as np
import Measurements
import fdNIRS
import matplotlib.pyplot as plt
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


measurementPath = Path("2023-11-20", "ARTERIAL", "7")
measurement1 = Measurements.DualSlopePhantomMeasurement(location=measurementPath,
                                                        common='detector')
measurement2 = Measurements.DualSlopePhantomMeasurement(location=measurementPath,
                                                        common='source')
timestamps = np.loadtxt(Path(measurementPath, 'times.txt'))
timestamps = timestamps - timestamps[0]
# print(measurement1.phase_slopes[:, 1].mean(axis=1)[23])
print('------------------------------')
# print(fdNIRS.compute_optical_parameters(measurement1.dual_amplitude_slopes_color1,
#                                         measurement1.dual_phase_slopes_color1,
#                                         measurement1.modulation_frequency)[0])
# print('------------------------------')
# print(fdNIRS.compute_optical_parameters(measurement1.dual_amplitude_slopes_color2,
#                                         measurement1.dual_phase_slopes_color2,
#                                         measurement1.modulation_frequency)[0])
absorption_coefficient_830 = fdNIRS.compute_optical_parameters(
    measurement1.dual_amplitude_slopes_color1,
    measurement1.dual_phase_slopes_color1,
    measurement1.modulation_frequency)[0].mean(axis=1)
absorption_coefficient_685 = -fdNIRS.compute_optical_parameters(
    measurement1.dual_amplitude_slopes_color2,
    measurement1.dual_phase_slopes_color2,
    measurement1.modulation_frequency)[0].mean(axis=1)
print('------------------------------')
# print(fdNIRS.compute_hemoglobins(absorption_coefficient_830,
#                                  absorption_coefficient_685))
# print('------------------------------')
s = absorption_coefficient_830.size // 10

absorption_coefficient_830 = fdNIRS.rolling_apply(np.mean, a=absorption_coefficient_830, w=15)
absorption_coefficient_685 = fdNIRS.rolling_apply(np.mean, a=absorption_coefficient_685, w=15)
HbO, HbR = fdNIRS.compute_hemoglobin_concentrations(absorption_coefficient_830,
                                                    absorption_coefficient_685)
# print(fdNIRS.compute_hemoglobin_concentrations(absorption_coefficient_830,
#                                                absorption_coefficient_685).mean(axis=1))
print('------------------------------')
print(HbO.size)
t = np.linspace(0, 300, HbO.size)
plt.figure()
plt.suptitle('Absorption coefficients')
plt.plot(t, absorption_coefficient_830, label='830nm', linewidth=3, alpha=0.75)
plt.plot(t, absorption_coefficient_685, label='685nm', linewidth=3, alpha=0.75)
plt.xlabel('Time(s)')
plt.ylabel(r'$mm^{-1}$')
# plt.tick_params(axis='x', bottom=False, labelbottom=False)
# plt.axvline(s*4)
# plt.axvline(s*7)
plt.axvspan(120, 180, alpha=0.15, color='blue', label='occlusion')
plt.legend()

plt.figure()
plt.suptitle('Hemoglobin concentrations')
plt.plot(t, 1e6 * fdNIRS.rolling_apply(np.mean, a=HbO, w=1), label='HbO', linewidth=3, alpha=0.75)
plt.plot(t, 1e6 * fdNIRS.rolling_apply(np.mean, a=HbR, w=1), label='HbR', linewidth=3, alpha=0.75)
plt.plot(t, 1e6 * fdNIRS.rolling_apply(np.mean, a=HbR+HbO, w=1), label='HbT', linewidth=3, alpha=0.75)
plt.xlabel('Time(s)')
plt.ylabel(r'$\mu$M $mm^{-2}$')
# plt.axvline(s*4)
# plt.axvline(s*7)
plt.axvspan(120, 180, alpha=0.15, color='blue', label='occlusion')
plt.legend()

plt.show()

import numpy as np

import fdNIRS
import matplotlib.pyplot as plt
from pathlib import Path


location = Path('2024-03-19', 'DUAL-SLOPE-BOTH', '1')
measurement = fdNIRS.DualSlopeMeasurement(location=location)
measurement.compute_hemoglobin_concentrations()
windowSize = 3
occlusionInterval = None
totalTime = 60*60

t = np.linspace(0, totalTime, measurement.oxy_hemoglobin_concentration.size)
measurement.plot_occlusion('Arterial',
                           total_time=totalTime,
                           occlusion_interval=occlusionInterval,
                           window_size=windowSize)
measurement.plot_raw(total_time=totalTime, occlusion_interval=occlusionInterval,
                     window_size=windowSize)
measurement.plot_slopes(total_time=totalTime,
                        occlusion_interval=occlusionInterval,
                        window_size=windowSize)
measurement.plot_absorption(total_time=totalTime,
                            occlusion_interval=occlusionInterval,
                            window_size=windowSize)
measurement.plot_scattering(total_time=totalTime,
                            occlusion_interval=occlusionInterval,
                            window_size=windowSize)

print(f'830nm ua: {measurement.absorption_coefficient_wavelength_1.mean()}, '
      f'error: {measurement.absorption_coefficient_wavelength_1.mean() / 0.0077 - 1}')
print(f'830nm us: {measurement.scattering_coefficient_wavelength_1.mean()}, '
      f'error: {measurement.scattering_coefficient_wavelength_1.mean() / 0.597 - 1}')

print(f'690nm ua: {measurement.absorption_coefficient_wavelength_2.mean()}, '
      f'error: {measurement.absorption_coefficient_wavelength_2.mean() / 0.0081 - 1}')
print(f'690nm us: {measurement.scattering_coefficient_wavelength_2.mean()}, '
      f'error: {measurement.scattering_coefficient_wavelength_2.mean() / 0.761 - 1}')

plt.show()

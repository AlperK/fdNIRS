import numpy as np

import fdNIRS
import matplotlib.pyplot as plt
from pathlib import Path


location = Path('2023-12-20', 'AO-5-3-2', '9')
measurement = fdNIRS.DualSlopeMeasurement(location=location)
measurement.compute_hemoglobin_concentrations()


t = np.linspace(0, 10*60, measurement.oxy_hemoglobin_concentration.size)
measurement.plot_occlusion('Arterial',
                           total_time=600,
                           occlusion_interval=(300, 480),
                           window_size=15)
measurement.plot_raw(total_time=600, occlusion_interval=(300, 480), window_size=10)
plt.show()

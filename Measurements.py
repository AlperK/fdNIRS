from Phantom import Phantom
import json
import numpy as np
from pathlib import Path


class FDNIRSMeasurement:

    def __init__(self, location):
        self.common = None
        self.measurement_location = location
        self.measurement_details = None
        self.measurement_type = None
        self._read_measurement_details()

    def _read_measurement_details(self):
        with open(Path.joinpath(self.measurement_location, "measurement settings.json")) as f:
            self.measurement_details = json.load(f)
        self.measurement_type = self.measurement_details["Measurement type"]
        self.measurement_method = self.measurement_details["Measurement method"]
        self.modulation_frequency = self.measurement_details["RF"] * 1e6 + self.measurement_details["IF"] * 1e3
        self.separations = np.array([self.measurement_details["Separations"]]).reshape(-1, 2)

        if self.measurement_type == "Phantom":
            self.phantom = Phantom(self.measurement_details["Phantom"])
        elif self.measurement_type == "Human":
            pass


class PhantomMeasurement(FDNIRSMeasurement):
    def __init__(self, location):
        super().__init__(location)


class DualSlopePhantomMeasurement(PhantomMeasurement):
    def __init__(self, location, common):
        super().__init__(location)

        self.common = common
        self.data = FDNIRSData(measurement=self)
        self.phase_slopes = self.get_phase_slopes()

    def load_measurement_data(self):
        """
        Loads the measurement data from amplitude.csv and phase.csv
        :return: Reshapes amplitudes and phases into (n, 2, 4) and returns
        """
        amplitudes = np.loadtxt(Path(self.measurement_location, "amplitude.csv"), delimiter=',')
        phases = np.loadtxt(Path(self.measurement_location, "phase.csv"), delimiter=',')

        return amplitudes.reshape((amplitudes.shape[0], 2, 4)), phases.reshape((phases.shape[0], 2, 4))

    def get_phase_slopes(self):
        return np.divide(np.diff(self.data.phases),
                         np.diff(self.separations))

    def _get_linearized_amplitudes(self):
        squared = (self.data.amplitudes *
                   self.separations**2)
        return np.log(squared)

    def get_amplitude_slopes(self):
        linearized_amplitudes = self._get_linearized_amplitudes()
        return np.divide(np.diff(linearized_amplitudes),
                         np.diff(self.separations))


class FDNIRSData:
    def __init__(self, measurement: FDNIRSMeasurement):
        self.measurement = measurement
        self.common = measurement.common
        self.amplitudes, self.phases = self.load_measurement_data()

        self.amplitudes_of_color1 = self.amplitudes[:, 0, :]
        self.phases_of_color1 = self.phases[:, 0, :, :]
        self.amplitudes_of_color2 = self.amplitudes[:, 1, :]
        self.phases_of_color2 = self.phases[:, 1, :, :]

    def load_measurement_data(self):
        """
        Loads the measurement data from amplitude.csv and phase.csv
        :return: Reshapes amplitudes and phases into (n, 2, 4) and returns
        """
        amplitudes = np.loadtxt(Path(self.measurement.measurement_location, "amplitude.csv"), delimiter=',')
        phases = np.loadtxt(Path(self.measurement.measurement_location, "phase.csv"), delimiter=',')

        amplitudes = self.reshape_data_into_pairs(amplitudes)
        phases = self.reshape_data_into_pairs(phases)
        return amplitudes, phases

    def reshape_data_into_pairs(self, data):
        if self.common == 'detector':
            return data.reshape(data.shape[0], 2, 2, 2)
        elif self.common == 'source':
            return data.reshape(data.shape[0], 2, 2, 2).swapaxes(2, 3)
            pass
            # TODO: Implement here

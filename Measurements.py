from Phantom import Phantom
import json
import numpy as np
from pathlib import Path

def asdasd:
    pass

class fdNIRSMeasurement:

    def __init__(self, location):
        self.measurement_location = location
        self.measurement_details = None
        self.measurement_type = None
        self._read_measurement_details()

    def _read_measurement_details(self):
        with open(Path.joinpath(self.measurement_location, "measurement settings.json")) as f:
            self.measurement_details = json.load(f)
        self.measurement_type = self.measurement_details["Measurement type"]
        self.measurement_method = self.measurement_details["Measurement method"]
        self.separations = self.measurement_details["Separations"]
        self.modulation_frequency = self.measurement_details["RF"] * 1e6 + self.measurement_details["IF"] * 1e3

        if self.measurement_type == "Phantom":
            self.phantom = Phantom(self.measurement_details["Phantom"])
        elif self.measurement_type == "Human":
            pass


class PhantomMeasurement(fdNIRSMeasurement):
    def __init__(self, location):
        super().__init__(location)


class DualSlopePhantomMeasurement(PhantomMeasurement):
    def __init__(self, location):
        super().__init__(location)

        self.data = fdNIRSData(measurement=self)
        self.amplitudes, self.phases = self.load_measurement_data()

        self.amplitudes_830 = self.data.amplitudes_830
        self.amplitudes_685 = self.data.amplitudes_685
        self.phases_830 = self.data.phases_830
        self.phases_685 = self.data.phases_685

    def load_measurement_data(self):
        """
        Loads the measurement data from amplitude.csv and phase.csv
        :return: Reshapes amplitudes and phases into (n, 2, 4) and returns
        """
        amplitudes = np.loadtxt(Path(self.measurement_location, "amplitude.csv"), delimiter=',')
        phases = np.loadtxt(Path(self.measurement_location, "phase.csv"), delimiter=',')

        return amplitudes.reshape((amplitudes.shape[0], 2, 4)), phases.reshape((phases.shape[0], 2, 4))


class fdNIRSData:
    def __init__(self, measurement: fdNIRSMeasurement):
        self.measurement = measurement

        self.amplitudes, self.phases = self.load_measurement_data()
        self.amplitudes_685 = self.amplitudes[:, 1, :]
        self.phases_685 = self.phases[:, 1, :]
        self.amplitudes_830 = self.amplitudes[:, 0, :]
        self.phases_830 = self.phases[:, 0, :]

    def load_measurement_data(self):
        """
        Loads the measurement data from amplitude.csv and phase.csv
        :return: Reshapes amplitudes and phases into (n, 2, 4) and returns
        """
        amplitudes = np.loadtxt(Path(self.measurement.measurement_location, "amplitude.csv"), delimiter=',')
        phases = np.loadtxt(Path(self.measurement.measurement_location, "phase.csv"), delimiter=',')

        return amplitudes.reshape((amplitudes.shape[0], 2, 4)), phases.reshape((phases.shape[0], 2, 4))

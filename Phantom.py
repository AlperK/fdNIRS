import numpy as np
import fdNIRS


class Phantom:
    def __init__(self, name):
        self.name = name
        self.absorption_685 = None
        self.scattering_685 = None
        self.absorption_830 = None
        self.scattering_830 = None

        self.Sac_685 = None
        self.Sph_685 = None
        self.Sac_830 = None
        self.Sph_830 = None

        self._assign_coefficients()

    def __repr__(self):
        return f'{self.name}'

    def _assign_coefficients(self):
        try:
            # assert
            if self.name.lower().startswith('fantini'):
                try:
                    assert int(self.name[-1]) in [1, 2, 3]
                    if self.name[-1] == "1":
                        self.absorption_685 = 0.0015
                        self.scattering_685 = 1.023
                        self.absorption_830 = 0.0013
                        self.scattering_830 = 0.785
                    elif self.name[-1] == "2":
                        self.absorption_685 = 0.0124
                        self.scattering_685 = 0.930
                        self.absorption_830 = 0.0106
                        self.scattering_830 = 0.763
                    elif self.name[-1] == "3":
                        self.absorption_685 = 0.0081
                        self.scattering_685 = 0.761
                        self.absorption_830 = 0.0077
                        self.scattering_830 = 0.597
                except AssertionError:
                    print('Phantom number was not recognized.\n'
                          'The phantom optical properties could not be set.')
            elif self.name.lower().startswith('pionirs'):
                try:
                    assert int(self.name[-1]) in [1, 2, 3]

                    if self.name[-1] == "1":
                        self.absorption_685 = 0.01
                        self.scattering_685 = 1.000
                        self.absorption_830 = 0.01
                        self.scattering_830 = 1.000
                except AssertionError:
                    print('Phantom number was not recognized.\n'
                          'The phantom optical properties could not be set.')
        except ValueError:
            print(f"The {self} Phantom did not found!")

    def _assign_slopes(self, slopes, wavelength):
        setattr(self, f"Sac_{wavelength}", slopes[0])
        setattr(self, f"Sph_{wavelength}", slopes[1])

    def get_slopes(self, frequency, wavelength):

        slopes = np.round(fdNIRS.get_slopes(absorption_coefficient=getattr(self, f"absorption_{wavelength}"),
                                            scattering_coefficient=getattr(self, f"scattering_{wavelength}"),
                                            frequency=frequency), 4)
        self._assign_slopes(slopes, wavelength)

        return slopes

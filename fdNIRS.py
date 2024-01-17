#TODO Phantom measurement
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from scipy.fft import fft, fftfreq
from scipy.signal.windows import blackman
from scipy import signal
from pathlib import Path


class fdNIRSData:
    def __init__(self, location):
        self.amplitude_slopes = None
        self.phase_slopes = None
        self.location = location
        self.amplitudes, self.phases = self.read_data_from_txt()

    def read_data_from_txt(self):
        amplitudes = np.loadtxt(Path(self.location, 'amplitude.csv'), delimiter=',')
        phases = np.loadtxt(Path(self.location, 'phase.csv'), delimiter=',')

        return amplitudes, phases


def absorption_coefficient_equation(s_ac, s_ph, f):
    """
    Returns the absorption coefficient from the given slopes and modulation frequency
    :param s_ac: Intensity slope (1/mm)
    :param s_ph: Phase slope (1/mm)
    :param f: Modulation frequency (Hz)
    :return: Absorption coefficient (1/mm)
    """
    w = f * 2 * np.pi
    v = 3e11 / 1.4

    return (s_ph / s_ac - s_ac / s_ph) * (w / (2 * v))


def scattering_coefficient_equation(s_ac, s_ph, absorption_coefficient):
    """
    Returns the scattering coefficient from the given slopes and modulation frequency
    :param s_ac: Intensity slope (1/mm)
    :param s_ph: Phase slope (1/mm)
    :param absorption_coefficient: Absorption coefficient (1/mm)
    :return: Scattering coefficient (1/mm)
    """
    return (s_ac ** 2 - s_ph ** 2) / (3 * absorption_coefficient)


def _slope_equations(s, *data):
    """
    Helper function to compute the slopes of a phantom from the known absorption and scattering coefficients
    at a certain modulation frequency
    :param s:
    :param data: Arguments for the absorption_coefficient_equation() and scattering_coefficient_equation().
    Should be 3 arguments exactly;
    [modulation frequency (Hz), absorption coefficient (1/mm), scattering coefficient (1/mm)]
    :return: The error between the
    """
    f = data[0]
    absorption = data[1]
    scattering = data[2]

    s_ac, s_ph = s[0], s[1]
    eq1 = absorption_coefficient_equation(s_ac, s_ph, f)
    eq2 = scattering_coefficient_equation(s_ac, s_ph, absorption)

    return eq1 - absorption, eq2 - scattering


def get_slopes_from_optical_parameters(absorption_coefficient, scattering_coefficient, frequency):
    slopes = fsolve(_slope_equations,
                    np.array([-0.1, 0.1]),
                    args=(frequency, absorption_coefficient, scattering_coefficient))

    return slopes


def compute_optical_parameters(
        amplitude_slopes,
        phase_slopes,
        modulation_frequency):
    w = 2 * np.pi * modulation_frequency
    n = 1.4
    c = 2.998e11

    absorption_coefficient = (
            (w / (2 * (c / n))) *
            (
                    np.divide(phase_slopes, amplitude_slopes) -
                    np.divide(amplitude_slopes, phase_slopes)
            )
    )

    scattering_coefficient = (
            (np.square(amplitude_slopes) - np.square(phase_slopes)) /
            (3 * absorption_coefficient)
    )
    return np.array([absorption_coefficient, scattering_coefficient])


def compute_hemoglobin_concentrations(
        absorption_coefficient_color_830,
        absorption_coefficient_color_685, ):
    oxy_hemoglobin_extinction_coefficient_830 = 974.0 / 10
    oxy_hemoglobin_extinction_coefficient_685 = 272.8 / 10

    deoxy_hemoglobin_extinction_coefficient_830 = 693.04 / 10
    deoxy_hemoglobin_extinction_coefficient_685 = 2188.24 / 10

    oxy_hemoglobin_concentration = (
            ((absorption_coefficient_color_830 * deoxy_hemoglobin_extinction_coefficient_685) -
             (absorption_coefficient_color_685 * deoxy_hemoglobin_extinction_coefficient_830)) /
            ((deoxy_hemoglobin_extinction_coefficient_685 * oxy_hemoglobin_extinction_coefficient_830) -
             (deoxy_hemoglobin_extinction_coefficient_830 * oxy_hemoglobin_extinction_coefficient_685))
    )

    deoxy_hemoglobin_concentration = (
            ((absorption_coefficient_color_685 * oxy_hemoglobin_extinction_coefficient_830) -
             (absorption_coefficient_color_830 * oxy_hemoglobin_extinction_coefficient_685)) /
            ((oxy_hemoglobin_extinction_coefficient_830 * deoxy_hemoglobin_extinction_coefficient_685) -
             (oxy_hemoglobin_extinction_coefficient_685 * deoxy_hemoglobin_extinction_coefficient_830))
    )

    return np.array([oxy_hemoglobin_concentration, deoxy_hemoglobin_concentration])


def rolling_apply(fun, a, w):
    if not w is None:
        r = np.empty(a.shape)
        r.fill(np.nan)

        for i in range(w - 1, a.shape[0]):
            if fun is not None:
                r[i] = fun(a[(i - w + 1):i + 1])
            else:
                r[i] = a[(i - w + 1):i + 1]
        return r
    else:
        return a


def apply_kalman_1d(initial_estimate,
                    initial_error_estimate,
                    error_in_measurement,
                    meas):
    initial_kalman_gain = initial_error_estimate / (initial_error_estimate + error_in_measurement)
    previous_estimate = initial_estimate
    previous_error_estimate = initial_error_estimate
    kalman_gain = initial_kalman_gain

    estimates = np.array([])

    for mea in meas:
        if np.isnan(mea) == False:
            # print(mea)
            current_estimate = previous_estimate + kalman_gain * (mea - previous_estimate)
            estimates = np.append(estimates, [current_estimate])
            # print(f'currentEstimate = {current_estimate}')
            # print(estimates)
            current_error_estimate = (1 - kalman_gain) * previous_error_estimate
            # print(f'currentErrorEstimate = {current_error_estimate}')

            kalman_gain = current_error_estimate / (current_error_estimate + error_in_measurement)
            # print(f'kalmanGain = {kalman_gain}')
            previous_estimate = current_estimate

    return estimates


def plot_fft(y):
    y = y - np.mean(y)
    N = np.size(y)
    T = 600 / N
    w = blackman(N)
    yf = fft(y*w)
    xf = fftfreq(N, T)[:N//2]

    plt.figure()
    plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
    plt.grid()


def apply_butterworth(y):
    plot_fft(y)
    m = np.mean(y)
    N = np.size(y)
    T = N / 600

    sos = signal.butter(14, [0.01, 0.02], fs=N/600, btype='bandstop', output='sos')
    filtered = signal.sosfilt(sos, y-m)
    # sos = signal.butter(14, [0.2, 0.25], fs=N/600, btype='bandstop', output='sos')
    # filtered = signal.sosfilt(sos, filtered)
    plot_fft(filtered)

    b, a = signal.butter(3, 0.015)
    filtered = signal.filtfilt(b, a, y-m)
    return filtered + m


class fdNIRS:

    def __init__(self, location):
        self.total_hemoglobin_concentration = None
        self.oxy_hemoglobin_concentration = None
        self.deoxy_hemoglobin_concentration = None
        self.location = location
        self.amplitudes, self.phases = self.read_data_from_txt()
        self.settings = self.read_measurement_settings()
        self.separations = np.array(self.settings['Separations'])
        self.modulation_frequency = float(self.settings['RF']) * 1e6

        self.absorption_coefficient_wavelength_1 = None
        self.absorption_coefficient_wavelength_2 = None
        self.scattering_coefficient_wavelength_1 = None
        self.scattering_coefficient_wavelength_2 = None

    def read_data_from_txt(self):
        amplitudes = np.loadtxt(Path(self.location, 'amplitude.csv'), delimiter=',')
        phases = np.loadtxt(Path(self.location, 'phase.csv'), delimiter=',')

        return amplitudes, phases

    def read_measurement_settings(self):
        with open(Path(self.location, 'measurement settings.json')) as f:
            settings = json.load(f)
        return settings

    def linearize_amplitudes(self, data: fdNIRSData):
        squared = data.amplitudes * self.separations ** 2
        return np.log(squared)

    def get_slopes(self, data: fdNIRSData):
        data.amplitude_slopes = np.divide(np.diff(self.linearize_amplitudes(data)),
                                          np.diff(self.separations))
        data.phase_slopes = np.divide(np.diff(np.deg2rad(data.phases)),
                                      np.diff(self.separations))

    def compute_optical_parameters(self, amplitude_slope, phase_slope):
        y = phase_slope
        y = np.ravel(y)
        amplitude_slope = np.ravel(amplitude_slope)
        m = np.mean(phase_slope)
        b, a = signal.butter(1, .045)
        # print(np.ravel(y))
        filtered = signal.filtfilt(b, a, y - m) + m
        phase_slope = filtered
        print(filtered.shape)

        w = 2 * np.pi * self.modulation_frequency
        n = 1.4
        c = 2.998e11

        absorption_coefficient = (
                (w / (2 * (c / n))) *
                (
                        np.divide(phase_slope, amplitude_slope) -
                        np.divide(amplitude_slope, phase_slope)
                )
        )
        print(absorption_coefficient.shape)
        scattering_coefficient = (
                (np.square(amplitude_slope) - np.square(phase_slope)) /
                (3 * absorption_coefficient)
        )
        return np.array([absorption_coefficient, scattering_coefficient])

    def compute_hemoglobin_concentrations(self):
        oxy_hemoglobin_extinction_coefficient_830 = 974.0 / 10
        oxy_hemoglobin_extinction_coefficient_685 = 272.8 / 10

        deoxy_hemoglobin_extinction_coefficient_830 = 693.04 / 10
        deoxy_hemoglobin_extinction_coefficient_685 = 2188.24 / 10

        self.oxy_hemoglobin_concentration = (
                                                    ((
                                                             self.absorption_coefficient_wavelength_1 * deoxy_hemoglobin_extinction_coefficient_685) -
                                                     (
                                                             self.absorption_coefficient_wavelength_2 * deoxy_hemoglobin_extinction_coefficient_830)) /
                                                    ((
                                                             deoxy_hemoglobin_extinction_coefficient_685 * oxy_hemoglobin_extinction_coefficient_830) -
                                                     (
                                                             deoxy_hemoglobin_extinction_coefficient_830 * oxy_hemoglobin_extinction_coefficient_685))
                                            ) * 1e6

        self.deoxy_hemoglobin_concentration = (
                                                      ((
                                                               self.absorption_coefficient_wavelength_2 * oxy_hemoglobin_extinction_coefficient_830) -
                                                       (
                                                               self.absorption_coefficient_wavelength_1 * oxy_hemoglobin_extinction_coefficient_685)) /
                                                      ((
                                                               oxy_hemoglobin_extinction_coefficient_830 * deoxy_hemoglobin_extinction_coefficient_685) -
                                                       (
                                                               oxy_hemoglobin_extinction_coefficient_685 * deoxy_hemoglobin_extinction_coefficient_830))
                                              ) * 1e6

        self.total_hemoglobin_concentration = self.oxy_hemoglobin_concentration + self.deoxy_hemoglobin_concentration

    def plot_occlusion(self, name: str, total_time, occlusion_interval, window_size=None):

        oxy = rolling_apply(np.mean, self.oxy_hemoglobin_concentration, window_size)
        deoxy = rolling_apply(np.mean, self.deoxy_hemoglobin_concentration, window_size)
        total = oxy + deoxy

        t = np.linspace(0, total_time, self.oxy_hemoglobin_concentration.size)
        plt.figure(f'{name} occlusion')
        plt.plot(t, oxy, color='red', alpha=0.35, linewidth=3, label='Oxy')
        plt.plot(t, deoxy, color='blue', alpha=0.35, linewidth=3, label='Deoxy')
        plt.plot(t, total, color='black', alpha=0.35, linewidth=3, label='Total')
        plt.axvspan(occlusion_interval[0], occlusion_interval[1], color='green', alpha=0.15, label=f'{name} occlusion')

        plt.legend()
        plt.title(f'{name} occlusion')
        plt.xlabel(f'Time(s)')
        plt.ylabel(r'$\mu M$')
        plt.tight_layout()

    def plot_raw(self, total_time, occlusion_interval=(0, 0), window_size=None):
        t = np.linspace(0, total_time, self.amplitudes[:, 1].size)
        fig = plt.figure('Raw Data', figsize=(12, 8), layout='constrained')
        fig.suptitle(f'Raw Amplitude and Phase data')
        subfig = fig.subfigures(1, 2)

        subfig[0].suptitle(f'For 830nm')

        axes = subfig[0].subplots(4, 2)
        for i in range(4):
            ax = axes[i][0]
            ax.plot(t, rolling_apply(np.mean, self.amplitudes[:, i], window_size))
            ax.axvspan(occlusion_interval[0], occlusion_interval[1], color='green', alpha=0.15)
            ax.set_title(f'{self.separations.ravel()[i]}mm separation, Pair {i//2 + 1}')
            ax.set_ylabel('mV')
            ax.set_xlabel('Time(s)')

            ax = axes[i][1]
            ax.plot(t, rolling_apply(np.mean, self.phases[:, i], window_size))
            ax.axvspan(occlusion_interval[0], occlusion_interval[1], color='green', alpha=0.15)
            ax.set_title(f'{self.separations.ravel()[i]}mm separation, Pair {i//2 + 1}')
            ax.set_ylabel('Degrees')
            ax.set_xlabel('Time(s)')

        subfig[1].suptitle(f'For 690nm')
        axes = subfig[1].subplots(4, 2)
        for i in range(4):
            ax = axes[i][0]
            ax.plot(t, rolling_apply(np.mean, self.amplitudes[:, i + 4], window_size))
            ax.axvspan(occlusion_interval[0], occlusion_interval[1], color='green', alpha=0.15)
            ax.set_title(f'{self.separations.ravel()[i + 4]}mm separation, Pair {i//2 + 1}')
            ax.set_ylabel('mV')
            ax.set_xlabel('Time(s)')

            ax = axes[i][1]
            # x = rolling_apply(np.mean, self.phases[:, i + 4], window_size)
            x = self.phases[:, i + 4]
            filtered = apply_butterworth(x)
            x = apply_kalman_1d(x[0], 1, 50, x)
            ax.plot(t, rolling_apply(np.mean, self.phases[:, i + 4], window_size))
            # ax.plot(t, x, color='red')
            ax.plot(t, filtered, color='red', alpha=0.5)
            ax.axvspan(occlusion_interval[0], occlusion_interval[1], color='green', alpha=0.15)
            ax.set_title(f'{self.separations.ravel()[i + 4]}mm separation, Pair {i//2 + 1}')
            ax.set_ylabel('Degrees')
            ax.set_xlabel('Time(s)')

        # plt.tight_layout()

    def plot_absorption(self, total_time, occlusion_interval=(0, 0), window_size=None):
        absoprtion_830 = rolling_apply(np.mean, self.absorption_coefficient_wavelength_1, window_size)
        absoprtion_690 = rolling_apply(np.mean, self.absorption_coefficient_wavelength_2, window_size)
        # total = oxy + deoxy

        t = np.linspace(0, total_time, absoprtion_830.size)
        plt.figure(f'Absorption coefficients')
        plt.plot(t, absoprtion_830, color='blue', alpha=0.35, linewidth=3, label='830nm')
        plt.plot(t, absoprtion_690, color='red', alpha=0.35, linewidth=3, label='690nm')
        plt.axvspan(occlusion_interval[0], occlusion_interval[1], color='green', alpha=0.15, label=f'Occlusion')

        plt.legend()
        plt.title(f'Absorption Coefficients')
        plt.xlabel(f'Time(s)')
        plt.ylabel(r'1/mm$')
        plt.tight_layout()


class DualSlopeMeasurement(fdNIRS):
    def __init__(self, common='detector', *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.common = common
        self.data = fdNIRSData(self.location)
        self.reshape_data()
        self.get_slopes(self.data)

        self.amplitude_slopes_wavelength_1 = self.data.amplitude_slopes[:, 0, :]
        self.amplitude_slopes_wavelength_2 = self.data.amplitude_slopes[:, 1, :]
        self.amplitude_dual_slope_wavelength_1 = np.mean(self.amplitude_slopes_wavelength_1, axis=1)
        self.amplitude_dual_slope_wavelength_2 = np.mean(self.amplitude_slopes_wavelength_2, axis=1)

        self.phase_slopes_wavelength_1 = self.data.phase_slopes[:, 0, :]
        self.phase_slopes_wavelength_2 = self.data.phase_slopes[:, 1, :]
        self.phase_dual_slope_wavelength_1 = np.mean(self.phase_slopes_wavelength_1, axis=1)
        self.phase_dual_slope_wavelength_2 = np.mean(self.phase_slopes_wavelength_2, axis=1)

        (self.absorption_coefficient_wavelength_1,
         self.scattering_coefficient_wavelength_1) = self.compute_optical_parameters(
            self.amplitude_dual_slope_wavelength_1,
            self.phase_dual_slope_wavelength_1)

        (self.absorption_coefficient_wavelength_2,
         self.scattering_coefficient_wavelength_2) = self.compute_optical_parameters(
            self.amplitude_dual_slope_wavelength_2,
            self.phase_dual_slope_wavelength_2)

    def reshape_data(self):
        if self.common == 'detector':
            self.data.amplitudes = (
                self.data.amplitudes.reshape(self.data.amplitudes.shape[0], 2, 2, 2))
            self.data.phases = (
                self.data.phases.reshape(self.data.phases.shape[0], 2, 2, 2))

        elif self.common == 'source':
            self.data.amplitudes = self.data.amplitudes.reshape(self.data.amplitudes.shape[0], 2, 2, 2).swapaxes(2, 3)
            self.data.phases = self.data.phases.reshape(self.data.phases.shape[0], 2, 2, 2).swapaxes(2, 3)

    def plot_slopes(self, total_time, occlusion_interval=None, window_size=None):

        t = np.linspace(0, total_time, self.data.amplitude_slopes[:, 0, 0].size)
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

        ax = axes[0][0]
        ax.set_title('Amplitude Slopes - 830nm')
        ax.plot(t, rolling_apply(np.mean, self.data.amplitude_slopes[:, 0, 0], window_size),
                linewidth=2, alpha=0.75, label='Pair 1', color='red')
        ax.plot(t, rolling_apply(np.mean, self.data.amplitude_slopes[:, 0, 1], window_size),
                linewidth=2, alpha=0.75, label='Pair 2', color='blue')
        ax.plot(t, rolling_apply(np.mean, self.amplitude_dual_slope_wavelength_1, window_size),
                linewidth=2, alpha=0.75, label='Dual', color='black')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('1/mm')

        ax = axes[1][0]
        ax.set_title('Phase Slopes - 830nm')
        ax.plot(t, np.rad2deg(rolling_apply(np.mean, self.data.phase_slopes[:, 0, 0], window_size)),
                linewidth=2, alpha=0.75, label='Pair 1', color='red')
        ax.plot(t, np.rad2deg(rolling_apply(np.mean, self.data.phase_slopes[:, 0, 1], window_size)),
                linewidth=2, alpha=0.75, label='Pair 2', color='blue')
        ax.plot(t, np.rad2deg(rolling_apply(np.mean, self.phase_dual_slope_wavelength_1, window_size)),
                linewidth=2, alpha=0.75, label='Dual', color='black')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Degrees/mm(°)')

        ax = axes[0][1]
        ax.set_title('Amplitude Slopes - 690nm')
        ax.plot(t, rolling_apply(np.mean, self.data.amplitude_slopes[:, 1, 0], window_size),
                linewidth=2, alpha=0.75, label='Pair 1', color='red')
        ax.plot(t, rolling_apply(np.mean, self.data.amplitude_slopes[:, 1, 1], window_size),
                linewidth=2, alpha=0.75, label='Pair 2', color='blue')
        ax.plot(t, rolling_apply(np.mean, self.amplitude_dual_slope_wavelength_2, window_size),
                linewidth=2, alpha=0.75, label='Dual', color='black')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('1/mm')

        ax = axes[1][1]
        ax.set_title('Phase Slopes - 690nm')
        ax.plot(t, np.rad2deg(rolling_apply(np.mean, self.data.phase_slopes[:, 1, 0], window_size)),
                linewidth=2, alpha=0.75, label='Pair 1', color='red')
        ax.plot(t, np.rad2deg(rolling_apply(np.mean, self.data.phase_slopes[:, 1, 1], window_size)),
                linewidth=2, alpha=0.75, label='Pair 2', color='blue')
        ax.plot(t, np.rad2deg(rolling_apply(np.mean, self.phase_dual_slope_wavelength_2, window_size)),
                linewidth=2, alpha=0.75, label='Dual', color='black')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Degrees/mm(°)')

        if occlusion_interval is not None:
            for i in range(2):
                for j in range(2):
                    axes[i][j].axvspan(occlusion_interval[0], occlusion_interval[1],
                                     alpha=0.15, label='Occlusion', color='green')

        plt.legend()
        plt.tight_layout()

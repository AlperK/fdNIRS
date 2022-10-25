import numpy as np
from scipy.optimize import fsolve


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


def get_slopes(absorption_coefficient, scattering_coefficient, frequency):
    slopes = fsolve(_slope_equations,
                    np.array([-0.1, 0.1]),
                    args=(frequency, absorption_coefficient, scattering_coefficient))

    return slopes

def asdasd:
    pass
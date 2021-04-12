"""Electrically large lossy circular waveguide

This suite of functions implements both the approximate and the exact solutions
for a lossy, air-filled circular waveguide embedded in a homogeneous medium.
The approximate one is only valid for electrically large tunnels, which can be
checked using this module as well. The exact solution works for any combination
of waveguide size and wavelength.

A major restriction currently is that the dielectric medium surrounding the
waveguide is homogeneous. These solutions also do not support partially-filled
waveguides. Unknown at this point whether they ever will.

The formulas for the approximate solution are based on the following paper,
and are pretty straightforward:

[1] E. A. J. Marcatili and R. A. Schmeltzer,
Hollow metallic and dielectric waveguides for long distance optical
transmission and lasers,
The Bell System Technical Journal, vol. 43, no. 4, pp. 1783–1809, Jul. 1964,
doi: 10.1002/j.1538-7305.1964.tb04108.x

The formulas for the exact solution are based on the following papers, and
unkile the approximate one are less straightforward, requiring numerical
solutions and finding approximate roots:

[2] D. G. Dudley and S. F. Mahmoud,
Linear source in a circular tunnel,
IEEE Transactions on Antennas and Propagation, vol. 54, no. 7, pp. 2034–2047,
Jul. 2006, doi: 10.1109/TAP.2006.877195

[3] C. L. Holloway, D. A. Hill, R. A. Dalke, and G. A. Hufford,
Radio wave propagation characteristics in lossy circular waveguides such as
tunnels, mine shafts, and boreholes,
IEEE Transactions on Antennas and Propagation, vol. 48, no. 9, pp. 1354–1366,
Sep. 2000, doi: 10.1109/8.898768

[4] J. I. Glaser,
Attenuation and Guidance of Modes on Hollow Dielectric Waveguides,
IEEE Transactions on Microwave Theory and Techniques, vol. 17, no. 3,
pp. 173–174, Mar. 1969, doi: 10.1109/TMTT.1969.1126923

The steps to follow are:
    1. Pick a mode
    2. Estimate roots for that particular mode using high-frequency
       approximation
    3. Find the exact roots for that particular mode using a numerical solver
    4. Derive the propagation constants from the exact roots

In the above, the term `root` refers to solutions to a particular modal
characteristic equation.
"""

import warnings
from typing import Tuple

import numpy as np
import scipy.special
from scipy.constants import speed_of_light
from scipy.optimize import newton

np.seterr(divide='raise')


def estimate_roots_te(freq: float, wvg_diameter: float,
                      real_permittivity: float, mode_m: int) -> complex:
    """Estimate root position for a TE0m mode

    Uses the high-frequency approximation to find an initial estimate for
    a root to the modal equation for a TE0m mode. This is used subsequently
    in the numerical solution for the propagation constant.

    Args:
        freq: A `float` with the frequency of interest. Units are GHz.
        wvg_diameter: A `float` with the diameter of the waveguide. Units are
                      metres.
        real_permittivity: A `float` with the real part of the complex relative
                           permittivity of the material surrounding the
                           waveguide
        mode_m: An `int` with the mode index of interest. Must be >= 1.

    Returns:
        A single `complex` number containing the approximate position of the
        root in the lambda-plane.

    Raises:
        ZeroDivisionError: In case the `freq` parameter is given as zero
    """

    wvg_radius = wvg_diameter / 2

    try:
        wavelength = speed_of_light / (freq * 1e9)
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0'). \
              with_traceback(error.__traceback__)

    wave_number = 2 * np.pi / wavelength

    imag_denom = wave_number * wvg_radius * np.sqrt(real_permittivity - 1)
    root_multiplier = 1 + 1j / imag_denom

    root_estimate = scipy.special.jn_zeros(1, mode_m)[-1]
    root_estimate = root_estimate * root_multiplier

    return root_estimate


def estimate_roots_tm(freq: float, wvg_diameter: float,
                      real_permittivity: float,
                      mode_m: int) -> Tuple[complex, complex, bool]:
    """Estimate root position for a TM0m mode

    Uses the high-frequency approximation to find an initial estimate for
    a root to the modal equation for a TM0m mode. This is used subsequently
    in the numerical solution for the propagation constant.

    Args:
        freq: A `float` with the frequency of interest. Units are GHz.
        wvg_diameter: A `float` with the diameter of the waveguide. Units are
                      metres.
        real_permittivity: A `float` with the real part of the complex relative
                           permittivity of the material surrounding the
                           waveguide
        mode_m: An `int` with the mode index of interest. Must be >= 1.

    Returns:
        Unlike the TE0m case, the TM0m root estimate returns two possible
        values plus a boolean check for the validity of the approximation.
        In practice, supplying both to the numerical solver improves
        performance.

    Raises:
        ZeroDivisionError: In case the `freq` parameter is given as zero
    """

    wvg_radius = wvg_diameter / 2

    try:
        wavelength = speed_of_light / (freq * 1e9)
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0'). \
              with_traceback(error.__traceback__)

    wave_number = 2 * np.pi / wavelength

    imag_denom = wave_number * wvg_radius * np.sqrt(real_permittivity - 1)

    root_multiplier = 1 + (1j * real_permittivity / imag_denom)

    root_estimate_1 = scipy.special.jn_zeros(1, mode_m)[-1]
    root_estimate_1 = root_estimate_1 * root_multiplier

    root_estimate_2 = scipy.special.jn_zeros(0, mode_m)[-1]

    root_limit = np.sqrt(imag_denom / real_permittivity)

    root_limit_check = root_estimate_2 > (10 * root_limit)

    root_estimate_2 += (
        (1j * imag_denom) / (real_permittivity * root_estimate_2)
    )

    return (root_estimate_1, root_estimate_2, root_limit_check)


def estimate_roots_hybrid(mode: str, mode_n: int, mode_m: int) -> float:
    """Estimate root position for a hybrid EH/HE mods

    Uses the high-frequency approximation to find an initial estimate for
    the real part of a root to the modal equation for an EH/HE mode.
    This is used subsequently in the numerical solution for the propagation
    constant.

    Args:
        mode_n: An `int` with the first mode index of interest. Must be >= 1.
        mode_m: An `int` with the second mode index of interest. Must be >= 1.

    Returns:
        A single `complex` number containing the approximate position of the
        root in the lambda-plane.

    Raises:
        Nothing
    """

    if 0 == mode_n or 0 == mode_m:
        raise ValueError("Hybrid modes cannot have a zero index")

    if "he" == mode.lower():
        root_estimate = scipy.special.jn_zeros(mode_n-1, mode_m)[-1]
    elif "eh" == mode.lower():
        root_estimate = scipy.special.jn_zeros(mode_n+1, mode_m)[-1]
    else:
        raise ValueError("Mode should be either EH or HE")

    return root_estimate


def root_to_prop_const(freq: float, mode_wavelength: complex) -> complex:
    """Calculates propagation constant from modal equation root

    This is used to move from the lambda (mode wavelength) to the gamma
    (mode propagation constant) plane.

    Args:
        freq: A `float` with the frequency of interest. Units are GHz.
        mode_wavelength: A `complex` number containing the mode
                         wavelength from the numerical solver.

    Returns:
        The mode propagation constant in the form `alpha + j * beta`, where
        `alpha` is the attenuation rate in Np/m and `beta` is the phase
        constant in `rad/m`.

    Raises:
        ZeroDivisionError: In case the `freq` parameter is given as zero.
    """

    try:
        wavelength = speed_of_light / (freq * 1e9)
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0'). \
              with_traceback(error.__traceback__)

    wave_number = 2 * np.pi / wavelength

    prop_const = np.power(wave_number, 2) - np.power(mode_wavelength, 2)
    prop_const = np.sqrt(prop_const)
    prop_const *= 1j

    return prop_const


def modal_equation_te_tm(root_estimate: complex, freq: float,
                         permittivity: complex, wvg_diameter: float,
                         mode: str, mode_n: int = 0) -> complex:
    """Modal equation for TE0m and TM0m modes

    Implements the modal equation for TE0m and TM0m modes. Also used for
    EH/HE modes. Used by `calc_propagation_constant_exact` to numerically
    find the roots, and therefore propagation constants, for a specified mode.

    Notes:
        1. The imaginary part of the complex relative permittivity should have
        a negative sign pre-applied before being passed as an argument to this
        function.

    Args:
        root_estimate: A `complex` number with an estimate for the root. Can
                       also be an arbitrary value if evaluating the equation
                       over a field of values.
        freq: A `float` with the frequency of interest. Units are GHz.
        permittivity: A `complex` number with the relative permittivity of the
                      material surrounding the waveguide.
        wvg_diameter: A `float` with the diameter of the waveguide. Units are
                      metres
        mode: A `str` with the mode of interest, i.e. TE or TM.
        mode_n: An `int` with the first mode index. It is 0 by default for the
                TE and TM modes, however can be 1 or greater if this is used
                for the hybrid HE/EH modes.

    Returns:
        A `complex` number with the value of the modal equation for the given
        root estimate.

    Raises:
        ValueError: If a mode different than TE or TM is specified
    """

    # * The only difference in the modal equation betweeh TE and TM modes is
    # * a single multiplier for one of the fractions
    if "tm" == mode.lower():
        mode_multiplier = permittivity
    elif "te" == mode.lower():
        mode_multiplier = 1
    else:
        raise ValueError("Mode must be TE or TM")

    wvg_radius = wvg_diameter / 2
    wave_number = (2 * np.pi * freq * 1e9) / speed_of_light

    lambda_1 = root_estimate / wvg_radius

    lambda_2 = np.power(wave_number, 2) * (permittivity - 1)
    lambda_2 += np.power(lambda_1, 2)
    lambda_2 = np.sqrt(lambda_2)

    lambda_ratio = lambda_1 / lambda_2

    bessel_argument = lambda_1 * wvg_radius
    hankel_argument = lambda_2 * wvg_radius

    bessel_fraction = scipy.special.jvp(mode_n, bessel_argument)
    bessel_fraction /= scipy.special.jv(mode_n, bessel_argument)

    hankel_fraction = scipy.special.h2vp(mode_n, hankel_argument)
    hankel_fraction /= scipy.special.hankel2(mode_n, hankel_argument)

    result = bessel_fraction
    result -= (mode_multiplier * lambda_ratio * hankel_fraction)

    return result


def modal_equation_hybrid(root_estimate: complex, freq: float,
                          permittivity: complex, wvg_diameter: float,
                          mode: str, mode_n: int = 1) -> complex:
    """Modal equation for HE and EH modes

    Implements the modal equation for HE and EH modes. Uses the function
    for TE and TM modes internally. Used by `calc_propagation_constant_exact`
    to numerically find the roots, and therefore propagation constants,
    for a specified mode.

    Notes:
        1. The imaginary part of the complex relative permittivity should have
        a negative sign pre-applied before being passed as an argument to this
        function.
        2. The notation `q1`, `q2`, and `q5` are taken from the paper on
        which this is based, [2].

    Args:
        root_estimate: A `complex` number with an estimate for the root. Can
                       also be an arbitrary value if evaluating the equation
                       over a field of values.
        freq: A `float` with the frequency of interest. Units are GHz.
        permittivity: A `complex` number with the relative permittivity of the
                      material surrounding the waveguide.
        wvg_diameter: A `float` with the diameter of the waveguide. Units are
                      metres
        mode: A `str` with the mode of interest, i.e. TE or TM.
        mode_n: An `int` with the first mode index. It is 1 by default, however
                in general it needs to be  1 or greater.

    Returns:
        A `complex` number with the value of the modal equation for the given
        root estimate.

    Raises:
        ValueError: If a mode different than HE or EH is specified
    """

    if "he" != mode.lower() and "eh" != mode.lower():
        raise ValueError("Mode must be HE or EH")

    wvg_radius = wvg_diameter / 2
    wave_number = (2 * np.pi * freq * 1e9) / speed_of_light

    lambda_1 = root_estimate / wvg_radius

    lambda_2 = np.power(wave_number, 2) * (permittivity - 1)
    lambda_2 += np.power(lambda_1, 2)
    lambda_2 = np.sqrt(lambda_2)

    q1_multiplier_1 = mode_n / root_estimate

    q1_multiplier_2 = root_estimate / (wave_number * wvg_radius)
    q1_multiplier_2 = np.sqrt(1 - np.power(q1_multiplier_2, 2))

    q1_multiplier_3 = 1 - np.power(lambda_1 / lambda_2, 2)

    q1 = q1_multiplier_1 * q1_multiplier_2 * q1_multiplier_3

    q2 = modal_equation_te_tm(
        root_estimate, freq, permittivity, wvg_diameter, 'te', mode_n
    )

    q5 = modal_equation_te_tm(
        root_estimate, freq, permittivity, wvg_diameter, 'tm', mode_n
    )

    result = np.power(q1, 2) - q2 * q5

    return result


def calc_propagation_constant_exact(freq: float, wvg_diameter: float,
                                    permittivity: complex, mode: str,
                                    mode_n: int, mode_m: int,
                                    imag_range_min: float = 0.2,
                                    imag_range_max: float = 2.0) -> complex:
    """Main function for electrically small lossy waveguides

    Combines other functions in this module to calculate the propagation
    constant for a particular mode of interest at a particular frequency of
    interest for a given combination of waveguide diameter and surrounding
    dielectric material. Uses a numerical solver to find the root.

    Notes:
        1. The imaginary part of the complex relative permittivity should have
        a negative sign pre-applied before being passed as an argument to this
        function.

    Args:
        freq: A `float` with the frequency of interest. Units are GHz.
        wvg_diameter: A `float` with the diameter of the waveguide. Units are
                      metres.
        permittivity: A `complex` number with the relative permittivity of the
                      material surrounding the waveguide
        mode: A `str` with the mode of interest. Must be TE, TM, HE, or EH.
        mode_n: An `int` with the first mode index. Ignored in the case of
                TE or TM modes.
        mode_m: An `int` with the second mode index. Must be >= 1.
        imag_range_min: Used in the search for approximate root positions for
                        the hybrid modes. Defaults to 0.2.
        imag_range_max: Used in the search for approximate root positions for
                        the hybrid modes. Defaults to 2.0.

    Raises:
        ValueError: If a mode different than TE, TM, HE, or EH is specified.

    Returns:
        The mode propagation constant for the specified mode in the form
        `alpha + j * beta`, where `alpha` is the attenuation rate in Np/m
        and `beta` is the phase constant in `rad/m`.
    """

    wvg_radius = wvg_diameter / 2.0

    if "te" == mode.lower():
        root_estimate = estimate_roots_te(freq, wvg_diameter,
                                          permittivity.real, mode_m)

        root_exact, root_result = newton(
            modal_equation_te_tm,
            x0=root_estimate,
            args=(freq, permittivity, wvg_diameter, mode),
            full_output=True
        )

    elif "tm" == mode.lower():
        root_estimate_1, root_estimate_2, _ = estimate_roots_tm(
            freq, wvg_diameter, permittivity.real, mode_m
        )

        # * We provide both estimates to the numerical solver to increase
        # * speed and accuracy.
        root_exact, root_result = newton(
            modal_equation_te_tm,
            x0=root_estimate_1, x1=root_estimate_2,
            args=(freq, permittivity, wvg_diameter, mode),
            full_output=True
        )

    elif "eh" == mode.lower() or "he" == mode.lower():
        root_estimate_real = estimate_roots_hybrid(mode, mode_n, mode_m)

        # ! The `estimate_roots_hybrid` function provides an estimate for the
        # ! real part only. We have to manually build a search field in a
        # ! complex plane, over which we evaluate the reciprocal of the modal
        # ! equation. The position of the maximum is our new root esimate,
        # ! which is passed to the numeric solver.
        root_estimate_real_range = np.linspace(
            root_estimate_real*0.9, root_estimate_real*1.1,
            200
        )

        root_estimate_imag_range = np.linspace(
            imag_range_min, imag_range_max, 200
        )

        root_search_field = np.zeros(
            (np.size(root_estimate_real_range),
             np.size(root_estimate_imag_range)),
            dtype=np.complex128
        )

        for real_idx, real_part in enumerate(root_estimate_real_range):
            for imag_idx, imag_part in enumerate(root_estimate_imag_range):
                root_estimate = complex(real_part, imag_part)
                func_value = modal_equation_hybrid(
                    root_estimate, freq, permittivity, wvg_diameter,
                    mode, mode_n
                )
                root_search_field[real_idx, imag_idx] = 1 / func_value

        pole_index = np.unravel_index(
            np.argmax(root_search_field), root_search_field.shape
        )

        root_estimate = complex(
            root_estimate_real_range[pole_index[0]],
            root_estimate_imag_range[pole_index[1]]
        )
        # ! The numeric solver does not work with the built-in `complex` type,
        # ! so we convert it to the one defined by NumPy.
        root_estimate = np.complex128(root_estimate)

        root_exact, root_result = newton(
            modal_equation_hybrid,
            x0=root_estimate,
            args=(freq, permittivity, wvg_diameter, mode, mode_n),
            full_output=True
        )

    else:
        raise ValueError("Unknown mode specified")

    if not root_result.converged:
        warnings.warn("Newton solver not converged",
                      category=RuntimeWarning)

    mode_wavelength = root_exact / wvg_radius

    prop_constant = root_to_prop_const(freq, mode_wavelength)

    return prop_constant


def check_electrical_size(freq: float, wvg_diameter: float,
                          permittivity: complex,
                          mode_n: int, mode_m: int,
                          largeness_factor: int = 10) -> float:
    """Electrical size check for a circular waveguide

    Uses the formula in Marcatilli and Schmeltzer's 1964 paper to determine
    if a given circular waveguide is electrically large at a given frequency
    and waveguide propagation mode combination. This is also dependent on the
    complex relative permittivity of the surrounding medium.

    Notes:
        1. The imaginary part of the complex permittivity should have a
        negative sign pre-applied before being passed as an argument to this
        function.

    Args:
        freq: A `float` with the frequency at which to perform the check.
              Units are GHz.
        wvg_diameter: A `float` with the diameter of the waveguide which
                      is being checked. Units are metres.
        permittivity: A `complex` value of the relative permittivity of
                      the material surrounding the circular waveguide.
        mode_n: The `n` index of the mode of interest
        mode_m: The `m` index of the mode of interest.
        largeness_factor: An `int` with a multiplication factor used to turn
                          the 'much greater than' inequality into a simple
                          'greater than or equal to'. Unitless.

    Returns:
        A single `float` value showing to what extend the waveguide is large
        electrically compared to the wavelength.

    Raises:
        ZeroDivisionError: In case the `freq` parameter is given as zero.
        ValueError: In case the `permittivity` is specified as a negative
                    real number.
    """

    freq *= 1e9

    try:
        wavelength = speed_of_light / freq
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0'). \
              with_traceback(error.__traceback__)

    check_value_1 = (np.pi * wvg_diameter) / wavelength

    refr_index = np.sqrt(permittivity)

    if np.isnan(refr_index):
        raise ValueError('Material permittivity cannot be real and < 0')

    # ! The `jn_zeros` function returns a list of length `mode_m`, but we are
    # ! only interested in the mth zero, i.e. index mode_m-1
    bessel_root = scipy.special.jn_zeros(mode_n-1, mode_m)[mode_m-1]

    check_value_2 = np.abs(refr_index) * bessel_root
    check_value_2 *= largeness_factor

    try:
        check_result = check_value_1 / check_value_2
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Largeness factor must be > 0'). \
              with_traceback(error.__traceback__)

    return check_result


def calc_mode_refr_index(permittivity: complex, mode: str) -> complex:
    """Helper function for propagation constant calculation

    This is a helper function that calculates a mode-specific refractive
    index, something necessary to determine the attenuation and phase
    constants of a particular mode in a lossy circular waveguide.

    Notes:
        1. The imaginary part of the complex permittivity should have a
        negative sign pre-applied before being passed as an argument to this
        function.

    Args:
        permittivity: A `complex` value of the relative permittivity of
                      the material surrounding the circular waveguide.
        mode: A `str` specifying what type of mode is propagating along the
              waveguide. Valid values are TE, TM, HE, or EH.

    Returns:
        A `complex` number with the mode-specific refractive index.

    Raises:
        ValueError: In case the `permittivity` is specified as a negative
                    real number.
        ZeroDivisionError: In case the `permittivity` is that of free space,
                           i.e. 1 + 0j.
    """

    common_root = permittivity - 1
    common_root = np.sqrt(common_root)

    if np.isnan(common_root):
        raise ValueError('Soil permittivity cannot be real and < 0')

    if 'tm' == mode.lower():
        numerator = permittivity
    elif 'te' == mode.lower():
        numerator = 1
    elif ('he' == mode.lower()) or ('eh' == mode.lower()):
        numerator = (permittivity + 1) / 2
    else:
        raise ValueError('Mode must be TE, TM, EH, or HE')

    try:
        mode_refr_index = numerator / common_root
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Surrounding material cannot be free space'). \
              with_traceback(error.__traceback__)

    return mode_refr_index


def calc_attenuation_constant(freq: float, wvg_diameter: float,
                              permittivity: complex,
                              mode: str, mode_n: int, mode_m: int,
                              largeness_factor: int = 10) -> float:
    """Calculate the attenuation constant of a lossy circular waveguide

    This function calculates the attenuation constant for a particular mode
    in an electrically large circular waveguide. It uses other functions
    from this module internally.

    Notes:
        1. This function does not do any error handling, this is done in the
        other functions in the module. In case of an error there this function
        will simply re-raise the error.
        2. The imaginary part of the complex permittivity should have a
        negative sign pre-applied before being passed as an argument to this
        function.

    Args:
        freq: A `float` with the frequency at which to perform the check.
              Units are GHz.
        wvg_diameter: A `float` with the diameter of the waveguide which
                      is being checked. Units are metres.
        permittivity: A `complex` value of the relative permittivity of
                      the material surrounding the circular waveguide.
        mode: A `str` specifying what type of mode is propagating along the
              waveguide. Valid values are TE, TM, HE, or EH.
        mode_n: The `n` index of the mode of interest
        mode_m: The `m` index of the mode of interest.
        largeness_factor: An `int` with a multiplication factor used to turn
                          the 'much greater than' inequality into a simple
                          'greater than or equal to'. Unitless.

    Returns:
        The attenuation rate in Np/m as a `float` number.

    Raises:
        RuntimeWarning: In case the waveguide is not electrically large.
    """

    elec_size_check = check_electrical_size(freq, wvg_diameter, permittivity,
                                            mode_n, mode_m, largeness_factor)
    if elec_size_check <= 1.0:
        warnings.warn('Waveguide is not electrically large',
                      category=RuntimeWarning)

    mode_refr_index = calc_mode_refr_index(permittivity, mode)
    mode_refr_index = mode_refr_index.real

    freq *= 1e9
    wavelength = speed_of_light / freq

    # ! The `jn_zeros` function returns a list of length `mode_m`, but we are
    # ! only interested in the mth zero, i.e. index mode_m-1
    bessel_root = scipy.special.jn_zeros(mode_n - 1, mode_m)[mode_m - 1]

    alpha_1 = np.float_power(bessel_root / (2 * np.pi), 2)

    alpha_2 = np.float_power(wavelength, 2)
    alpha_2 /= np.float_power(wvg_diameter / 2, 3)

    alpha = alpha_1 * alpha_2 * mode_refr_index

    return alpha


def calc_phase_constant(freq: float, wvg_diameter: float,
                        permittivity: complex,
                        mode: str, mode_n: int, mode_m: int,
                        largeness_factor: int = 10) -> float:
    """Calculate the phase constant of a lossy circular waveguide

    This function calculates the phase constant for a particular mode
    in an electrically large circular waveguide. It uses other functions
    from this module internally.

    Notes:
        1. This function does not do any error handling, this is done in the
        other functions in the module. In case of an error there this function
        will simply re-raise the error.
        2. The imaginary part of the complex permittivity should have a
        negative sign pre-applied before being passed as an argument to this
        function.

    Args:
        freq: A `float` with the frequency at which to perform the check.
              Units are GHz.
        wvg_diameter: A `float` with the diameter of the waveguide which
                      is being checked. Units are metres.
        permittivity: A `complex` value of the relative permittivity of
                      the material surrounding the circular waveguide.
        mode: A `str` specifying what type of mode is propagating along the
              waveguide. Valid values are TE, TM, HE, or EH.
        mode_n: The `n` index of the mode of interest
        mode_m: The `m` index of the mode of interest.
        largeness_factor: An `int` with a multiplication factor used to turn
                          the 'much greater than' inequality into a simple
                          'greater than or equal to'. Unitless.

    Returns:
        The phase constant in rad/m as a `float` number.

    Raises:
        RuntimeWarning: In case the waveguide is not electrically large.
    """

    elec_size_check = check_electrical_size(freq, wvg_diameter, permittivity,
                                            mode_n, mode_m, largeness_factor)
    if elec_size_check <= 1.0:
        warnings.warn('Waveguide is not electrically large',
                      category=RuntimeWarning)

    mode_refr_index = calc_mode_refr_index(permittivity, mode)

    freq *= 1e9
    wavelength = speed_of_light / freq

    # ! The `jn_zeros` function returns a list of length `mode_m`, but we are
    # ! only interested in the mth zero, i.e. index mode_m-1
    bessel_root = scipy.special.jn_zeros(mode_n - 1, mode_m)[mode_m - 1]

    beta_0 = 2 * np.pi / wavelength

    beta_1 = mode_refr_index * wavelength / (np.pi * wvg_diameter / 2)
    beta_1 = beta_1.imag + 1

    beta_2 = bessel_root * wavelength / (np.pi * wvg_diameter)
    beta_2 = (np.float_power(beta_2, 2)) / 2

    beta = beta_0 * (1 - beta_2 * beta_1)

    return beta

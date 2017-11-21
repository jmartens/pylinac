"""
The NCS-18 module contains a number of helper functions and classes that can calculate parameters for performing the
NCS-18 absolute linac dose calibration.
Functions include all relevant calculations for NCS-18. Where values/equations are used they are specified in the
documentation.

Classes include photon, using cylindical and electron calibrations using plan-parallel chambers. Pass all the relevant
raw measurements and the class will compute all corrections and corrected readings and dose at 10cm and dmax/dref.
"""

from pylinac import Q_
import numpy
from warnings import warn


"""Quadratic fit coefficients for pulsed radiation as a function of the voltage ratio U1/U2). Data are taken from 
Weinhous and Meli
Source: NCS 18, Table 4 (p. 36)
"""
RECOMBINATION_COHORTS = {
    2.0: [2.299, -3.636, 2.337],
    2.5: [1.114, -1.587, 1.474],
    3.0: [0.677, -0.875, 1.198],
    3.5: [0.463, -0.542, 1.080],
    4.0: [0.341, -0.363, 1.022],
    5.0: [0.2135, -0.1875, 0.9745],
    6.0: [0.1495, -0.1075, 0.9584],
    8.0: [0.08750, -0.03732, 0.9502],
    10.0: [0.05909, -0.01041, 0.9516],
}

"""
List of the recommended graphite-walled cylindrical chambers together with the parameters of the sigmoid fit
Source: NCS 18, Table 10 (p. 54)
"""
CHAMBERS_PHOTONS = {
    # PTW
    '30012': {'X0': 0.9198, 'C': 11.67, 'A': 0.80},

    # IBA
    'FC65-G': {'X0': 0.9198, 'C': 11.67, 'A': 0.80},

    # Other
    '2561': {'X0': 0.8971, 'C': 15.15, 'A': 0.80},
    '2571': {'X0': 0.9198, 'C': 11.67, 'A': 0.80},
    '2611A': {'X0': 0.8971, 'C': 15.15, 'A': 0.80},
}

"""
Quadratic fit coefficients for pulsed radiation as a function of the voltage ratio U1/U2). Data are taken from Weinhous 
and Meli [1]_
Source: NCS-18, Table 4 (p. 36)

[1]_ Weinhous M.S., Meli A.M., Determinating Pion, the correction factor for recombination losses in an ionisation 
chamber, Med. Phys. 11 846-849, 1984.
"""
CHAMBERS_ELECTRONS = {
    # PTW
    '30010': {'A': 0.9345, 'B': 0.0057, 'C': 0.7733},

    # IBA
    'FC65G': {'A': 0.9345, 'B': 0.0057, 'C': 0.7733},

    # Other
    '2571': {'A': 0.9345, 'B': 0.0057, 'C': 0.7733},
    'NACP02': {'A': 1.1955, 'B': 0.2274, 'C': 0.1479},
    'Roos': {'A': 1.1376, 'B': 0.1700, 'C': 0.1835},
}


def k_tp(temp=Q_(20, 'celsius'), press=Q_(101.325, 'kPa')):
    """Calculate the temperature and pressure correction according to NCS-18 A.2 (p. 33) as defined in equation 17

    Parameters
    ----------
    temp : pint.Quantity
        The temperature as pint.Quantity or as float.
    press : pint.Quantity
        The pressure as pint.Quantity or as float.
    """
    if type(temp) != Q_:
        t = Q_(temp, 'celsius')
    else:
        # TODO: throw warning that we assume degrees Celsius
        t = temp

    if type(press) != Q_:
        # TODO: throw warning that we assume mbar
        p = Q_(press, 'mbar')
    else:
        p = press

    # Reference temperature and pressure as defined in NCS-18, A.2 (p. 33)
    t0 = Q_(20, 'celsius')
    p0 = Q_(101.325, 'kPa')

    # Calculate k_tp using base units to allow for quantities in different units
    return t.to_base_units()/t0.to_base_units() * p0.to_base_units()/p.to_base_units()


def k_pol(m_reference=(1, 2), m_negative=(-3, -4), m_positive=(5, 6)):
    """
    Calculate the polarity correction according to NCS-18 A.2 (p. 33) as defined in equation 18

    Parameters
    ----------
    m_reference : iterable
        The readings of the ion chamber at the reference polarity and voltage.
    m_negative : iterable
        The readings of the ion chamber at the negative polarity.
        This value should be of the opposite sign of the M positive value.
        If it's not, its sign will automatically be flipped.
    m_positive : iterable
        The readings of the ion chamber at the positive polarity.
        This value should be of the opposite sign of the M reference value.
    """
    m_negative_avg = numpy.mean(m_negative)
    m_positive_avg = numpy.mean(m_positive)
    m_reference_avg = numpy.mean(m_reference)
    # Technically, they're opposite charges, but some pass positive values for both, if same sign given, flip one.
    if numpy.sign(m_negative_avg) == numpy.sign(m_positive_avg):
        m_negative_avg = -m_negative_avg

    return (abs(m_positive_avg) + abs(m_negative_avg)) / (2 * m_reference_avg)


def r_50(i_50):
    """Calculate the R50 depth of an electron beam based on the I50 depth.

    Parameters
    ----------
    i_50 : float
        The value of I50 in cm.
    """
    if i_50 <= 10:
        r50 = 1.029 * i_50 - 0.06
    else:
        r50 = 1.59 * i_50 - 0.37
    return r50


def k_s(volt_normal=300, volt_low=100, m_normal=(1, 2), m_low=(3, 4)):
    """Calculate the recombination correction according to NCS-18 A.2 (p. 36) as defined in equation 21

    Parameters
    ----------
    volt_normal : int
        The "high" voltage; same as the NCS 18 measurement voltage.
    volt_low : int
        The "low" voltage; usually a third or less of the high voltage.
    m_normal : float, iterable
        The readings of the ion chamber at the "high" voltage.
    m_low : float, iterable
        The readings of the ion chamber at the "low" voltage.
    """

    ratio = round(volt_normal / volt_low, 1)
    if ratio < 3:
        warn('A voltage ratio < 3 is not recommended, see for details NCS 18, A. 2 (p. 36)', RuntimeWarning)

    if (ratio < min(RECOMBINATION_COHORTS)) | (ratio > max(RECOMBINATION_COHORTS)):
        raise ValueError('Unsupported ratio of voltages')

    try:
        fit_coefficients = RECOMBINATION_COHORTS[ratio]
        poly = numpy.poly1d(fit_coefficients)
        return poly(numpy.mean(m_normal) / numpy.mean(m_low))
    except KeyError:
        raise ValueError('Unsupported ratio of voltages')


def k_q(model='30012', tpr=None, r_50=None):
    """Calculate beam quality correction as described in NCS-18 A.2 (p. 50) as defined in equation 30 and fit parameters
    in Table 10 (p. 54)

    Parameters
    ----------
    model : str
        The model of the chamber. Valid values are those listed in
        Table III of Muir and Rodgers and Table I of the TG-51 Addendum.
    tpr : {>0.623, <0.805}
        The TPR ratio of the 20cm measurement divided by the 10cm measurement.
    r_50 : float
        The R50 value in cm of an electron beam.

    .. warning::
        Only 1 of  ``tpr`` or ``r_50`` can be defined.
    """
    # TODO: check if this range (from TG51) applies to NCS 18
    TPR_LOW = 0.623
    TPR_HIGH = 0.805

    # error checking
    if not any((tpr, r_50)):
        raise ValueError("At least one of the parameters tpr or r_50 must be defined.")
    if tpr and r_50 is not None:
        raise ValueError("Cannot define both a photon component (TPR) and an electron component (R50)")

    if tpr is not None:
        if tpr > TPR_HIGH or tpr < TPR_LOW:
            raise ValueError("Measured TPR is out of range; must be between {:2.2} and {:2.2}.".format(TPR_LOW, TPR_HIGH))
        else:
            ch = CHAMBERS_PHOTONS[model]
            # NCS-18 A.4, p. 40, eq 30
            return ch['A'] + (1 - ch['A']) * (1 + numpy.exp(ch['C'] * (0.57 - ch['X0'])))/(1 + numpy.exp(ch['C'] * (tpr - ch['X0'])))

    if r_50 is not None:
        ch = CHAMBERS_ELECTRONS[model]
        # Farmer type chambers, eq 37
        return ch['A'] - ch['B'] * r_50**ch['C']


def m_corrected(k_tp=1.0, k_pol=1.0, k_s=1.0, m_raw=(1.1, 2.2)):
    """Calculate M_corrected, the ion chamber reading with all corrections applied.

    Parameters
    ----------
    k_tp : float
        The temperature & pressure correction.
    k_pol : float
        The polarity correction.
    k_s : float
        The recombination correction.
    m_raw : float, iterable
        The raw ion chamber readings.

    Returns
    -------
    float
    """
    return k_tp * k_pol * k_s * numpy.mean(m_raw)

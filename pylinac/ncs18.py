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


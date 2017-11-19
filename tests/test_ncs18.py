from unittest import TestCase

from pylinac import ncs18
from pylinac import Q_


class TestFunctions(TestCase):

    def test_k_tp(self):
        temperatures = (Q_(20, 'celsius'), Q_(20, 'celsius').to('fahrenheit'), 25, 15)
        pressures = (Q_(1013.25, 'mbar'), Q_(1013.25, 'mbar').to('mmHg'), 1023.25, 1003.25)
        expected_k_tps = (1.0, 1.0, 1.007, 0.9927)
        for temperature, pressure, exp in zip(temperatures, pressures, expected_k_tps):
            self.assertAlmostEqual(ncs18.k_tp(temperature, pressure), exp, delta=0.001)

    def test_k_pol(self):
        m_negatives = (300, -202, 198)
        m_positives = (-300, 198, 201)
        m_references = (300, 198, -201)
        expected_k_pols = (1.0, 1.010, -0.993)
        for reference, negative, positive, expected_k_pol in zip(m_references, m_negatives, m_positives, expected_k_pols):
            self.assertAlmostEqual(ncs18.k_pol(reference, negative, positive), expected_k_pol, delta=0.001)

    def test_k_s(self):
        low_vals = (20, 20.05)
        high_vals = (20, 20.1)
        expected_pion = (1.0, 1.0012)
        for low, high, exp in zip(low_vals, high_vals, expected_pion):
            self.assertAlmostEqual(ncs18.k_s(300, 100, high, low), exp, delta=0.001)

    def test_k_s_cohort_lbound(self):
        with self.assertRaises(ValueError):
            ncs18.k_s(100, 100, 1, 2)

    def test_k_s_cohort_ubound(self):
        with self.assertRaises(ValueError):
            ncs18.k_s(1100, 100, 1, 2)

    def test_k_s_cohort_ratio_too_low(self):
        with self.assertWarns(RuntimeWarning):
            ncs18.k_s(200, 100, 1, 2)

    def test_k_s_cohort_ratio_not_found(self):
        with self.assertRaises(ValueError):
            ncs18.k_s(360, 100, 1, 2)

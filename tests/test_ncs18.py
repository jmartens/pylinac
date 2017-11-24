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
        normal_vals = (20, 20.1)
        expected_pion = (1.0, 1.0012)
        for low, normal, exp in zip(low_vals, normal_vals, expected_pion):
            self.assertAlmostEqual(ncs18.k_s(300, 100, normal, low), exp, delta=0.001)

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

    def test_k_q_r50(self):
        self.assertAlmostEqual(ncs18.k_q(model='Roos', r_50=ncs18.r_50(6.16)), 0.899, delta=0.001)

    def test_r_50_le_10(self):
        self.assertAlmostEqual(ncs18.r_50(9.9), 10.127, delta=0.001)

    def test_r_50_eq_10(self):
        self.assertAlmostEqual(ncs18.r_50(10), 10.229, delta=0.001)

    def test_r_50_gt_10(self):
        self.assertAlmostEqual(ncs18.r_50(10.1), 15.689, delta=0.001)

    def test_k_q_no_tpr_and_r50(self):
        with self.assertRaises(ValueError):
            ncs18.k_q()

    def test_k_q_tpr_and_r50(self):
        with self.assertRaises(ValueError):
            ncs18.k_q(model='30012', tpr=0.0667, r_50=ncs18.r_50(0.67))

    def test_k_q_tpr_lbound(self):
        with self.assertRaises(ValueError):
            ncs18.k_q(model='30012', tpr=0)

    def test_k_q_tpr_ubound(self):
        with self.assertRaises(ValueError):
            ncs18.k_q(model='30012', tpr=1)

    def test_k_q_tpr(self):
        self.assertAlmostEqual(ncs18.k_q(model='30012', tpr=0.667), 0.993, delta=0.001)

    def test_m_corrected(self):
        self.assertAlmostEqual(ncs18.m_corrected(k_tp=1.0, k_pol=1.0, k_s=1.0, m_raw=(1.1, 2.2)), 1.65, delta=0.001)


class TestNCS18Base(TestCase):
    temperature = Q_(20, 'celsius')
    pressure = Q_(1013.25, 'mbar')
    model = '30012'
    volt_high = -400
    volt_low = -133
    m_raw = (20, 20, 20)
    m_opp = (20, 20, 20)
    m_low = (20, 20, 20)

    def setUp(self):
        self.ncs18 = ncs18.NCS18Base(
            temp=self.temperature,
            press=self.pressure,
            model=self.model,
            volt_reference=self.volt_high,
            volt_low=self.volt_low,
            m_raw=self.m_raw,
            m_reference=self.m_raw,
            m_opposite=self.m_opp,
            m_low=self.m_low,
            adjusted_m_raw=self.m_raw)

    def test_k_tp(self):
        self.assertAlmostEqual(self.ncs18.k_tp, 1, delta=0.0005)

    def test_k_pol(self):
        self.assertAlmostEqual(self.ncs18.k_pol, 1, delta=0.0005)

    def test_k_s(self):
        self.assertAlmostEqual(self.ncs18.k_s, 1, delta=0.0005)

    def test_m_corrected(self):
        self.assertAlmostEqual(self.ncs18.m_corrected, 20, delta=0.0005)

    def test_adjusted_m_corrected(self):
        self.assertAlmostEqual(self.ncs18.adjusted_m_corrected, 20, delta=0.0005)

    def test_output_was_adjusted(self):
        self.assertEqual(self.ncs18.output_was_adjusted, True)

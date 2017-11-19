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

"""
The NCS-18 module contains a number of helper functions and classes that can calculate parameters for performing the
NCS-18 absolute linac dose calibration.
Functions include all relevant calculations for NCS-18. Where values/equations are used they are specified in the
documentation.

Classes include photon, using cylindical and electron calibrations using plan-parallel chambers. Pass all the relevant
raw measurements and the class will compute all corrections and corrected readings and dose at 10cm and dmax/dref.
"""
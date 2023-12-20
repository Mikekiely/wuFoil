import logging

import numpy as np

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

class FlightConditions:
    """
    Calculates atmospheric Conditions at given altitude and mach number

    Cal

    Parameters
    -   h: Height (ft)
    -   mach: Mach number
    -   L: Characteristic Length
    -   Units: 'ft' or 'm'
    """

    p = 0           # Pressure (lb/ft^2)
    t = 0           # Temperature (Rankine)
    rho = 0         # Density (slug/ft^3)
    mu = 0          # Dynamic viscosity (Whatever stupid unit the imperial system uses for viscosity)
    a = 0           # Speed of sound (ft/s)
    v = 0           # Velocity (ft/s)
    mach = 0        # Mach Number
    q = 0           # Dynamic Pressure
    re = 0          # Reynolds Number
    altitude = 0    # Altitude (ft)
    gamma = 1.4     # Gas Constant

    def __init__(self, h: float, mach: float, length: float, aoa: float = None, cl: float = None, input_units: str = 'ft'):
        # Calculates everything in imperial then converts it back to SI

        # TODO: Make flight conditions calculate in metric, this is stupid
        if cl or cl == 0:
            self.aoa = aoa
            self.cl = None
        elif aoa or aoa == 0:
            self.aoa = 1
            self.cl = cl
        else:
            logging.warning('Please set either angle of attack or Cl before running any analysis')

        if input_units == 'm':
            length *= .3048
            h *= .3048

        # Constants
        self.length = length
        g = 32.2  # gravity
        t0 = 519  # Rankine
        p0 = 2116  # lb/ft^2
        rho0 = .0023769
        a0 = 1116  # ft/s
        L = 0.003575  # temp lapse rate (R/ft)
        R = 1716.5  # Gas Constant
        gamma = 1.4

        t = t0 - L * h
        theta = t / t0
        delta = theta ** (g / (L * R))
        p = delta * p0
        sigma = delta / theta

        # Stratosphere Correction
        if delta < 0.224:
            h0 = 36000  # Troposphere Ceiling
            beta = 20790
            dh = h - 36000
            theta = .752
            delta = theta ** (g / (L * R)) * np.exp(-dh / beta)
            sigma = delta / theta
            t = t0 * theta

        # Sutherland's law (calculates in metric, converts back to imperial
        c1 = 1.458e-6
        S = 110.4
        tref = t * 5 / 9  # Rankine to Kelvin
        self.mu = c1 * (tref ** 1.5) / (tref + S)
        self.mu *= .0208854  # metric to imperial

        self.t = t
        self.p = delta * p0
        self.rho = sigma * rho0
        self.a = np.sqrt(gamma * R * t)

        self.v = mach * self.a
        self.mach = mach
        self.q = .5 * self.rho * self.v ** 2
        self.altitude = h
        self.re = self.rho * self.v * self.length / self.mu
        self._to_SI()

    def _to_SI(self):
        """
        Converts to metric units for input into SU2
        Really all that matters here is freestream pressure and temperature
        """
        ft2meters = .3048
        self.p *= 47.8803       # Pascals
        self.rho *= 515.379     # Kg/m^3
        self.v *= ft2meters     # m/s
        self.t *= 5/9           # Kelvin
        self.mu *= 1/.0208854
        self.a *= ft2meters
        self.q = .5 * self.rho * self.v ** 2
        self.altitude *= ft2meters
        self.length *= ft2meters
        self.re = self.rho * self.v * self.length / self.mu

import numpy as np


class FlightConditions:
    """
    Calculates atmospheric Conditions at given altitude and mach number

    Parameters
    -   h: Height (ft)
    -   mach: Mach number
    -   L: Characteristic Length
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

    def __init__(self, h: float, mach: float, length: float, aoa: float = None, cl: float = None, units='imperial'):
        # Set aoa, set it to default value and set cl if fixed cl mode
        if aoa or aoa == 0:
            self.aoa = aoa
            self.cl = None
        elif cl or cl == 0:
            self.aoa = 1
            self.cl = cl

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

        if units.lower() == 'si':
            ft2meters = 0.3048

            self.length = length * ft2meters
            g = g / ft2meters
            

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

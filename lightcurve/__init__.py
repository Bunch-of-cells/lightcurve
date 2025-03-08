from scipy.optimize import fsolve
import numpy as np

V_SUN = 4.83
AU_BY_R_SUN = 215.1


class Orbit:
    def __init__(self, a: float, e: float, i: float, Omega: float, omega: float):
        """
        # Parameters
        a: semi-major axis (AU)
        e: eccentricity
        i: inclination
        Omega: longitude of the ascending node
        omega: argument of periapsis
        """
        self.a = a
        self.e = e
        self.i = i
        self.Omega = Omega
        self.omega = omega
        self.init_conversion()

    def init_conversion(self):
        m3 = np.array(
            [
                [np.cos(self.Omega), np.sin(self.Omega), 0],
                [-np.sin(self.Omega), np.cos(self.Omega), 0],
                [0, 0, 1],
            ]
        )
        m2 = np.array(
            [
                [1, 0, 0],
                [0, np.cos(self.i), np.sin(self.i)],
                [0, -np.sin(self.i), np.cos(self.i)],
            ]
        )
        m1 = np.array(
            [
                [np.cos(self.omega), np.sin(self.omega), 0],
                [-np.sin(self.omega), np.cos(self.omega), 0],
                [0, 0, 1],
            ]
        )
        m = (m1 @ m2 @ m3).T
        self.conversion = m

    def init_primaries(
        self, M1: float, M2: float, R1: float, R2: float, L1: float, L2: float
    ):
        """
        Initialize the properties of the two primaries

        # Parameters
        M1: mass of the first primary (solar mass)
        M2: mass of the second primary (solar mass)
        R1: radius of the first primary (solar radius)
        R2: radius of the second primary (solar radius)
        L1: luminosity of the first primary (solar luminosity)
        L2: luminosity of the second primary (solar luminosity)
        """
        self.M1 = M1
        self.M2 = M2
        self.R1 = R1
        self.R2 = R2
        self.L1 = L1
        self.L2 = L2
    
    def init_distance(self, D: float):
        """
        Set the distance between the system and the observer
        
        # Parameters
        dist: distance between the system and the observer (pc)
        """
        self.D = D

    def posiiton_sys(self, t: float):
        """
        Calculate the distance between two primaries in the system's coordinate system.
        Origin (0, 0, 0) is at the foci of the ellipse.
        periapsis = (a(1-e), 0, 0)

        # Parameters
        t: time since periapsis (years)

        # Returns
        distance between the primaries (AU)
        """
        M = 2 * np.pi * t / np.sqrt(self.a**3 / (self.M1 + self.M2))
        E = fsolve(lambda x: x - self.e * np.sin(x) - M, M)[-1]
        r = (
            self.a * (np.cos(E) - self.e),
            self.a * np.sqrt(1 - self.e**2) * np.sin(E),
            0.0,
        )
        return r

    def position_ref(self, t: float):
        """
        Calculate the distance between two primaries in the reference coordinate system (x-y plane is the plane of the sky)

        # Parameters
        t: time since periapsis (years)

        # Returns
        distance between the primaries (AU)
        """
        return self.conversion.dot(self.posiiton_sys(t))

    def luminosity(self, t: float) -> float:
        """
        Calculate the luminosity of the system, as seen in the reference frame (x-y plane is the plane of the sky)

        # Parameters
        t: time since periapsis (years)

        # Returns
        luminosity of the system (solar luminosity)
        """
        r = self.position_ref(t)
        d = np.sqrt(r[0] ** 2 + r[1] ** 2)
        blocked = area_blocked(self.R1, self.R2, d * AU_BY_R_SUN)
        b = self.L1 + self.L2
        if r[2] < 0:
            b -= blocked * self.L2 / np.pi / self.R2**2
        else:
            b -= blocked * self.L1 / np.pi / self.R1**2
        return b

    def abs_mag(self, t: float) -> float:
        """
        Calculate the absolute magnitude of the system, as seen in the reference frame (x-y plane is the plane of the sky)

        # Parameters
        t: time since periapsis (years)

        # Returns
        absolute magnitude of the system
        """
        return V_SUN - 2.5 * np.log10(self.luminosity(t))
    
    def app_mag(self, t: float) -> float:
        """
        Calculate the absolute magnitude of the system, as seen in the reference frame (x-y plane is the plane of the sky)

        # Parameters
        t: time since periapsis (years)

        # Returns
        absolute magnitude of the system
        """
        return self.abs_mag(t) - 5 + 5 * np.log10(self.D)


def area_blocked(r1: float, r2: float, d: float) -> float:
    """
    Calculate the blocked area by intersection of 2 disks.

    # Parameters
    r1: radius of the first disk
    r2: radius of the second disk
    d: distance between the centers of the disks

    # Returns
    visible area of the disk
    """
    if r1 + r2 <= d:
        return 0
    if r1 >= r2 + d:
        return np.pi * r2**2
    if r2 >= r1 + d:
        return np.pi * r1**2
    a1 = np.arccos((r1**2 + d**2 - r2**2) / (2 * r1 * d))
    a2 = np.arccos((r2**2 + d**2 - r1**2) / (2 * r2 * d))
    return r1**2 * a1 + r2**2 * a2 - r1 * d * np.sin(a1)

from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

def posiiton(e, M):
    """
    Calculate the position of a mass in an elliptical orbit.
    Origin (0, 0, 0) is at the foci of the ellipse.
    semi-major axis = 1
    periapsis = (1-e, 0, 0)
    
    # Parameters
    e: eccentricity
    M: mean anomaly
    
    # Returns
    position of the mass
    """ 
    E = fsolve(lambda x: x - e*np.sin(x) - M, M)[-1]
    r = (np.cos(E) - e, np.sqrt(1 - e**2)*np.sin(E), 0.0)
    return r

def projected_coordinates(r, i, Omega, omega):
    """
    Calculate the projected coordinates of a mass in an elliptical orbit.
    
    # Parameters
    r: position of the mass
    i: inclination
    Omega: longitude of the ascending node
    omega: argument of periapsis
    
    # Returns
    projected position of the mass
    """
    m3 = np.array([[np.cos(Omega), np.sin(Omega), 0], [-np.sin(Omega), np.cos(Omega), 0], [0, 0, 1]])
    m2 = np.array([[1, 0, 0], [0, np.cos(i), np.sin(i)], [0, -np.sin(i), np.cos(i)]])
    m1 = np.array([[np.cos(omega), np.sin(omega), 0], [-np.sin(omega), np.cos(omega), 0], [0, 0, 1]])
    m = (m1 @ m2 @ m3).T
    return m.dot(r)

def disk_visible(r):
    """
    Calculate the visible area of a disk.
    unocculted area of disk = 1
    
    # Parameters
    r: distance of origin from the center of the disk
    
    # Returns
    visible area of the disk
    """
    if r >= 1:
        return 1
    return (1 - (np.arccos(r) - r * np.sqrt(1 - r**2))/np.pi)


def luminosity(a, e, i, Omega, omega, R, M):
    """
    Calculate the luminosity (visible disk area) of a mass in an elliptical orbit.
    
    # Parameters
    e: eccentricity
    i: inclination
    Omega: longitude of the ascending node
    omega: argument of periapsis
    M: mean anomaly
    a: semi-major axis
    R: radius of the disk
    
    # Returns
    lightcurve of the mass
    """
    r = posiiton(e, M)
    r = projected_coordinates(r, i, Omega, omega)
    r = np.sqrt(r[0]**2 + r[1]**2)
    return disk_visible(r*a/R)


def main():
    a = 3
    e = 0.6
    i = 0.45*np.pi
    Omega = 0
    omega = 0.4*np.pi
    
    t = np.arange(-0.5, 2.5, 0.01)
    s = [np.log10(luminosity(a, e, i, Omega, omega, 1, 2 * np.pi * k)) for k in t]

    fig, ax = plt.subplots()
    ax.plot(t, s)

    ax.set(xlabel='time (s)', ylabel='magnitude',
        title='Lightcurve')
    ax.grid()

    fig.savefig("test.png")
    plt.show()


if __name__ == "__main__":
    main()

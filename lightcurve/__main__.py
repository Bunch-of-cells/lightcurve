import matplotlib.pyplot as plt
import numpy as np

from __init__ import Orbit

def main():
    a = 3
    e = 0.5
    i = 0.6*np.pi
    Omega = 1
    omega = 0.8
    orbit = Orbit(a, e, i, Omega, omega)
    orbit.init_primaries(1, 0.01, 1, 2, 10, 1)

    t = np.arange(-0.5, 2.5, 0.01)
    s = [orbit.magnitude(2 * np.pi * k) for k in t]

    fig, ax = plt.subplots()
    ax.plot(t, s)

    ax.set(xlabel='time (orbital period)', ylabel='magnitude',
        title='Lightcurve')
    ax.invert_yaxis()
    ax.grid()

    fig.savefig("test.png")
    plt.show()


if __name__ == "__main__":
    main()

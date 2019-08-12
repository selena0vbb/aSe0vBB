# Script to calculate best fit for the W_ehp(E) function
import numpy as np
import scipy.optimize as scp
import scipy.special
import matplotlib.pyplot as plt

# field = np.arange(10, 24, 2)
# wehp = np.array([62.5, 56, 51, 47.5, 43.5, 40.5, 37.8]) / 1000

field = np.arange(15, 45, 5)
wehp = 122.0 / np.array([425, 520, 610, 695, 780, 865])
# wehp = [35.44107434, 28.86732946, 24.60788534, 21.63541832, 19.44924246, 17.77693871]


def we(field, *params):
    """ functional form of wehp. w(E) = param[0] + param[1]/E """
    return params[0] + params[1] / field


def columnWe(field, *params):
    e = 1.6e-19
    kb = 1.38e-23
    T = 293
    eps0 = 8.85e-12
    epsr = 6
    w0 = 0.0055
    sin2thet = 0.55

    x = params[1] ** 2 * field ** 2 * sin2thet * e ** 2 / (kb ** 2 * T ** 2)

    wehp = w0 * (
        1
        + e ** 2
        / (4 * np.pi * eps0 * epsr * kb * T)
        * params[0]
        * np.exp(x)
        * scipy.special.kn(0, x)
    )
    return wehp


def fitWehp():
    """ Fit the wehp(E) to an 1/F function """

    # Fit the data
    initialParams = [4.5e8, 2.0e-9]
    popt, pcov = scp.curve_fit(
        columnWe, field * 1.0e6, wehp, initialParams, maxfev=10000
    )
    print (popt)
    eFieldAxis = np.linspace(10, 40, 200)

    # plot the result
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.scatter(field, wehp, marker="X", label="Data")
    ax.plot(
        eFieldAxis,
        columnWe(eFieldAxis * 1.0e6, *popt),
        "r",
        linewidth=2,
        label="Fit: %.3f + %.2f/E" % (popt[0], popt[1]),
    )
    ax.set_xlabel(r"$\vec{E}$ [$\rm V\ \mu m^{-1}$]", fontsize=18)
    ax.set_ylabel(r"$\rm W_{ehp} \ [keV]$", fontsize=18)
    ax.legend(fontsize=16)
    # ax.tick_params()
    plt.show()


if __name__ == "__main__":
    fitWehp()
    print (columnWe(field * 1e6, 4.5e8, 2.0e-9))

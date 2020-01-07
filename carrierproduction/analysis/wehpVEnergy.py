# Script to calculate best fit for the W_ehp(E) function
import numpy as np
import scipy.optimize as scp
import scipy.special
import matplotlib.pyplot as plt

# field = np.arange(10, 24, 2)
# wehp = np.array([62.5, 56, 51, 47.5, 43.5, 40.5, 37.8]) / 1000

field = np.arange(15, 45, 5)
wehp = 122.0 / np.array([425, 520, 610, 695, 780, 865]) * 0.2

# From Bubon et al papers
# field = np.array([10., 15., 20., 25., 30., 40., 50.,])
# wehp = np.array([42., 32.5, 27., 23.5, 20., 17., 14 ]) / 1000


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
    sin2thet = 0.6

    x = params[1] ** 2 * field ** 2 * sin2thet * e ** 2 / (4 * kb ** 2 * T ** 2)
    y = e ** 2 / (4 * np.pi * eps0 * epsr * kb * T) * params[0]
    # print(x)
    # print(y)
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
    eFieldAxis = np.linspace(10, field[-1], 200)

    # plot the result
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.scatter(field, wehp, marker="X", label="Data")
    ax.plot(
        eFieldAxis,
        columnWe(eFieldAxis * 1.0e6, *popt),
        "r",
        linewidth=2,
        label="Fit: $N_0$=%.3e,  b=%.2e" % (abs(popt[0]), abs(popt[1])),
    )
    ax.set_xlabel(r"$\vec{E}$ [$\rm V\ \mu m^{-1}$]", fontsize=18)
    ax.set_ylabel(r"$\rm W_{ehp} \ [keV / ehp]$", fontsize=18)
    ax.legend(fontsize=16)
    # ax.tick_params()
    plt.show()
    return popt


if __name__ == "__main__":
    par = fitWehp()
    # par = [5.46e8, 1.27e-9]

    testfield = 20
    print (columnWe(testfield * 1e6, *par))

    # Get parameters from the settings file
    b = par[1] * 1e3
    N0 = par[0] * 1e-3
    w0 = 0.0055

    # Compute diffusion and and recombination coefficients
    voltConvCoeff = 1.0e-6  # Conversion to convert volts from SI to mm and us
    mu = 29.02e-6 + 0.661e-6
    kb = 1.38e-29  # in mm^2 kg s^-2 K^-1
    T = 293
    q = 1.6e-19
    D = mu * kb * T / q

    eps0 = 8.85e-15 / voltConvCoeff  # in C V^-1 mm^-1
    epsr = 6
    alpha = (q * mu) / (eps0 * epsr)

    sinThetaSquare = 0.55
    testfield = 20 * 1e3

    # Compute x, defined in paper
    x = (
        mu ** 2
        * sinThetaSquare
        * b ** 2
        * (testfield * voltConvCoeff) ** 2
        / (4 * D ** 2)
    )

    # print(x)
    # print((alpha * N0 / (4 * np.pi * D)))

    # N0 is the linear charge density. Assume that energy is lost linearly along the length of the track.
    # N0 = event.energy / (w0 * trackLength)
    wehp = w0 * (
        1 + (alpha * N0 / (4 * np.pi * D)) * np.exp(x) * scipy.special.kn(0, x)
    )

    print (wehp)

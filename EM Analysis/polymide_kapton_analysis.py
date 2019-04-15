# Polymide Kapton Thickness analysis

from plot import *
import matplotlib.pyplot as plt
import numpy as np
import os
import brewer2mpl


def effectiveAreaOfDifferentKaptonWidth():

    # Data directory and get all file names
    fileLoc = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer"
    filelist = os.listdir(fileLoc)
    bias = np.arange(100, 301, 20)
    kaptonThickness = {
        "500nm": 0.5,
        "1um": 1,
        "2um": 2,
        "5um": 5,
        "10um": 10,
        "15um": 15,
        "3um": 3,
        "4um": 4,
    }
    legend = []
    effectiveArea = np.zeros((len(filelist), bias.size))

    # set the boundaries for a positive induced current on electrode A
    nElectrodes = 4
    electrodeABoundary = (np.arange(0, 150, 20) - 5).clip(min=0)
    for findx, file in enumerate(filelist):

        # Get data from file
        _, y, z, data = readComsolFileGrid(fileLoc + "\\" + file)
        y, z = y * 1e6, z * 1e6  # convert to um

        # get the kapton thickness from the filename
        fileIdentifier = file.split("_")[3]
        thc = kaptonThickness[fileIdentifier]
        legend.append(fileIdentifier)

        # Remove NaN
        data = np.nan_to_num(data, copy=False)

        # Split data into phi, Ey, and Ez
        phi, Ey, Ez = (
            data[:, 0 : np.amin(data.shape) : 3],
            data[:, 1 : np.amin(data.shape) : 3],
            data[:, 2 : np.amin(data.shape) : 3],
        )

        # for value in the electrostatic variable, we will find the efefctive area of the collection
        ystart = np.linspace(10, 140, 250)
        zstart = 39 + thc  # start near the top

        vDriftHoles = 0.19e2  # um^2/(V*us)
        dt = 0.002  # us

        boundaryLimits = [0, 150, thc + 0.5, 40 + thc]  # xmin, xmax, ymin, ymax limits

        # iterate over all the different bias potentials
        for i in range(np.amin(phi.shape)):
            qSign = np.zeros(ystart.size)
            E = [
                y,
                z,
                Ey[:, i] / 1e6,
                Ez[:, i] / 1e6,
            ]  # building the E object and converting e fields to V/um

            # Iterate over different y values
            for j in range(ystart.size):
                path = findMotion(
                    (ystart[j], zstart),
                    E,
                    vDriftHoles,
                    dt,
                    q=1.6e-19,
                    limits=boundaryLimits,
                )

                # check if last value falls in limits for electrode A
                for k in range(nElectrodes):
                    if (
                        path[-1, 0] > electrodeABoundary[2 * k]
                        and path[-1, 0] < electrodeABoundary[2 * k + 1]
                    ):
                        qSign[j] = 1

            effectiveArea[findx, i] = sum(qSign) / len(qSign)
            print(effectiveArea)

    plt.rc("font", family="serif")
    fig, ax = plt.subplots()
    bmap = brewer2mpl.get_map("Set2", "Qualitative", 8).mpl_colors
    ax.plot(bias, effectiveArea.T, linewidth=2)
    ax.set_title(
        "Effective Area vs Bias Potential for Different Width Kapton", fontsize=16
    )
    ax.set_xlabel("Bias Potential (V)", fontsize=14)
    ax.set_ylabel("Effective Area", fontsize=14)

    for idx, line in enumerate(ax.get_lines()):
        line.set_color(bmap[idx])

    ax.legend(legend)
    plt.show()

    return None


def darkCurrentCalculation():
    # Uses the resistivity of polymide kapton and efield from comsol to get current density
    # Finds the current by adding the current density through many small areas

    fileloc = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_noHV"
    filelist = os.listdir(fileloc)
    # filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_noHV\kapton_layer_analysis_8um_spacing_nohv.txt'
    kaptonThickness = {
        "500nm": 0.5,
        "1um": 1,
        "2um": 2,
        "5um": 5,
        "10um": 10,
        "15um": 15,
        "3um": 3,
        "8um": 8,
    }
    resistivity = 1.4e21  # ohm*um. From Dupont data sheet
    thck = 5.0
    dy = 1.0  # 1 um dx
    ydim = 150.0
    xdim = 100.0
    nArea = ydim / dy  # number of grids
    bias = np.arange(100, 301, 20)

    # actual detector parameters
    detectorLength = 2e3  # um. Corresponds to 2 mm
    detectorArea = detectorLength ** 2

    print(filelist)

    # Create positions of small area
    ycent = np.arange(0, ydim, dy) + dy / 2  # center position of all the area elements

    current = np.zeros((len(filelist), bias.size))
    legend = []
    for findx, file in enumerate(filelist):
        print(findx)
        # Read data
        _, y, z, data = readComsolFileGrid(fileloc + "\\" + file)
        y, z = y * 1e6, z * 1e6  # convert to um

        # Remove NaN
        data = np.nan_to_num(data, copy=False)

        # Split data into phi, Ey, and Ez
        phi, Ey, Ez = (
            data[:, 0 : np.amin(data.shape) : 3],
            data[:, 1 : np.amin(data.shape) : 3],
            data[:, 2 : np.amin(data.shape) : 3],
        )

        # get the kapton thickness from the filename
        fileIdentifier = file.split("_")[3]
        thck = kaptonThickness[fileIdentifier]
        legend.append(fileIdentifier)
        zcent = np.ones(ycent.size) * thck / 2

        for i in range(bias.size):
            biasCurrent = 0
            E = [y, z, Ey[:, i] / 1e6, Ez[:, i] / 1e6]
            einterp = interpEField2D(ycent, zcent, E)

            for j in range(ycent.size):
                biasCurrent += einterp[j, 1] / resistivity * dy

            current[findx, i] = biasCurrent

    meanCurrentPerUm = current / ydim
    darkCurrent = meanCurrentPerUm * detectorArea
    print(darkCurrent)

    # Making plots of the dark current related data
    bmap = brewer2mpl.get_map("Set2", "Qualitative", 8).mpl_colors
    plt.rc("font", family="serif")
    fig, ax = plt.subplots()

    ax.plot(bias, darkCurrent.T * 1e12, linewidth=3)
    ax.set_title("Dark Current Through the Polymide Kapton", fontsize=16)
    ax.set_xlabel("Bias Voltage (V)", fontsize=16)
    ax.set_ylabel("Current (pA)", fontsize=16)

    # necessary bias voltage index
    biasVoltageIndex = [8, 10, 3, 3, 4, 5, 6]

    for i, line in enumerate(ax.get_lines()):
        line.set_color(bmap[i])
        ax.plot(
            bias[biasVoltageIndex[i]],
            darkCurrent[i, biasVoltageIndex[i]] * 1e12,
            "o",
            markersize=8,
            markeredgecolor=bmap[i],
            markeredgewidth=3,
            markerfacecolor="None",
        )

    ax.legend(legend)
    plt.show()

    return None


if __name__ == "__main__":
    # List of calculations to do. Comment out unnecessary ones to save perfomance

    # effectiveAreaOfDifferentKaptonWidth()
    darkCurrentCalculation()

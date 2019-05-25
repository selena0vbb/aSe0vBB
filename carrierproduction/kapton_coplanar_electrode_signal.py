import particledata as pd
import seleniumconfig as sc
import numpy as np
import matplotlib.pyplot as plt
import os
import brewer2mpl


def signalAnalysis():
    filename = r"C:\Users\alexp\Documents\UW\Research\Selenium\aSe0vBB\particle\selenium-build\output\122_keV_testTuple.root"
    emfilename = [
        r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_1um_spacing_fullsize.txt",
        r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_3um_spacing_fullsize.txt",
        r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_5um_spacing_fullsize.txt",
        r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_8um_spacing_fullsize.txt",
    ]
    configfilename = r"./config.txt"

    print ("Read settings")
    settings = sc.readConfigFile(configfilename)

    print ("Read Geant4 Particle data")
    newCollection = pd.gEventCollection(filename)
    event = newCollection.collection[140]
    creatorProc = event.GetHits()[1]["creatorProcess"].split("\x00")[0]

    print ("Calculating induced charge signal")
    for file in emfilename:
        biasVoltIndex = [0, 1, 2, 3, 4]
        plt.rc("font", family="serif")
        fig, ax = plt.subplots()
        bmap = brewer2mpl.get_map("Set1", "Qualitative", 5).mpl_colors
        biasString = ["120 V", "160 V", "200 V", "240 V", "280 V"]
        for indx in biasVoltIndex:
            print (indx)
            settings["VOLTAGE_SWEEP_INDEX"] = indx
            time, q = pd.computeChargeSignal(event, file, **settings)
            ax.plot(time, -q, linewidth=2, color=bmap[indx], label=biasString[indx])

        ax.set_xlabel(r"Time ($\mu s$)", fontsize=14)
        ax.set_ylabel(r"Induced Charge (C)", fontsize=14)
        ax.set_title(
            "Induced Charge Signal at Amplifier for %s Kapton layer"
            % file.split("\\")[-1].split("_")[3]
        )
        ax.legend()

    return fig, ax


def testScaleFactor():

    filename = r"C:\Users\alexp\Documents\UW\Research\Selenium\aSe0vBB\particle\selenium-build\output\122_keV_testTuple.root"

    emfilename = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_5um_spacing_fullsize.txt"

    gCollection = pd.gEventCollection(filename)
    event = gCollection.collection[146]
    configfilename = r"./config.txt"

    allhits = event.GetHits()
    etot = 0
    for i in allhits:
        etot += i["energy"]

    print (etot)
    print (allhits[1]["y"])

    event.plotH2()

    settings = sc.readConfigFile(configfilename)

    settings["SCALE_WEIGHTED_PHI"] = 1

    t, q = pd.computeChargeSignal(event, emfilename, **settings)

    settings["SCALE_WEIGHTED_PHI"] = 0

    tt, qq = pd.computeChargeSignal(event, emfilename, **settings)

    # fig, ax = plt.subplots()
    # ax.plot(t, -q, 'k', linewidth=3)
    # ax.plot(tt, -qq, '-b', linewidth=3)
    # ax.set_xlabel(r'Time ($\mu s$)', fontsize=14)
    # ax.set_ylabel(r'Induced Charge (C)', fontsize=14)
    # ax.legend(['Scaled', 'Not Scaled'])

    # return fig, ax
    return etot


if __name__ == "__main__":

    # fig, ax = signalAnalysis()

    fig, ax = testScaleFactor()

    plt.show()

import sys
import getopt
import argparse
import numpy as np

sys.path.append("./carrierproduction")
import particledata as pd
import processedPulse as proc
import seleniumconfig as sc


def main(argv):
    """ Main program that acts as an entry point to the simulation """

    # Set default values
    configfile = "./carrierproduction/config.txt"
    nEvents = 0

    parser = argparse.ArgumentParser(
        description="Process command line arguments for Selena simulation"
    )

    # Add arguments to the parsers
    parser.add_argument(
        "-c", "--configfile", nargs="?", default="./carrierproduction/config.txt"
    )
    parser.add_argument("-n", "--number", nargs="?", default=0, type=int)
    parser.add_argument("-v", "--configvar")
    parser.add_argument("configvarValue", nargs="*", default=["None"])

    commandArgs = parser.parse_args()

    # Assign command line variables
    configfile = commandArgs.configfile
    nEvents = commandArgs.number

    # Setup and run the simulation
    settings = sc.readConfigFile(configfile)
    particlefilename = settings["PARTICLE_FILENAME"]
    emfilename = settings["EM_FILENAME"]

    if nEvents <= 0:
        nEvents = pd.numberOfEvents(particlefilename)

    # Print simulation information
    print ("Config file is %s" % configfile)
    print ("Number of events is %i\n" % nEvents)

    simObj = pd.CarrierSimulation(emfilename=emfilename, configfile=configfile)

    # Iterate over the values passed to the config var value
    for val in commandArgs.configvarValue:

        filesize = settings["NEVENTS_PER_FILE"]
        outdir = settings["OUTPUT_DIR"]
        outfilename = settings["OUTPUT_FILE"]

        if val != "None":
            print ("\nConfig variable: " + commandArgs.configvar + " = " + val)
            outfilename = (
                settings["OUTPUT_FILE"].split(".")[0]
                + "_"
                + commandArgs.configvar
                + "_"
                + val
                + "."
                + settings["OUTPUT_FILE"].split(".")[1]
            )

            # Try to parse config value into double if that is the form of the command line arguments
            try:
                val = float(val)
            except ValueError:
                pass
            except Exception as e:
                print (str(e))
                return

            # Change the configvar in the settings and apply
            simObj.settings[commandArgs.configvar] = val
            simObj.settings["OUTPUT_FILE"] = outfilename
            simObj.applySettings()

        # Chunk event collection into smaller pieces if needed
        for j in range(int(np.ceil(float(nEvents) / filesize))):
            indx = []
            # Create new event collect with the smaller chunck size
            newEventCollection = pd.gEventCollection(
                particlefilename,
                eventCounterRange=[j * filesize, (j + 1) * filesize - 1],
            )
            simObj.newEventCollection(newEventCollection)
            for i in range(len(newEventCollection.collection)):
                event = newEventCollection.collection[i]
                zr = np.max(np.abs(event.z))
                yr = np.max(np.abs(event.y))
                xr = np.max(np.abs(event.x))
                eevent = np.sum(event.energy)

                zmean = np.sum(np.array(event.z) * np.array(event.energy)) / np.sum(
                    event.energy
                )

                if (
                    xr < settings["X_MAX"]
                    and yr < settings["Y_MAX"]
                    and zr < settings["Z_MAX"]
                    # and zmean < 0.05 and zmean > 0.01
                    # and eevent > 115
                ):
                    indx.append(i)

            simOutput = simObj.processMultipleEvents(
                indx, processes=int(settings["NPROCESSORS"])
            )
            # simOutput = simObj.processMultipleEvents(
            #     [indx[50]], processes=int(settings["NPROCESSORS"])
            # )

            simOutput.setGitInfo(settings["GIT_DIR"])
            simObj.outputfile = outfilename % j
            simObj.outputdir = outdir
            simObj.savedata(simOutput)

    return


def printhelp():
    print ("python selena.py -c <path to config file> -n <number to simulate>")


if __name__ == "__main__":
    main(sys.argv[1:])

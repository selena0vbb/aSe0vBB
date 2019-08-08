# Module containing classes/functions for handling the Geant4 particle data
# First step for electron hole pair production


try:
    import matplotlib.pyplot as plt
except:
    pass
import numpy as np
import sys
import os
import warnings
import seleniumconfig as sc
from multiprocessing import Pool
from pathlib import Path
import tqdm
import cexprtk
import processedPulse as pout
import scipy.special

# add the EM plot module to the path and import it
sys.path.append(os.path.join(os.path.dirname(sys.path[0]), "EM Analysis"))
sys.path.append(os.path.join(sys.path[0], "EM Analysis"))
from plot import (
    readComsolFileGrid,
    readComsolFileGrid3d,
    findMotion,
    inducedCharge,
    inducedChargeSingle,
    interpEField2D,
    interpEField3D,
)

# Libraries to import
import ROOT as rt


class gEvent(object):
    """
    gEvent is an object to handle G4 Events. It takes the information from a sensitive detector
    and stores it as a series of arrays
    """

    def __init__(self, gEventID=-1, hits=None, geometry=None):
        """ Initializes the event class
        Inputs:
            gEventID - int geant4 event ID number
            hits - list of hits. Each hit is a dictionairy containing the information provided by geant4 for hits (position, energy, track and parent ID, particle type, creator process name)
        """
        # Initialize variables
        self.gEventID = gEventID
        if hits:
            self.trackID = hits["trackID"]
            self.parentID = hits["parentID"]
            self.x = hits["x"]
            self.y = hits["y"]
            self.z = hits["z"]
            self.xi = hits["xi"]
            self.yi = hits["yi"]
            self.zi = hits["zi"]
            self.energy = hits["energy"]
            self.particle = hits["particle"]
            self.creatorProcess = hits["creatorProcess"]
            self.secondaryTrackID = hits["secondaryTrackID"]
        else:
            self.trackID = []
            self.parentID = []
            self.x = []
            self.y = []
            self.z = []
            self.xi = []
            self.yi = []
            self.zi = []
            self.energy = []
            self.particle = []
            self.creatorProcess = []
            self.secondaryTrackID = []

        if geometry:
            self.geomtry = geometry
        else:
            self.geometry = {"x": (-5, 5), "y": (-2, 2), "z": (-0.1, 0.1)}

    def __str__(self):
        return "Geant4 Event. Event ID: %i, Total Energy: %.1f" % (
            self.gEventID,
            np.sum(self.energy),
        )

    def GetSpecificHit(self, index):
        """ Tries to return a single dictionairy hit at index """
        try:
            hit = {}
            hit["trackID"] = self.trackID[index]
            hit["parentID"] = self.parentID[index]
            hit["x"] = self.x[index]
            hit["y"] = self.y[index]
            hit["z"] = self.z[index]
            hit["xi"] = self.xi[index]
            hit["yi"] = self.yi[index]
            hit["zi"] = self.zi[index]
            hit["energy"] = self.energy[index]
            hit["particle"] = self.particle[index]
            hit["creatorProcess"] = self.creatorProcess[index]
            hit["secondaryTrackID"] = self.secondaryTrackID[index]
            return self.hit

        except Exception as e:
            raise e

    # Member functions
    def AddHit(self, hit):
        """ Adds a hit to the collection of hits. Function checks to make sure all entries are present and data types are correct
        Input:
            hit - dictionairy containing hit data
        """

        # Check to make sure all types of data in hit are what is expected
        try:
            boolTypeCheck = []
            boolTypeCheck.append(type(hit["trackID"]) is int)
            boolTypeCheck.append(type(hit["parentID"]) is int)
            boolTypeCheck.append(type(hit["x"]) is float)
            boolTypeCheck.append(type(hit["y"]) is float)
            boolTypeCheck.append(type(hit["z"]) is float)
            boolTypeCheck.append(type(hit["xi"]) is float)
            boolTypeCheck.append(type(hit["yi"]) is float)
            boolTypeCheck.append(type(hit["zi"]) is float)
            boolTypeCheck.append(type(hit["energy"]) is float)
            boolTypeCheck.append(type(hit["particle"]) is str)
            boolTypeCheck.append(type(hit["creatorProcess"]) is str)
            boolTypeCheck.append(type(hit["secondaryTrackID"]) is int)

            if all(boolTypeCheck):
                self.trackID.append(hit["trackID"])
                self.parentID.append(hit["parentID"])
                self.r.append([hit["x"], hit["y"], hit["z"]])
                self.ri.append([hit["xi"], hit["yi"], hit["zi"]])
                self.energy.append(hit["energy"])
                self.particle.append(hit["particle"])
                self.creatorProcess.append(hit["creatorProcess"])
                self.secondaryTrackID.append(hit["secondaryTrackID"])
            else:
                print ("One or more data type in the hit is invalid")

        except Exception as e:
            raise e

    def plotH1(self, x="z", y="energy", nbins=200, figH=None):
        """
        Plots histogram of two compatable variables. Default are energy vs z axis
        Inputs:

        Outputs:
            val - 1d array of the values of each histogram bin
            bins - 1d array of nbins+1 values of the edge of each bin
            ax - handle to the axis
            fig - handle to the plotted figure
        """
        if not figH:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            fig = figH[0]
            ax = figH[1]

        histX = getattr(self, x)
        histY = getattr(self, y)

        val, bins, _ = ax.hist(
            histX,
            bins=nbins,
            range=self.geometry[x],
            weights=histY,
            log=True,
            histtype="step",
        )

        # Set some default axis
        ax.set_xlabel(x + " (mm)", fontsize=14)
        ax.set_ylabel(y + " (keV)", fontsize=14)
        ax.set_title(
            y + " vs " + x + " for Event ID %i " % self.GetEventID(), fontsize=16
        )

        return val, bins, ax, fig

    def plotH2(
        self, x="y", y="z", z="energy", nbins=[200, 200], figH=None, delete=False
    ):

        histX = getattr(self, x)
        histY = getattr(self, y)
        histZ = getattr(self, z)

        # setting up the range
        histRange = [self.geometry[x], self.geometry[y]]

        if delete:
            val, binx, biny = np.histogram2d(
                histX, histY, range=histRange, bins=nbins, weights=histZ
            )
            fig, ax = None, None

        else:
            if not figH:
                fig = plt.figure()
                ax = fig.add_subplot(111)
            else:
                fig = figH[0]
                ax = figH[1]

            val, binx, biny, cax = ax.hist2d(
                histX, histY, range=histRange, bins=nbins, weights=histZ
            )

            # Set some default axis
            ax.set_xlabel(x + " (mm)", fontsize=14)
            ax.set_ylabel(y + " (mm)", fontsize=14)
            ax.set_title(
                z + " vs " + x + " and " + y + " for Event ID %i" % self.GetEventID(),
                fontsize=16,
            )

            fig.colorbar(cax)

        return val, [binx, biny], ax, fig

    def plotHN(self, varname=("x", "y", "z"), fname="energy", testbins=[200, 200, 200]):
        """
        Bins a function into N different axis. By default, this function bins energy into x,y,z axis
        """

        # Extract the variables and function from the flatten event data
        var = [getattr(self, v) for v in varname]
        func = getattr(self, fname)

        # Convert the list of variables into a NxD array (N number of entries, D dimensions)
        var = np.array(var).T

        # histogram the data
        val, edges = np.histogramdd(var, bins=testbins, weights=func)

        # Compute bin centers from the edges
        binCenters = [edge[:-1] + np.diff(edge) / 2 for edge in edges]

        return val, binCenters

    def createCarriersBin(self, wehpFunction=None, efield=None, **kwargs):
        """
        Function that creates carriers from the energy distribution of the incoming gamma rays. Standard is to create in 3D
        Inputs:
            Wehp - work function of amorphous selenium
            kwargs - key word arguments from the config file of the settings file
        Outputs
            nehp - the nonzero number of electron hole pairs created in spacce associated with the xv, yv, zv indicess
            nehpFluctuations - nonzero nehp fluctuations
            binCenters - list of the bincenters along each (x,y,z) axis
            [xv, yv, zv] - indexes in the nehp matrix where the values are non zero
        """

        # Create the bins for the histogramming of energy events
        binnames = ["X", "Y", "Z"]
        bins = []
        for name in binnames:
            axisMinValue = kwargs["%s_MIN" % name]
            axisMaxValue = kwargs["%s_MAX" % name]
            dx = kwargs["D%s" % name]
            bins.append(np.arange(axisMinValue, axisMaxValue + dx, dx))

        val, binCenters = self.plotHN(
            varname=("x", "y", "z"), fname="energy", testbins=bins
        )

        # Find the non zero bin entries
        if not kwargs["THREE_DIMENSIONS"]:
            # If working in 2D projection, collapse the x axis
            val = np.sum(val, axis=0, keepdims=True)

        # x, y, z arrays of the index of the non zero elements of the nehpf array
        xv, yv, zv = np.nonzero(val)

        # Find the work function with E Field dependence
        if wehpFunction and efield:
            wehp = 0.05 * np.ones(val.shape)
            # Interpolate E field on non zero bin centers
            if not kwargs["THREE_DIMENSIONS"]:
                Eyz = interpEField2D(binCenters[1][yv], binCenters[2][zv], efield)
                EMag = np.sqrt(np.sum(Eyz ** 2, axis=1))
            else:
                Exyz = interpEField3D(
                    binCenters[0][xv], binCenters[1][yv], binCenters[2][zv], efield
                )
                EMag = np.sqrt(np.sum(Exyz ** 2, axis=1))

            # Iterate over the points with non zero nehp and find the work function at that position
            for i in range(xv.size):
                wehpFunction.symbol_table.variables["E"] = EMag[i]
                wehp[xv[i], yv[i], zv[i]] = wehpFunction()
        else:
            # Use standard value of 50 eV if no E field information is passed\
            wehp = 0.05 * np.ones(val.shape)

        # histogrammed data are in a (nbinx x nbiny x nbinz) array. Simply divide by Wehp to get the energy as a distribution of position
        nehp = val / wehp
        if kwargs["CARRIER_GENERATION_POISSON"]:
            nehpFluctuation = np.random.poisson(nehp)
        else:
            nehpFluctuation = np.round(nehp)

        nehp = np.round(nehp)

        # Return
        return (
            nehp[xv, yv, zv],
            nehpFluctuation[xv, yv, zv],
            binCenters,
            xv,
            yv,
            zv,
            wehp[xv, yv, zv],
            val[xv, yv, zv],
        )

    def createCarriersTrack(
        self, useColumnarRecomb, wehpFunction=None, efield=None, **settings
    ):

        # Compute the work function with the specified technique (columnar or functional)
        if settings["COLUMNAR_RECOMBINATION"]:
            wehp, _ = CarrierSimulation.columnarRecombinationWorkFunction(
                self, efield, **settings
            )
        else:
            wehp = 0.05 * np.ones(len(self.energy))

            # If have efield and wehp function compute the functional form
            if wehpFunction and efield:
                if settings["THREE_DIMENSIONS"]:
                    Exyz = interpEField3D(self.x, self.y, self.z, efield)
                    EMag = np.sqrt(np.sum(Exyz ** 2, axis=1))
                else:
                    Eyz = interpEField2D(self.y, self.z, efield)
                    EMag = np.sqrt(np.sum(Eyz ** 2, axis=1))

                for i in range(len(self.energy)):
                    wehpFunction.symbol_table.variables["E"] = EMag[i]
                    wehp[i] = wehpFunction()

        # Convert energy to charge using wehp. Include statistical fluctuations
        nehp = self.energy / wehp

        # Apply fluctuations
        if settings["CARRIER_GENERATION_POISSON"]:
            nehpFluctuation = np.random.poisson(nehp)
        else:
            nehpFluctuation = np.round(nehp)

        nehp = np.round(nehp)

        return (nehp, nehpFluctuation, self.x, self.y, self.z, wehp, self.energy)


class gEventCollection(object):
    """Container for reading, parsing and storing many Geant4 events"""

    def __init__(self, rootFilename=None, eventCounterRange=None, **kwargs):

        self.rootFilename = rootFilename
        self.collection = []

        if rootFilename:
            # Read the data from the file
            f = rt.TFile(rootFilename)
            # Gets a tree with the name aSeData. Current name from Geant4 simulation.
            tree = f.Get("aSeData")
            eventID = -1
            eventCounter = 0

            # Iterate over all of the entries in the tree, extracting tuple i
            for entry in tree:
                if eventID != entry.EventID:
                    if eventID != -1:
                        if eventCounterRange != None:
                            if eventCounter >= eventCounterRange[0]:
                                self.collection.append(
                                    gEvent(gEventID=eventID, hits=hitsInfo)
                                )
                        else:
                            self.collection.append(
                                gEvent(gEventID=eventID, hits=hitsInfo)
                            )
                        eventCounter += 1

                    hitsInfo = {}
                    eventID = entry.EventID
                    hitsInfo["trackID"] = [entry.TrackID]
                    hitsInfo["parentID"] = [entry.ParentID]
                    hitsInfo["x"] = [entry.x]
                    hitsInfo["y"] = [entry.y]
                    hitsInfo["z"] = [entry.z]
                    hitsInfo["energy"] = [entry.energy * 1e3]
                    hitsInfo["particle"] = [entry.ParticleType]
                    hitsInfo["creatorProcess"] = [entry.ProcessName]

                    # Add additional parameters to the hit. For backwords compatibility
                    try:
                        hitsInfo["secondaryTrackID"] = [entry.SecondaryTrackID]
                        hitsInfo["xi"] = [entry.xi]
                        hitsInfo["yi"] = [entry.yi]
                        hitsInfo["zi"] = [entry.zi]
                    except AttributeError:
                        hitsInfo["secondaryTrackID"] = [-1]
                        hitsInfo["xi"] = [0]
                        hitsInfo["yi"] = [0]
                        hitsInfo["zi"] = [0]  # Default values. No information
                    except:
                        pass

                else:
                    hitsInfo["trackID"].append(entry.TrackID)
                    hitsInfo["parentID"].append(entry.ParentID)
                    hitsInfo["x"].append(entry.x)
                    hitsInfo["y"].append(entry.y)
                    hitsInfo["z"].append(entry.z)
                    hitsInfo["energy"].append(entry.energy * 1e3)
                    hitsInfo["particle"].append(entry.ParticleType)
                    hitsInfo["creatorProcess"].append(entry.ProcessName)

                    # Add additional parameters to the hit. For backwords compatibility
                    try:
                        hitsInfo["secondaryTrackID"].append(entry.SecondaryTrackID)
                        hitsInfo["xi"].append(entry.xi)
                        hitsInfo["yi"].append(entry.yi)
                        hitsInfo["zi"].append(entry.zi)
                    except AttributeError:
                        hitsInfo["secondaryTrackID"].append(-1)
                        hitsInfo["xi"].append(0)
                        hitsInfo["yi"].append(0)
                        hitsInfo["zi"].append(0)  # Default values. No information
                    except:
                        pass

                if eventCounterRange != None and eventCounter >= eventCounterRange[1]:
                    break

            self.totalNEvents = len(self.collection)
        else:
            self.totalNEvents = 0

    def __str__(self):
        return "Geant4 Event Collection. N=%i" % len(self.collection)

    def printInfo(self):
        """ Prints information regarding the collection of events"""
        print ("Particle data comes from " + self.rootFilename)
        print ("There are %i events stored in this collection" % (len(self.collection)))

    def findEvent(self, eventID):
        """Searches the collection for an event with eventID. Returns if finds it"""
        foundEvent = None

        for event in self.collection:
            if event.gEventID == eventID:
                return event

        if not foundEvent:
            print ("Could not find a event matching the ID %i" % eventID)


class CarrierSimulation(object):
    """Wrapper class for all the functions and data surrounding the charge carrier simulation. All info is housed in one class for efficiency purposes"""

    def __init__(self, emfilename=None, eventCollection=None, configfile=None):
        """
        Initialize the class
        """
        self.emfilname = emfilename
        self.configfile = configfile

        # Default work function of 0.05 keV/ehp
        self.symbolTable = cexprtk.Symbol_Table({"E": 10}, add_constants=True)
        self.wehpExpression = cexprtk.Expression("0.05", self.symbolTable)

        # Read and store settings
        print ("Read config")
        if configfile:
            self.settings = sc.readConfigFile(configfile)
            self.applySettings()
        else:
            self.settings = None
            warnings.warn(
                "No Settings File Read. Please read (or create) class settings and run applySettings()."
            )
            return

        # Read and store particle data
        print ("Read Geant4 Simulation Data")
        if type(eventCollection) is str:
            self.eventCollection = gEventCollection(eventCollection)
        else:
            self.eventCollection = eventCollection

        # Read and store emdata
        print ("Read Comsol EM Data")
        if emfilename:
            self.newEmFile(emfilename)
        else:
            self.x = self.y = self.z = self.data = []

        if self.settings["SCALE_WEIGHTED_PHI"]:
            self.computeScaleFactor()
        else:
            self.scale = 1

    def newEmFile(self, emfilename):
        """
        Replace the em data with what is in the new file
        Inputs:
            emfilename - string of the path to the comsol simulation file
        """
        if self.use3D:
            _, pos, self.data = readComsolFileGrid3d(emfilename)
            self.x, self.y, self.z = pos[0], pos[1], pos[2]
        else:
            _, self.x, self.y, self.data = readComsolFileGrid(emfilename)
            self.z = np.array([])
        self.createFields()

    def newEventCollection(self, eventCollection):
        """
        Replaces the current event collection with the new eventCollection
        Inputs:
            eventCollection - either a string to the file location or an gEventCollection ojbect
        """
        if type(eventCollection) is str:
            self.eventCollection = gEventCollection(eventCollection)
        else:
            self.eventCollection = eventCollection

    def newSettings(self, configfile):

        self.settings = sc.readConfigFile(configfile)

    def applySettings(self):
        if self.settings:
            print ("Apply Relevant Settings")

            self.outputdir = self.settings["OUTPUT_DIR"]
            self.outputfile = self.settings["OUTPUT_FILE"]
            self.use3D = self.settings["THREE_DIMENSIONS"]

            # Create the function to calculate nehp
            self.symbolTable = cexprtk.Symbol_Table({"E": 10}, add_constants=True)
            self.wehpExpression = cexprtk.Expression(
                self.settings["WORK_FUNCTION"], self.symbolTable
            )

        else:
            warnings.warn("No settings to apply. Returning without applying settings.")

    def computeScaleFactor(self, depth=[0.5, 1], width=[0.35, 0.65]):
        """
        Computes the scale factor between weighted potentials for tiered coplanar design. Save it as an attribute of the class so it only needs to be done once
        """

        # Get dimensions and sets the limits that we want to take
        xdim = self.x.size
        ydim = self.y.size
        xmin, xmax = int(xdim * width[0]), int(xdim * width[1])
        ymin, ymax = int(ydim * depth[0]), int(ydim * depth[1])

        phi = self.data[:, 0 : np.amin(self.data.shape) : 3]
        wPhiA = np.reshape(phi[:, -2], (ydim, xdim))
        wPhiB = np.reshape(phi[:, -1], (ydim, xdim))

        phiDiffA = np.diff(wPhiA[ymin : ymax - 1, xmin : xmax - 1], axis=0)
        phiDiffB = np.diff(wPhiB[ymin : ymax - 1, xmin : xmax - 1], axis=0)

        self.scale = np.mean(phiDiffA / phiDiffB)

        return self.scale

    def createFields(self):
        """
        Takes the raw comsol data and parses it into the form used by findmotion and inducedCharge
        """

        # Gets the EM information
        x, y, z, data = self.x, self.y, self.z, self.data
        convFactor = self.settings["UNIT_CONVERSION_FACTOR"]
        x *= convFactor
        y *= convFactor
        z *= convFactor  # convert m to mm
        fieldScale = self.settings["FIELD_SCALE"]

        # Parses the data depending on the mode. Parse 2D and 3D simulations differently.
        if self.use3D:
            phiAll, ExAll, EyAll, EzAll = (
                data[:, 0 : np.amin(data.shape) : 4],
                data[:, 1 : np.amin(data.shape) : 4],
                data[:, 2 : np.amin(data.shape) : 4],
                data[:, 3 : np.amin(data.shape) : 4],
            )
            if self.settings["CHARGE_DIFFERENCE"]:
                wPhiA = [x, y, z, phiAll[:, -2]]
                wPhiB = [x, y, z, phiAll[:, -1]]
                wPhi = [wPhiA, wPhiB]
            else:
                wPhi = [x, y, z, phiAll[:, -1]]

            Ex, Ey, Ez = (
                ExAll[:, self.settings["VOLTAGE_SWEEP_INDEX"]] / convFactor,
                EyAll[:, self.settings["VOLTAGE_SWEEP_INDEX"]] / convFactor,
                EzAll[:, self.settings["VOLTAGE_SWEEP_INDEX"]] / convFactor,
            )  # in V/mm

            # Creates an instance variable of the E field and wPhi
            self.Etot = [x, y, z, Ex * fieldScale, Ey * fieldScale, Ez * fieldScale]
            self.wPhi = wPhi

        else:
            phiAll, ExAll, EyAll = (
                data[:, 0 : np.amin(data.shape) : 3],
                data[:, 1 : np.amin(data.shape) : 3],
                data[:, 2 : np.amin(data.shape) : 3],
            )
            if self.settings["CHARGE_DIFFERENCE"]:
                wPhiA = [x, y, phiAll[:, -2]]
                wPhiB = [x, y, phiAll[:, -1]]
                wPhi = [wPhiA, wPhiB]
            else:
                wPhi = [x, y, phiAll[:, -1]]

            Ex, Ey = (
                ExAll[:, self.settings["VOLTAGE_SWEEP_INDEX"]] / convFactor,
                EyAll[:, self.settings["VOLTAGE_SWEEP_INDEX"]] / convFactor,
            )  # in V/mm

            # Creates an instance variable of the E field and wPhi
            self.Etot = [x, y, Ex * fieldScale, Ey * fieldScale]
            self.wPhi = wPhi

        return None

    def chargeArray(self, nehp, nHoleSteps, nElectronSteps):
        """
        Generates the charge arrays used in the induced charge calculations
        """

        # Get data from the settings
        tauHole = self.settings["TAU_HOLES"]
        tauElectron = self.settings["TAU_ELECTRONS"]
        dtHole = self.settings["DT_HOLES"]
        dtElectron = self.settings["DT_ELECTRONS"]
        elementaryCharge = self.settings["ELEMENTARY_CHARGE"]

        qHole = elementaryCharge * nehp
        qElectron = -qHole

        # Create array for number of charges particles at each point of the path
        qHoleArray = qHole * np.ones(nHoleSteps)
        qElectronArray = qElectron * np.ones(nElectronSteps)

        # If carrier lifetime noise is considered, do a geometric distribution to find the number of steps until a success. Iterate over that an substract
        if self.settings["CARRIER_LIFETIME_GEOMETRIC"]:
            holeTrapTimeIndx = np.random.geometric(dtHole / tauHole, nehp)
            electronTrapTimeIndx = np.random.geometric(dtElectron / tauElectron, nehp)

            # Eliminate values that are greater than the size of the length of the path
            holeTrapTimeIndx = holeTrapTimeIndx[holeTrapTimeIndx < qHoleArray.size]
            electronTrapTimeIndx = electronTrapTimeIndx[
                electronTrapTimeIndx < qElectronArray.size
            ]

            # iterate over all the trap location. Subtract a charge from that location in the charge array
            for k in range(holeTrapTimeIndx.size):
                qHoleArray[holeTrapTimeIndx[k] :] -= elementaryCharge
            for k in range(electronTrapTimeIndx.size):
                qElectronArray[electronTrapTimeIndx[k] :] += elementaryCharge

        return qHoleArray, qElectronArray

    @staticmethod
    def columnarRecombinationWorkFunction(event, Efield, **settings):
        """
        Computes the charge yield after applying a columnar recombination model on the data following Ausman and Mclean. https://doi.org/10.1063/1.88104

        Inputs:
            event - gEvent object
            energy - Nx1 numpy array of the energy in each spatial bin assosciated with the hits
        Outputs:
            nehp - Nx1 numpy array of the effective number of electron hole pairs produced in each binned energy location (N bins)
        """

        # Get parameters from the settings file
        b = settings["TRACK_RADIUS"]
        N0 = settings["LINEAR_IONIZATION_DENSITY"]
        w0 = settings["SATURATION_WF"]

        # Compute diffusion and and recombination coefficients
        voltConvCoeff = 1.0e-6  # Conversion to convert volts from SI to mm and us
        mu = settings["MU_HOLES"] + settings["MU_ELECTRONS"]
        kb = 1.38e-29  # in mm^2 kg s^-2 K^-1
        T = 293
        q = settings["ELEMENTARY_CHARGE"]
        D = mu * kb * T / q

        eps0 = 8.85e-15 / voltConvCoeff  # in C V^-1 mm^-1
        epsr = 6
        alpha = (q * mu) / (eps0 * epsr)

        # Calculate each step of the track
        pos = ["x", "y", "z"]
        drTrack = np.array(
            [
                np.array(getattr(event, p)) - np.array(getattr(event, p + "i"))
                for p in pos
            ]
        ).T

        # Calculate the length of each track segment
        trackLength = np.sqrt( np.sum( drTrack**2, axis=1))

        # Find the mean E-field over the extent of the track
        rAve = [
            np.sum(np.array(getattr(event, i)) * np.array(event.energy))
            / np.sum(event.energy)
            for i in pos
        ]

        if settings["THREE_DIMENSIONS"]:
            efieldAve = interpEField3D(rAve[0], rAve[1], rAve[2], Efield)
            dotProduct = np.dot(drTrack, np.squeeze(efieldAve))
            drTrackMag = np.sqrt(np.sum(drTrack ** 2, axis=1))
        else:
            efieldAve = interpEField2D(rAve[1], rAve[2], Efield)
            dotProduct = np.dot(drTrack[:, 1:], np.squeeze(efieldAve))
            drTrackMag = np.sqrt(np.sum(drTrack[:, 1:] ** 2, axis=1))

        sinThetaSquare = (
            np.sin(
                np.arccos(dotProduct / (drTrackMag * np.sqrt(np.sum(efieldAve ** 2))))
            )
            ** 2
        )

        # Compute x, defined in paper
        x = (
            mu ** 2
            * sinThetaSquare
            * b ** 2
            * np.sum(
                (efieldAve * voltConvCoeff) ** 2
            )  # Convert E field units to be in mm and us
            / (4 * D ** 2)
        )

        # N0 is the linear charge density. Assume that energy is lost linearly along the length of the track.
        N0 = event.energy / (w0 * trackLength)
        wehp = w0 * (
            1 + (alpha * N0 / (4 * np.pi * D)) * np.exp(x) * scipy.special.kn(0, x)
        )

        return wehp, sinThetaSquare

    def computeChargeSignal(self, eventIdx):
        """
        Uses the number of electron hole pairs spatial distribution from the event and the electrostatic simulation from COMSOL to determine the induced charge signal on the detector plates

        Inputs:
            event - filled gEvent object
            emFilename - string to the path containing the comsol electromagnetic data for the detector geometry
            **kwargs - a variety of keyword arguments that govern simulation settings. See the config.txt for more information about what each one does
        Outputs:
            time - Nx1 numpy array where N is the number of time steps for the last charge to reach the electrodes
            qInduced - Nx1 numpy array of the induced charge signal at the electrodes. For single electrode outputs just induced charge, for coplanar design (id 'CHARGE_DIFFERENCE' == True) outputs the difference in induced charge

        """

        # Creates a simulated pulse object to store information
        simPulseObj = pout.SimulatedPulse(**self.settings)

        # Gets the carrier information
        event = self.eventCollection.collection[eventIdx]
        # Build charge carriers using the binning technique
        if self.settings["CARRIER_SCHEMA"] == "bin":
            nehp, nehpf, bins, xv, yv, zv, wehp, energy = event.createCarriersBin(
                self.wehpExpression, self.Etot, **self.settings
            )
            xpos, ypos, zpos = bins[0][xv], bins[1][yv], bins[2][zv]
        # Build charge carriers using the track segment technique
        elif self.settings["CARRIER_SCHEMA"] == "track":
            useColumnarRecomb = self.settings
            nehp, nehpf, xarr, yarr, zarr, wehp, energy = event.createCarriersTrack(
                useColumnarRecomb, self.wehpExpression, self.Etot, **self.settings
            )
            xpos, ypos, zpos = np.array(xarr), np.array(yarr), np.array(zarr)
        else:
            return 0

        simPulseObj.setG4Event(event)
        simPulseObj.setBinnedPosition(xpos, ypos, zpos)
        simPulseObj.setBinnedEnergy(energy)
        simPulseObj.setNehp(nehp, nehpf)
        simPulseObj.setWehp(wehp)
        simPulseObj.setScaleFactor(self.scale)

        # Set up the boundaries of the em simulation
        if self.settings["USE_BOUNDARY_LIMITS"]:
            if self.use3D:
                limits = [
                    self.settings["X_MIN"],
                    self.settings["X_MAX"],
                    self.settings["Y_MIN"],
                    self.settings["Y_MAX"],
                    self.settings["Z_MIN"],
                    self.settings["Z_MAX"],
                ]
            else:
                limits = [
                    self.settings["Y_MIN"],
                    self.settings["Y_MAX"],
                    self.settings["Z_MIN"],
                    self.settings["Z_MAX"],
                ]
        else:
            limits = []

        # Get parameters from the settings necessary for the simulations
        dtHole = self.settings["DT_HOLES"]  # us
        dtElectron = self.settings["DT_ELECTRONS"]  # us
        muHole, muElectron = (
            self.settings["MU_HOLES"],
            self.settings["MU_ELECTRONS"],
        )  # mm^2/(V*us)
        totalTimeHoles, totalTimeElectrons = (
            self.settings["TOTAL_TIME_HOLES"],
            self.settings["TOTAL_TIME_ELECTRONS"],
        )  # us
        maxtime = 0
        allInduced = []
        elecAInduced = []
        elecBInduced = []

        # Find the number of points in the nehp array that are nonzero. Iterate over those
        for i in range(xpos.size):

            # Find trajectory of this group of particles
            # Assign the sign of the charge of holes and electrons
            qHoles = 1
            qElectrons = -qHoles

            # Find path of the electrons and holes given the location by (xi), yi, zi
            if self.use3D:
                pathHoles = findMotion(
                    (xpos[i], ypos[i], zpos[i]),
                    self.Etot,
                    muHole,
                    dtHole,
                    totalTime=totalTimeHoles,
                    q=qHoles,
                    limits=limits,
                )
                pathElectrons = findMotion(
                    (xpos[i], ypos[i], zpos[i]),
                    self.Etot,
                    muElectron,
                    dtElectron,
                    totalTime=totalTimeElectrons,
                    q=qElectrons,
                    limits=limits,
                )
            else:
                pathHoles = findMotion(
                    (ypos[i], zpos[i]),
                    self.Etot,
                    muHole,
                    dtHole,
                    totalTime=totalTimeHoles,
                    q=qHoles,
                    limits=limits,
                )
                pathElectrons = findMotion(
                    (ypos[i], zpos[i]),
                    self.Etot,
                    muElectron,
                    dtElectron,
                    totalTime=totalTimeElectrons,
                    q=qElectrons,
                    limits=limits,
                )

            qHoleArray, qElectronArray = self.chargeArray(
                nehpf[i], pathHoles.shape[0], pathElectrons.shape[0]
            )

            # Keep the longest time for reference to create the total pulse at the end
            try:
                if np.max([pathHoles[-1, -1], pathElectrons[-1, -1]]) > maxtime:
                    maxtime = np.max([pathHoles[-1, -1], pathElectrons[-1, -1]])
            except:
                pass

            if self.settings["CHARGE_DIFFERENCE"]:
                # Compute difference in induced charge, qA-qB, for both electrons and holes. Incorporates scale factor (which is =1 for no scaling)
                qAHoles, qBHoles, _ = inducedCharge(
                    self.wPhi[0],
                    self.wPhi[1],
                    pathHoles,
                    q=qHoleArray,
                    roundFinalVal=self.settings["ROUND_FINAL_WEIGHTED_PHI"],
                )
                qAElectrons, qBElectrons, _ = inducedCharge(
                    self.wPhi[0],
                    self.wPhi[1],
                    pathElectrons,
                    q=qElectronArray,
                    roundFinalVal=self.settings["ROUND_FINAL_WEIGHTED_PHI"],
                )
                indChargeHoles = qAHoles - self.scale * qBHoles
                indChargeElectrons = qAElectrons - self.scale * qBElectrons

                # append to  list of induced charge
                allInduced.append(indChargeHoles)
                allInduced.append(indChargeElectrons)
                if self.settings["SAVE_PARTIAL_PULSES"]:
                    elecAInduced.append(qAHoles)
                    elecAInduced.append(qAElectrons)
                    elecBInduced.append(self.scale * qBHoles)
                    elecBInduced.append(self.scale * qBElectrons)

            else:
                # Compute induced charge on a single electrode
                indChargeHoles = inducedChargeSingle(self.wPhi, pathHoles, q=qHoleArray)
                indChargeElectrons = inducedChargeSingle(
                    self.wPhi, pathElectrons, q=qElectronArray
                )
                allInduced.append(indChargeHoles)
                allInduced.append(indChargeElectrons)

        timeHoles = np.arange(0, maxtime + dtHole, dtHole)
        timeElectrons = np.arange(0, maxtime + dtElectron, dtElectron)
        qIndHole = np.zeros(timeHoles.shape)
        qIndElectron = np.zeros(timeElectrons.shape)

        # Iterate over the induced charge by each grid contribution and add them to the total charge induced on the electrode
        for indx, charge in enumerate(allInduced):
            if indx % 2 == 0:
                qIndHole[: len(charge)] += charge
                # fill the rest of the time values with the last value so that the charge doesn't "disapear"
                qIndHole[len(charge) :] += charge[-1]
            else:
                qIndElectron[: len(charge)] += charge
                qIndElectron[len(charge) :] += charge[-1]

        # Interpolate the Electron induced charge to match the holes so we can add them together
        qInduced = qIndHole + np.interp(timeHoles, timeElectrons, qIndElectron)

        # Create total induced charge for electrodes A and B (if setting is selected)
        if self.settings["SAVE_PARTIAL_PULSES"] and self.settings["CHARGE_DIFFERENCE"]:
            qIndHoleA = np.zeros(timeHoles.shape)
            qIndElectronA = np.zeros(timeElectrons.shape)
            qIndHoleB = np.zeros(timeHoles.shape)
            qIndElectronB = np.zeros(timeElectrons.shape)

            for partialIndx in range(len(elecAInduced)):
                if partialIndx % 2 == 0:
                    qIndHoleA[: len(elecAInduced[partialIndx])] += elecAInduced[
                        partialIndx
                    ]
                    qIndHoleA[len(elecAInduced[partialIndx]) :] += elecAInduced[
                        partialIndx
                    ][-1]
                    qIndHoleB[: len(elecBInduced[partialIndx])] += elecBInduced[
                        partialIndx
                    ]
                    qIndHoleB[len(elecBInduced[partialIndx]) :] += elecBInduced[
                        partialIndx
                    ][-1]
                else:
                    qIndElectronA[: len(elecAInduced[partialIndx])] += elecAInduced[
                        partialIndx
                    ]
                    qIndElectronA[len(elecAInduced[partialIndx]) :] += elecAInduced[
                        partialIndx
                    ][-1]
                    qIndElectronB[: len(elecBInduced[partialIndx])] += elecBInduced[
                        partialIndx
                    ]
                    qIndElectronB[len(elecBInduced[partialIndx]) :] += elecBInduced[
                        partialIndx
                    ][-1]

            qInducedA = qIndHoleA + np.interp(timeHoles, timeElectrons, qIndElectronA)
            qInducedB = qIndHoleB + np.interp(timeHoles, timeElectrons, qIndElectronB)

            if self.settings["FAST"]:
                signal = np.array([qInduced, qInducedA, qInducedB]).T
                simPulseObj.setSignalMaximum(
                    np.max(np.abs(signal), axis=0) / self.settings["ELEMENTARY_CHARGE"]
                )
            else:
                simPulseObj.setTimeSeries(
                    timeHoles, np.array([qInduced, qInducedA, qInducedB]).T, dim=3
                )
            return simPulseObj

        if self.settings["FAST"]:
            simPulseObj.setSignalMaximum(
                np.max(np.abs(signal)) / self.settings["ELEMENTARY_CHARGE"]
            )
        else:
            simPulseObj.setTimeSeries(timeHoles, qInduced)
        return simPulseObj

    def processMultipleEvents(self, eventIdxs, processes=1, chunksize=8):
        """
        Tool for computing multiple Geant for events with one function call. Hass option for parallel computing as well

        Inputs:
            eventIdxs - a list of event IDs to iterate over
            processes - number of processors to uses. Only relevant when PARALLEL_PROCESS setting is true
            chunksize - number of chunks to break the eventIdx into and send to each worker at a time.

        Outputs:
            simOutputObj - simulation output object
        """

        simOutputObj = pout.SimulatedOutputFile(settings=self.settings)
        if self.settings["PARALLEL_PROCESS"]:
            pool = Pool(processes=processes)
            for simPulse in tqdm.tqdm(
                pool.imap_unordered(
                    worker, ((self, event) for event in eventIdxs), chunksize=chunksize
                ),
                total=len(eventIdxs),
            ):
                simOutputObj.addPulses(simPulse)
            pool.close()
            pool.join()
        else:
            for event in tqdm.tqdm(eventIdxs):
                simPulse = self.computeChargeSignal(event)
                simOutputObj.addPulses(simPulse)

        return simOutputObj

    def savedata(self, obj):
        outdir = self.outputdir
        outfile = self.outputfile

        outpath = os.path.join(outdir, outfile)

        # Serialize obj/data using pick
        pout.savePickleObject(obj, outpath)

        return


def worker(args):
    """
    Defines a worker function for parallel processing since the Pool.map() function takes only one argument. args is a tuple of different parameters that are needed to run the simulation
    """
    obj, indx = args
    return obj.computeChargeSignal(indx)


def computeDarkCurrentNoise(E, binx, biny, rho, dt):
    """
    Calculates the dark current through a specific area of the detector (deliminated by the binx, biny range), converts that dark current to number of charge particles (probabalistically) and finds the induced charge from them at a time step
    """
    return None


def scaleWeightPhi(wPhiA, wPhiB, dimSize, depth=0.5, xrange=[0.35, 0.65]):
    """
    Finds the scaling factor between weighted potential required to make the the different in induced charge 0 in the bulk

    Inputs:
        wPhiA - NxM np array of the weighted potential of electrode A (collection)
        wPhiB - NxM array of weighted potential for B
        dimSize = array of [nXelements, nYelements]
        depth - what fraction of the weighted potential to look at. Must be small enough to avoid non-linearities of weighted potential near coplanar electrode
        xrange - fraction of x values to take. Typically the edges of detector have distortions from linear potential, therefore we want to exclude
    Outpus:
        scaleFactor - float factor between the linear components of the weighted potential
    """

    # Reshape weighted potential into nx x ny array
    wPhiA = np.reshape(wPhiA, (dimSize[1], dimSize[0]))
    wPhiB = np.reshape(wPhiB, (dimSize[1], dimSize[0]))
    zdim = wPhiA.shape[0]
    zindx = int(depth * zdim)
    xmin, xmax = int(xrange[0] * wPhiA.shape[1]), int(xrange[1] * wPhiA.shape[1])

    # Calculates the slope of the weighted potential in the linear region. Then finds the ratio at every point and averages them
    phiDiffA = np.diff(wPhiA[zindx:, xmin:xmax], axis=0)
    phiDiffB = np.diff(wPhiB[zindx:, xmin:xmax], axis=0)

    return np.mean(phiDiffA / phiDiffB)


def readFreqResponse(freqResponseFile):
    """
    Reads the LTSpice frequency response. Outputs frequency and gain in dB
    """
    f = open(freqResponseFile)

    fcircuit = []
    gaindB = []

    # Clear header
    f.readline()

    for line in f:
        # Split line at tabs to get all the elements in the line
        lineElement = line.split("\t")

        # Append frequency and gain
        fcircuit.append(float(lineElement[0]))
        gaindB.append(float(lineElement[1].split("dB")[0][1:]))

    fcircuit = np.array(fcircuit)
    gaindB = np.array(gaindB)

    f.close()

    return fcircuit, gaindB


def filterSignal(signal, time, freqResponseFile):
    """
    Filters the signal given the frequency response of the circuit provided in a .txt file
    Inputs:
        signal - Nx1 dimensional numpy array containing the time domain singal
        time - Nx1 dimensional numpy array of time associated with the signale
        freqResponseFile - spice frequency response file name for the amplification circuit
    """

    fcircuit, gain = readFreqResponse(freqResponseFile)

    dt = time[1] - time[0]
    N = signal.size

    # Compute fft of charge time domain signal
    yf = np.fft.fft(signal)
    xf = np.fft.fftfreq(N, dt)

    # interpolate gain
    gainInterpolate = np.interp(xf[: N // 2], fcircuit, 10 ** (gain / 10))

    # Finds the filtered signal in the frequency domain. Multiplication in frequency domain
    signalFilterFreq = yf * gainInterpolate

    # Inverse fourrier transform
    signalFilter = np.fft.ifft(signalFilterFreq)

    return signalFilter, time, signalFilterFreq, xf, yf


def numberOfEvents(rootFile, key="EneDepSe"):
    """ Reads the root file and calculates the number of entries in the file"""
    f = rt.TFile(rootFile)
    obj = f.Get(key)

    return obj.GetEntries()


if __name__ == "__main__":
    # used for debuggin information. If the particledata.py file is run this segment of the code will run
    filename = r"C:\Users\alexp\Documents\UW\Research\Selenium\aSe0vBB\particle\selenium-build\output\pixel_sio2_122kev_0degangle_200k.root"
    filenameAngle = r"C:\Users\alexp\Documents\UW\Research\Selenium\aSe0vBB\particle\selenium-build\output\pixel_sio2_122kev_75degangle_200k.root"
    emfilename = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_5um_spacing_fullsize.txt"
    configfilename = r"./config.txt"

    settings = sc.readConfigFile(configfilename)

    # Create new collection
    newCollection = gEventCollection(filename, [0, 150])
    newCollection.printInfo()

    newCollection.plotAngularDistribution1D()
    # Create simulation object
    # simObject = CarrierSimulation(emfilename=emfilename, eventCollection=newCollection, configfile=configfilename)

    # print('running simulations')
    # # Plot results
    # fig, ax = plt.subplots()
    # indx = [125, 129, 130, 131, 132]
    # results = simObject.processMultipleEvents(indx, processes=1)
    # for i in indx:
    #   print(i)
    #   t, q = simObject.computeChargeSignal(i)
    #   ax.plot(t, -q/1.6e-19*0.05, linewidth=3, label='%i'%i)
    # ax.set_xlabel('Time ($\mu) s$', fontsize=14)
    # ax.set_ylabel('Charge (keV)', fontsize=14)
    # ax.legend()
    # plt.show()

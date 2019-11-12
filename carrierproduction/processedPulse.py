# stores the output of the simulation in a more class containing more

import numpy as np
import scipy as scp
import cPickle as pickle
import os
import git


class SimulatedPulse(object):
    """docstring for ProcessedPulse"""

    statlist = ["meanPosition", "meanWeightedPosition", "sigmaPosition", "totalEnergy"]
    q = 1.6e-19

    def __init__(self, **kwargs):
        super(SimulatedPulse, self).__init__()

        # Geant4 Event data
        self.g4Event = []
        self.binnedPosition = []
        self.binnedEnergy = []

        # Charge yield
        self.wehp = []
        self.nehp = []
        self.nehpf = []
        self.sinTheta = 0

        # Pulse information
        self.timeSeries = []
        self.signalMaximum = []

        # Set the relevant carrier settings
        self.setSimSettings(**kwargs)

        # Define the statistic values we will compute
        for stat in self.statlist:
            setattr(self, stat, None)

    def __str__(self):
        outstr = ""
        # Raw pulse info
        outstr += str(self.g4Event)
        outstr += "\nBins with energy information: "
        outstr += "\n\tx: "
        outstr += str(self.binnedPosition[0])
        outstr += "\n\ty: "
        outstr += str(self.binnedPosition[1])
        outstr += "\n\tz: "
        outstr += str(self.binnedPosition[2])
        outstr += "\n\tenergy: "
        outstr += str(self.binnedEnergy)
        outstr += "\nCharge Yield information: "
        outstr += "\n\tWork function: "
        outstr += str(self.wehp)
        outstr += "\n\tTheoretical Nehp: "
        outstr += str(self.nehp)
        outstr += "\n\tFluctuating Nehp: "
        outstr += str(self.nehpf)
        outstr += "\n\tTrack-Field Mean Sin Theta: "
        outstr += str(self.sinTheta)
        outstr += "\nPulse information: "
        outstr += "\n\tNumber of time steps: "
        outstr += str(self.timeSeries[:, 0].size)
        outstr += "\n\tMax time: "
        outstr += str(self.timeSeries[-1, 0])
        outstr += "\n\tNumber of Signals: "
        outstr += str(self.timeSeries.shape[1] - 1)
        outstr += "\n\tMax Pulse Height (e): "
        outstr += str(self.signalMaximum)
        outstr += "\n\tScale Factor (for Coplanar): "
        outstr += str(self.scaleFactor)

        # Carrier Dynamic setting information
        settingDictionary = self.getSimSettings()
        outstr += "\nCarrier Dynamics settings: "
        for key, val in sorted(settingDictionary.iteritems()):
            outstr += str("\n\t%s: " % key)
            outstr += str(val)

        # Statistic information
        outstr += "\nComputed Statistics on Event: "
        for stat in self.statlist:
            outstr += str("\n\t%s: " % stat)
            outstr += str(getattr(self, stat))

        return outstr

    # Define setters
    def setSimSettings(self, **kwargs):
        """ 
		Take the dictionary of keyword arguments of settings and try to set the relevant sim settings 
		Sets the parameters relevant to the dynamics of the simulation. Other settings are stored in the
		container object	
		"""

        try:
            self.muHoles = kwargs["MU_HOLES"]
            self.muEletrons = kwargs["MU_ELECTRONS"]
            self.tauHoles = kwargs["TAU_HOLES"]
            self.tauElectron = kwargs["TAU_ELECTRONS"]
            self.dtHoles = kwargs["DT_HOLES"]
            self.dtElectrons = kwargs["DT_ELECTRONS"]
            self.timeTotalHoles = kwargs["TOTAL_TIME_HOLES"]
            self.timeTotalElectrons = kwargs["TOTAL_TIME_ELECTRONS"]
            self.workFunction = kwargs["WORK_FUNCTION"]
        except KeyError:
            print ("Invalid keys in kwargs")
            return

    def setG4Event(self, g4Event):
        self.g4Event = g4Event

    def setBinnedPosition(self, x, y, z):
        """ Given a list of x,y,z pairs, convert into a numpy array of position """
        self.binnedPosition = np.array([x, y, z])

    def setBinnedEnergy(self, energy):
        self.binnedEnergy = energy

    def setWehp(self, wehp):
        self.wehp = wehp

    def setNehp(self, nehp, nehpf=None):
        self.nehp = nehp
        self.nehpf = nehpf

    def setTimeSeries(self, t, f, dim=1):
        """ t is an Nx1 array and f is an NxM array, where M is the different number of pulses
		Result is an Nx(M+1) array """
        data = np.zeros((t.size, dim + 1))
        data[:, 0] = t
        if f.ndim == 1:
            data[:, 1] = f
            self.setSignalMaximum(np.max(np.abs(f)) / self.q)
        else:
            data[:, 1:] = f
            self.setSignalMaximum(np.max(np.abs(f), axis=0) / self.q)
        self.timeSeries = data

    def setSignalMaximum(self, signalMax):
        self.signalMaximum = signalMax

    def setScaleFactor(self, scaleFactor):
        self.scaleFactor = scaleFactor

    def setSinThetaSquare(self, sinTheta):
        self.sinTheta = sinTheta

    # Define getters
    def getSimSettings(self):
        """ Returns a dictionary of the simulation settings """

        simSettings = {
            "muHoles": self.muHoles,
            "muElectrons": self.muEletrons,
            "tauHoles": self.tauHoles,
            "tauElectron": self.tauElectron,
            "dtHoles": self.dtHoles,
            "dtElectrons": self.dtElectrons,
            "timeTotalHoles": self.timeTotalElectrons,
            "timeTotalElectrons": self.timeTotalElectrons,
            "workFunction": self.workFunction,
        }

        return simSettings

    def getG4Event(self):
        return self.g4Event

    def getBinnedPosition(self):
        return self.binnedPosition

    def getBinnedEnergy(self):
        return self.binnedEnergy

    def getWehp(self):
        return self.wehp

    def getNehp(self):
        return self.nehp, self.nehpf

    def getTimeSeries(self):
        return self.timeSeries

    def getScaleFactor(self):
        return self.scaleFactor

    def getSinThetaSquare(self):
        return self.sinTheta

    # Computing characteristics of the pulse
    def computeStats(self):
        """ Compute a set of statistics about the event """
        self.meanPosition = self.computeMeanPosition(
            self.g4Event.x, self.g4Event.y, self.g4Event.z
        )
        self.meanWeightedPosition = self.computeMeanWeightedPosition(
            self.g4Event.x, self.g4Event.y, self.g4Event.z
        )
        self.sigmaPosition = self.computeSigmaPosition(
            self.g4Event.x, self.g4Event.y, self.g4Event.z
        )
        self.totalEnergy = self.computeTotalEnergy(self.g4Event.energy)

    @staticmethod
    def computeMeanPosition(x, y, z):
        return np.array([np.average(x), np.average(y), np.average(z)])

    @staticmethod
    def computeMeanWeightedPosition(x, y, z, w):
        return np.array(
            [
                np.average(x, weights=w),
                np.average(y, weights=w),
                np.average(z, weights=w),
            ]
        )

    @staticmethod
    def computeSigmaPosition(x, y, z):
        return np.array([np.std(x), np.std(y), np.std(z)])

    @staticmethod
    def computeTotalEnergy(energy):
        return np.sum(energy)


class SimulatedOutputFile(object):
    """Simulated"""

    def __init__(self, settings=None, outputfile=None):
        super(SimulatedOutputFile, self).__init__()

        if type(settings) is not dict:
            print ("Must pass a dictionary of settings.")
            return

        self.settings = settings

        try:
            pass
        except Exception as e:
            raise

        if outputfile:
            self.outputfile = outputfile
        else:
            try:
                self.outputfile = os.path.join(
                    self.settings["OUTPUT_DIR"], self.settings["OUTPUT_FILE"]
                )
            except KeyError:
                print ("No valid output file dir/name keys in settings.")
                return

        self.pulses = []

    def __str__(self):
        """ Return class information as string """
        return "Carrier Simulation Object with %i simulated Events" % len(self.pulses)

    def printInfo(self):

        infostr = ""

        # Add class info
        infostr += str(self)
        infostr += "\nOutput filename: "
        infostr += str(self.outputfile)

        # Add settings info
        for key, val in self.settings.iteritems():
            infostr += "\n%s: " % key
            infostr += str(val)

        # Add git info
        infostr += self.getGitInfoStr()

        print (infostr)

    def addPulses(self, simulatedPulses):
        try:
            self.pulses.extend(simulatedPulses)
        except TypeError:
            self.pulses.append(simulatedPulses)
        except:
            print ("Error Adding pulse to object")

    def getPulses(self):
        return self.pulses

    def getPulse(self, index):
        try:
            return self.pulses[index]
        except IndexError:
            print ("Invalid index for simulated pulses\n")

    def setGitInfo(self, repopath=None):
        """ Scrapes information about the git status and adds to the class """
        gitInfo = {}

        if type(repopath) is not str:
            repopath = "./"

        # Get the repo
        repo = git.Repo(repopath)

        # Create commit dictionary
        commit = repo.active_branch.commit
        commitDict = {
            "commit": str(commit.hexsha),
            "message": str(commit.summary),
            "author": str(commit.author.name) + " (" + str(commit.author.email) + ")",
            "date": str(commit.authored_datetime),
        }

        gitInfo["repoLocation"] = str(repo.remote().url)
        gitInfo["repoPath"] = str(repo.git_dir)
        gitInfo["activeBranch"] = str(repo.active_branch)
        gitInfo["commit"] = commitDict
        gitInfo["dirty"] = repo.is_dirty()
        gitInfo["dirtyFiles"] = [str(file.a_path) for file in repo.index.diff(None)]

        self.git = gitInfo

    def getGitInfo(self):
        """ returns the git information """
        return self.git

    def getGitInfoStr(self):
        """ returns the git information in human readable string """

        gitInfoStr = ""
        gitInfoStr += "\nGit Remote Location: " + self.git["repoLocation"]
        gitInfoStr += "\nGit Local Location: " + self.git["repoPath"]
        gitInfoStr += "\nCurrent Branch: " + self.git["activeBranch"]
        gitInfoStr += "\nCurrent Commit: "
        gitInfoStr += "\n\tHash: " + self.git["commit"]["commit"]
        gitInfoStr += "\n\tMessage: " + self.git["commit"]["message"]
        gitInfoStr += "\n\tAuthor: " + self.git["commit"]["author"]
        gitInfoStr += "\n\tDate: " + self.git["commit"]["date"]
        gitInfoStr += "\nIs Dirty: " + str(self.git["dirty"])

        # Print the dirty files
        gitInfoStr += "\nDirty Files: "
        for file in self.git["dirtyFiles"]:
            gitInfoStr += "\n\t" + file

        return gitInfoStr

    def computeResolution(self, bins=140):
        """
        Computes the resolution (sigma/mean) of the main peak
        """
        amplitude = []
        for pulse in self.pulses:
            amplitude.append(pulse.signalMaximum)

        # Histogram the amplitude of each pulse
        f, bins = np.histogram(amplitude, bins=bins)
        x = bins[1:] - np.diff(bins)[0]

        # Define the fit range (only looking for the main peak)
        fr = [2.0 / 3 * bins[-1], bins[-1]]

        # Perform the fit
        param, covariance = fitGaussian(x, f, fr)

        return param[2] / param[1]


def savePickleObject(obj, filename):
    """ Pickles and object and saves it to disk """

    with open(filename, "wb") as outfile:
        pickle.dump(obj, outfile, pickle.HIGHEST_PROTOCOL)


def loadPickleObject(filename):
    """ Load object and unpickles it """

    with open(filename, "rb") as infile:
        obj = pickle.load(infile)

    return obj


def fitGaussian(x, y, fitRange=[]):
    """ 
    Perform a gaussian fit to data (x, y) in the fitRange. If fitRange is empty, do fit over entire range
    Outputs the paramters that minimize the least squares and the covariance matrix
    """

    xslice, yslice = [], []

    if fitRange:
        for i in range(len(x)):
            if x[i] >= fitRange[0] and x[i] <= fitRange[1]:
                xslice.append(x[i])
                yslice.append(y[i])
        xslice = np.array(xslice)
        yslice = np.array(yslice)
    else:
        xslice = x
        yslice = y

    mui = np.sum(xslice * yslice) / np.sum(yslice)
    sigmai = np.sqrt(np.sum(yslice * (xslice - mui) ** 2) / np.sum(yslice))
    pi = [1 / (np.sqrt(2 * np.pi) * sigmai), mui, sigmai]

    param, cov = scp.optimize.curve_fit(gaus, xslice, yslice, p0=pi)

    return param, cov


def gaus(x, *par):
    """
    Function definition of gaussian.
    Inputs:
        x - numpy array of x values to evaluate the value of the gaussian at
        par[0] - N, amplitude
        par[1] - mu, mean
        par[2] - sigma, standard deviation
    Ouptus:
        Gaussian with the given parameters for each x values
    """
    N, mu, sigma = par
    return N * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


if __name__ == "__main__":
    # Testing the operation of these objects

    import seleniumconfig as sc
    import particledata as pd

    # Create the data for the simulated pulse objects

    # Load test particle data
    eventCol = pd.gEventCollection(
        "/home/apiers/mnt/rocks/selena/data/particle/pixel_sio2_122kev_12degBeam_100k.root",
        eventCounterRange=[0, 5],
    )

    # Load config file
    config = sc.readConfigFile(
        "/home/apiers/mnt/rocks/aSe0vBB/carrierproduction/config.txt"
    )

    # Create fake data
    n = 100
    t = np.linspace(0, 10, n)
    s1 = t[0 : int(np.random.uniform(1, n))] * np.random.uniform(0, 1)
    s2 = t[0 : int(np.random.uniform(1, n))] * np.random.uniform(0, 1)
    signal = np.concatenate((s1, s1[-1] * np.ones(n - s1.size))) + np.concatenate(
        (s2, s2[-1] * np.ones(n - s2.size))
    )
    x, y, z = [1, 1, 1], [1, 2, 3], [5, 7, 7]
    energy = np.array([12, 10, 100])
    wehp = np.array([0.04, 0.041, 0.042])
    nehp = energy / wehp
    nehpf = np.random.poisson(nehp)

    # Create and test functionality of simulated pulse object
    simpulse = SimulatedPulse(**config)
    simpulse.setG4Event(eventCol.collection[1])
    simpulse.setBinnedPosition(x, y, z)
    simpulse.setBinnedEnergy(energy)
    simpulse.setWehp(wehp)
    simpulse.setNehp(nehp, nehpf)
    simpulse.setTimeSeries(t, signal)
    simpulse.computeStats()
    print (simpulse)
    oldclass = str(simpulse)

    # Create simulated output file object
    simfile = SimulatedOutputFile(settings=config, outputfile="./test.npy")
    simfile.addPulses(simpulse)
    simfile.setGitInfo("/home/apiers/mnt/rocks/aSe0vBB")
    print ("\n")
    simfile.printInfo()
    print ("\n")
    print (simfile.getGitInfoStr())

    # Save file
    savePickleObject(simfile, simfile.outputfile)

    # Load file and print results to compare
    loadSimfile = loadPickleObject(simfile.outputfile)
    newclass = str(loadSimfile.pulses[0])
    assert oldclass == newclass

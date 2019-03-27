# stores the output of the simulation in a more class containing more

import numpy as np
import scipy as scp
import cPickle as pickle
import os

class SimulatedPulse(object):
	"""docstring for ProcessedPulse"""
	def __init__(self, **kwargs):
		super(ProcessedPulse, self).__init__()
		
		# Geant4 Event data
		self.g4Event = []
		self.binnedPosition = []
		self.binnedEnergy = []

		# Charge yield
		self.wehp = []
		self.nehp = []
		self.nehpf = []

		# Pulse information
		self.timeSeries = []
		self.path = []

		# Set the relevant carrier settings
		self.setSimSettings(**kwargs)

	# Define setters
	def setSimSettings(self, **kwargs):
		""" 
		Take the dictionary of keyword arguments of settings and try to set the relevant sim settings 
		Sets the parameters relevant to the dynamics of the simulation. Other settings are stored in the
		container object	
		"""

		try:
			self.muHoles 			= kwargs["MU_HOLES"]
			self.muEletrons 		= kwargs["MU_ELECTRONS"]
			self.tauHoles 			= kwargs["TAU_HOLES"]
			self.tauElectron 		= kwargs["TAU_ELECTRONS"]
			self.dtHoles 			= kwargs["DT_HOLES"]
			self.dtElectrons 		= kwargs["DT_ELECTRONS"]
			self.timeTotalHoles 	= kwargs["TOTAL_TIME_HOLES"]
			self.timeTotalElectrons = kwargs["TOTAL_TIME_ELECTRONS"]
			self.workFunction 		= kwargs["WORK_FUNCTION"]
		except KeyError:
			print("Invalid keys in kwargs")
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
		data = np.zeros( (t.size, dim+1) )
		data[:, 0]  = t
		data[:, 1:] = f
		self.timeSeries = data

	# Define getters
	def getSimSettings(self):
		""" Returns a dictionary of the simulation settings """

		simSettings = {
			"muHoles"				: self.muHoles,
			"muElectrons"			: self.muEletrons,
			"tauHoles"				: self.tauHoles,
			"tauElectron"			: self.tauElectron,
			"dtHoles"				: self.dtHoles,
			"dtElectrons"			: self.dtElectrons,
			"timeTotalHoles"		: self.timeTotalElectrons,
			"timeTotalElectrons"	: self.timeTotalElectrons,
			"workFunction"			: self.workFunction,
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

	# Computing characteristics of the pulse
	def computeStats(self):
		""" Compute a set of statistics about the event """
		g4Data = self.g4Event.flattenEvent();
		self.meanPosition = computeMeanPosition(g4Data["x"], g4Data["y"], g4Data["z"])
		self.meanWeightedPosition = computeMeanWeightedPosition(g4Data["x"], g4Data["y"], g4Data["z"], g4Data["energy"])
		self.sigmaPosition = computerSigmaPosition(g4Data["x"], g4Data["y"], g4Data["z"])
		self.totalEnergy = computeTotalEnergy(g4Data["energy"])

	def computeMeanPosition(x, y, z):
		return np.array([np.mean(x), np.mean(y), np.mean(z)]) 

	def computeMeanWeightedPosition(x, y, z, w):
		return np.array([np.mean(x*w), np.mean(y*w), np.mean(z*w)]) / (np.sum(w))

	def computeSigmaPosition(x, y, z):
		return np.array([np.std(x), np.std(y), np.std(z)])

	def computeTotalEnergy(energy):
		return np.sum(energy)

class SimulatedOutputFile(object):
	"""Simulated"""
	def __init__(self, settings=None, outputfile=None):
		super(SimulatedOutputFile, self).__init__()
		
		if type(settings) is not dict:
			print("Must pass a dictionary of settings.")
			return

		self.settings=settings

		try:
			pass
		except Exception as e:
			raise

		if outputfile:
			self.outputfile = outputfile
		else:
			try:
				self.outputfile = os.path.join(self.settings["OUTPUT_DIR"], self.settings["OUTPUT_FILE"])
			except KeyError:
				print("No valid output file dir/name keys in settings.")
				return

		self.pulses = []



def savePickleObject(obj, filename):
	""" Pickles and object and saves it to disk """

	with open(filename, "wb") as outfile:
		pickle.dump(obj, outfile, pickle.HIGHEST_PROTOCOL)


def loadPickleObject(filename):
	""" Load object and unpickles it """

	with open(filename, "rb") as infile:
		obj = pickle.load(infile)

	return obj
# Module containing classes/functions for handling the Geant4 particle data
# First step for electron hole pair production

# Libraries to import
import ROOT as rt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import sys

# add the EM plot module to the path and import it
sys.path.append(r"..\EM Analysis")
from plot import *

import timeit


class gEvent(object):
	"""docsring for gEvent"""
	def __init__(self, gEventID=-1, hits=[], geometry=None):
		""" Initializes the event class
		Inputs:
			gEventID - int geant4 event ID number
			hits - list of hits. Each hit is a dictionairy containing the information provided by geant4 for hits (position, energy, track and parent ID, particle type, creator process name)
		"""
		self.gEventID = gEventID
		self.hits = hits
		if geometry:
			self.geomtry = geometry
		else:
			self.geometry = {'x':(-5, 5), 'y':(-5, 5), 'z':(-0.1, 0.1)}


	# Define Setters
	def SetEventID(eventID):
		""" Sets the geant4 event ID """
		self.gEventID = eventID

	def SetHits(hits):
		""" Sets the list of hits """
		self.hits = hits

	# Define Getters
	def GetEventID(self):
		""" Returns the event ID """
		return self.gEventID

	def GetHits(self):
		""" Returns the list of hits """
		return self.hits

	def GetSpecificHit(self, index):
		""" Tries to return a single dictionairy hit at index """
		try:
			return self.hits[index]
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
			boolTypeCheck.append(type(hit["energy"]) is float)
			boolTypeCheck.append(type(hit["particle"]) is str)
			boolTypeCheck.append(type(hit["creatorProcess"]) is str)

			if all(boolTypeCheck):
				self.hits.append(hit)
			else:
				print("One or more data type in the hit is invalid")

		except Exception as e:
			raise e


	def flattenEvent(self):
		""" Transforms the hit data from a list of dictionary into a dictionary of lists--therefore calling the dict member gets you list of all values in that event. Useful for plotting and data manipulation. Documentation assumes that there are N different hits.
		Outputs:
			eventData = {'trackID' : [], 'parentID' : [], 'x' : [], 'y' : [], 'z' : [], 'energy' : [], 'particle' : [], 'creatorProcess' : []}

		"""

		flattenedData = {'trackID' : [], 'parentID' : [], 'x' : [], 'y' : [], 'z' : [], 'energy' : [], 'particle' : [], 'creatorProcess' : []}
		for hit in self.GetHits():
			flattenedData['trackID'].append(hit['trackID'])
			flattenedData['parentID'].append(hit['parentID'])
			flattenedData['x'].append(hit['x'])
			flattenedData['y'].append(hit['y'])
			flattenedData['z'].append(hit['z'])
			flattenedData['energy'].append(hit['energy'])
			flattenedData['particle'].append(hit['particle'])
			flattenedData['creatorProcess'].append(hit['creatorProcess'])
		
		return flattenedData 


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

		flatData = self.flattenEvent()
		histX = flatData[x]
		histY = flatData[y]

		val, bins, _ = ax.hist(histX, bins=nbins, range=self.geometry[x], weights=histY, log=True, histtype='step')

		# Set some default axis
		ax.set_xlabel(x + " (mm)", fontsize=14)
		ax.set_ylabel(y + " (keV)", fontsize=14)
		ax.set_title(y + " vs " + x + " for Event ID %i " % self.GetEventID(), fontsize=16)

		return val, bins, ax, fig 


	def plotH2(self, x="y", y="z", z="energy", nbins=[200, 200], figH=None, delete=False):

		if not figH:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		else:
			fig = figH[0]
			ax = figH[1]

		flatData = self.flattenEvent()
		histX = flatData[x]
		histY = flatData[y]
		histZ = flatData[z]

		# setting up the range
		histRange = [self.geometry[x], self.geometry[y]]
		val, binx, biny, cax = ax.hist2d(histX, histY, range=histRange, bins=nbins, weights=histZ)

		# Set some default axis
		ax.set_xlabel(x + " (mm)", fontsize=14)
		ax.set_ylabel(y + " (mm)", fontsize=14)
		ax.set_title(z + " vs " + x + " and " + y + " for Event ID %i" % self.GetEventID(), fontsize=16)

		fig.colorbar(cax)

		if delete:
			plt.close(fig)

		return val, [binx, biny], ax, fig 


	def createCarriers(self, modelOptions=None):
		""" 
		Function that creates carriers from the energy distribution of the incoming gamma rays
		"""
		Wehp = 0.05 # keV. 50 eV per electron hole pair

		# create histogram to bin the data in 1 um increments. Pass this into EM model
		# number of z and y bins depend on the geometry
		nbinx, nbiny = int(np.diff(self.geometry["y"])*1000), int(np.diff(self.geometry["z"])*1000)
		val, bins, ax, fig = self.plotH2(nbins=[nbinx, nbiny], delete=True)

		# get center points of all bins
		binxCenter = (np.array(bins[0][0:-1]) + np.array(bins[0][1:]))/2.
		binyCenter = (np.array(bins[1][0:-1]) + np.array(bins[1][1:]))/2.

		# histogrammed data are in a (nbinx x nbiny) array. Simply divide by Wehp to get the energy as a distribution of position
		nehp = np.round(val/Wehp)
		nehpFluctuation = np.random.poisson(nehp)

		# Add noise to number of electron hole pair
		return nehp, nehpFluctuation, binxCenter, binyCenter

class gEventCollection(object):
	"""docstring for gEventCollection"""
	def __init__(self, rootFilename):
		
		self.rootFilename = rootFilename
		self.collection = []

		# Read the data from the file
		f = rt.TFile(rootFilename)
		# Gets a tree with the name aSeData. Current name from Geant4 simulation. 
		tree = f.Get("aSeData")		
		eventID = -1
		hitsList = []

		# Iterate over all of the entries in the tree, extracting tuple i
		for entry in tree:
			if eventID != entry.EventID:
				if eventID != -1:
					self.collection.append(gEvent(gEventID=eventID, hits=hitsList))
				hitsList = []
				eventID = entry.EventID
				hit = {'trackID' : entry.TrackID, 'parentID' : entry.ParentID, 'x' : entry.x, 'y' : entry.y, 'z' : entry.z, 'energy' : entry.energy*1e3, 'particle' : entry.ParticleType, 'creatorProcess' : entry.ProcessName}
				hitsList.append(hit)

			else:
				hit = {'trackID' : entry.TrackID, 'parentID' : entry.ParentID, 'x' : entry.x, 'y' : entry.y, 'z' : entry.z, 'energy' : entry.energy*1e3, 'particle' : entry.ParticleType, 'creatorProcess' : entry.ProcessName}
				hitsList.append(hit)


	def printInfo(self):
		""" Prints information regarding the collection of events"""
		print("Particle data comes from " + self.rootFilename)
		print("There are %i events stored in this collection" % (len(self.collection)))


	def findEvent(self, eventID):
		"""Searches the collection for an event with eventID. Returns if finds it"""
		foundEvent = None

		for event in self.collection:
			if event.GetEventID() == eventID:
				foundEvent = event

		if not foundEvent:
			print("Could not find a event matching the ID %i" % eventID)

		return foundEvent

def computeChargeSignal(event, emFilename):
	"""
	Uses the number of electron hole pairs spatial distribution from the event and the electrostatic simulation from COMSOL to determine the induced charge signal on the detector plates
	"""

	# Gets the carrier information
	nehp, nehpf, binx, biny = event.createCarriers()

	# Gets the EM information
	_, x, y, data = readComsolFileGrid(emFilename)
	Ex = data[:, 1] / 1e3 # Convert to V/mm
	Ey = data[:, 2] / 1e3 # Convert to V/mm
	Etot = [x*1e3, y*1e3, Ex, Ey] # Construct the E field. Put lengths in mm to match Comsol
	wPhi = [x*1e3, y*1e3, data[:, 3]]
	dt = 0.01 # in us. so 10 ns
	muHole, muElectron = 19e-6, 1e-6 # mm^2/(V*us)
	maxtime = 0
	allInduced = [] 

	# iterate over all the grid points of electron hole pairs
	for i in range(binx.shape[0]):
		for j in range(biny.shape[0]):
			if nehpf[i,j] == 0:
				pass
			else:
				# Find trajectory of this group of particles
				qHoles = 1.6e-19*nehpf[i,j]
				qElectrons = -qHoles

				pathHoles = findMotion((binx[i], biny[j]), Etot, muHole, dt, q=qHoles)
				pathElectrons = findMotion((binx[i], biny[j]), Etot, muElectron, dt, q=qElectrons)

				# Keep the longest time for reference in 
				if np.max([pathHoles[-1,2], pathElectrons[-1,2]]) > maxtime:
					maxtime = np.max([pathHoles[-1,2], pathElectrons[-1,2]])

				indChargeHoles = inducedChargeSingle(wPhi, pathHoles, q=qHoles)
				indChargeElectrons = inducedChargeSingle(wPhi, pathElectrons, q=qElectrons)
				allInduced.append(indChargeHoles)
				allInduced.append(indChargeElectrons)

	time = np.arange(0, maxtime + dt, dt)
	qInduced = np.zeros(time.shape)

	# Iterate over the induced charge by each grid contribution and add them to the total charge induced on the electrode
	for charge in allInduced:
		qInduced[:len(charge)] += charge
		# fill the rest of the time values with the last value so that the charge doesn't "disapear"
		qInduced[len(charge):] += charge[-1] 

	return time, qInduced




	return poissonDistMatrix
if __name__ == '__main__':
	# used for debuggin information. If the particledata.py file is run this segment of the code will run
	filename = r"C:\Users\alexp\Documents\UW\Research\Selenium\aSe0vBB\particle\selenium-build\output\122_keV_testTuple.root" 

	emfilename = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\real_electrode.txt"

	newCollection = gEventCollection(filename)
	newCollection.printInfo()

	event = newCollection.collection[124]
	cProc = event.GetHits()[1]['creatorProcess'].split('\x00')[0]
	# event.createCarriers()

	t, q = computeChargeSignal(event, emfilename)
	fig, ax= plt.subplots()

	ax.plot(t,-1*q, linewidth=3)
	ax.set_xlabel(r'Time ($\mu$s)', fontsize=14)
	ax.set_ylabel(r'Induced Charge (C)', fontsize=14)
	ax.set_title('Q(t) at a Unipolar electrode for event %i. Initial electron creator process: %s' % (event.GetEventID(), cProc) , fontsize=14)
	event.plotH1(x="z", y="energy", nbins=200)
	event.plotH2()


	# iterate many times to watch fluctuation
	fig1, ax1 = plt.subplots()

	iterations = 50
	for i in range(iterations):
		print(i)

		t, q = computeChargeSignal(event, emfilename)
		ax1.plot(t,-1*q,color='grey', alpha=0.5)
		if i == 0:
			qtot = q
		else:
			qtot += q

	ax1.plot(t,-1*qtot/iterations, color='black', linewidth=3)
	ax1.set_xlabel(r'Time ($\mu$s)', fontsize=14)
	ax1.set_ylabel(r'Induced Charge (C)', fontsize=14)
	ax1.set_title('Q(t) at a Unipolar electrode with Poisson Fluctuation of Carriers for event %i' % (event.GetEventID()) , fontsize=16)

	plt.show()



	plt.show()

	print('stopping')
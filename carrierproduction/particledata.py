# Module containing classes/functions for handling the Geant4 particle data
# First step for electron hole pair production

# Libraries to import
import ROOT as rt
import matplotlib.pyplot as plt
import numpy as np
import sys
import seleniumconfig as sc

# add the EM plot module to the path and import it
sys.path.append(r"..\EM Analysis")
from plot import readComsolFileGrid, findMotion, inducedCharge, inducedChargeSingle


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
	def SetEventID(self, eventID):
		""" Sets the geant4 event ID """
		self.gEventID = eventID

	def SetHits(self, hits):
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


		flatData = self.flattenEvent()
		histX = flatData[x]
		histY = flatData[y]
		histZ = flatData[z]

		# setting up the range
		histRange = [self.geometry[x], self.geometry[y]]

		if delete:
			val, binx, biny = np.histogram2d(histX, histY, range=histRange, bins=nbins, weights=histZ)
			fig, ax = None, None

		else:
			if not figH:
				fig = plt.figure()
				ax = fig.add_subplot(111)
			else:
				fig = figH[0]
				ax = figH[1]

			val, binx, biny, cax = ax.hist2d(histX, histY, range=histRange, bins=nbins, weights=histZ)

			# Set some default axis
			ax.set_xlabel(x + " (mm)", fontsize=14)
			ax.set_ylabel(y + " (mm)", fontsize=14)
			ax.set_title(z + " vs " + x + " and " + y + " for Event ID %i" % self.GetEventID(), fontsize=16)

			fig.colorbar(cax)


		return val, [binx, biny], ax, fig


	def createCarriers(self, **kwargs):
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
		if kwargs['CARRIER_GENERATION_POISSON']:
			nehpFluctuation = np.random.poisson(nehp) - nehp
		else:
			nehpFluctuation = np.zeros(nehp.shape)

		# Add noise to number of electron hole pair
		return nehp, nehpFluctuation, binxCenter, binyCenter

class gEventCollection(object):
	"""docstring for gEventCollection"""
	def __init__(self, rootFilename, **kwargs):

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

class CarrierSimulation(object):
	"""Wrapper class for all the functions and data surrounding the charge carrier simulation. All info is housed in one class for efficiency purposes"""
	def __init__(self, emfilename=None, eventCollection=None, configfile=None):
		"""
		Initialize the class
		"""
		self.emfilname = emfilename
		self.configfile = configfile

		print('Read Geant4 Simulation Data')
		if type(eventCollection) is str:
			self.eventCollection = gEventCollection(eventCollection)
		else:
			self.eventCollection = eventCollection

		# Read the emfile and store data in the class
		print('Read Comsol EM Data')
		if emfilename:
			_, self.x, self.y, self.data = readComsolFileGrid(emfilename)
		else:
			self.x = self.y = self.data = []

		# Read and store settings
		print('Read config')
		if configfile:
			self.settings = sc.readConfigFile(configfile)
		else:
			self.settings = None

		# Perform initialization actions based on the settings
		print('Apply Relevant')
		if settings['SCALE_WEIGHTED_PHI']:
			self.computeScaleFactor()
		else:
			self.scale = 1

	def newEmFile(self, emfilename):
		"""
		Replace the em data with what is in the new file
		Inputs:
			emfilename - string of the path to the comsol simulation file
		"""
		_, self.x, self.y, self.data = readComsolFileGrid(emfilename)

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

	def computeScaleFactor(self, depth=[0.5, 1], width=[0.35, 0.65]):
		"""
		Computes the scale factor between weighted potentials for tiered coplanar design. Save it as an attribute of the class so it only needs to be done once
		"""

		# Get dimensions and sets the limits that we want to take
		xdim = self.x.size
		ydim = self.y.size
		xmin, xmax = int(xdim*width[0]), int(xdim*width[1])
		ymin, ymax = int(ydim*depth[0]), int(ydim*depth[1])

		phi = self.data[:,0:np.amin(self.data.shape):3]
		wPhiA = np.reshape(phi[:,-2], (ydim, xdim))
		wPhiB = np.reshape(phi[:,-1], (ydim, xdim))

		phiDiffA = np.diff(wPhiA[ymin:ymax-1, xmin:xmax-1], axis=0)
		phiDiffB = np.diff(wPhiB[ymin:ymax-1, xmin:xmax-1], axis=0)

		self.scale = np.mean(phiDiffA/phiDiffB)

		return self.scale

	def createFields(self):
		"""
		Takes the raw comsol data and parses it into the form used by findmotion and inducedCharge
		"""

		# Gets the EM information
		x, y, data = self.x, self.y, self.data
		convFactor = self.settings['UNIT_CONVERSION_FACTOR']
		x *= convFactor
		y *= convFactor # convert m to mm

		# Parses the data depending on the mode
		if self.settings['CHARGE_DIFFERENCE']:
			phiAll, ExAll, EyAll = data[:,0:np.amin(data.shape):3], data[:,1:np.amin(data.shape):3], data[:,2:np.amin(data.shape):3]
			wPhiA = [x, y, phiAll[:,-2]]
			wPhiB = [x, y, phiAll[:,-1]]
			wPhi = [wPhiA, wPhiB]
			Ex, Ey = ExAll[:,self.settings['VOLTAGE_SWEEP_INDEX']]/convFactor, EyAll[:,self.settings['VOLTAGE_SWEEP_INDEX']]/convFactor # in V/mm
		else:
			phiAll, ExAll, EyAll = data[:,0:np.amin(data.shape):3], data[:,1:np.amin(data.shape):3], data[:,2:np.amin(data.shape):3]
			wPhi = [x, y, phiAll[:,-1]]
			Ex, Ey = ExAll[:,self.settings['VOLTAGE_SWEEP_INDEX']]/convFactor, EyAll[:,self.settings['VOLTAGE_SWEEP_INDEX']]/convFactor # in V/mm

		Etot = [x, y, Ex, Ey] # Construct the E field. Put lengths in mm to match Comsol

		return Etot, wPhi

	def chargeArray(self, nehp, nHoleSteps, nElectronSteps):
		"""
		Generates the charge arrays used in the induced charge calculations
		"""

		# Get data from the settings
		tauHole = self.settings['TAU_HOLES']
		tauElectron = self.settings['TAU_ELECTRONS']
		dtHole = self.settings['DT_HOLES']
		dtElectron = self.settings['DT_ELECTRONS']
		elementaryCharge = self.settings['ELEMENTARY_CHARGE']

		qHole = elementaryCharge*nehp
		qElectron = -qHole

		# Create array for number of charges particles at each point of the path
		qHoleArray = qHole*np.ones(nHoleSteps)
		qElectronArray = qElectron*np.ones(nElectronSteps)

		# If carrier lifetime noise is considered, do a geometric distribution to find the number of steps until a success. Iterate over that an substract
		if self.settings['CARRIER_LIFETIME_GEOMETRIC']:
			holeTrapTimeIndx = np.random.geometric(dtHole/tauHole, nehp)
			electronTrapTimeIndx = np.random.geometric(dtElectron/tauElectron, nehp)

			# Eliminate values that are greater than the size of the length of the path
			holeTrapTimeIndx = holeTrapTimeIndx[holeTrapTimeIndx < qHoleArray.size]
			electronTrapTimeIndx = electronTrapTimeIndx[electronTrapTimeIndx < qElectronArray.size]

			# iterate over all the trap location. Subtract a charge from that location in the charge array
			for k in range(holeTrapTimeIndx.size):
				qHoleArray[holeTrapTimeIndx[k]:] -= elementaryCharge
			for k in range(electronTrapTimeIndx.size):
				qElectronArray[electronTrapTimeIndx[k]:] += elementaryCharge

		return qHoleArray, qElectronArray

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

		# Gets the carrier information
		event = self.eventCollection.collection[eventIdx]
		nehp, nehpNoise, binx, biny = event.createCarriers(**self.settings)
		nehpf = nehp + nehpNoise	# number of electron hole pairs w/ fluctuations is the average plus the noise noise determined by the specific model
		nehpf = nehpf.astype(int)

		Etot, wPhi = self.createFields()

		# Set up the boundaries of the em simulation
		if self.settings['USE_BOUNDARY_LIMITS']:
			limits = [self.settings['X_MIN'], self.settings['X_MAX'], self.settings['Y_MIN'], self.settings['Y_MAX']]
		else:
			limits=[]

		# Get parameters from the settings necessary for the simulations
		dtHole = self.settings['DT_HOLES']
		dtElectron = self.settings['DT_ELECTRONS']
		muHole, muElectron = self.settings['MU_HOLES'], self.settings['MU_ELECTRONS']# mm^2/(V*us)
		maxtime = 0
		allInduced = []

		# iterate over all the grid points of electron hole pairs
		for i in range(binx.shape[0]):
			for j in range(biny.shape[0]):
				if nehpf[i,j] == 0:
					pass
				else:
					# Find trajectory of this group of particles
					# Assign the sign of the charge of holes and electrons
					qHoles = 1
					qElectrons = -qHoles

					# Find path of the electrons and holes at grid point i,j
					pathHoles = findMotion((binx[i], biny[j]), Etot, muHole, dtHole, q=qHoles, limits=limits)
					pathElectrons = findMotion((binx[i], biny[j]), Etot, muElectron, dtElectron, q=qElectrons, limits=limits)

					qHoleArray, qElectronArray = self.chargeArray(nehpf[i,j], pathHoles.shape[0], pathElectrons.shape[0])

					# Keep the longest time for reference in
					if np.max([pathHoles[-1,2], pathElectrons[-1,2]]) > maxtime:
						maxtime = np.max([pathHoles[-1,2], pathElectrons[-1,2]])

					if self.settings['CHARGE_DIFFERENCE']:
						# Compute difference in induced charge, qA-qB, for both electrons and holes. Incorporates scale factor (which is =1 for no scaling)
						qAHoles, qBHoles, _ = inducedCharge(wPhi[0], wPhi[1], pathHoles, q=qHoleArray)
						qAElectrons, qBElectrons, _ = inducedCharge(wPhi[0], wPhi[1], pathElectrons, q=qElectronArray)
						indChargeHoles = qAHoles - self.scale * qBHoles
						indChargeElectrons = qAElectrons - self.scale * qBElectrons

						# append to  list of induced charge
						allInduced.append(indChargeHoles)
						allInduced.append(indChargeElectrons)

					else:
						# Compute induced charge on a single electrode
						indChargeHoles = inducedChargeSingle(wPhi, pathHoles, q=qHoleArray)
						indChargeElectrons = inducedChargeSingle(wPhi, pathElectrons, q=qElectronArray)
						allInduced.append(indChargeHoles)
						allInduced.append(indChargeElectrons)

		timeHoles = np.arange(0, maxtime + dtHole, dtHole)
		timeElectrons = np.arange(0, maxtime + dtElectron, dtElectron)
		qIndHole = np.zeros(timeHoles.shape)
		qIndElectron = np.zeros(timeElectrons.shape)

		# Iterate over the induced charge by each grid contribution and add them to the total charge induced on the electrode
		for indx, charge in enumerate(allInduced):
			if indx % 2 == 0:
				qIndHole[:len(charge)] += charge
				# fill the rest of the time values with the last value so that the charge doesn't "disapear"
				qIndHole[len(charge):] += charge[-1]
			else:
				qIndElectron[:len(charge)] += charge
				qIndElectron[len(charge):] += charge[-1]

		# Interpolate the Electron induced charge to match the holes so we can add them together

		qInduced = qIndHole + np.interp(timeHoles, timeElectrons, qIndElectron)

		return timeHoles, qInduced

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
	zindx = int(depth*zdim)
	xmin, xmax = int(xrange[0]*wPhiA.shape[1]), int(xrange[1]*wPhiA.shape[1])

	# Calculates the slope of the weighted potential in the linear region. Then finds the ratio at every point and averages them
	phiDiffA = np.diff(wPhiA[zindx:, xmin:xmax], axis=0)
	phiDiffB = np.diff(wPhiB[zindx:, xmin:xmax], axis=0)

	return np.mean(phiDiffA/phiDiffB)

if __name__ == '__main__':
	# used for debuggin information. If the particledata.py file is run this segment of the code will run
	filename = r"C:\Users\alexp\Documents\UW\Research\Selenium\aSe0vBB\particle\selenium-build\output\122_keV_testTuple.root"
	emfilename = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_5um_spacing_fullsize.txt"
	configfilename = r"./config.txt"

	settings = sc.readConfigFile(configfilename)

	# Create new collection
	newCollection = gEventCollection(filename)
	newCollection.printInfo()

	# Create simulation object
	simObject = CarrierSimulation(emfilename=emfilename, eventCollection=newCollection, configfile=configfilename)

	t, q = simObject.computeChargeSignal(125)

	# Plot results
	fig, ax = plt.subplots()
	ax.plot(t, -q/1.6e-19*0.05, linewidth=3)
	ax.set_xlabel('Time ($\mu s$', fontsize=14)
	ax.set_ylabel('Charge (keV)', fontsize=14)
	plt.show()




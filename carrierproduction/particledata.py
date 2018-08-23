# Module containing classes/functions for handling the Geant4 particle data
# First step for electron hole pair production

# Libraries to import
import ROOT as rt
try:
	import matplotlib.pyplot as plt
except:
	pass
import numpy as np
import sys
import seleniumconfig as sc
from pathlib import Path

# add the EM plot module to the path and import it
sys.path.append(Path('../EM Analysis'))
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

def computeChargeSignal(event, emFilename, **kwargs):
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
	nehp, nehpNoise, binx, biny = event.createCarriers(**kwargs)
	nehpf = nehp + nehpNoise 	# number of electron hole pairs w/ fluctuations is the average plus the noise
								# noise determined by the specific model

	# Gets the EM information
	_, x, y, data = readComsolFileGrid(emFilename)
	convFactor = 1.e3 # for converting m to mm
	x *= convFactor
	y *= convFactor # convert m to mm


	if kwargs['CHARGE_DIFFERENCE']:
		phiAll, ExAll, EyAll = data[:,0:np.amin(data.shape):3], data[:,1:np.amin(data.shape):3], data[:,2:np.amin(data.shape):3]
		wPhiA = [x, y, phiAll[:,-2]]
		wPhiB = [x, y, phiAll[:,-1]]
		Ex, Ey = ExAll[:,kwargs['VOLTAGE_SWEEP_INDEX']]/convFactor, EyAll[:,kwargs['VOLTAGE_SWEEP_INDEX']]/convFactor # in V/mm
	else:
		phiAll, ExAll, EyAll = data[:,0:np.amin(data.shape):3], data[:,1:np.amin(data.shape):3], data[:,2:np.amin(data.shape):3]
		wPhi = [x, y, phiAll[:,-1]]
		Ex, Ey = ExAll[:,kwargs['VOLTAGE_SWEEP_INDEX']]/convFactor, EyAll[:,kwargs['VOLTAGE_SWEEP_INDEX']]/convFactor # in V/mm

	# Set up the boundaries of the em simulation
	if kwargs['USE_BOUNDARY_LIMITS']:
		limits = [kwargs['X_MIN'], kwargs['X_MAX'], kwargs['Y_MIN'], kwargs['Y_MAX']]
	else:
		limits=[]

	Etot = [x, y, Ex, Ey] # Construct the E field. Put lengths in mm to match Comsol
	dt = 0.02 # in us. so 10 ns
	muHole, muElectron = kwargs['MU_HOLES'], kwargs['MU_ELECTRONS']# mm^2/(V*us)
	tauHole, tauElectron = kwargs['TAU_HOLES'], kwargs['TAU_ELECTRONS']
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


				# Find path of the electrons and holes at grid point i,j
				pathHoles = findMotion((binx[i], biny[j]), Etot, muHole, dt, q=qHoles, limits=limits)
				pathElectrons = findMotion((binx[i], biny[j]), Etot, muElectron, dt, q=qElectrons, limits=limits)

				# Create array for number of charges particles at each point of the path
				qHoleArray = qHoles*np.ones(pathHoles.shape[0])
				qElectronArray = qElectrons*np.ones(pathElectrons.shape[0])

				# If carrier lifetime noise is considered, do a geometric distribution to find the number of steps until a success. Iterate over that an substract
				if kwargs['CARRIER_LIFETIME_GEOMETRIC']:
					holeTrapTimeIndx = np.random.geometric(dt/tauHole, int(nehpf[i,j]))
					electronTrapTimeIndx = np.random.geometric(dt/tauElectron, int(nehpf[i,j]))

					# Eliminate values that are greater than the size of the length of the path
					holeTrapTimeIndx = holeTrapTimeIndx[holeTrapTimeIndx < qHoleArray.size]
					electronTrapTimeIndx = electronTrapTimeIndx[electronTrapTimeIndx < qElectronArray.size]

					# iterate over all the trap location. Subtract a charge from that location in the charge array
					for k in range(holeTrapTimeIndx.size):
						qHoleArray[holeTrapTimeIndx[k]:] -= 1.6e-19
					for k in range(electronTrapTimeIndx.size):
						qElectronArray[electronTrapTimeIndx[k]:] += 1.6e-19

					# iterate

				# Keep the longest time for reference in
				if np.max([pathHoles[-1,2], pathElectrons[-1,2]]) > maxtime:
					maxtime = np.max([pathHoles[-1,2], pathElectrons[-1,2]])

				if kwargs['CHARGE_DIFFERENCE']:

					# Implement scaled weighted potential if necessary. Finds teh ration between phiA and phiB and then scales qB by the same amount
					if kwargs['SCALE_WEIGHTED_PHI']:
						scale = scaleWeightPhi(wPhiA[2], wPhiB[2], [x.size, y.size])
						qAHoles, qBHoles, _ = inducedCharge(wPhiA, wPhiB, pathHoles, q=qHoleArray, roundFinalVal=kwargs['ROUND_FINAL_WEIGHTED_PHI'])
						qAElectrons, qBElectrons, _ = inducedCharge(wPhiA, wPhiB, pathElectrons, q=qElectronArray, roundFinalVal=kwargs['ROUND_FINAL_WEIGHTED_PHI'])
						indChargeHoles = qAHoles - scale * qBHoles
						indChargeElectrons = qAElectrons - scale * qBElectrons
					else:
						_, _, indChargeHoles = inducedCharge(wPhiA, wPhiB, pathHoles, q=qHoleArray, roundFinalVal=kwargs['ROUND_FINAL_WEIGHTED_PHI'])
						_, _, indChargeElectrons = inducedCharge(wPhiA, wPhiB, pathElectrons, q=qElectronArray, roundFinalVal=kwargs['ROUND_FINAL_WEIGHTED_PHI'])
					allInduced.append(indChargeHoles)
					allInduced.append(indChargeElectrons)
				else:
					indChargeHoles = inducedChargeSingle(wPhi, pathHoles, q=qHoleArray, roundFinalVal=kwargs['ROUND_FINAL_WEIGHTED_PHI'])
					indChargeElectrons = inducedChargeSingle(wPhi, pathElectrons, q=qElectronArray, roundFinalVal=kwargs['ROUND_FINAL_WEIGHTED_PHI'])
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

	newCollection = gEventCollection(filename)
	newCollection.printInfo()

	event = newCollection.collection[126]
	cProc = event.GetHits()[1]['creatorProcess'].split('\x00')[0]
	# event.createCarriers()

	t, q = computeChargeSignal(event, emfilename, **settings)
	fig, ax= plt.subplots()

	ax.plot(t,-1*q, linewidth=3)
	ax.set_xlabel(r'Time ($\mu$s)', fontsize=14)
	ax.set_ylabel(r'Induced Charge (C)', fontsize=14)
	ax.set_title('Q(t) at a Unipolar electrode for event %i. Initial electron creator process: %s' % (event.GetEventID(), cProc) , fontsize=14)
	event.plotH1(x="z", y="energy", nbins=200)
	event.plotH2()


	# # iterate many times to watch fluctuation
	fig1, ax1 = plt.subplots()
	settings['CARRIER_GENERATION_POISSON'] = 1
	settings['CHARGE_DIFFERENCE'] = 1
	settings['VOLTAGE_SWEEP_INDEX'] = 4
	settings['SCALE_WEIGHTED_PHI'] = 0

	iterations = 50
	for i in range(iterations):
		print(i)

		t, q = computeChargeSignal(event, emfilename, **settings)
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

	print('stopping')

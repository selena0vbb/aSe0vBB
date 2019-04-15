import sys
import getopt
import numpy as np

sys.path.append("./carrierproduction")
print(sys.path)
import particledata as pd
import processedPulse as proc
import seleniumconfig as sc

def main(argv):
	""" Main program that acts as an entry point to the simulation """

	# Set default values
	configfile = "./carrierproduction/config.txt"
	nEvents = 0

	# Process command line arguments
	try:
		opts, args = getopt.getopt(argv, "hc:n:")
	except getopt.GetoptError:
		printhelp()
		sys.exit(2)

	for opt, arg in opts:
		if opt == "-h":
			printhelp()
			sys.exit()
		elif opt == "-c":
			configfile = arg
		elif opt == "-n":
			nEvents = int(arg)

	print("Config file is %s" % configfile)

	# Setup and run the simulation
	settings = sc.readConfigFile(configfile)
	particlefilename=settings['PARTICLE_FILENAME']
	emfilename=settings['EM_FILENAME']

	if nEvents <= 0:
		nEvents = pd.numberOfEvents(particlefilename)

	simObj = pd.CarrierSimulation(emfilename=emfilename, configfile=configfile)

	filesize=settings['NEVENTS_PER_FILE']
	outdir=settings['OUTPUT_DIR']
	outfilename=settings['OUTPUT_FILE']


    # Chunk event collection into smaller pieces if needed
	for j in range(int(np.ceil(float(nEvents)/filesize))):
		indx = []
		# Create new event collect with the smaller chunck size
		newEventCollection = pd.gEventCollection(particlefilename, eventCounterRange=[j*filesize, (j+1)*filesize-1])
		simObj.newEventCollection(newEventCollection)
		for i in range(len(newEventCollection.collection)):
			event = newEventCollection.collection[i]
			flat = event.flattenEvent()
			zr = np.max(np.abs(flat['z']))
			yr = np.max(np.abs(flat['y']))
			xr = np.max(np.abs(flat['x']))
			eevent = np.sum(flat['energy'])

			if xr < settings['X_MAX'] and yr < settings['Y_MAX'] and zr < settings['Z_MAX']:
				indx.append(i)	

		simOutput = simObj.processMultipleEvents(indx, processes=int(settings['NPROCESSORS']))
		simOutput.setGitInfo("./")
		simObj.outputfile = outfilename%j	
		simObj.outputdir = outdir
		simObj.savedata(simOutput)

	return

def printhelp():
	print("python selena.py -c <path to config file> -n <number to simulate>")

if __name__ == '__main__':
	main(sys.argv[1:])
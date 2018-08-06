import particledata as pd
import seleniumconfig as sc
import numpy as np
import matplotlib.pyplot as plt
import os
import brewer2mpl
# required files
filename = r"C:\Users\alexp\Documents\UW\Research\Selenium\aSe0vBB\particle\selenium-build\output\122_keV_testTuple.root"
emfilename = [r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_1um_spacing_fullsize.txt",
				r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_3um_spacing_fullsize.txt",
				r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_5um_spacing_fullsize.txt",
				r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_8um_spacing_fullsize.txt"]
configfilename = r"./config.txt"

print('Read settings')
settings = sc.readConfigFile(configfilename)

print('Read Geant4 Particle data')
newCollection = pd.gEventCollection(filename)
event = newCollection.collection[140]
creatorProc = event.GetHits()[1]['creatorProcess'].split('\x00')[0]

print('Calculating induced charge signal')
for file in emfilename:
	biasVoltIndex = [0,1,2,3,4]
	plt.rc('font', family='serif')
	fig, ax = plt.subplots()
	bmap = brewer2mpl.get_map('Set1', 'Qualitative', 5).mpl_colors
	biasString = ['120 V', '160 V', '200 V', '240 V', '280 V']
	for indx in biasVoltIndex:
		print(indx)
		settings['VOLTAGE_SWEEP_INDEX'] = indx
		time, q = pd.computeChargeSignal(event, file, **settings)
		ax.plot(time, -1*q, linewidth=2, color=bmap[indx], label=biasString[indx])

	ax.set_xlabel(r'Time ($\mu s$)', fontsize=14)
	ax.set_ylabel(r'Induced Charge (C)', fontsize=14)
	ax.set_title('Induced Charge Signal at Amplifier for %s Kapton layer' % file.split('\\')[-1].split('_')[3])
	ax.legend()


# Plot results


plt.show()

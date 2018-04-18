import plot
import numpy as np
import matplotlib.pyplot as plt
import os
import re


# Data analysis of varied width of coplanar electrode

# Find all files in a location
pathLocation = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data'

N = 100
label = []
yCoord = {'5':943, '10':965, '25':912, '50':1028, '100':847}
z = np.linspace(20, 220, N)

# Predefines array for all induced charge data
allInduced = np.zeros((N, len(os.listdir(pathLocation))))

# iterates over all files in the correct directory
for i,  file in enumerate(os.listdir(pathLocation)):
	coplanarSpacing = re.findall('\d+',file)
	label.append(coplanarSpacing[0] + r' $\mu m$') # Defining legend strings

	# import data from file
	header, data = plot.readComsolFile(os.path.join(pathLocation, file))

	y = np.repeat(yCoord[coplanarSpacing[0]], N)
	pos = np.concatenate((y,z)).reshape(2,N).T
	qA, qB, qDiff = plot.inducedCharge(data[:,(1,2,6)], data[:,(1,2,9)], pos, method='linear')

	# Putting the induced charge into 1 array for stored data
	allInduced[:,i] = qDiff[:,0]


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(z-20, allInduced, linewidth=3)
ax.set_title(r'Induced Charge vs Depth for a 2000x2000x2000 $\mu m$ Coplanar Electrode with Different Finger Spacing', fontsize=16)
ax.set_xlabel(r'Depth ($\mu m$)', fontsize=14)
ax.set_ylabel(r'Induced Charge at Coplanar Electrode (C)', fontsize=14)
ax.legend(label)
plt.show()
	
	
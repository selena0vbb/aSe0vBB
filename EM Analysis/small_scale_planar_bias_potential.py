# Small scale voltage bias analysis

import numpy as np 
import matplotlib.pyplot as plt 
from plot import *
import os

filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\bias_voltage\bias_voltage_small_scale_10umspace_0_360_range_unequal_electrode_size.txt'
bias = np.arange(0, 400, 40)
header, y, z, data = readComsolFileGrid(filename)
y, z = y*1e6, z*1e6
test = data[:,np.amin(data.shape)-2]
print(test.shape)
# Remove NaNs
locNan = np.argwhere(np.isnan(data))
locNanXy = np.argwhere(np.isnan(data[:,0]))
print(np.isnan(data[:,0]).shape)
data = np.nan_to_num(data, copy=False)
# data = np.delete(data, locNan)
# y, z = np.delete(y, locNanXy), np.delete(z, locNanXy)
phi, Ey, Ez = data[:,0:np.amin(data.shape):3], data[:,1:np.amin(data.shape):3], data[:,2:np.amin(data.shape):3] 
print(data.shape)
test = data[:,1:4]
ax = []


for i in range(np.amin(phi.shape)-2):
	# locNan = np.argwhere(np.isnan(phi[:,i]))
	# _, tempAxis = plotPhi((np.delete(y, locNan), np.delete(z, locNan)), np.delete(phi[:,i], locNan), i, type='contour')
	print(i)
	_, tempAxis = plotPhi((y, z), phi[:,i], i, type='contour')
	tempAxis.set_title(r'Coplanar Potential for $\Phi_A=0$ and $\Phi_B=$'+str(bias[i]), fontsize=16)
	tempAxis.set_xlabel(r'x ($\mu m$)', fontsize=14)
	tempAxis.set_ylabel(r'Depth ($\mu m$)', fontsize=14)
	ax.append(tempAxis)




# Looking at where the field lines end. Start a particle in a bunch of locations and look at where it ends by looking at the difference in induced charge
# If the particle ends at electrode A, the qDiff will be one sign, and electrode B will give qDiff of the other sign
ystart = np.linspace(15, 125, 100)
zstart = 29

vDriftHoles = 0.19e2 # um^2/(V*us)
dt = 0.001 # us
wPhiA = [y, z, phi[:,-2]]
wPhiB = [y, z, phi[:,-1]]

effectiveArea = []
for j in range(bias.size):
	qSign = []
	E = [y, z, Ey[:,j]/1e6, Ez[:,j]/1e6]
	for i in range(ystart.size):
		path = findMotion((ystart[i], zstart), E, vDriftHoles, dt, q=1.6e-19, limits=[10, 130, 1.5, 29.5])
		_, _, qdiff = inducedCharge(wPhiA, wPhiB, path[:,(0,1)], q=1.6e-19)
		if qdiff[-1] > 0:
			qSign.append(1)
		else:
			qSign.append(0)
	effectiveArea.append(1-sum(qSign)/len(qSign))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(bias, np.array(effectiveArea), linewidth=3)
ax.set_ylabel('Normalized Effective Area', fontsize=14)
ax.set_xlabel('Coplanar Bias Voltage (V)', fontsize=14)
ax.set_title('Effective Area of the Collecting Electrode as a Function of Bias Voltage', fontsize=16)

plt.show()
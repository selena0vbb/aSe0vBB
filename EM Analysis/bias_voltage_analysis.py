# Bias voltage between electrodes A and B on coplanar electrode analysis

import numpy as np 
import matplotlib.pyplot as plt 
from plot import *
import os
from scipy.ndimage.interpolation import zoom

# Import data
filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\bias_voltage\bias_voltage_10umspace_0_200_range.txt'

header, y, z, phi = readComsolFileGrid(filename)

# Split potential into weighted potentials and potentials under real electrode biases
weightPhiA, weightPhiB = phi[:,-3], phi[:,-6]
data = phi[:,np.arange(0,np.amin(phi.shape)-6)]


ax = plotPhi((y,z), weightPhiA, 0, type='contour')[1]


# Populating a more dense grid and doing the plot
# zoomFactor = 4
# xi = np.linspace(0, 1000, y.size*zoomFactor)
# yi = np.linspace(20, 220, z.size*zoomFactor)
# Phii = zoom(np.reshape(weightPhiA,(z.size, y.size)), zoomFactor, order=3)
# ax1 = plotPhi((xi,yi), Phii.flatten(), 0, type='contour')[1]


# Solving the equation of motion for a give potential
E = [y, z, data[:,1]/1e6, data[:,2]/1e6] # Setting up the E field input to finding the motion. Convert E to V/um

# Parameters needed for motion
vDriftHoles = 0.19e2 # um^2/(V*us)
dt = 0.0005 # us
yrange = (100, 500) # um
zrange = (20, 220) # um
# initialPos = (np.random.uniform(yrange[0], yrange[1]), np.random.uniform(zrange[0], 0.3*zrange[1]))
initialPos = (244, 50)
particleMotion = findMotion(initialPos, E, vDriftHoles, dt, q=1.6e-19)


fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.scatter(particleMotion[:,0], particleMotion[:,1], linewidth=2)
ax2.set_ylabel(r'Depth ($\mu m$', fontsize=14)
ax2.set_xlabel(r'Width ($\mu m$', fontsize=14)
ax2.set_title(r'Position of Charged Particle. Va=0, Vb=0', fontsize=16)
ax2.set_ylim([zrange[0], zrange[1]])
ax2.set_xlim([yrange[0], yrange[1]])

# Compute induced charge based on position
qa, qb, qdiff = inducedCharge([y, z, weightPhiA], [y, z, weightPhiB], particleMotion[:,(0,1)], q=1.6e-19)
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.plot(particleMotion[:,1], qdiff, linewidth=2)
ax3.set_xlabel(r'Depth ($\mu m$', fontsize=14)
ax3.set_ylabel(r'Induced Charge (C)', fontsize=14)
ax3.set_title(r'Induced Charge Difference. Va=0, Vb=0', fontsize=16)

x = np.linspace(75,375,300)
zStart = 50
electrode = []
for i in range(x.size):
	xzPath = findMotion((x[i], zStart), E, vDriftHoles, dt, q=1.6e-19)
	_, _, qdiff = inducedCharge([y, z, weightPhiA], [y, z, weightPhiB], xzPath[:,(0,1)], q=1.6e-19)
	if qdiff[-1] > 0:
		electrode.append(1)
	else:
		electrode.append(0)

fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
ax6.plot(x, np.array(electrode), linewidth=2)


# Monte Carlo Simulation. Finding the induced charge for a variety of different initial conditions
N = 50
qDiffAll = []
fig4 = plt.figure()
ax4 = fig4.add_subplot(111)

# Finding the induced charge for a variety of intial positions
for i in range(N):
	initialPos = (np.random.uniform(1.1*yrange[0], 0.9*yrange[1]), np.random.uniform(zrange[0], 0.5*zrange[1]))
	particleMotion = findMotion(initialPos, E, vDriftHoles, dt, q=1.6e-19)
	try:
		qa, qb, qdiff = inducedCharge([y, z, weightPhiA], [y, z, weightPhiB], particleMotion[:,(0,1)], q=1.6e-19)
	except:
		print(particleMotion.shape)
	ax4.plot(particleMotion[:,1], qdiff, 'k', linewidth=0.5)

ax4.set_title('Difference in Induced Charged for N='+str(N)+' Cases', fontsize=16)
ax4.set_ylabel(r'Induced Charge (C)', fontsize=14)
ax4.set_xlabel(r'Depth ($\mu m$', fontsize=14)

# Finding which electrode each charged particle
# Monte carlo with M different intial positions
M = 100
fracAtElectrode = []
biasVoltage =  [0, 40, 80, 120, 160, 200]

for j, volt in enumerate(biasVoltage):
	countHoles =  0 
	for i in range(M):
		initialPos = (np.random.uniform(1.1*yrange[0], 0.9*yrange[1]), np.random.uniform(zrange[0], 0.5*zrange[1]))
		particleMotion = findMotion(initialPos, [y, z, data[:,1+3*j]/1e6, data[:,2+3*j]/1e6], vDriftHoles, dt, q=1.6e-19)
		try:
			qa, qb, qdiff = inducedCharge([y, z, weightPhiA], [y, z, weightPhiB], particleMotion[:,(0,1)], q=1.6e-19)
		except:
			print(particleMotion.shape)
		print(qdiff[-1])
		if qdiff[-1] > 0:
			countHoles = countHoles + 1

	# Calculates the fraction of that hit one of the electrodes
	fracAtElectrode.append(countHoles/M)

fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
ax5.plot(np.array(biasVoltage), np.array(fracAtElectrode), linewidth=3)
ax5.set_title('Fraction of Holes Arriving at Electrode A', fontsize=16)
ax5.set_ylabel(r'Fraction', fontsize=14)
ax5.set_xlabel(r'Voltage Bias (V)', fontsize=14)


plt.show()


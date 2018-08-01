# Weighted potential of the tiered design with Polymide kapton spacing

from plot import *
import matplotlib.pyplot as plt
import numpy as np


folder = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data'
filelist = ['kapton_layer_analysis_1um_spacing_fullsize.txt',
			'kapton_layer_analysis_3um_spacing_fullsize.txt',
			'kapton_layer_analysis_5um_spacing_fullsize.txt',
			'kapton_layer_analysis_8um_spacing_fullsize.txt']

# filelist = ['kapton_layer_analysis_1um_spacing_fullsize.txt']

figList = []
axList = [] 
for findx, file in enumerate(filelist):
	_, y, z, data = readComsolFileGrid(folder + r'/' + file)
	y, z = y*1e6, z*1e6
	phi = data[:,0:np.amin(data.shape):3]
	wPhiA = phi[:,-2]
	wPhiB = phi[:,-1] 

	plt.rc('font', family='serif')

	# Regular zoom
	fig, (ax1, ax2) = plt.subplots(1,2, sharey=True)
	fig, ax1, cax1 = plotPhi([y,z], wPhiA, 1, figH=(fig, ax1))
	fig, ax2, cax2 = plotPhi([y,z], wPhiB, 1, figH=(fig, ax2))

	filename = file.split('_')
	spacingText = filename[3] + ' ' + filename[4]

	fig.suptitle('Weighted Potential for ' + spacingText + ' Kapton', fontsize=16)
	ax1.set_title(r'Collection Electrode', fontsize=14)
	ax1.set_xlabel(r'y ($\mu m$)', fontsize=14)
	ax1.set_ylabel(r'Depth ($\mu m$)', fontsize=14)

	ax2.set_title(r'Control Electrode', fontsize=14)
	ax2.set_xlabel(r'y ($\mu m$)', fontsize=14)

	fig.colorbar(cax1, ax=ax1)
	fig.colorbar(cax2, ax=ax2)

	# zoomed in
	fig2, (ax3, ax4) = plt.subplots(1,2, sharey=True)
	fig, ax3, cax3 = plotPhi([y,z], wPhiA, 1, figH=(fig2, ax3))
	fig, ax4, cax4 = plotPhi([y,z], wPhiB, 1, figH=(fig2, ax4))

	fig2.suptitle('Weighted Potential for ' + spacingText + ' Kapton', fontsize=16)
	ax3.set_title(r'Collection Electrode', fontsize=14)
	ax3.set_xlabel(r'y ($\mu m$)', fontsize=14)
	ax3.set_ylabel(r'Depth ($\mu m$)', fontsize=14)
	ax3.set_ylim([np.amin(z), -80])
	ax3.set_xlim([-125, 125])

	ax4.set_title(r'Control Electrode', fontsize=14)
	ax4.set_xlabel(r'y ($\mu m$)', fontsize=14)
	ax4.set_ylim([np.amin(z), -75])
	ax4.set_xlim([-100, 100])
	fig2.colorbar(cax3, ax=ax3)
	fig2.colorbar(cax4, ax=ax4)

	# Plotting a 1D slice of weighted potential
	wPhiAReshape = np.reshape(wPhiA,(z.size, y.size))
	wPhiBReshape = np.reshape(wPhiB,(z.size, y.size))

	# Get nearest index to 5 um on the y axis
	if findx == 0:
		ySlice = 25
	else:
		ySlice = 5
	yIndx = (np.abs(y-ySlice).argmin())
	fig3, ax5 = plt.subplots()
	phiASlice = wPhiAReshape[:,yIndx]
	phiBSlice = wPhiBReshape[:,yIndx]


	slopeA = (phiASlice[-1]-phiASlice[100])/(z[-1]-z[100])
	slopeB = (phiBSlice[-1]-phiBSlice[100])/(z[-1]-z[100])

	print(slopeA/slopeB)

	ax5.plot(-z, wPhiAReshape[:,yIndx], 'k:', -z, wPhiBReshape[:,yIndx], 'k:', -z, wPhiAReshape[:,yIndx]-wPhiBReshape[:,yIndx], 'k', linewidth=3)
	ax5.plot(-z, phiBSlice*slopeA/slopeB, 'r:', -z, phiASlice - phiBSlice*slopeA/slopeB, 'r', linewidth=3)
	ax5.set_xlabel(r'Depth ($\mu m$)', fontsize=14)
	ax5.set_ylabel(r'Weighted Potentail', fontsize=14)
	ax5.set_title('Y Slice of Weighted Potential and Difference for ' + spacingText, fontsize=16)
	ax5.set_xlim([-100,100])
	ax5.set_ylim([0,1])
	lines = ax5.get_lines()
	ax5.legend((lines[2], lines[4]), ['Not Scaled', 'Scaled'])
	ax5.text(-75, 0.3, 'Scaling Coefficient = %0.2f' % (slopeA/slopeB))

	# figList.extend((fig, fig2, fig3))
	# axList.extend((ax1, ax2, ax3, ax4, ax5))

plt.show()
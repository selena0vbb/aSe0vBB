import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as scp
import re
from mpl_toolkits.mplot3d import axes3d

def readComsolFile(filename):
	"""
	Reads E&M simulation files into data structures for plotting

	Inputs:
		filename - path to Comsol data file
	Outputs:
		header - list of all the Comsol header information such as the model used, dimension, units, functions exported, etc.
		data - an NxM numpy array of the data contained in the Comsol output file. N is the number of points in the model and M is 
		the number of dimensions + the number of functions
	"""
	header, data = ([],[],)
	file = open(filename,'r')
	for lines in file:

		# Seperate into header and data based of leading % sign. Cleanup of uncessesary symbols
		if lines[0] == '%':
			header.append(re.sub('[\n%]','',lines.replace(' ', '')))	
		else:
			# Split the line into individual values and convert to floats
			splitString = re.sub('[\n]','',lines).split() 
			data.append([float(i)for i in splitString]) # Rounds the input data to avoid some of Comsol variations
		
	file.close()
	
	return header, np.array(data)


def plotSlice(data, sliceIdx, sliceVal, funcIdx, eps=1e-12, gridSize=100, type='contour'):
	"""
	Plots of 2D slice of 3D Comsol Data 

	Inputs:
		data - NxM numpy array of Comsol data organized in columns of axis position and then columns of the function values
		sliceIdx - axis index that you want to take a slice. Numerical index that should correspond to x, y, or z axis
		sliceVal - the value of the index you want a slice upon. Will return an empty data set and raise and error if no
			data found with that value
		funcIdx - index of the function to be plotted in the data set (column number)
		eps - parameter to select similar values in a slice The way that Comsol sets up the grid means there are very
			small variations even within what we would consider one slice. The eps makes it so we include everything
			within the range of sliceVal +/- eps
	Outputs:
		figHandle - handle to the python plot created by the function

	"""

	# Find all values of the sliceVal within the sliceAxis. Use the indices that satisfy this condition for the rows
	# and the remaining axis (x, y, or z) plus the function for the index of the columns
	sliceAxis = data[:,sliceIdx]
	columns = [ax for ax in range(0,3) if ax != sliceIdx]
	columns.append(funcIdx)
	sliceIndices = np.nonzero((sliceAxis < sliceVal + eps) & (sliceAxis > sliceVal - eps))

	# Create a grid from the indices calculated and apply it to the data to get the correct slice
	idxR, idxC = np.meshgrid(sliceIndices, columns)
	slicedData = data[idxR, idxC]

	# Define a regular grid that we will interpolate the data to
	xi = np.linspace(np.amin(slicedData[0,:]), np.amax(slicedData[0,:]), gridSize)
	yi = np.linspace(np.amin(slicedData[1,:]), np.amax(slicedData[1,:]), gridSize)
	val = scp.griddata((slicedData[0,:], slicedData[1,:]), slicedData[2,:], (xi[None,:], yi[:,None]), method='linear')

	# Plot the interpolated values
	if type is 'contour':
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.contour(xi, yi, val, 10, linewidth=0.5, colors='k')
		cax = ax.contourf(xi, yi, val, 10, cmap=plt.cm.jet)
		fig.colorbar(cax)
		return fig, ax

	elif type is 'mesh':
		X, Y = np.meshgrid(xi,yi)
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1, projection='3d')
		ax.plot_wireframe(X, Y, val)
		return fig, ax

	else:
		print('Error: a type of %s is not a valid input', type)
		return None

def findMotion(xi, E, vDrift, dt):
	"""
	Computes the position as a function of time (and therefore can be used with the weighted potential) using the drift velocity and E fields of the charge particle
	The units must be consistent. For example the units of xi and the distance in vDrift must be the same. Same goes with xi and electric field.

	Inputs
		xi - initial position. Tuple of 3 points, x,y,z
		E - the electric field at all point in the model. Nx6 matrix where N is the number of different grid points and 6 corresponds to x,y,z,Ex,Ey,Ez to fully describe the vector field
		vDrift - drift velocity
		dt = time step
	Outputs
		xt - the position of the function as a function of time. Nx4 matrix where N is the number of time steps and columns are x,y,z,t
	"""
	return None

def inducedCharge(wPotentialA, wPotentialB, path, q=1.6e-19):
	"""
	Finds the induced charge at each electrode given a path of the the charged particle

	Inputs:
		wPotentialA - the weighted potential (V) from the A electrode. Nx4 numpy array
		wPotentialB - the weighted potential (V) from the B electrode. Nx4 numpy array
		path - the position to compute the weighted potential at. Nx1 numpy array
		q - charge of the particle
	Outputs:
		qA - charge induced at electrode A
		qB - charge induced at electrode B
		qDiff - difference in the charge induced at electrode A and B
	"""
	return None

if __name__ == '__main__':
	filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\test_weighted_potential.txt'
	testHeader, testData = readComsolFile(filename)
	
	# Creating contour and wireframe plot
	print('Test plot function\n')
	figC, axC = plotSlice(testData, 0, 1000, 4, gridSize=1000, type='contour')
	figM, axM = plotSlice(testData, 0, 1000, 4, gridSize=1000, type='mesh')
	# Contour plot
	axC.set_title('Contour Plot of the Weighted Potential for a Coplanar Selenium Detector', fontsize=16)
	axC.set_xlabel(r"Width ($\mu m$)", fontsize=14)
	axC.set_ylabel(r"Depth ($\mu m$)", fontsize=14)
	# Wireframe Plot
	axM.set_title('Wireframe Plot of the Weighted Potential for a Coplanar Selenium Detector', fontsize=16)
	axM.set_xlabel("\n"+r"Width ($\mu m$)", fontsize=14)
	axM.set_ylabel("\n"+r"Depth ($\mu m$)", fontsize=14)
	axM.set_zlabel(r"Weighted Potential (V)", fontsize=14)
	plt.show()
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as scp
import re
from mpl_toolkits.mplot3d import axes3d
import brewer2mpl

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

def readComsolFileGrid(filename):
	"""
	Reads E&M simulation files from Comsol that are output onto a regular grid. Returns data in numpy array
	The outputs from this function are in a format such that interp2D accepts them correctly

	Inputs:
		filename - path to Comsol data file
	Outputs:
		header - list of all the Comsol header information
		x - x positions in numpy array
		y - y poisition in numpy array
		data - Comsol data. NxM where N is len(x)*len(y) and M is the number of functions in the output file
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
	data = np.array(data)
	
	func = data[:,np.arange(2,len(splitString))] # Seperates imported data into position and function
	x, y = np.unique(data[:,0]), np.unique(data[:,1])
	
	return header, x, y, func

def takeSlice(data, sliceIdx, sliceVal, funcIdx, eps=1e-12):
	"""
	Performs the same action as Plot Slice except returns the data instead of the figure
	
	"""

	# Find all values of the sliceVal within the sliceAxis. Use the indices that satisfy this condition for the rows
	# and the remaining axis (x, y, or z) plus the function for the index of the columns
	sliceAxis = data[:,sliceIdx]
	columns = [ax for ax in range(0,3) if ax != sliceIdx]
	columns.extend(funcIdx)
	sliceIndices = np.nonzero((sliceAxis < sliceVal + eps) & (sliceAxis > sliceVal - eps))

	# Create a grid from the indices calculated and apply it to the data to get the correct slice
	idxR, idxC = np.meshgrid(sliceIndices, columns)
	slicedData = data[idxR, idxC]

	return slicedData


def plotPhi(pos, data, funcIdx, type='contour', figH=None, color='Reds'):
	"""
	Plots the potential (phi) of a 2D slice of comsol simulation

	Inputs:
		data - NxM numpy array of Comsol data organized in columns of axis position and then columns of the function values
		funcIdx - index of the function to be plotted in the data set (column number)
		eps - parameter to select similar values in a slice The way that Comsol sets up the grid means there are very
			small variations even within what we would consider one slice. The eps makes it so we include everything
			within the range of sliceVal +/- eps
		figH - tuple with (fig, ax)
	Outputs:
		figHandle - handle to the python plot created by the function

	"""

	# Take the specific column of the funcIdx
	if data.ndim > 1:
		phi = np.reshape(data[:,funcIdx],(pos[1].size, pos[0].size))
	else:
		phi = np.reshape(data,(pos[1].size, pos[0].size))

	if not figH:
		fig = plt.figure()
		ax = fig.add_subplot(111)
	else:
		fig = figH[0]
		ax = figH[1]

	# Plot the interpolated values
	if type is 'contour':
		ax.contour(pos[0], pos[1], phi, 8, linewidth=0.5, colors='k')
		cax = ax.contourf(pos[0], pos[1], phi, 8, cmap=brewer2mpl.get_map(color, 'sequential', 8).mpl_colormap, vmin=np.amin(phi), vmax=np.amax(phi))
		return fig, ax, cax

	elif type is 'mesh':
		X, Y = np.meshgrid(pos[0],pos[1])
		# ax = fig.add_subplot(1,1,1, projection='3d')
		ax.plot_wireframe(X, Y, phi)
		return fig, ax

	else:
		print('Error: a type of %s is not a valid input', type)
		return None

def findMotion(xi, E, vDrift, dt, method='linear', q=-1.6e-19, limits=[]):
	"""
	Computes the position as a function of time (and therefore can be used with the weighted potential) using the drift velocity and E fields of the charge particle
	The units must be consistent. For example the units of xi and the distance in vDrift must be the same. Same goes with xi and electric field.
	Computation done on a 2D slice (assumes translational invariance in one direction)

	Inputs
		xi - initial position. Tuple of 2 points, x,y
		E - the electric field at all point in the model. Nx4 matrix where N is the number of different grid points and 4 corresponds to x,y,Ex,Ey to fully describe the vector field
		vDrift - drift velocity (length^2/(V*time))
		dt = time step
	Outputs
		xt - the position of the function as a function of time. Nx3 matrix where N is the number of time steps and columns are x,y,t
	"""
	
	# Defining coordinates and finding the max and min potentials values
	t = 0
	x = xi[0]
	y = xi[1]
	if not limits:
		xmin, xmax = np.amin(E[0]), np.amax(E[0])
		ymin, ymax = np.amin(E[1]), np.amax(E[1])
	else:
		xmin, xmax = limits[0], limits[1]
		ymin, ymax = limits[2], limits[3]

	xt = [] 
	
	# Create interpolating functions for the E fields
	ExInter = scp.interp2d(E[0], E[1], np.reshape(E[2],(E[1].size, E[0].size)), kind=method) 
	EyInter = scp.interp2d(E[0], E[1], np.reshape(E[3],(E[1].size, E[0].size)), kind=method)

	# While the charge carrier is in the selenium, keep finding the position
	while(x < xmax and x > xmin and y < ymax and y > ymin):
		xt.append([x,y,t])

		# Interpolate for values of Ex and Ey at the specific position
		Ex = ExInter(x, y)
		Ey = EyInter(x, y)
	
		# Solve equation of motion 
		xNext = x + vDrift*Ex*np.sign(q)*dt
		yNext = y + vDrift*Ey*np.sign(q)*dt

		# Assign the new version of x, y, and t
		x = xNext
		y = yNext
		t = t + dt

	# Add value at the limit
	if x > xmax:
		xt.append([xmax, y, t])
	elif x < xmin:
		xt.append([xmin, y, t])
	elif y > ymax:
		xt.append([x, ymax, t])
	else:
		xt.append([x, ymin, t])
		
	# convert to numpy array and return the data
	return np.array(xt)

def plotEField(Efield):
	"""
	Makes a vector field plot of the electric field using matplotlib 

	Inputs:
		Efield - list in the form [x, y, Ex, Ey] just like the other functions.
				x - Nx1 numpy array of x data points
				y - Mx1 numpy array of y data points
				Ex - (M*N)x1 numpy array corresponding to an Ex at every point on the grid
				Ey - (M*N)x1 numpy array corresponding to an Ey at every point on the grid
	Outputs:
		ax - handle to the axis of the created figure
	"""

	# Paceholder code
	return None


def interpEField2D(x, y, E, method='linear'):
	"""
	Wrapper function for interpolating points of the Efield over the Comsol grid
	
	Inputs:
		x - single value or list of x positions
		y - single value or list of y positions
		E - list containing x, y grid data and E field data
	Outputs:
		EInterp - interpolated E field. Nx2 numpy array where N is the number of xy coordinate pairs
	"""

	
	# Create interpolating functions for the E fields
	ExInter = scp.interp2d(E[0], E[1], np.reshape(E[2],(E[1].size, E[0].size)), kind=method) 
	EyInter = scp.interp2d(E[0], E[1], np.reshape(E[3],(E[1].size, E[0].size)), kind=method)

	# Create output array of zeros
	EInterp = np.zeros((x.size, 2))

	for i in range(x.size):
		EInterp[i, 0] = ExInter(x[i], y[i])[0]
		EInterp[i, 1] = EyInter(x[i], y[i])[0]

	return EInterp

def inducedChargeSingle(wPotential, path, q=1.6e-19, method='linear'):

	qi = []

	wPInter = scp.interp2d(wPotential[0], wPotential[1], np.reshape(wPotential[2],(wPotential[1].size, wPotential[0].size)), kind=method)
	
	for i in range(path.shape[0]):

		wP = wPInter(path[i,0], path[i,1])
		qi.append(-q*wP[0])

	return np.array(qi) 

def inducedCharge(wPotentialA, wPotentialB, path, q=-1.6e-19, method='linear'):
	"""
	Finds the induced charge at each electrode given a path of the the charged particle

	Inputs:
		wPotentialA - list including [x, y, Phi]. The x,y position pairs and the potential occuring at this. For electrode A
		wPotentialB - list including [x, y, Phi]. The x,y position pairs and the potential occuring at this. For electrode B
		path - the position to compute the weighted potential at. Nx2 (x,y position at different time steps) numpy array
		q - charge of the particle
	Outputs:
		qA - charge induced at electrode A
		qB - charge induced at electrode B
		qDiff - difference in the charge induced at electrode A and B
	"""
	qA = []
	qB = []	

	# Definte interplation functions
	VaInter = scp.interp2d(wPotentialA[0], wPotentialA[1], np.reshape(wPotentialA[2],(wPotentialA[1].size, wPotentialA[0].size)), kind=method)
	VbInter = scp.interp2d(wPotentialB[0], wPotentialB[1], np.reshape(wPotentialB[2],(wPotentialA[1].size, wPotentialA[0].size)), kind=method)


	# Iterated over all the positions in the path
	for i in range(path.shape[0]):
		Va = VaInter(path[i,0], path[i,1])
		Vb = VbInter(path[i,0], path[i,1])
		# print(Va, path[i,0], path[i,1])
		# Find the q induced via the Shokley-Ramo Theorem
		qA.append(-q*Va[0])
		qB.append(-q*Vb[0])

	qA, qB = np.array(qA), np.array(qB)
	
	return qA, qB, qA-qB


if __name__ == '__main__':
	filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\test_weighted_potential.txt'
	testHeader, testData = readComsolFile(filename)
	
	# # Creating contour and wireframe plot
	# print('Test plot function\n')
	# figC, axC = plotSlice(testData, 0, 1000, [6], gridSize=1000, type='contour')
	# figM, axM = plotSlice(testData, 0, 1000, [6], gridSize=1000, type='mesh')
	# # Contour plot
	# axC.set_title('Contour Plot of the Weighted Potential for a Coplanar Selenium Detector', fontsize=16)
	# axC.set_xlabel(r"Width ($\mu m$)", fontsize=14)
	# axC.set_ylabel(r"Depth ($\mu m$)", fontsize=14)
	# # Wireframe Plot
	# axM.set_title('Wireframe Plot of the Weighted Potential for a Coplanar Selenium Detector', fontsize=16)
	# axM.set_xlabel("\n"+r"Width ($\mu m$)", fontsize=14)
	# axM.set_ylabel("\n"+r"Depth ($\mu m$)", fontsize=14)
	# axM.set_zlabel(r"Weighted Potential (V)", fontsize=14)
	# # plt.show()

	initialPos = (1000, 100) # in um
	vDriftHoles = 0.19e6 # cm^2/(V*s)

	eField = takeSlice(testData, 0, 1000, [4,5])/(1e6) # Converting V/m to V/um
	# xt = findMotion(initialPos, eField.T,  vDriftHoles, 0.001)

	# For finding the current induced (vs the signal generated), we can more easily just set y to a fixed value
	# and vary z across the range of depth. If we choose y to be 
	N = 500
	z = np.linspace(21,219,N)
	y = np.repeat(443,N)
	pos = np.concatenate((y,z)).reshape(2,N).T

	qA, qB, qDiff = inducedCharge(testData[:,(1,2,6)], testData[:,(1,2,9)], pos)
	
	

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(pos[:,1]-20,qA,'--', pos[:,1]-20,qB,'--', pos[:,1]-20, qDiff, linewidth=3)
	ax.set_title(r'Induced Charge at a Coplanar Electrode. 100 $\mu m$ Spacing between fingers.', fontsize=16)
	ax.set_xlabel(r'Depth ($\mu m$)', fontsize=14)
	ax.set_ylabel(r'Induced Charge (C)', fontsize=14)
	ax.legend(['qA','qB','qDiff'])
	
	
	gridFile = r'C:\Users\alexp\Documents\UW\Research\Selenium\test_export.txt'
	header, x, y, V = readComsolFileGrid(gridFile)
	print(V.shape)
	xx, yy = np.meshgrid(x,y)
	z = np.reshape(V, (x.size, y.size))
	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111)
	ax2.contourf(x, y, z)
	plt.show()


	
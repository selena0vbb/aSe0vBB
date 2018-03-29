import numpy as np
import matplotlib.pyplot as plt
import re


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
			data.append([float(i) for i in splitString])
		
	file.close()
	
	return header, np.array(data)


def plotSlice(data, sliceIdx, sliceVal, funcIdx, eps=1e-15):
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


	# 2D interpolation of a regular grid against the known Comsol data so that we can perform a colormap

def 2DInterpolate(rawData, x, y):
	"""
	Takes in the regular grid we want interpolation points of (x, y) and interpolates them against the known results.
	Returns a regular grid of the function

	"""

if __name__ == '__main__':
	filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\test_export.txt'
	testHeader, testData = readComsolFile(filename)
	
	print(len(testHeader))
	print(len(testData))
	

	print('Test plot function\n')
	plotSlice(testData, 0, 0, 3)
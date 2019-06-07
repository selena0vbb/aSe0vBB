import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import palettable
except:
    pass
import scipy.interpolate as scp
import re
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
    header, data = ([], [])
    file = open(filename, "r")
    for lines in file:

        # Seperate into header and data based of leading % sign. Cleanup of uncessesary symbols
        if lines[0] == "%":
            header.append(re.sub("[\n%]", "", lines.replace(" ", "")))
        else:
            # Split the line into individual values and convert to floats
            splitString = re.sub("[\n]", "", lines).split()
            data.append(
                [float(i) for i in splitString]
            )  # Rounds the input data to avoid some of Comsol variations

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
    header, data = ([], [])
    file = open(filename, "r")
    for lines in file:

        # Seperate into header and data based of leading % sign. Cleanup of uncessesary symbols
        if lines[0] == "%":
            header.append(re.sub("[\n%]", "", lines.replace(" ", "")))
        else:
            # Split the line into individual values and convert to floats
            splitString = re.sub("[\n]", "", lines).split()
            data.append(
                [float(i) for i in splitString]
            )  # Rounds the input data to avoid some of Comsol variations
    file.close()
    data = np.array(data)

    func = data[
        :, np.arange(2, len(splitString))
    ]  # Seperates imported data into position and function
    x, y = np.unique(data[:, 0]), np.unique(data[:, 1])

    return header, x, y, func


def readComsolFileGrid3d(filename):
    header, data = [], []
    file = open(filename, "r")

    # read each line and populate data array
    for lines in file:

        # Separate into header and data based of leading % sign.
        if lines[0] == "%":
            header.append(re.sub("[\n%]", "", lines.replace(" ", "")))
        else:
            # Split the line into individual values and convert to floats
            splitString = re.sub("[\n]", "", lines).split()
            data.append([float(i) for i in splitString])

    file.close()
    data = np.array(data)
    func = data[:, np.arange(3, len(splitString))]
    x, y, z = np.unique(data[:, 0]), np.unique(data[:, 1]), np.unique(data[:, 2])

    return header, (x, y, z), func


def takeSlice(data, sliceIdx, sliceVal, funcIdx, eps=1e-12):
    """
	Performs the same action as Plot Slice except returns the data instead of the figure

	"""

    # Find all values of the sliceVal within the sliceAxis. Use the indices that satisfy this condition for the rows
    # and the remaining axis (x, y, or z) plus the function for the index of the columns
    sliceAxis = data[:, sliceIdx]
    columns = [ax for ax in range(0, 3) if ax != sliceIdx]
    columns.extend(funcIdx)
    sliceIndices = np.nonzero(
        (sliceAxis < sliceVal + eps) & (sliceAxis > sliceVal - eps)
    )

    # Create a grid from the indices calculated and apply it to the data to get the correct slice
    idxR, idxC = np.meshgrid(sliceIndices, columns)
    slicedData = data[idxR, idxC]

    return slicedData


def plotPhi(pos, data, funcIdx, type="contour", figH=None, color=None, ncontour=9):
    """
	Plots the potential (phi) of a 2D slice of comsol simulation

	Inputs:
		pos - (x,y) data
		data - comsol data in regular grid format. If multiple columns supplied, selcts one with funcIdx
		funcIdx - index of the function to be plotted in the data set (column number)
		figH - tuple with (fig, ax)
		color - brewer2mpl color to use
		colortype - brewer2mpl color scheme to use. sequential, diverging, or qualitative
	Outputs:
		figHandle - handle to the python plot created by the function

	"""

    # Take the specific column of the funcIdx
    if data.ndim > 1:
        phi = np.reshape(data[:, funcIdx], (pos[1].size, pos[0].size))
    else:
        phi = np.reshape(data, (pos[1].size, pos[0].size))

    if not figH:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = figH[0]
        ax = figH[1]

    # Plot the interpolated values

    if color == None:
        color = brewer2mpl.get_map("Reds", "sequential", ncontour).mpl_colormap
    if type is "contour":
        ax.contour(pos[0], pos[1], phi, ncontour, linewidth=0.5, colors="k")
        cax = ax.contourf(
            pos[0],
            pos[1],
            phi,
            ncontour,
            cmap=color,
            vmin=np.amin(phi),
            vmax=np.amax(phi),
        )
        return fig, ax, cax

    elif type is "mesh":
        X, Y = np.meshgrid(pos[0], pos[1])
        # ax = fig.add_subplot(1,1,1, projection='3d')
        ax.plot_wireframe(X, Y, phi)
        return fig, ax

    else:
        print ("Error: a type of %s is not a valid input", type)
        return None


def findMotion(
    xi, E, vDrift, dt, totalTime=100, method="linear", q=-1.6e-19, limits=[]
):
    """
	Computes the position as a function of time (and therefore can be used with the weighted potential) using the drift velocity and E fields of the charge particle
	The units must be consistent. For example the units of xi and the distance in vDrift must be the same. Same goes with xi and electric field.
	Computation done on a 2D slice (assumes translational invariance in one direction)

	Inputs
		xi - initial position. Tuple of 2 or 3 points, x,y,z
		E - the electric field at all point in the model. Nx4 or Nx6 matrix where N is the number of different grid points and 4 (6) corresponds to the 2 (3) positions and their corresponding E fields
		vDrift - drift velocity (length^2/(V*time))
		dt - time step of the motion
		totalTime - maximum time allowed to calculated the motion over. Prevents code getting stuck when E field is 0
		method - interpolation method
		q - fundamental unit of charge. The sign is used to determine what direction the charges go
		limits - physical (or self imposed) limits of the detector geometry. Once it goes out of bounds, the code stops and the motion is returned
	Outputs
		xt - the position of the function as a function of time. Nx3 matrix where N is the number of time steps and columns are x,y,t
	"""

    # Defining coordinates and finding the max and min potentials values
    t = 0
    x = xi[0]
    y = xi[1]
    if len(xi) == 3:
        threeDim = True
        z = xi[2]
    else:
        threeDim = False

    if not limits:
        xmin, xmax = np.amin(E[0]), np.amax(E[0])
        ymin, ymax = np.amin(E[1]), np.amax(E[1])
        if threeDim:
            zmin, zmax = np.amin(E[2]), np.amax(E[2])
    else:
        xmin, xmax = limits[0], limits[1]
        ymin, ymax = limits[2], limits[3]
        if threeDim:
            zmin, zmax = limits[4], limits[5]

    # Check to ensure the initial points is within the limits. Otherwise return None
    if x > xmax or x < xmin or y > ymax or y < ymin:
        return np.array([])

    xt = []

    # Create interpolating functions for the E fields
    if threeDim:
        ExInter = scp.RegularGridInterpolator(
            (E[0], E[1], E[2]),
            np.reshape(E[3], (E[0].size, E[1].size, E[2].size), order="F"),
            bounds_error=False,
            fill_value=None,
        )
        EyInter = scp.RegularGridInterpolator(
            (E[0], E[1], E[2]),
            np.reshape(E[4], (E[0].size, E[1].size, E[2].size), order="F"),
            bounds_error=False,
            fill_value=None,
        )
        EzInter = scp.RegularGridInterpolator(
            (E[0], E[1], E[2]),
            np.reshape(E[5], (E[0].size, E[1].size, E[2].size), order="F"),
            bounds_error=False,
            fill_value=None,
        )
    else:
        ExInter = scp.RegularGridInterpolator(
            (E[0], E[1]),
            np.reshape(E[2], (E[0].size, E[1].size), order="F"),
            bounds_error=False,
            fill_value=None,
        )
        EyInter = scp.RegularGridInterpolator(
            (E[0], E[1]),
            np.reshape(E[3], (E[0].size, E[1].size), order="F"),
            bounds_error=False,
            fill_value=None,
        )

    # While the charge carrier is in the selenium, keep finding the position
    inBound = True
    stepTotal = float(totalTime) / dt
    stepCounter = 0
    while inBound and stepCounter < stepTotal:
        stepCounter += 1
        if threeDim:
            xt.append([x, y, z, t])

            # Interpolate for values of Ex and Ey at the specific position
            Ex = ExInter([x, y, z])
            Ey = EyInter([x, y, z])
            Ez = EzInter([x, y, z])

            # Solve equation of motion
            xNext = x + vDrift * Ex[0] * np.sign(q) * dt
            yNext = y + vDrift * Ey[0] * np.sign(q) * dt
            zNext = z + vDrift * Ez[0] * np.sign(q) * dt

            # Assign the new version of x, y, and t
            x = xNext
            y = yNext
            z = zNext

            inBound = (
                x < xmax
                and x > xmin
                and y < ymax
                and y > ymin
                and z < zmax
                and z > zmin
            )

        else:
            xt.append([x, y, t])

            # Interpolate for values of Ex and Ey at the specific position
            Ex = ExInter([x, y])
            Ey = EyInter([x, y])

            # Solve equation of motion
            xNext = x + vDrift * Ex[0] * np.sign(q) * dt
            yNext = y + vDrift * Ey[0] * np.sign(q) * dt

            # Assign the new version of x, y, and t
            x = xNext
            y = yNext

            inBound = x < xmax and x > xmin and y < ymax and y > ymin

        t = t + dt

    if np.isnan(x):
        x = xt[-1][0]
    if np.isnan(y):
        y = xt[-1][1]

    # Add value at the limit
    if threeDim:
        if x > xmax:
            xt.append([xmax, y, z, t])
        elif x < xmin:
            xt.append([xmin, y, z, t])
        elif y > ymax:
            xt.append([x, ymax, z, t])
        elif y < ymin:
            xt.append([x, ymin, z, t])
        elif z > zmax:
            xt.append([x, y, zmax, t])
        elif z < zmin:
            xt.append([x, y, zmin, t])
    else:
        if x > xmax:
            xt.append([xmax, y, t])
        elif x < xmin:
            xt.append([xmin, y, t])
        elif y > ymax:
            xt.append([x, ymax, t])
        elif y < ymin:
            xt.append([x, ymin, t])

    # convert to numpy array and return the data
    return np.array(xt)


def plotEField(field, cmap=None, type="mag", figH=None, norm="log"):
    """
	Plots the vector field or magnitude of the given electrid field

	Inputs:
		field - list in the form [x, y, Ex, Ey] just like the other functions. x and y here are conventions for the 
				x - Nx1 numpy array of x data points
				y - Mx1 numpy array of y data points
				Ex - (N*M)x1 numpy array corresponding to an Ex at every point on the grid
				Ey - (N*M)x1 numpy array corresponding to an Ey at every point on the grid
		cmap - matplotlib colormap
		type - "mag" for cont
	Outputs:
		ax - handle to the axis of the created figure
	"""

    # Define default colormap
    if not cmap:
        cmap = palettable.cmocean.sequential.Amp_20.mpl_colormap

    # Define normalization
    if norm == "log":
        norm = LogNorm()
    else:
        norm = None

    if figH:
        fig, ax = figH[0], figH[1]
    else:
        fig, ax = plt.subplots(1, 1, figsize=(14, 9))

    # Generate x-y pairs
    x = field[0]
    y = field[1]

    Ex = np.reshape(field[2], (x.size, y.size))
    Ey = np.reshape(field[3], (x.size, y.size))
    Ez = np.reshape(field[4], (x.size, y.size))

    # Plot results
    if type == "mag":

        # Compute graph value and axis
        Emag = np.sqrt(Ex ** 2 + Ey ** 2 + Ez ** 2)
        extent = [np.min(x), np.max(x), np.min(y), np.max(y)]
        im = ax.imshow(
            Emag.T,
            origin="lower",
            interpolation="bilinear",
            aspect="auto",
            extent=extent,
            cmap=cmap,
            norm=norm,
        )

        # Return the figure
        return fig, ax, im

    elif type == "vector":
        Print("Not yet implemented")
        return None
    else:
        Print("Type: %s not supported" % type)
        return None


def interpEField2D(x, y, E, method="linear"):
    """
	Wrapper function for interpolating points of the Efield over the Comsol grid

	Inputs:
		x - single value or a 1D array of x positions
		y - single value or a 1D of y positions
		E - list containing x, y grid data and E field data
	Outputs:
		EInterp - interpolated E field. Nx2 numpy array where N is the number of xy coordinate pairs
	"""

    # Create interpolating functions for the E fields
    Ex = scp.RegularGridInterpolator(
        (E[0], E[1]),
        np.reshape(E[2], (E[0].size, E[1].size), order="F"),
        bounds_error=False,
        method=method,
        fill_value=None,
    )
    Ey = scp.RegularGridInterpolator(
        (E[0], E[1]),
        np.reshape(E[3], (E[0].size, E[1].size), order="F"),
        bounds_error=False,
        method=method,
        fill_value=None,
    )

    # Create nd array of position
    pos = np.array([x, y]).T

    # Create output array of zeros
    EInterp = np.zeros(pos.shape)

    # Interpolate positions points
    EInterp[:, 0] = Ex(pos)
    EInterp[:, 1] = Ey(pos)

    return EInterp


def interpEField3D(x, y, z, E, method="linear"):
    """ Wrapper function for interpolating 3D field. Same function schema as interpEField2D but everything extended to 3D """

    # Create interpolating functions for the E fields
    Ex = scp.RegularGridInterpolator(
        (E[0], E[1], E[2]),
        np.reshape(E[3], (E[0].size, E[1].size, E[2].size), order="F"),
        bounds_error=False,
        method=method,
        fill_value=None,
    )
    Ey = scp.RegularGridInterpolator(
        (E[0], E[1], E[2]),
        np.reshape(E[4], (E[0].size, E[1].size, E[2].size), order="F"),
        bounds_error=False,
        method=method,
        fill_value=None,
    )
    Ez = scp.RegularGridInterpolator(
        (E[0], E[1], E[2]),
        np.reshape(E[5], (E[0].size, E[1].size, E[2].size), order="F"),
        bounds_error=False,
        method=method,
        fill_value=None,
    )

    # Create nd array of position
    pos = np.array([x, y, z]).T

    # Create output array of zeros
    EInterp = np.zeros(pos.shape)

    # Interpolate positions points
    EInterp[:, 0] = Ex(pos)
    EInterp[:, 1] = Ey(pos)
    EInterp[:, 2] = Ez(pos)

    return EInterp


def inducedChargeSingle(
    wPotential, path, q=1.6e-19, method="linear", roundFinalVal=False, stepback=0.0005
):

    qi = []

    if len(path) == 0:
        return np.zeros(1)

    if type(q) != np.ndarray:
        q = q * np.ones(path.shape[0])

    # Check if we are in 2 or 3D
    if len(wPotential) == 4:
        threeDim = True
    else:
        threeDim = False

    # Define potential interpolation function
    if threeDim:
        wPInter = scp.RegularGridInterpolator(
            (wPotential[0], wPotential[1], wPotential[2]),
            np.reshape(
                wPotential[3],
                (wPotential[0].size, wPotential[1].size, wPotential[2].size),
                order="F",
            ),
            bounds_error=False,
            method=method,
            fill_value=None,
        )
    else:
        wPInter = scp.RegularGridInterpolator(
            (wPotential[0], wPotential[1]),
            np.reshape(
                wPotential[2], (wPotential[0].size, wPotential[1].size), order="F"
            ),
            bounds_error=False,
            method=method,
            fill_value=None,
        )

    # Find charge induced at each point in the path
    for i in range(path.shape[0]):
        if threeDim:
            wP = wPInter([path[i, 0], path[i, 1], path[i, 2]])
        else:
            wP = wPInter([path[i, 0], path[i, 1]])

        if i == path.shape[0] - 1 and roundFinalVal:
            qi.append(-q[i] * round(wP[0]))
        else:
            qi.append(-q[i] * wP[0])

    # ensure the last value of qi is not nan
    if threeDim:
        zFinal = path[-1, 2]
        sign = np.sign(np.mean(np.diff(path[:, 2])))
    else:
        zFinal = path[-1, 1]
        sign = np.sign(np.mean(np.diff(path[:, 1])))

    count = 0
    maxCount = 10
    while np.isnan(qi[-1]):
        zFinal -= sign * stepback
        if count == maxCount:
            return np.zeros(1)

        # Recompute the last value with the step backed z value
        if threeDim:
            if roundFinalVal:
                qi[-1] = -q[-1] * round(wPInter([path[-1, 0], path[-1, 1], zFinal])[0])
            else:
                qi[-1] = -q[-1] * wPInter([path[-1, 0], path[-1, 1], zFinal])[0]
        else:
            if roundFinalVal:
                qi[-1] = -q[-1] * round(wPInter([path[-1, 0], zFinal])[0])
            else:
                qi[-1] = -q[-1] * wPInter([path[-1, 0], zFinal])[0]

        count += 1

    return np.array(qi)


def inducedCharge(
    wPotentialA,
    wPotentialB,
    path,
    q=-1.6e-19,
    method="linear",
    roundFinalVal=False,
    stepback=0.0005,
):
    """
	Finds the induced charge at each electrode given a path of the the charged particle

	Inputs:
		wPotentialA - list including [x, y, (z), Phi]. The x,y(z) position pairs and the potential occuring at this. For electrode A
		wPotentialB - list including [x, y, (z), Phi]. The x,y,(z) position pairs and the potential occuring at this. For electrode B
		path - the position to compute the weighted potential at. Nx2 (x,y position at different time steps) numpy array
		q - charge of the particle
		roundFinalVal - boolean value indicating whether or not this is an electron/hole moving through an object. If true, rounds the value of the last weighted potential
		stepback - double of how much to back off in the z direction if we encounter a NaN final value
	Outputs:
		qA - charge induced at electrode A
		qB - charge induced at electrode B
		qDiff - difference in the charge induced at electrode A and B
	"""
    qA = []
    qB = []

    if len(path) == 0:
        return np.zeros(1), np.zeros(1), np.zeros(1)

    # Check if we are in 2 or 3D
    if len(wPotentialA) == 4:
        threeDim = True
    else:
        threeDim = False

    # Definte interplation functions
    if threeDim:
        VaInter = scp.RegularGridInterpolator(
            (wPotentialA[0], wPotentialA[1], wPotentialA[2]),
            np.reshape(
                wPotentialA[3],
                (wPotentialA[0].size, wPotentialA[1].size, wPotentialA[2].size),
                order="F",
            ),
            bounds_error=False,
            method=method,
            fill_value=None,
        )
        VbInter = scp.RegularGridInterpolator(
            (wPotentialB[0], wPotentialB[1], wPotentialB[2]),
            np.reshape(
                wPotentialB[3],
                (wPotentialB[0].size, wPotentialB[1].size, wPotentialB[2].size),
                order="F",
            ),
            bounds_error=False,
            method=method,
            fill_value=None,
        )
    else:
        VaInter = scp.RegularGridInterpolator(
            (wPotentialA[0], wPotentialA[1]),
            np.reshape(
                wPotentialA[2], (wPotentialA[0].size, wPotentialA[1].size), order="F"
            ),
            bounds_error=False,
            method=method,
            fill_value=None,
        )
        VbInter = scp.RegularGridInterpolator(
            (wPotentialB[0], wPotentialB[1]),
            np.reshape(
                wPotentialB[2], (wPotentialB[0].size, wPotentialB[1].size), order="F"
            ),
            bounds_error=False,
            method=method,
            fill_value=None,
        )

    if type(q) != np.ndarray:
        q = q * np.ones(path.shape[0])

    # Iterated over all the positions in the path
    for i in range(path.shape[0]):
        if threeDim:
            Va = VaInter([path[i, 0], path[i, 1], path[i, 2]])
            Vb = VbInter([path[i, 0], path[i, 1], path[i, 2]])
        else:
            Va = VaInter([path[i, 0], path[i, 1]])
            Vb = VbInter([path[i, 0], path[i, 1]])

        # Find the q induced via the Shokley-Ramo Theorem
        if i == path.shape[0] - 1 and roundFinalVal:
            qA.append(-q[i] * round(Va[0]))
            qB.append(-q[i] * round(Vb[0]))
        else:
            qA.append(-q[i] * Va[0])
            qB.append(-q[i] * Vb[0])

    # ensure the last value of qi is not nan
    if threeDim:
        zFinal = path[-1, 2]
        sign = np.sign(np.mean(np.diff(path[:, 2])))
    else:
        zFinal = path[-1, 1]
        sign = np.sign(np.mean(np.diff(path[:, 1])))

    count = 0
    maxCount = 10
    while np.isnan(qA[-1]) or np.isnan(qB[-1]):
        zFinal -= sign * stepback
        if count == maxCount:
            return np.zeros(1), np.zeros(1), np.zeros(1)

        if threeDim:
            if roundFinalVal:
                qA[-1] = -q[-1] * round(VaInter([path[-1, 0], path[-1, 1], zFinal])[0])
                qB[-1] = -q[-1] * round(VbInter([path[-1, 0], path[-1, 1], zFinal])[0])
            else:
                qA[-1] = -q[-1] * VaInter([path[-1, 0], path[-1, 1], zFinal])[0]
                qB[-1] = -q[-1] * VbInter([path[-1, 0], path[-1, 1], zFinal])[0]
        else:
            if roundFinalVal:
                qA[-1] = -q[-1] * round(VaInter([path[-1, 0], zFinal])[0])
                qB[-1] = -q[-1] * round(VbInter([path[-1, 0], zFinal])[0])
            else:
                qA[-1] = -q[-1] * VaInter([path[-1, 0], zFinal])[0]
                qB[-1] = -q[-1] * VbInter([path[-1, 0], zFinal])[0]
        count += 1

    qA, qB = np.array(qA), np.array(qB)

    return qA, qB, qA - qB


def computeEffectiveArea(
    E, goodRange, xbounds, dx, ypos=0, vDrift=14.0e-6, dt=0.1, limits=[]
):
    """
    Computes the "effective area" of a multi electrode detector. Effective area is the ratio of number of initial charge particls
    that at arrive at the electrode of interest to the total number of charge particles simulated. 
    Inputs:
        E - list of [pos, Efield]. Either 4 or 6 elements depending on ndims in the detector
        wphiA - weighted potential of interested electrode. list of [pos weightedPotential]
        wphiB - weighted potentail at the other electrode. 
        xbounds - range in the x direction we scan over. [xmin, xmax]
        dx - x step size
        ypos - y position to set the inital charge deposition at
        vDrift - drift velocity in mm^2/ (V*us)
        dt - time in us
    Outputs
        effectiveArea - effective area of the detector
    """

    x = np.arange(xbounds[0], xbounds[1] + dx, dx)
    carrierAtElectrodeCount = 0

    for i in range(x.size):
        ri = [x[i], ypos]
        path = findMotion(ri, E, vDrift, dt, q=1, limits=limits)
        # _, _, wpDiff = inducedCharge(wphiA, wphiB, path, q=1)

        # if wpDiff[-1] < 0:
        if np.any(
            np.logical_and(path[-1, 0] > goodRange[:, 0], path[-1, 0] < goodRange[:, 1])
        ):
            carrierAtElectrodeCount += 1

    return float(carrierAtElectrodeCount) / x.size


if __name__ == "__main__":
    # (?# filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\test_weighted_potential.txt')
    # testHeader, testData = readComsolFile(filename)

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

    # initialPos = (1000, 100) # in um
    # vDriftHoles = 0.19e6 # cm^2/(V*s)

    # eField = takeSlice(testData, 0, 1000, [4,5])/(1e6) # Converting V/m to V/um
    # # xt = findMotion(initialPos, eField.T,  vDriftHoles, 0.001)

    # # For finding the current induced (vs the signal generated), we can more easily just set y to a fixed value
    # # and vary z across the range of depth. If we choose y to be
    # N = 500
    # z = np.linspace(21,219,N)
    # y = np.repeat(443,N)
    # pos = np.concatenate((y,z)).reshape(2,N).T

    # qA, qB, qDiff = inducedCharge(testData[:,(1,2,6)], testData[:,(1,2,9)], pos)

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(pos[:,1]-20,qA,'--', pos[:,1]-20,qB,'--', pos[:,1]-20, qDiff, linewidth=3)
    # ax.set_title(r'Induced Charge at a Coplanar Electrode. 100 $\mu m$ Spacing between fingers.', fontsize=16)
    # ax.set_xlabel(r'Depth ($\mu m$)', fontsize=14)
    # ax.set_ylabel(r'Induced Charge (C)', fontsize=14)
    # ax.legend(['qA','qB','qDiff'])

    # gridFile = r'C:\Users\alexp\Documents\UW\Research\Selenium\test_export.txt')
    # header, x, y, V = readComsolFileGrid(gridFile)
    # print(V.shape)
    # xx, yy = np.meshgrid(x,y)
    # z = np.reshape(V, (x.size, y.size))
    # fig2 = plt.figure()
    # ax2 = fig2.add_subplot(111)
    # ax2.contourf(x, y, z)
    # plt.show()

    # Test 3D weighted potentail
    filename = "/home/apiers/mnt/rocks/selena/data/em/pixel_detector_2mmdiam_200umAse_3d_small.txt"
    _, pos, field = readComsolFileGrid3d(filename)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 9))
    wPotential = [pos[0] * 1.0e6, pos[1] * 1.0e6, pos[2] * 1.0e6, field[:, -4]]
    x = np.linspace(-1800, 1800, 400)
    y = np.linspace(-1800, 1800, 400)
    z = np.linspace(-99, 99, 400)
    ax1.plot(
        x,
        inducedChargeSingle(
            wPotential, np.array([x, np.zeros(400), np.zeros(400)]).T, q=1
        ),
        linewidth=2,
    )
    ax2.plot(
        y,
        inducedChargeSingle(
            wPotential, np.array([np.zeros(400), y, np.zeros(400)]).T, q=1
        ),
        linewidth=2,
    )
    ax3.plot(
        z,
        inducedChargeSingle(
            wPotential, np.array([np.zeros(400), np.zeros(400), z]).T, q=1
        ),
        linewidth=2,
    )
    plt.show()

# Plot potential of for real detector
import numpy as np
import matplotlib.pyplot as plt
from plot import *
import os

filename = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\real_electrode.txt"
header, y, z, data = readComsolFileGrid(filename)
print (data.shape)
print (y.shape)
print (z.shape)

# Convert units from meters to mm
y, z = y * 1000, z * 1000

# Get the weighted Potentials
wPhiBot = data[:, 3]
wPhiTop = data[:, 6]

# Get real E fields
phi, Ey, Ez = data[:, 0], data[:, 1], data[:, 2]

fig, ax, cax = plotPhi((y, z), phi, 1)
ax.set_title(r"Potential from Experiment Detector", fontsize=16)
ax.set_xlabel(r"x (mm)", fontsize=14)
ax.set_ylabel(r"Depth (mm)", fontsize=14)
fig.colorbar(cax, ax=ax)

plt.show()

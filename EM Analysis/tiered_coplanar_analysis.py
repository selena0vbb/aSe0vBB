from plot import *
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.interpolate as scp

# Code to analyze the tiered coplanar electrode simultation

filename = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\tier_electrode\tiered_coplanar_electrode_15and5um_bias_voltage_sap.txt"
filename2 = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\bias_voltage\bias_voltage_small_scale_10umspace_0_360_range_unequal_electrode_size_fine.txt"

# Read in both files (to compare weighted potentials and whatnot from different)
header, yTier, zTier, dataTier = readComsolFileGrid(filename)
header, y, z, data = readComsolFileGrid(filename2)


# Replace NaNs with zeros
data = np.nan_to_num(data, copy=False)
dataTier = np.nan_to_num(dataTier, copy=False)

wPhiATier, wPhiBTier = dataTier[:, -6], dataTier[:, -3]
wPhiA, wPhiB = data[:, -6], data[:, -3]
phi, Ey, Ez = (
    data[:, 0 : np.amin(dataTier.shape) : 3] / 1e6,
    data[:, 1 : np.amin(dataTier.shape) : 3] / 1e6,
    data[:, 2 : np.amin(dataTier.shape) : 3],
)

yTier, zTier, y, z = yTier * 1e6, zTier * 1e6, y * 1e6, z * 1e6

# Creating contour plots of the weighted potential
subFig, (ax1, ax2) = plt.subplots(1, 2)
_, ax, cax = plotPhi((y, z), wPhiA, 1, figH=(subFig, ax1))
_, axT, caxT = plotPhi((yTier, zTier), wPhiATier, 1, figH=(subFig, ax2))

subFig.colorbar(cax, ax=ax1)
subFig.colorbar(caxT, ax=ax2)


# Axes and Titles
ax1.set_ylabel(r"Depth ($\mu m$)", fontsize=14)
ax1.set_xlabel(r"Width ($\mu m$)", fontsize=14)
ax1.set_title("Weighted Potential for Flat Coplanar", fontsize=14)
ax2.set_ylabel(r"Depth ($\mu m$)", fontsize=14)
ax2.set_xlabel(r"Width ($\mu m$)", fontsize=14)
ax2.set_title("Weighted Potential for Tiered Coplanar", fontsize=14)

# Plotting different comparisons between flat and tiered

# Comparing wPhiA, wPhiB, and wPhiDiff for both cases
# Tiered axis first
fig1, (axA1, axB1, axD1) = plt.subplots(1, 3)
_, _, caxA1 = plotPhi((y, z), wPhiATier, 1, figH=(fig1, axA1))
_, _, caxA2 = plotPhi((y, z), wPhiBTier, 1, figH=(fig1, axB1))
_, _, caxA3 = plotPhi((y, z), wPhiATier - wPhiBTier, 1, figH=(fig1, axD1))
fig1.colorbar(caxA1, ax=axA1)
fig1.colorbar(caxA2, ax=axB1)
fig1.colorbar(caxA3, ax=axD1)
axA1.set_ylabel(r"Depth ($\mu m$)", fontsize=14)
axA1.set_xlabel(r"Width ($\mu m$)", fontsize=14)
axB1.set_xlabel(r"Width ($\mu m$)", fontsize=14)
axD1.set_xlabel(r"Width ($\mu m$)", fontsize=14)
axA1.set_title("Weighted Potential Electrode A", fontsize=14)
axB1.set_title("Weighted Potential Electrode B", fontsize=14)
axD1.set_title("Weighted Potential Difference", fontsize=14)


# Regular coplanar electrode
fig2, (axA2, axB2, axD2) = plt.subplots(1, 3)
_, _, caxB1 = plotPhi((y, z), wPhiA, 1, figH=(fig2, axA2))
_, _, caxB2 = plotPhi((y, z), wPhiB, 1, figH=(fig2, axB2))
_, _, caxB3 = plotPhi((y, z), wPhiA - wPhiB, 1, figH=(fig2, axD2))
fig2.colorbar(caxA1, ax=axA2)
fig2.colorbar(caxA2, ax=axB2)
fig2.colorbar(caxA3, ax=axD2)
axA2.set_ylabel(r"Depth ($\mu m$)", fontsize=14)
axA2.set_xlabel(r"Width ($\mu m$)", fontsize=14)
axB2.set_xlabel(r"Width ($\mu m$)", fontsize=14)
axD2.set_xlabel(r"Width ($\mu m$)", fontsize=14)
axA2.set_title("Weighted Potential Electrode A", fontsize=14)
axB2.set_title("Weighted Potential Electrode B", fontsize=14)
axD2.set_title("Weighted Potential Difference", fontsize=14)

# Comparison of the difference in weighted potential between different types of electrodes
fig3, (axD3, axD4) = plt.subplots(1, 2)
_, _, caxC1 = plotPhi((y, z), wPhiA - wPhiB, 1, figH=(fig3, axD3))
_, _, caxC2 = plotPhi((y, z), wPhiATier - wPhiBTier, 1, figH=(fig3, axD4))
fig2.colorbar(caxC1, ax=axD3)
fig2.colorbar(caxC2, ax=axD4)
axD3.set_ylabel(r"Depth ($\mu m$)", fontsize=14)
axD3.set_xlabel(r"Width ($\mu m$)", fontsize=14)
axD4.set_xlabel(r"Width ($\mu m$)", fontsize=14)
axD3.set_title(r"$\Phi_{Wa}-\Phi_{Wb}$ Flat Electrode", fontsize=14)
axD4.set_title(r"$\Phi_{Wa}-\Phi_{Wb}$ Tiered Electrode", fontsize=14)


# Plot of the WPhi Difference for the flat electrode. Requested by ALvaro in May/31/18 email
fig4, ax4 = plt.subplots(1, 1)
_, _, cax4 = plotPhi((y, z), wPhiA - wPhiB, 1, figH=(fig4, ax4))
fig4.colorbar(cax4, ax=ax4)
ax4.set_ylabel(r"Depth ($\mu m$)", fontsize=14)
ax4.set_xlabel(r"Width ($\mu m$)", fontsize=14)
ax4.set_title(r"$\Phi_{Wa}-\Phi_{Wb}$ Flat Electrode", fontsize=14)

plt.show()

# Computing the max Efield near the collection electrodes in the different case

N = 200
zi = 1  # z position, near the electrode
yi = np.arange(1, 195, N)

Emax = []
for i in range(np.amin(Ey.shape)):
    # Create interpolation functions for E

    EyInter = scp.interp2d(
        yTier, zTier, np.reshape(Ey[:, i], (zTier.size, yTier.size)), kind="linear"
    )
    EzInter = scp.interp2d(
        yTier, zTier, np.reshape(Ez[:, i], (zTier.size, yTier.size)), kind="linear"
    )

    Enorm = []
    for j in range(yi.size):

        Enorm.append(
            np.sqrt(EyInter(yi[j], zi) ** 2 + EzInter(yi[j], zi) ** 2)
        )  # Keep track of the magnitude of the efield at different points

    Emax.append(np.amax(np.array(Enorm)))

print(Emax)

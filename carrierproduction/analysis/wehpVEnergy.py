# Script to calculate best fit for the W_ehp(E) function
import numpy as np
import scipy.optimize as scp
import matplotlib.pyplot as plt

field = np.arange(10, 24, 2)
wehp = np.array([62.5, 56, 51, 47.5, 43.5, 40.5, 37.8])/1000

def we(field, *params):
	""" functional form of wehp. w(E) = param[0] + param[1]/E """
	return params[0] + params[1]/field

def fitWehp():
	""" Fit the wehp(E) to an 1/F function """

	# Fit the data
	initialParams = [6, 400]
	popt, pcov = scp.curve_fit(we, field, wehp, initialParams)
	print(popt)
	eFieldAxis = np.linspace(10, 25, 100)

	# plot the result
	fig, ax = plt.subplots(1,1, figsize=(12,8))
	ax.scatter(field, wehp, marker="X", label="Data")
	ax.plot(eFieldAxis, we(eFieldAxis, *popt), "r", linewidth=2, label="Fit: %.3f + %.2f/E"%(popt[0], popt[1]))
	ax.set_xlabel(r"$\vec{E}$ [$\rm V\ \mu m^{-1}$]", fontsize=18)
	ax.set_ylabel(r"$\rm W_{ehp} \ [keV]$", fontsize=18)
	ax.legend(fontsize=16)
	# ax.tick_params()
	plt.show()

if __name__ == '__main__':
	fitWehp()
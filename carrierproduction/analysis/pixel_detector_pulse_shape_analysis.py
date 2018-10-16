import numpy as np
import matplotlib.pyplot as plt
import particledata as pd
import brewer2mpl
import os
from scipy.signal import butter, lfilter, freqz

def plotPulsesRaw(datafile, nPerGraph, nGraph, titlestring='', show=True):
    """
    plots a series of raw pulses vs time. Plots nPerGraph pulses on the same graph, and creates nGraphs
    """

    data = np.load(datafile)
    bmap = brewer2mpl.get_map('Paired','Qualitative',min([nPerGraph, 12])).mpl_colors
    fig = []
    ax = []
    totalcount = 0
    savecount = 0
    w = 0.05
    q = 1.6e-19
    outputdir = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\charge collection\pixel_detector'


    for i in range(nGraph):
        figt, axt = plt.subplots(figsize=(16,9))
        for j in range(nPerGraph):
            t = data[totalcount][0]
            sig = data[totalcount][1]
            axt.plot(t, -sig*w/q, linewidth=2, color=bmap[j])
            totalcount += 1

        axt.set_xlabel('Time ($\mu s$)', fontsize=16)
        axt.set_ylabel('Charge Induced (keV)', fontsize=16)
        axt.set_title('Pulse Shape (Inverted) for %i Different Events. %s' % (nPerGraph, titlestring), fontsize=18)
        fileloc = os.path.join(outputdir, 'pixel_detector_pulses%i.png'%savecount)
        # figt.savefig(fileloc, bbox_inches='tight')
        # figt.savefig(fileloc)
        savecount += 1

        fig.append(figt)
        ax.append(axt)

    if show:
        plt.show()

    return None

def plotPulsesFiltered(datafile, nPerGraph, nGraph, show=True):

    data = np.load(datafile)
    bmap = brewer2mpl.get_map('Paired','Qualitative',min([nPerGraph, 12])).mpl_colors
    fig = []
    ax = []
    totalcount = 0
    savecount = 0
    w = 0.05
    q = 1.6e-19
    outputdir = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\charge collection\pixel_detector'

    # Define filter parameters
    flow, fhigh = 5.30e3, 0.53e6


    for i in range(nGraph):
        figt, axt = plt.subplots(figsize=(16,9))
        for j in range(nPerGraph):
            t = data[totalcount][0]
            dt = np.diff(t)[0]*1.e-6
            sig = -data[totalcount][1]
            axt.plot(t, sig*w/q, linewidth=2, color=bmap[j])
            axt.plot(t, butterBandpassFilter(sig, flow, fhigh, dt, order=1)*w/q, linewidth=2, color=bmap[j], linestyle=':')
            totalcount += 1

        axt.set_xlabel('Time ($\mu s$)', fontsize=16)
        axt.set_ylabel('Charge Induced (keV)', fontsize=16)
        axt.set_title('Pulse Shape (Inverted) and Filtered Signal for %i Different Events.' % (nPerGraph), fontsize=18)
        fileloc = os.path.join(outputdir, 'pixel_detector_pulses%ifilter.png'%savecount)
        figt.savefig(fileloc, bbox_inches='tight')
        # figt.savefig(fileloc)
        savecount += 1

        fig.append(figt)
        ax.append(axt)

    if show:
        plt.show()

    return None

def energyHistogram():

    basedir = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector'

    filenames = ['pixel_particle_event.npy', 'pixel_particle_event_poisson.npy', 'pixel_particle_event_poisson_trap.npy', 'pixel_particle_event_trap.npy']
    legstring = ['Pixel Particle Event', 'Poisson Fluctuation', 'Poisson + Trapping', 'Trapping']
    # filenames = ['pixel_particle_event.npy']

    fig, ax = plt.subplots()

    bmap = brewer2mpl.get_map('Paired','Qualitative',2*max([len(filenames), 3])).mpl_colors
    w = 0.05
    q = 1.6e-19
    flow, fhigh = 5.30e3, 0.53e6
    dt = 1.e-8
    bins=np.arange(0,150)
    for j, file in enumerate(filenames):
        print(file)
        finalEnergy = []
        finalEnergyFilter = []
        data = np.load(os.path.join(basedir, file))
        for i in range(data.shape[0]):
            finalEnergy.append(-np.min(data[i][1])*w/q)
            finalEnergyFilter.append(np.max(butterBandpassFilter(-data[i][1], flow, fhigh, dt, order=1)*w/q))
        ax.hist(finalEnergy, bins=bins, linewidth=2, histtype='step', label=legstring[j] + ' $\mu$=%.2f, $\sigma$=%.2f keV'%(np.mean(finalEnergy), np.std(finalEnergy)), color=bmap[2*j])
        # ax.hist(finalEnergyFilter, bins=bins, linewidth=2, histtype='step', label=legstring[j]+' Filter. $\mu$=%.2f, $\sigma$=%.2f keV'%(np.mean(finalEnergyFilter), np.std(finalEnergyFilter)), color=bmap[2*j+1], linestyle=':')

    ax.set_yscale('log')
    ax.legend(loc='upper left')
    ax.set_xlabel('Energy (keV)', fontsize=16)
    ax.set_ylabel('Number of Events', fontsize=16)
    ax.set_title('Energy Resolution Histogram for Raw and Filtered Signal', fontsize=18)

    plt.show()

    return None

def filterParam(f, rHigh, cHigh, rLow, cLow):
    """
    Creates the frequency response of the filter parameters
    """
    tauHigh = rHigh*cHigh
    tauLow = rLow*cLow
    N = 1e5

    # f = np.linspace(-50e6, 50e6, int(N))

    Hjw = (2*np.pi*abs(f)*tauHigh*1j) / ((1+2*np.pi*abs(f)*tauHigh*1j ) * (1+2*np.pi*abs(f)*tauLow*1j))

    return Hjw

def testScpFilter():

    f = np.logspace(2, 7.7, 1000)
    Hjw =filterParam(f, 1.e6, 30.e-12, 1.e4, 30.e-12)

    flow, fhigh = 5.30e3, 0.53e6
    dt = 10.e-9

    bb, aa  = butterBandpass(flow, fhigh, dt, order=1)
    ww, hh = freqz(bb, aa , worN=2**15)

    fig, ax = plt.subplots()
    ax.plot(0.5/(np.pi*dt)*ww, abs(hh), label='Butterworth', linewidth=2, color='k')
    ax.plot(f, abs(Hjw), label='Analytic Filter', linewidth=2, color='blue')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=16)
    ax.set_ylabel('Gain', fontsize=16)
    ax.set_title('Frequency Response of Different Filter Implementations', fontsize=18)
    ax.legend()

    data = np.concatenate((np.zeros(100), np.ones(10000)))

    y = butterBandpassFilter(data, flow, fhigh, dt, order=1)

    time = np.arange(0, dt*y.size, dt)

    fig2, ax2 = plt.subplots()
    ax2.plot(time*1e6, y, 'b', linewidth=2, label='Unit Step Response')
    ax2.plot(time*1e6, data, 'k', linewidth=2, label='Filtered Response')
    ax2.set_xlabel('Time ($\mu s$)', fontsize=16)
    ax2.set_title('Operation of Butterworth Filter', fontsize=18)
    ax2.legend()

    plt.show()

def butterBandpass(flow, fhigh, dt, order=5):
    nyqf = 0.5 / dt
    low = flow / nyqf
    high = fhigh / nyqf
    b, a  = butter(order, [low, high], btype='band')
    return b, a

def butterBandpassFilter(data, flow, fhigh, dt, order=5):
    b, a = butterBandpass(flow, fhigh, dt, order=order)
    y = lfilter(b, a, data)
    return y

def spectrumWidthFilter(filename):

    data = np.load(filename)
    flow = np.logspace(1, 4, 30)
    fhigh =[0.1e6, 0.53e6, 1.e6]
    w=0.05
    q=1.6e-19

    sigma = []

    for i, high in enumerate(fhigh):
        sigmaIntermidiate = []
        for j, low in enumerate(flow):
            datamax = []
            for k in range(data.shape[0]):

                datamax.append(np.max(butterBandpassFilter(-data[k][1]*w/q, low, high, 1.e-8, order=1)))

            sigmaIntermidiate.append(np.std(datamax))

        sigma.append(sigmaIntermidiate)

    fig, ax = plt.subplots(figsize=(16,9))
    sigma = np.array(sigma)
    ax.plot(flow, sigma.T, linewidth=2)
    ax.plot([5.3e3, 5.3e3], [0, np.max(sigma)], linewidth=2, color='black', linestyle='--', label='Current Lowpass Frequency')
    legstring = ['$f_{high}$ = %.2f MHz' % (x/1e6) for x in fhigh]
    legstring.append('Current Lowpass Frequency')
    ax.legend(legstring)
    ax.set_xlabel('Lowpass Filter Frequency (Hz)', fontsize=16)
    ax.set_xscale('log')
    ax.set_ylabel('Sigma of Energy Distribution (keV)', fontsize=16)
    ax.set_title('Spread of Energy Spectrum Due to Lowpass Filter Frequency', fontsize=18)


if __name__ == '__main__':

    filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\pixel_particle_event.npy'

    # plotPulsesRaw(filename, 12, 3, show=False)
    energyHistogram()
    # testFreqResponse()
    # filterParam(1.e6, 30.e-12, 1.e4, 30.e-12)
    # testScpFilter()
    # plotPulsesFiltered(filename, 8, 4, show=False)
    # energyHistogram()
    # spectrumWidthFilter(filename)
    plt.show()
    # checkData(filename)

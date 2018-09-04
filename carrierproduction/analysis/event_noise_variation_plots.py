# Noise variations

import numpy as np
import matplotlib.pyplot as plt
import brewer2mpl

def make_histograms():

    filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\charge collection\noise_histograms\data\planar_particle_event_trap.npy'


    data = np.load(filename)
    chartTitle = 'Planar Detector, Trapping, N=%i 122 keV Geant4 Events'%len(data)
    fig, ax = plt.subplots()
    ax.set_title(chartTitle, fontsize=16)
    ax.set_xlabel('Charge Induced in terms of Energy (keV)', fontsize=14)
    nbin = 140
    histrange = (0, 140)
    times = [-1, 0.75, 1.5]
    histdata = []

    e = 1.6e-19
    wehp = 0.05
    bmap = brewer2mpl.get_map('Dark2', 'Qualitative', max(len(times), 3)).mpl_colors

    for i, time in enumerate(times):
        chargeenergy = []
        for event in data:


            if time == -1:
                tindx = -1
                timelabel = 'final time'
                labelstring = '%s  Mean = %0.2f  rms = %0.2f'%(timelabel, np.mean(chargeenergy), np.std(chargeenergy))
            else:
                try:
                    tindx = np.where(event[0] > time)[0][0]
                    timelabel = str(time)
                except(IndexError):
                    tindx = -1
                    timelabel = str(time)
                labelstring = '%s $\mu s$  Mean = %0.2f  rms = %0.2f'%(timelabel, np.mean(chargeenergy), np.std(chargeenergy))

            chargeenergy.append(abs(event[1][tindx]/1.6e-19*wehp))
        ax.hist(chargeenergy, bins=nbin, range=histrange, color=bmap[i], linewidth=3, histtype='step', label=labelstring)
        histdata.append(chargeenergy)

    # create the histogram

    ax.legend()


    return fig, ax

def width_of_signal():
    """Looks at the delta t of the rise in the signal"""

    filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\charge collection\noise_histograms\data\coplanar_scaling_particle_event.npy'

    data = np.load(filename)

    tThresh = 2
    etot = 122
    e = 1.6e-19
    wehp = 0.05
    signalLength =  []
    index = []
    for i, event in enumerate(data):
        t1 = np.where(-event[1]/e * wehp > tThresh)[0][0]
        t2 = np.where(-event[1]/e * wehp > etot - tThresh)[0][0]
        signalLength.append(event[0][t2]-event[0][t1])
        if signalLength[i] > 2.5:
            index.append(i)

    fig, ax = plt.subplots()
    ax.set_xlabel(r'Time ($\mu s$)', fontsize=16)
    ax.set_title('Distribution of Rise Time of Charge Signal', fontsize=18)
    ax.hist(signalLength, bins=100, histtype='step', linewidth=2, label='Mean=%0.2f  rms=%0.2f'%(np.mean(signalLength), np.std(signalLength)))
    ax.legend()
    print(len(index))
    return None

if __name__ == '__main__':
    # fig, ax = make_histograms()
    width_of_signal()
    plt.show()

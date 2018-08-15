import particledata as pd
import seleniumconfig as sc
import matplotlib.pyplot as plt
import brewer2mpl
import plot
import numpy as np

def geometric_carrier_trapping(eventIDs):
    filename = r"C:\Users\alexp\Documents\UW\Research\Selenium\aSe0vBB\particle\selenium-build\output\122_keV_testTuple.root"
    emfilename = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_5um_spacing_fullsize.txt"
    configfilename = r"./config.txt"

    settings = sc.readConfigFile(configfilename)

    newCollection = pd.gEventCollection(filename)
    newCollection.printInfo()
    bmap = brewer2mpl.get_map('Set2', 'Qualitative', max(len(eventIDs)+1,3)).mpl_colors

    fig, ax = plt.subplots()

    for colorj, j in enumerate(eventIDs):
        settings['CARRIER_LIFETIME_GEOMETRIC'] = 1
        event = newCollection.collection[j]

        # Iterate over many scenarios to watch fluctuations
        N = 30
        qtot = 0

        for i in range(N):
            print(i)

            t, q = pd.computeChargeSignal(event, emfilename, **settings)
            ax.plot(t, -q/1.6e-19*0.05, color=bmap[colorj], alpha=0.3)
            qtot += q

        ax.plot(t, -qtot/N/1.6e-19*0.05, color=bmap[colorj], linewidth=2, label='Noise, Event ID %i' % event.GetEventID())
        ax.set_xlabel(r'Time ($\mu$s)', fontsize=14)
        ax.set_ylabel(r'Induced Charge (keV)', fontsize=14)
        ax.set_title('Q(t), expressed as energy, on a Coplanar Detector (5 $\mu m$ thickness kapton, 280 V bias) \n with Carrier Trapping Fluctuation. Scaling' , fontsize=16)


        settings['CARRIER_LIFETIME_GEOMETRIC'] = 0
        t, q = pd.computeChargeSignal(event, emfilename, **settings)
        ax.plot(t, -q/1.6e-19*0.05, color='black', linewidth=1, label='No Noise')
    ax.legend()
    return fig, ax

def noise_histogram():
    filename = r"C:\Users\alexp\Documents\UW\Research\Selenium\aSe0vBB\particle\selenium-build\output\122_keV_testTuple.root"
    emfilename = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_5um_spacing_fullsize.txt"
    configfilename = r"./config.txt"

    settings = sc.readConfigFile(configfilename)

    # For planar detector
    settings['CHARGE_DIFFERENCE'] = 1
    settings['VOLTAGE_SWEEP_INDEX'] = 4
    settings['SCALE_WEIGHTED_PHI'] = 1

    # color scheme
    bmap = brewer2mpl.get_map('Dark2', 'Qualitative',5).mpl_colors

    # Read event data
    collection = pd.gEventCollection(filename)
    event = collection.collection[125]

    N = 50
    wehp = 0.05 # keV
    e = 1.6e-19

    energy = []
    charge = []
    charge1usNoTrapping = []
    charge20usNoTrapping = []
    charge1usTrapping = []
    charge20usTrapping = []
    allfig = []
    allax = []
    flatData = event.flattenEvent()


    for i in range(N):

        print(i)

        energy.append(sum(flatData['energy']))

        # nehp, nehpFluctuations, _, _ = event.createCarriers(**settings)
        # totalCharge = np.sum((nehp + nehpFluctuations) * wehp)
        # charge.append(totalCharge)


        # No Trapping
        settings['CARRIER_LIFETIME_GEOMETRIC'] = 0
        t, qNo = pd.computeChargeSignal(event, emfilename, **settings)
        indx1us = np.where(t > 1)[0][0]
        indx20us = np.where(t > 20)[0][0]
        charge1usNoTrapping.append(-qNo[indx1us]/e * wehp)
        charge20usNoTrapping.append(-qNo[indx20us]/e * wehp)

        # Trapping
        settings['CARRIER_LIFETIME_GEOMETRIC'] = 1
        t, qTrap = pd.computeChargeSignal(event, emfilename, **settings)
        indx1us = np.where(t > 0.5)[0][0]
        indx20us = np.where(t > 22)[0][0]
        charge1usTrapping.append(-qTrap[indx1us]/e * wehp)
        charge20usTrapping.append(-qTrap[indx20us]/e * wehp)

    fig, ax = plt.subplots()
    ax.plot(t, -qTrap)
    # Energy histogram
    fig, ax = plt.subplots()
    ax.hist(energy, bins=100, range=(100,140),  histtype='step', linewidth=2, color=bmap[0])
    ax.set_title('Energy Deposited in Event %i' % event.GetEventID(), fontsize=16)
    ax.set_xlabel('Energy (keV)', fontsize=14)
    allfig.append(fig), allax.append(ax)

    # Charge creation histogram
    # fig, ax = plt.subplots()
    # ax.hist(charge, bins=100, range=(100,140),  histtype='step', linewidth=2, color=bmap[0])
    # ax.set_title('Charge Created with Poisson Fluctuations in Event %i' % event.GetEventID(), fontsize=16)
    # ax.set_xlabel('Energy (keV)', fontsize=14)

    # Induced charge, trapping
    fig, ax = plt.subplots()
    ax.hist(charge1usTrapping, bins=122, range=(0,122), histtype='step', linewidth=2, color=bmap[0], label='1 $\mu s$')
    ax.hist(charge20usTrapping, bins=122, range=(0,122), histtype='step', linewidth=2, color=bmap[1], label='20 $\mu s$')
    ax.set_title('Charge Induced at Planar Detector (Poisson and Trapping Fluctuations) in Event %i' % event.GetEventID(), fontsize=16)
    ax.set_xlabel('Energy (keV)', fontsize=14)
    ax.legend()
    allfig.append(fig), allax.append(ax)

    # Induced charge, no trapping
    fig, ax = plt.subplots()
    ax.hist(charge1usNoTrapping, bins=122, range=(0,122), histtype='step', linewidth=2, color=bmap[0], label='1 $\mu s$')
    ax.hist(charge20usNoTrapping, bins=122, range=(0,122), histtype='step', linewidth=2, color=bmap[1], label='20 $\mu s$')
    ax.set_title('Charge Induced at Planar Detector (Poisson and No Trapping Fluctuations) in Event %i' % event.GetEventID(), fontsize=16)
    ax.set_xlabel('Energy (keV)', fontsize=14)
    ax.legend()
    allfig.append(fig), allax.append(ax)

    allfig.append(fig), allax.append(ax)

    return allfig, allax

if __name__ == '__main__':

    # functions that perform different analysis. Comment out to remove from computation

    # fig, ax = geometric_carrier_trapping([125, 126, 129, 143])
    # fig, ax = geometric_carrier_trapping([126])

    noise_histogram()
    plt.show()

    # filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\kapton_layer_analysis_5um_spacing_fullsize.txt'
    # _, x, y, data = plot.readComsolFileGrid(filename)
    # x, y = x*1e6, y*1e6

    # phi = data[:,0:np.amin(data.shape):3]
    # wPhiA = np.reshape(phi[:,-2], (y.size, x.size))
    # wPhiB = np.reshape(phi[:,-1], (y.size, x.size))

    # # scaleFactor = pd.scaleWeightPhi(wPhiA, wPhiB)

    # Ex, Ey = data[:,10], data[:,11]
    # Etot = [x, y, Ex/1e6, Ey/1e6]

    # xi = np.arange( -200, 200, 2)
    # yi = 90
    # limits = [-5000, 5000, -100, 100]

    # fig, ax = plt.subplots()

    # for x in xi:
    #     print(x)
    #     path = plot.findMotion((x, yi), Etot, 19, 0.03, limits=limits, q=1)
    #     ax.plot(path[:,0], path[:,1], 'k')

    # plt.show()





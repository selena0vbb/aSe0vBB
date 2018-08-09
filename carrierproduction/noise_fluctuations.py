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
        N = 25
        qtot = 0

        for i in range(N):
            print(i)

            t, q = pd.computeChargeSignal(event, emfilename, **settings)
            ax.plot(t, -q/1.6e-19*0.05, color=bmap[colorj], alpha=0.3)
            qtot += q

        ax.plot(t, -qtot/N/1.6e-19*0.05, color=bmap[colorj], linewidth=2, label='Noise, Event ID %i' % event.GetEventID())
        ax.set_xlabel(r'Time ($\mu$s)', fontsize=14)
        ax.set_ylabel(r'Induced Charge (keV)', fontsize=14)
        ax.set_title(r'Q(t), expressed as energy, on a Coplanar Detector (5 $\mu m$ thickness kapton, 280 V bias) with Carrier Trapping Fluctuation. No Scaling' , fontsize=16)


        settings['CARRIER_LIFETIME_GEOMETRIC'] = 0
        t, q = pd.computeChargeSignal(event, emfilename, **settings)
        ax.plot(t, -q/1.6e-19*0.05, color='black', linewidth=1, label='No Noise')
    ax.legend()
    return fig, ax

if __name__ == '__main__':

    # functions that perform different analysis. Comment out to remove from computation

    fig, ax = geometric_carrier_trapping([125, 126, 129, 143])
    # fig, ax = geometric_carrier_trapping([126])
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





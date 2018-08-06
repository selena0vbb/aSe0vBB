import particledata as pd
import seleniumconfig as sc
import matplotlib.pyplot as plt
import brewer2mpl

def geometric_carrier_trapping():

    filename = r"C:\Users\alexp\Documents\UW\Research\Selenium\aSe0vBB\particle\selenium-build\output\122_keV_testTuple.root"
    emfilename = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\real_electrode.txt"
    configfilename = r"./config.txt"

    settings = sc.readConfigFile(configfilename)
    settings['CARRIER_LIFETIME_GEOMETRIC'] = 1

    newCollection = pd.gEventCollection(filename)
    newCollection.printInfo()

    event = newCollection.collection[124]

    # Iterate over many scenarios to watch fluctuations
    N = 60
    qtot = 0
    fig, ax = plt.subplots()
    for i in range(N):
        print(i)

        t, q = pd.computeChargeSignal(event, emfilename, **settings)
        ax.plot(t,q, color='grey', alpha=0.5)
        qtot += q

    ax.plot(t, qtot/N, color='black', linewidth=2, label='Charge Trap Noise')
    ax.set_xlabel(r'Time ($\mu$s)', fontsize=14)
    ax.set_ylabel(r'Induced Charge (C)', fontsize=14)
    ax.set_title('Q(t) at a Unipolar electrode with Carrier Trapping Fluctuation for event %i' % (event.GetEventID()) , fontsize=16)


    settings['CARRIER_LIFETIME_GEOMETRIC'] = 0
    t, q = pd.computeChargeSignal(event, emfilename, **settings)
    ax.plot(t, q, color='blue', linewidth=2, label='No Noise')
    ax.legend()
    return fig, ax

if __name__ == '__main__':

    # functions that perform different analysis. Comment out to remove from computation

    fig, ax = geometric_carrier_trapping()
    plt.show()




import numpy as np
import matplotlib.pyplot as plt
import brewer2mpl
import particledata as pd


def test_filter_algorithm():
    filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\charge collection\noise_histograms\data\ccoplanar__scaling_particle_event.npy'

    filterfilename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Circuits\Draft22.txt'

    data = np.load(filename)

    t = data[0][0] * 1e-6
    sig = data[0][1]

    signalFilt, signalFiltFreq, xf, yf = pd.filterSignal(sig, t, filterfilename )

    return None




if __name__ == '__main__':

    test_filter_algorithm()

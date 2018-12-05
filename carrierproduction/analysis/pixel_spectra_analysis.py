# Analysis of pixel detector energy spectra

import numpy as np
import pulseAnalysis as pa
import matplotlib.pyplot as plt
import brewer2mpl

intensity122 = 0.856
intensity136 = 0.1068


def compare_partial_particle_tracks():

    # Partial particle vs standard for 122 keV
    outdirNoPartial = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector'
    outdirPartial = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\simple'
    partialPattern = ['pixel_122kev\d+.npy', 'pixel_122kev_400k_pt\d+.npy']
    noPartialPattern = ['pixel_particle_event_all.npy', 'pixel_particle_event_all_poisson_trap.npy']

    ePartial = pa.getEnergySpectrum(pa.readBatchFiles(outdirPartial, partialPattern[0]))
    eNoPartial = pa.getEnergySpectrum(pa.readBatchFiles(outdirNoPartial, noPartialPattern[0]))

    ePartialPT = pa.getEnergySpectrum(pa.readBatchFiles(outdirPartial, partialPattern[1]))
    eNoPartialPT = pa.getEnergySpectrum(pa.readBatchFiles(outdirNoPartial, noPartialPattern[1]))

    bmap = brewer2mpl.get_map("Set1", "Qualitative", 3).mpl_colors
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    bins = np.arange(0, 140)
    ax1.hist(ePartial, bins=bins, histtype='step', linewidth=3, color=bmap[0], label='Partial Tracks Included', density=True)
    ax1.hist(eNoPartial, bins=bins, histtype='step', linewidth=3, color=bmap[1], label='Partial Tracks Excluded', density=True)
    ax1.set_yscale('log')
    ax1.legend(loc='upper left', fontsize=12)
    ax1.set_title('Raw 122 keV Spectrum', fontsize=14)

    ax2.hist(ePartialPT, bins=bins, histtype='step', linewidth=3, color=bmap[0], label='Partial Tracks Included', density=True)
    ax2.hist(eNoPartialPT, bins=bins, histtype='step', linewidth=3, color=bmap[1], label='Partial Tracks Excluded', density=True)
    ax2.set_yscale('log')
    ax2.legend(loc='upper left', fontsize=12)
    ax2.set_xlabel('Energy (keV)', fontsize=14)
    ax2.set_title('Poisson and Trapping 122 keV Spectrum', fontsize=14)

    fig.suptitle('Comparison of Spectrum when including/excluding Partial Tracks', fontsize=18)


def higher_statistics_spectra():

    # 122 kev line
    outdir =r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\simple'
    pattern = ['pixel_122kev\d+.npy', 'pixel_122kev_400k_pt\d+.npy']

    e122 = pa.getEnergySpectrum(pa.readBatchFiles(outdir, pattern[0]))
    e122pt = pa.getEnergySpectrum(pa.readBatchFiles(outdir, pattern[1]))

    bmap = brewer2mpl.get_map("Set1", "Qualitative", 3).mpl_colors

    fig, ax = plt.subplots()
    bins = np.arange(0, 140)
    ax.hist(e122, bins=bins, histtype='step', linewidth=3, color=bmap[0], label='Raw')
    ax.hist(e122pt, bins=bins, histtype='step', linewidth=3, color=bmap[1], label='Poisson + Trapping')
    ax.set_yscale('log')
    ax.legend(loc='upper left', fontsize=16)
    ax.set_xlabel('Energy (keV)', fontsize=16)
    ax.set_title('Higher Statistics (N=%i) Energy Spectra 122 keV'%len(e122), fontsize=18)

    # Co57 spectrum
    outdir136 = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\136keV\simple'
    pattern136 = [r'pixel_136kev_400k\d+.npy', r'pixel_136kev_400k_pt\d+.npy']

    e136 = pa.getEnergySpectrum(pa.readBatchFiles(outdir136, pattern136[0]))
    e136pt = pa.getEnergySpectrum(pa.readBatchFiles(outdir136, pattern136[1]))

    coRaw, bins = pa.histogramCo57(e122, e136)
    coPT, bins = pa.histogramCo57(e122pt, e136pt)

    fig, ax = plt.subplots()
    ax.hist(bins[:-1], bins=bins, weights=coRaw, histtype='step', linewidth=3, color=bmap[0], label='Raw')
    ax.hist(bins[:-1], bins=bins, weights=coPT, histtype='step', linewidth=3, color=bmap[1], label='Poisson + Trapping')
    ax.set_yscale('log')
    ax.legend(loc='upper left', fontsize=16)
    ax.set_xlabel('Energy (keV)', fontsize=16)
    ax.set_title('Higher Statistics Normalized Energy Spectra Co-57', fontsize=18)


def new_pixel_particle_spectra():

    # 122 kev line
    outdir = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\sio2'
    pattern = ['pixel_122kev_sio2_5M\d+.npy','pixel_122kev_sio2_5M_pt\d+.npy']

    e122pt = pa.getEnergySpectrum(pa.readBatchFiles(outdir, pattern[1]))
    e122ptSimple = pa.getEnergySpectrum(pa.readBatchFiles(r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\simple', 'pixel_122kev_400k_pt\d+.npy'))

    # Co57 spectrum
    outdir136 = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\136keV\sio2'
    pattern136 = ['pixel_136kev_sio2_5M\d+.npy', 'pixel_136kev_sio2_5M_pt\d+.npy']

    e136pt = pa.getEnergySpectrum(pa.readBatchFiles(outdir136, pattern136[1]))
    e136ptSimple = pa.getEnergySpectrum(pa.readBatchFiles(r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\136keV\simple', 'pixel_136kev_400k_pt\d+.npy'))

    # 136 lines
    bmap = brewer2mpl.get_map("Set1", "Qualitative", 3).mpl_colors

    fig, ax = plt.subplots()
    bins = np.arange(0, 140)
    ax.hist(e136pt, bins=bins, histtype='step', linewidth=3, color=bmap[0], label='SiO2', density=True)
    ax.hist(e136ptSimple, bins=bins, histtype='step', linewidth=3, color=bmap[1], label='Simple', density=True)
    ax.set_yscale('log')
    ax.legend(loc='upper left', fontsize=16)
    ax.set_xlabel('Energy (keV)', fontsize=16)
    ax.set_title('Comparison of Spectra for Diff Geometry. 136 keV', fontsize=18)

    coSimple, bins = pa.histogramCo57(e122ptSimple, e136ptSimple)
    coSio2, bins = pa.histogramCo57(e122pt, e136pt)

    fig, ax = plt.subplots()
    ax.hist(bins[:-1], bins=bins, weights=coSimple, histtype='step', linewidth=3, color=bmap[1], label='Simple Detector')
    ax.hist(bins[:-1], bins=bins, weights=coSio2, histtype='step', linewidth=3, color=bmap[0], label='SiO2')
    # ax.set_yscale('log')
    ax.legend(loc='upper left', fontsize=16)
    ax.set_xlabel('Energy (keV)', fontsize=16)
    ax.set_title('Comparison of Spectra for Diff Geometry. Co-57', fontsize=18)


def add_noise_to_pulses(outdir, pattern):

    adcConversion = 0.05 / 1.6e-19 * 750./122.

    # 122 kev line
    outdir =r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\sio2'
    pattern = 'pixel_122kev_sio2_5M_pt\d+.npy'

    data = pa.readBatchFiles(outdir, pattern)

    # Define noise parameters
    mu, sigma = 0., 18.
    dt = np.diff(data[0][0])[0]
    appendTime = np.zeros(int(4/dt))

    for i in range(data.shape[0]):
        pulseLength = len(data[i][1])
        noise = np.random.normal(mu, sigma, pulseLength + 2*len(appendTime))

        # append necessary parts to the pulse
        totLength = len(data[i][0]) + 2*len(appendTime)
        t = np.arange(0, dt*totLength, dt)
        data[i][0] = t
        # scale the pulses (based on Xinran adc) and add noise
        data[i][1] = np.concatenate([appendTime, data[i][1], np.min(data[i][1])*np.ones(len(appendTime))])
        data[i][1] *= adcConversion
        data[i][1] += noise

        # Check to make sure they are the same size
        if len(data[i][1]) != len(data[i][0]):
            if len(data[i][1]) > len(data[i][0]):
                data[i][1] = data[i][1][0:len(data[i][0])]
            else:
                data[i][0] = data[i][0][0:len(data[i][1])]

    return data


def noise_spectrum():

    e122NoNoise = pa.readBatchFiles(r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\sio2', 'pixel_122kev_sio2_5M_pt\d+.npy')
    e136NoNoise = pa.readBatchFiles(r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\136keV\sio2', 'pixel_136kev_sio2_5M_pt\d+.npy')
    # Get energy spectrum with no noise
    e122NoSpec = pa.getEnergySpectrum(e122NoNoise)
    e136NoSpec = pa.getEnergySpectrum(e136NoNoise)

    e122 = add_noise_to_pulses(r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\sio2', 'pixel_122kev_sio2_5M_pt\d+.npy')
    e136 = add_noise_to_pulses(r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\136keV\sio2', 'pixel_136kev_sio2_5M_pt\d+.npy')

    energyConversion = 0.05 / 1.6e-19
    adcConversion = energyConversion * 750/122
    adc2energy = 122./750.
    # 122 kev data
    # Perform fits and evaluate energy
    e122Spectrum = []
    e136Spectrum = []
    for i in range(1000):
        print(i)
        # reject events below like 25 kev
        if -np.min(e122[i][1])*adc2energy < 30:
            pass

        # Fit lines to
        fitobj, result = pa.fitLinearPiecewise(e122[i][0], e122[i][1], 4)

        e122Spectrum.append(-fitobj.predict([result[-2]])[0]*adc2energy)

    fig, ax = plt.subplots()
    ax.hist(e122Spectrum, bins=np.arange(0, 150), color='k', label='Noise', histtype='step', linewidth=3, density=True)
    ax.hist(e122NoSpec, bins=np.arange(0, 150), color='r', label='No Noise', histtype='step', linewidth=3, density=True)
    ax.legend(fontsize=16, loc='upper left')
    ax.set_xlabel('Energy (keV)', fontsize=16)
    ax.set_title('Energy Spectrum Comparison with Noise', fontsize=16)
    plt.show()

    for i in range(1000):
        print(i)
        # reject events below like 25 kev
        if -np.min(e136[i][1])*adc2energy < 30:
            pass

        # Fit lines to
        fitobj, result = pa.fitLinearPiecewise(e136[i][0], e136[i][1], 4)

        e136Spectrum.append(-fitobj.predict([result[-2]])[0]*adc2energy)

    coNoNoise, bins = pa.histogramCo57(e122NoSpec, e136NoSpec)
    coNoise, _ = pa.histogramCo57(e122Spectrum, e136Spectrum)

    fig, ax = plt.subplots()
    ax.hist(coNoise, bins=bins, color='k', label='Noise', histtype='step', linewidth=3, density=True)
    ax.hist(coNoNoise, bins=bins, color='r', label='No Noise', histtype='step', linewidth=3, density=True)
    ax.legend(fontsize=16, loc='upper left')
    ax.set_xlabel('Energy (keV)', fontsize=16)
    ax.set_title('Energy Spectrum Comparison (Co-57)', fontsize=16)

    plt.show()


def compare_angled_spectrum(logplot=False):
    """ Analysis comparing the energy spectrum between normally incident particle source and at an oblique angle
    """
    normaldata = pa.readBatchFiles(r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\sio2\Efield', 'pixel_122kev_sio2_4000V_5M_pt\d+.npy')
    obliquedata = pa.readBatchFiles(r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\sio2\Efield', 'pixel_122kev_sio2_4000V_75deg_200k_pt\d+.npy')

    # get spectrums
    normalSpectrum = pa.getEnergySpectrum(normaldata)
    obliqueSpectrum = pa.getEnergySpectrum(obliquedata)
    bins = np.arange(0,140)

    # Plot spectrum with redsidue
    f, (ax0, ax1) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

    # Plot spectra
    normhist, binx, _ = ax0.hist(np.array(normalSpectrum), bins=bins, density=True, histtype='step', label='Normal Incidence', linewidth=2, color='blue')
    obliquehist, _, _ =ax0.hist(np.array(obliqueSpectrum), bins=bins, density=True, histtype='step', label='Oblique Incidence', linewidth=2, color='red', log=logplot)

    ax0.set_ylabel('Normalized Spectra')
    ax0.set_xlabel('Energy (keV)', fontsize=12)
    ax0.set_title('Pixel Detector Spectra with Different Angle of Incidence', fontsize=14)

    # Calculating error bars
    xaxis = bins[0:-1]+np.diff(bins)
    errNormal = np.sqrt(normhist / len(normalSpectrum))
    errOblique = np.sqrt(obliquehist / len(obliqueSpectrum))
    ax0.fill_between(xaxis, normhist-errNormal, normhist+errNormal, color='blue', alpha=0.5, label='Normal $\sigma$')
    ax0.fill_between(xaxis, obliquehist-errOblique, obliquehist+errOblique, color='red', alpha=0.5, label='Oblique $\sigma$')

    ax0.legend(loc='upper left')

    # Plot residual
    ax1.plot(xaxis, normhist-obliquehist, linewidth=3, color='black',label='residual')
    ax1.fill_between(xaxis, normhist-obliquehist-(errNormal+errOblique), normhist-obliquehist+(errNormal+errOblique), color='black', alpha=0.5, label='Residual $\sigma$')
    ax1.set_title('Difference between Spectra', fontsize=14)
    ax1.set_xlabel('Energy (keV)', fontsize=12)
    ax1.set_ylabel('Residual')
    ax1.legend(loc='lower left')
    f.tight_layout()
    plt.show()

def plot_noise_vs_no_noise():

    data = pa.readBatchFiles(r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\sio2', 'pixel_122kev_sio2_5M_pt\d+.npy')
    dataNoise = add_noise_to_pulses(r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\sio2', 'pixel_122kev_sio2_5M_pt\d+.npy')
    energyConversion = 0.05 / 1.6e-19
    adcConversion = energyConversion * 750/122
    fig, (ax1, ax2) = plt.subplots(1,2, sharey=True)

    for i in range(5):
        start = 11
        ax1.plot(data[start+i][0], data[start+i][1]*adcConversion, linewidth=2)
        ax2.plot(dataNoise[start+i][0], dataNoise[start+i][1], linewidth=2)

    ax1.set_xlabel('Time ($\mu s$)', fontsize=14)
    ax1.set_ylabel('ADC', fontsize=14)
    ax1.set_title('Pulses No Noise', fontsize=14)


    ax2.set_xlabel('Time ($\mu s$)', fontsize=14)
    ax2.set_title('Pulses w/ Noise', fontsize=14)


if __name__ == '__main__':
    # higher_statistics_spectra()
    # new_pixel_particle_spectra()
    # compare_partial_particle_tracks()
    # data = add_noise_to_pulses('test', 'test')
    # noise_spectrum()
    # plot_noise_vs_no_noise()
    compare_angled_spectrum(False)
    plt.show()

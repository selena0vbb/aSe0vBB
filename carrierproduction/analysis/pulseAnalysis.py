import numpy as np
import os
import regex as re
import pwlf
import matplotlib.pyplot as plt
import ROOT as rt
import scipy.signal as scp
import processedPulse as proc


def readBatchFiles(filepath, pattern="*.npy"):
    """
    Reads all .npy data files from a folder that match the pattern.

    Match patter is regex. For example to match the files 122_keV_pixel%i.npy you would use the pattern string '122_keV_pixel\d+.npy'
    """

    # Get all files
    files = [
        f for f in os.listdir(filepath) if os.path.isfile(os.path.join(filepath, f))
    ]

    # Filter by pattern
    r = re.compile(pattern)

    # filters all files in directory with the regex expression
    files2Read = list(filter(r.match, files))

    # for each of the files, reads the data and concatenates into a single array
    for i, fr in enumerate(files2Read):
        if i == 0:
            data = np.load(os.path.join(filepath, fr))
        else:
            dataTemp = np.load(os.path.join(filepath, fr))
            data = np.concatenate((data, dataTemp), axis=0)

    return data


def readBatchFilesCso(filepath, pattern="*.cso"):
    """
    Reads all .cso data files from a folder that match the pattern. Combines (if they have the same settings) into a single carrier simulation object.

    Match patter is regex. For example to match the files 122_keV_pixel%i.cso you would use the pattern string '122_keV_pixel\d+.cso'
    """

    # Get all files
    files = [
        f for f in os.listdir(filepath) if os.path.isfile(os.path.join(filepath, f))
    ]

    # Filter by pattern
    r = re.compile(pattern)

    # filters all files in directory with the regex expression
    files2Read = list(filter(r.match, files))

    for i, fr in enumerate(files2Read):
        if i == 0:
            simOutFile = proc.loadPickleObject(fr)
        else:
            simOutFileTemp = procc.loadPickleObject(fr)
            # Check if settings are the same before combining
            if simOutFile.settings == simOutFileTemp.settings:
                simOutpFile.addPulses(simOutFileTemp.getPulses())
            else:
                Print("Incompatible settings. Files not concatenated.")
                return


def getEnergySpectrum(data, w=0.05):
    """
    From the pulse shapes reconstructs the energy specturm from the given work function
    """

    energy = []
    q = 1.6e-19
    for i in range(max(data.shape)):
        energy.append(abs(data[i][1][-1] / q * w))

    return energy


def histogramCo57(energy122, energy136, inten122=0.866, inten136=0.1068):
    """
    Scale the energy histograms appropriately for Co-57 spectrum
    """
    bins = np.arange(0, 150)

    hist122kev, edges = np.histogram(energy122, bins=bins, density=True)
    hist136kev, _ = np.histogram(energy136, bins=bins, density=True)

    histCo57 = inten122 * hist122kev + inten136 * hist136kev

    return histCo57, bins


def fitLinearPiecewise(t, signal, nPieces=4):

    # Using pwlf module, load the data and fit
    fitObj = pwlf.PiecewiseLinFit(t, signal)
    fitResult = fitObj.fit(nPieces)

    # xx = np.linspace(0, t[-1], 1000)
    # yy = fitObj.predict(xx)

    # fig, ax = plt.subplots()
    # ax.plot(t, signal, linewidth=2, color='k')
    # ax.plot(xx, yy, linewidth=2, color='r')

    # ax.legend(['Raw Data', 'Fit'])
    # ax.set_xlabel('Time ($\mu s$)', fontsize=16)
    # ax.set_ylabel('ADC', fontsize=16)
    # plt.show()

    return fitObj, fitResult


def writeBinFiles(
    outfile,
    infilepath,
    inpattern="*.npy",
    Ns=50000,
    conversion=800 * 0.05 / (122 * 1.6e-19),
    fs=250.0,
    fcut=[5300.0 * 1e-6, 530000.0 * 1e-6],
    addNoise=False,
    filterPulse=False,
    minEnergy=20,
    returnPulse=False,
    returnPulseIdx=0,
    addTwice=False,
):

    data = readBatchFiles(infilepath, inpattern)
    Ns = 50000
    nNoise = 3000
    noise = np.zeros(Ns)
    q = 1.6e-19
    w = 0.05
    N = data.shape[0]

    # parameters for chunking data into smaller files
    numberOfPulsesPerFile = 2500
    maxNumberOfSamples = numberOfPulsesPerFile * Ns
    totalNumberOfSamples = N * Ns
    numberOfIterations = int(np.ceil(float(N) / numberOfPulsesPerFile))
    filesWritten = []

    if addNoise:
        rootfile = r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\Se_2000V_Co57_newfilter_Oct31_1_FFT_500mu_e.root"  # hardcoding rootfile for now
        f = rt.TFile(rootfile)
        setree = f.Get("SeTree")
        noiseLeaf = setree.GetLeaf("wf")
        entryCounter = 0
        nEntries = setree.GetEntries()

    # Iterate over the number of chunked files we need to make
    for j in range(numberOfIterations):
        splitfilename = outfile.split(".")
        chunkOutfile = splitfilename[0] + "%i." % j + splitfilename[1]
        filesWritten.append(chunkOutfile)
        print(filesWritten)
        # Create the numpy array to store all pulses in (save time over constantly appending a list)
        pulses = np.zeros(maxNumberOfSamples, "uint16")
        pulseCounter = 0

        # Iterate over each pulse in the data, convert to ADC unsigned 16bit int
        print(min(numberOfPulsesPerFile, N - numberOfPulsesPerFile * j))
        for i in range(min(numberOfPulsesPerFile, N - numberOfPulsesPerFile * j)):
            cnt = numberOfPulsesPerFile * j + i
            print(i)
            try:
                energyBoolean = (
                    abs(data[cnt][1][np.logical_not(np.isnan(data[cnt][1]))][-1])
                    > minEnergy * q / w
                )
            except:
                energyBoolean = False
            if energyBoolean:
                print(
                    abs(data[cnt][1][np.logical_not(np.isnan(data[cnt][1]))][-1])
                    * w
                    / q
                )
                singlePulse = convertPulse2ADC(
                    data[cnt][0], data[cnt][1], fs=fs, Ns=Ns, conversion=conversion
                )

                if data.shape[1] == 4:
                    singlePulseA = convertPulse2ADC(
                        data[cnt][0], data[cnt][2], fs=fs, Ns=Ns, conversion=conversion
                    )
                    singlePulseB = convertPulse2ADC(
                        data[cnt][0], data[cnt][3], fs=fs, Ns=Ns, conversion=conversion
                    )

                if filterPulse:
                    # Filter the pulse with the given RC parameters
                    singlePulse = butter_bandpass_electronic_filter(
                        singlePulse, lowcut=fcut[0], highcut=fcut[1], fs=fs
                    )

                    if data.shape[1] == 4:
                        singlePulseA = butter_bandpass_electronic_filter(
                            singlePulseA, lowcut=fcut[0], highcut=fcut[1], fs=fs
                        )
                        singlePulseB = butter_bandpass_electronic_filter(
                            singlePulseB, lowcut=fcut[0], highcut=fcut[1], fs=fs
                        )

                if addNoise:
                    k = 0

                    while k < Ns:
                        if k % nNoise == 0:
                            setree.GetEntry(entryCounter)
                            entryCounter = (entryCounter + 1) % nEntries
                            print(noise[0])

                            # matches baseline noise between two events
                            if k != 0:
                                nMean = 30
                                endPreviousSegmentMean = np.mean(noise[k - nMean : k])
                                beginCurrentSegMean = 0
                                for m in range(nMean):
                                    beginCurrentSegMean += noiseLeaf.GetValue(
                                        m % nNoise
                                    )
                                beginCurrentSegMean /= nMean
                                noiseOffset = (
                                    endPreviousSegmentMean - beginCurrentSegMean
                                )
                                print(noiseOffset)
                            else:
                                noiseOffset = 0

                        if addTwice:
                            noise[k] = (
                                noiseLeaf.GetValue(k % nNoise) * 2 - 3600 + noiseOffset
                            )
                        else:
                            noise[k] = noiseLeaf.GetValue(k % nNoise) + noiseOffset

                        # Increment the counter
                        k += 1

                    # add noise to pulse
                    try:
                        singlePulse += noise
                        pulses[
                            pulseCounter * Ns : (pulseCounter + 1) * Ns
                        ] = singlePulse.astype("uint16")
                        if data.shape[1] == 4:
                            singlePulseA += noise
                            singlePulseA.astype("uint16")
                            singlePulseB += noise
                            singlePulseB.astype("uint16")

                        pulseCounter += 1
                    except ValueError:
                        pass

            if returnPulse:
                if i == returnPulseIdx:
                    if data.shape[1] == 4:
                        return (singlePulse, singlePulseA, singlePulseB)
                    else:
                        return singlePulse
        pulses[0 : pulseCounter * Ns].tofile(chunkOutfile)

    # read the binary files using file operations and keep writing to the same file
    # combineDatFiles(filesWritten, outfile)


# conversion factor calculation 122 kev: 1.6e-19 * 122/0.05 = 500 ADC
# Therefore 1e = 500 ADC * 0.05 energy per charge / 122


def combineDatFiles(filelist, outfile):

    finalFile = open(outfile, "w")

    for file in filelist:
        partialFile = open(file, "rb")
        data = partialFile.read()
        finalFile.write(data)
        partialFile.close()

    finalFile.close()


def generateRootFile(path, pattern, outputfile):

    data = readBatchFiles(path, pattern)

    # Iterate over data to get depth
    depth = []
    energy = []
    q, w = (1.6e-19, 0.05)
    print(data.shape[0])
    for i in range(5):  # data.shape[0]):
        print(i)
        time = data[i][0]
        pulse = data[i][1]
        if data[i][0].shape == (1,):
            pass
        else:
            fitobj, result = fitLinearPiecewise(time, pulse, 2)
            depth.append(fitobj.predict([result[1]])[0])
            energy.append(abs(data[i][1][-1] / q * w))

    h = rt.TH1D("eneSim", "Simulation Spectra", 140, 0, 140)
    rt.gROOT.ProcessLine(
        "struct simdata{ \
             Int_t depthT; \
             Int_t energyT; \
             }; "
    )
    simdata = rt.simdata()
    tree = rt.TTree("trueData", "trueData")
    tree.Branch("depthT", rt.AddressOf(simdata, "depthT"))
    tree.Branch("energyT", rt.AddressOf(simdata, "energyT"))

    for i in range(len(energy)):
        h.Fill(energy[i])
        simdata.depthT = depth[i]
        simdata.energyT = energy[i]

        tree.Fill()

    tree.Print()
    tf = rt.TFile(outputfile, "RECREATE")
    h.Write()
    tree.Write()
    tf.Close()


def convertPulse2ADC(
    t, signal, fs=250, Ns=50000, prePulseTime=100, conversion=1, datatype="int16"
):

    # Check to make sure there is data in the pulse
    if t.size <= 1:
        return np.array([]).astype(datatype)
    sampleRate = 1.0 / fs
    tTotal = Ns * sampleRate
    dt = np.diff(t)[0]

    # Create the pulse to pre and post append
    prePulse = np.zeros(int(prePulseTime / dt))
    postPulse = np.min(signal) * np.ones(
        int(round((tTotal - t.size * dt - prePulseTime) / dt))
    )

    # Pre and post append the pulse
    newSignal = np.concatenate([prePulse, signal, postPulse])
    newTime = np.arange(0, tTotal, dt)

    # Convert the signal to adc values
    newSignal *= conversion

    # Create new time axis corresponding to adc sampling rate
    # Interpolate the signal with the new rate
    tADC = np.arange(0, sampleRate * Ns, sampleRate)
    signalADC = np.interp(tADC, newTime, newSignal)

    # Return the ADC sampled pulses
    return (signalADC).astype(datatype)


def butter_bandpass_electronic(lowcut, highcut, fs, order=1):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = scp.butter(order, [low, high], btype="band")
    return b, a


def butter_bandpass_electronic_filter(data, lowcut, highcut, fs, order=1):
    """ Applies a butterworth bandpass of order=order and lowcut and highcut cutoff frequencies to the daata
    """

    b, a = butter_bandpass_electronic(lowcut, highcut, fs, order=order)
    # zi = scp.lfilter_zi(b, a)
    y = scp.lfilter(b, a, data)  # , zi=zi*data[0])
    # y = scp.filtfilt(b, a, data)
    return y


if __name__ == "__main__":
    testdata = readBatchFiles(
        r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\sio2\Efield",
        "pixel_122kev_sio2_4000V_5M_pt\d+.npy",
    )
    # energy = getEnergySpectrum(testdata)
    # plt.hist(energy, bins=np.arange(0,140))
    # plt.show()
    # fitLinearPiecewise(testdata[4][0], testdata[4][1], 2)
    # readBinFiles(r'C:\Users\alexp\Documents\UW\Research\Selenium\realData\Se_2500V_Co57_newfilter_Oct31_1.dat')
    # writeBinFiles(r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\co57\co57_4000V_75deg_bin.dat',r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\co57', 'pixel_Co57_sio2_4000V_75deg_200k_pt\d+.npy', addNoise=True, filterPulse=True)
    # writeBinFiles(r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\sio2\Efield\122_4000V_bin.dat', r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\sio2\Efield', 'pixel_122kev_sio2_4000V_5M_pt\d+.npy', addNoise=True, filterPulse=True)
    generateRootFile(
        r"C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\sim_data\pixel_detector\122keV\sio2\Efield",
        r"pixel_122kev_sio2_4000V_5M_pt\d+.npy",
        "./pixel_122kev_sio2_4000V_sim_results.root",
    )

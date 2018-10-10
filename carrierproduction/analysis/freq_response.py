import numpy as np
import matplotlib.pyplot as plt
import brewer2mpl
# import particledata as pd

def signalFFT(signal, time, n=None):

    dt = np.diff(time)[0]
    print(dt)


    if n:
        if n > signal.size:
            signal = np.concatenate((signal, signal[-1]*np.ones(n-signal.size)))


    # Compute fourier tranform
    yf = np.fft.fft(signal, n)
    xf = np.fft.fftfreq(n, dt)

    return xf, yf





def readNoiseFile(filename):
    count = 0
    freq = []
    ampl = []

    f = open(filename)

    for line in f:
        if count < 1:
            count += 1
        else:
            count += 1
            freq.append(np.array(line.split(',')).astype(float)[0])
            ampl.append(np.array(line.split(',')).astype(float)[1])

    f.close()

    return np.array(freq), np.array(ampl), count

def testReadNoiseFile():
    filename = './Noise_withCurrent.txt'
    f, amplitude, size = readNoiseFile(filename)
    print(size)
    fig, ax = plt.subplots()
    ax.plot(f,amplitude)
    ax.set_xscale('log')
    plt.show()


def test_filter_algorithm():
    filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\charge collection\noise_histograms\data\coplanar_scaling_particle_event.npy'

    filterfilename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Circuits\Draft22.txt'

    data = np.load(filename)

    t = data[0][0] * 1e-6
    sig = data[0][1]
    omega = 40e6


    testsig = np.sin(omega*t)
    plt.plot(t, testsig)
    plt.show()
    # signalFilt, signalFiltFreq, xf, yf = pd.filterSignal(sig, t, filterfilename )

    return None

def FourierTransformTestSignals():

        filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\charge collection\noise_histograms\data\coplanar_scaling_particle_event.npy'

        data = np.load(filename)

        # Plot all time signals
        fig, (ax1, ax2, ax3) = plt.subplots(3,1)


        for run in data[0:20]:
            t = run[0]
            dt = np.diff(t)[0]
            signal = run[1]
            dsignaldt = np.diff(signal)/dt*1e6


            # tindx = np.nonzero(np.abs(signal) > 0.01*np.max(abs(signal)))[0][0]
            tindx = 1
            timeLength = 99

            ax1.plot(t[:timeLength], -signal[:timeLength], color='grey', alpha=0.3)
            ax1.set_xlabel('t ($\mu s)')
            ax1.set_ylabel('q(t) - (C)', fontsize=14)

            ax2.plot(t[:timeLength-tindx], -dsignaldt[tindx:timeLength]*1e9, color='red', alpha=0.3)
            ax2.set_ylabel('i(t) - (nA)', fontsize=14)
            ax2.set_xlabel('t ($\mu s)')

            # Compute Fourier transform of all Signals
            n = 100000
            xf, yf = signalFFT(dsignaldt, t*1e-6, n=n)


            ax3.plot(xf[:len(xf)//2], abs(yf[:len(yf)//2]), color='blue', alpha=0.3)
            ax3.set_ylabel('I($j\omega$)', fontsize=14)
            ax3.set_xlabel('f (MHz)')


        plt.show()


def getNoise(psFilename, size, dt):

    # Generate white gaussian noise
    nt = np.random.normal(size=size)

    # Read power spectruem
    f, ps, _ = readNoiseFile(psFilename)

    # Get frequency domain of time signal
    nfy = np.fft.fft(nt)
    nfx = np.fft.fftfreq(nt.size, dt)

    # Apply powerspectrum to noise
    nfyPowerSpectrum = nfy * np.interp(nfx, f, 10**(ps/10))

    # inverse fourier tranform to get time
    ntPowerSpectrum = np.fft.ifft(nfyPowerSpectrum)




    # fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)

    # # # plot spectrum of noise
    # t=np.arange(0, size*dt, dt)*1e6

    # ax1.plot(t, nt, linewidth=2)
    # ax1.set_ylabel('White Noise Amplitude', fontsize=14)
    # ax2.plot(t, np.real(ntPowerSpectrum), 'k', linewidth=2)
    # ax2.plot(t, np.imag(ntPowerSpectrum), 'r', linewidth=2)
    # ax2.set_ylabel('Filtered Noise Amplitude', fontsize=14)
    # ax2.set_xlabel('Time ($\mu s$)', fontsize=14)
    # ax2.legend(['Real', 'Imaginary'])
    # ax1.hist(nt, bins=100, color='b', histtype='step')
    # ax2.hist(abs(ntPowerSpectrum), bins=100, color='r', histtype='step')
    # ax.plot(f, ps, linewidth=2, color='r')
    # ax.plot(nfx[:nfx.size//2], 10*np.log10(abs(nfy[:nfy.size//2])), linewidth=2, color='b')
    # ax.set_xscale('log')

    plt.show()
    return ntPowerSpectrum


def applyNoisetoSignal():

    filename = r'C:\Users\alexp\Documents\UW\Research\Selenium\Coplanar Detector\charge collection\noise_histograms\data\coplanar_scaling_particle_event.npy'

    data = np.load(filename)

    psfilename = './Noise_withCurrent.txt'

    t = data[1][0]
    sig = data[1][1]
    dt = np.diff(t)[0]*1e-6
    current = np.diff(sig)/dt

    noise = getNoise(psfilename, 3000, dt)


    totalSig = current + np.real(noise[:current.size])

    fig, ax = plt.subplots()
    ax.plot(t[:-1], -current*1e9, linewidth=2, color='red')
    ax.plot(t[:-1], -totalSig*1e9, linewidth=2, color='black')
    ax.set_xlabel('Time ($\mu s$)', fontsize=16)
    ax.set_ylabel('Current (nA)', fontsize=16)
    ax.set_title('Current Simulation Signal w/ Noise', fontsize=18)
    # ax.plot(t[:-1], np.array([current, totalSig]).T, linewidth=2)

    ax.legend(['No Noise', 'Noise'])
    plt.show()




if __name__ == '__main__':

    # FourierTransformTestSignals()
    # testReadNoiseFile()
    # getNoise('./Noise_withCurrent.txt', 1000, 2.e-8)
    applyNoisetoSignal()

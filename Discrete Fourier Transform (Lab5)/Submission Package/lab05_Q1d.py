import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft, rfftfreq
import seaborn as sns

# for formatting the graph
sns.set()


def program(instrument_name):
    # load the file
    waveform = np.loadtxt(instrument_name + '.txt')
    
    # plot the signal
    plt.figure(figsize=(8,5))
    plt.plot( np.arange(len(waveform)), waveform)
    plt.title("Waveform of {}".format(instrument_name), fontsize=16)    
    plt.ylabel('Waveform value', fontsize=14)
    plt.xlabel("Signal index", fontsize=14)
    plt.tight_layout()
    plt.savefig(instrument_name + 'a.pdf')
    
    # apply fourier transform
    ck = rfft(waveform)
    # plot for the first 10,000 coeff
    plt.figure(figsize=(8,5))
    plt.plot(np.arange(10000), np.abs(ck[:10000]))
    plt.title("Magnitude of DFT coefficients for {}".format(instrument_name), fontsize=14) 
    plt.ylabel("Magnitude of the coefficient ($|C_K|$)", fontsize=14)
    plt.xlabel("$K$", fontsize=14)
    plt.tight_layout()
    plt.savefig(instrument_name + 'b.pdf')
    
    ########## Part b (find the note of the instrument) ##########
    # sampling rate for audio signal
    sample_rate = 44100 # Hz
    # Get the frequency using the inbuilt function
    freq = rfftfreq(len(waveform), 1/sample_rate)
    # get power spectrum
    ck_squared = np.abs(ck)**2
    
    # plot |ck|^2 vs freq
    plt.figure(figsize=(8,5))
    plt.plot(freq, ck_squared)
    # plt.tight_layout()
    plt.title(" $|c_k|^2$ vs frequency for {}".format(instrument_name), fontsize=14) 
    plt.ylabel("Magnitude squared ($|C_K|^2$)", fontsize=14)
    plt.xlabel("frequency", fontsize=14)
    plt.savefig(instrument_name + 'c.pdf')
    
    # The first peak we are interested in is greater than 0.2 * 10^16 for both 
    # the instrunments.
    peak_threshold =  0.2 * 10**16
    # Find the location of peaks greater than this threshold
    peak_index = np.where(ck_squared > peak_threshold)
    # the index of the first peak. Take min of the first element as np.where 
    # returns a tuple
    first_index = min(peak_index[0])
    # corresponding frequency
    freq_first_peak = freq[first_index]

    # print the frequency of the first peak
    print("The frequency of the first peak magnitude |C_K|^2 for {} is {}."\
          .format(instrument_name, freq_first_peak))
        
    # find the note corresponding to this frequency using the formula given in 
    # the handout
    note = round(12 * np.log2(freq_first_peak / 440))
    
    print("The note number being played by the {} is {}".format(instrument_name,
                                                                   note))
    
    
# run the above program for piano and trumpet
program('piano')
program('trumpet')

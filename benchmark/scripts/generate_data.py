import numpy as np
from matplotlib import pyplot as plt

class GaussianClass:
    def __init__(self):
        self.position = 0.0
        self.sigma = 1.0
        self.amplitude = 1.0
    
    def value(self, freq):
        return np.exp(-0.5*(freq - self.position)**2/self.sigma**2)/(self.sigma*(2*np.pi)**0.5)

if "__main__":

    ###############################################################
    ##          Generate synthetic spectral functions
    ###############################################################

    # class and params for artificial spectral function
    # two Gaussian peeks included
    peek1, peek2 = GaussianClass(), GaussianClass()
    peek1.position, peek2.position = -2.0, 2.0
    peek1.sigma, peek2.sigma = 1.0, 1.0
    peek1.amplitude, peek2.amplitude = 1.0, 1.0

    # grids params
    freq_min, freq_max = -10.0, 10.0
    freq_num = int(1e3)

    spectrum = np.empty(shape=[2, freq_num], dtype=float) 
    spectrum[0,:] = np.linspace(freq_min, freq_max, freq_num, endpoint=False)
    spectrum[1,:] = 0.5 * (peek1.value(spectrum[0,:]) + peek2.value(spectrum[0,:]))

    # write spectral function to file
    out_path_spec = "data/exact_spec.dat"
    with open(out_path_spec, mode='w') as outfile:
        for idf in range(freq_num):
            outfile.write("{:>15d}{:>20.8f}{:>20.8f}\n".format(idf,spectrum[0,idf],spectrum[1,idf]))
    outfile.close()


    ###############################################################
    ##          Generate biased correlation functions
    ###############################################################
    
    # calculate correlation functions from synthetic spectral function
    beta = 8.0
    lt = 160
    tau_grids = np.linspace(0.0, beta, lt, endpoint=False)

    # generate kernel matrix
    freq_grids = spectrum[0,:]
    kernel = np.exp(-freq_grids * np.array(np.mat(tau_grids).transpose())) / (1.0+np.exp(-beta*freq_grids))
    correlation = (freq_max-freq_min)/freq_num * kernel.dot(spectrum[1,:])
    
    # add Gaussian-like noise 
    # generate noisy data in bins
    nbin = 1000
    sigma = 0.01
    bootstrap_sample = np.empty(shape=[nbin, lt], dtype=float)
    for idb in range(nbin):
        bootstrap_sample[idb,:] = correlation + np.random.normal(loc=0.0, scale=sigma, size=lt)

    # output of biased corrlation functions
    out_path_corr = "../input/benchmark/cor.dat"
    with open(out_path_corr, mode='w') as outfile:
        outfile.write("{:>20d}{:>20d}\n".format(nbin, lt))
        for idb in range(nbin):
            for idt in range(lt):
                outfile.write("{:>20d}{:>20d}{:>20.10f}\n".format(idb, idt, bootstrap_sample[idb,idt]))
    outfile.close()

    # output of imaginary-time grids
    out_path_tau = "../input/benchmark/tau.dat"
    with open(out_path_tau, mode='w') as outfile:
        outfile.write("{:>20d}{:>20.5f}\n".format(lt, beta))
        for idt in range(lt):
            outfile.write("{:>20d}{:>20.10f}\n".format(idt, tau_grids[idt]))
    outfile.close()

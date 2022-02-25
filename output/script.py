import numpy as np
from matplotlib import pyplot as plt


def read_data(filename):
    x = []
    y = []

    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                data = line.split()
                if ( abs(float(data[1])) < 10.0):
                    x.append(float(data[1]))
                    y.append(float(data[2]))
    file.close()
    return x, y


def plot(x, y, label):
    plt.figure()
    plt.grid(linestyle='-.')
    plt.plot(x, y, 's', ms=4, label=label, ls=":")
    plt.xlabel('${\omega}$\n', fontsize = 13)
    plt.ylabel(label, fontsize = 13)
    plt.legend(fontsize = 12)
    plt.savefig('L4b4U4k1.01.0/fig.pdf')
    plt.show()


def analyse(filename, u):
    # read spectrum from data
    freq, spec  = read_data(filename=filename)
    
    # measure freqency in uint of interaction strength U
    freq = [f/u for f in freq]
    # normalization of spectrum
    sum_of_spec = sum(spec)
    spec = [s/sum_of_spec for s in spec]
    print(sum_of_spec)

    return freq, spec

if __name__ == "__main__":

    uint = 1.0
    spec_interval = 1e-2

    # freq1, spec1 = analyse(filename="L8b4U-4k0.500.50/spec.dat", u=uint)
    # freq2, spec2 = analyse(filename="L8b4U-4k0.750.75/spec.dat", u=uint)
    # freq3, spec3 = analyse(filename="L8b4U-4k1.001.00/spec.dat", u=uint)
    # freq4, spec4 = analyse(filename="L8b4U-4k0.250.25/spec.dat", u=uint)

    freq, spec = analyse(filename="benchmark/spec.dat", u=uint)

    freq_base, spec_base = [], []
    with open("../input/benchmark/spec_base.dat", mode='r') as file:
        for line in file:
            if len(line) != 0:
                data = line.split()
                freq_base.append(float(data[0]))
                spec_base.append(float(data[1]))
    file.close()
    spec_base = [ spec_interval*s for s in spec_base ]

    plt.figure()
    plt.grid(linestyle='-.')

    # plt.plot(freq1, spec1, linewidth=3.0, label="${k=(\pi/2, \pi/2)}$")
    # plt.plot(freq2, spec2, linewidth=3.0, label="${k=(3\pi/4, 3\pi/4)}$")
    # plt.plot(freq3, spec3, linewidth=3.0, label="${k=(\pi, \pi)}$")
    # plt.plot(freq4, spec4, linewidth=3.0, label="${k=(\pi/4, \pi/4)}$")

    plt.plot(freq_base, spec_base, linewidth=3.0, linestyle="-.", label="Synthetic")
    plt.plot(freq, spec, linewidth=2.0, label="SAC")
    
    plt.xlabel('${\omega}$\n', fontsize = 13)
    plt.ylabel("${A(\omega)}$", fontsize = 13)
    plt.legend(fontsize = 12)
    plt.title("Benchmark", fontsize=13)
    plt.xlim(xmin=-6.0, xmax=10.0)
    plt.savefig('benchmark/fig.pdf')
    plt.show()
    
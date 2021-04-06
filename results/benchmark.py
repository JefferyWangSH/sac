# Script to generate standard fermion spectrum functions.

# For benchmark use.

import random
import numpy as np
from matplotlib import pyplot as plt

def benchmarkA(omega_list, peak_position, peak_width):
    # omega: frequency
    # peak_position: position of Gaussian peak
    # peak_width: half-height width of the peak
    return [1/2 * abs(peak_width)*(2*np.pi)**(-0.5) * np.exp(-(omega-peak_position)**2 / (2*peak_width**2)) for omega in omega_list]


def randomError(errorMin, errorMax):
    return (-1)**random.randint(0, 1) * (errorMin + random.random() * (errorMax - errorMin))


def BenchmarkG(tau_list, omega_list, A_list, beta):
    G_list = []
    err_list = []

    for t in range(len(tau_list)):
        G_temp = 0.0
        for i in range(len(A_list)):
            G_temp += A_list[i] * (np.exp(-tau_list[t]*omega_list[i]) + np.exp(-(beta-tau_list[t])*omega_list[i])) / np.pi
        # random error of order 10**-3
        error = randomError(1.0*G_temp*10**(-2), 5.0*G_temp*10**(-2))
        G_temp += error
        G_list.append(G_temp)
        err_list.append(abs(error))
    return G_list, err_list


def writeBenchmarkData(filename, x_list, y_list, z_list=[]):
    with open(filename, mode='w') as file:
        assert (len(x_list) == len(y_list))
        assert (len(z_list) == 0 or len(z_list) == len(x_list))

        for i in range(len(x_list)):
            str_temp = (str(x_list[i]) + 20*" ")[0:23:1]
            file.write(str_temp)
            if len(z_list) == 0:
                file.write(str(y_list[i]) + "\n")
            else:
                file.write(str(y_list[i]) + 20*" " + str(z_list[i]) + "\n")
    file.close()


def plotBenchmark(xlist, ylist, xlabel, ylabel):
    plt.figure()
    plt.grid(linestyle='-.')

    plt.plot(xlist, ylist, 's', ms=4, label=ylabel, ls=":")

    plt.xlabel(xlabel + '\n')
    plt.ylabel(ylabel)

    plt.legend()
    plt.show()


# the main program
if __name__ == "__main__" :

    beta = 4.0
    ll = 80

    omega_list = [0.1 * i for i in range(1, 51)]
    tau_list = [ i*beta/ll for i in range(ll)]

    A_list = benchmarkA(omega_list, peak_position=2.5, peak_width=0.4)
    G_list, err_list = BenchmarkG(tau_list, omega_list, A_list, beta=4.0)

    plotBenchmark(omega_list, A_list, xlabel="${\\omega}$", ylabel="A(${\\omega}$)")

    plotBenchmark(tau_list, G_list, xlabel="${\\tau}$", ylabel="G(${\\tau}$)")

    writeBenchmarkData("benchmark_a.txt", omega_list, A_list)

    writeBenchmarkData("benchmark_g.txt", tau_list, G_list, err_list)

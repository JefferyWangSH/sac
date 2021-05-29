import numpy as np
from matplotlib import pyplot as plt


# function to read hamiltonian data from input file
def read_hamiltonian_from_file(filename):

    p = []
    lnU = []
    errlnU = []

    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split()
                p.append(float(strList[0]))
                lnU.append(float(strList[1]))
                errlnU.append(float(strList[2]))
    file.close()

    return p, lnU, errlnU


# function to read configs data n(x) from input file
def read_configs_from_file(filename):
    
    # key:   alpha in p
    # value: [x, n]
    data = {}

    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split()

                if float(strList[0]) in data.keys():
                    data[float(strList[0])][0].append(float(strList[1]))
                    data[float(strList[0])][1].append(float(strList[2]))
                else:
                    data[float(strList[0])] = [[float(strList[1])], [float(strList[2])]]
                
    file.close()

    p = list(data.keys())
    sorted(p, reverse=False)

    return p, data


if __name__ == "__main__":

    # read data from file
    p, lnU, errlnU = read_hamiltonian_from_file(filename = "h-alpha.txt")
    p, data = read_configs_from_file(filename="configs-alpha.txt")

    # model params
    nalpha = len(p)
    nconfig = len(data[0][0])

    omega_min = -10.0
    omega_max =  10.0

    # critical alpha where phase transition occurs
    p_c = 50

    # recover fermion spectrum from averaged configs n(x)
    x = data[0][0]
    omega = [ omega_min + (omega_max - omega_min) * _x for _x in x]
    n = []
    n_alpha = [data[alpha][1] for alpha in p]
    
    for i in range(len(x)):
        n_averaged = 0.0
        for _p in range(p_c, nalpha - 1):
            n_averaged = n_averaged + ( np.exp(lnU[_p]) - np.exp(lnU[_p+1]) ) * n_alpha[_p][i]
        n_averaged = n_averaged / ( np.exp(lnU[p_c]) - np.exp(lnU[nalpha-1]) )
        n.append(n_averaged)

    # plot figures
    plt.figure()
    plt.grid(linestyle='-.')

    plt.errorbar(p, lnU, errlnU, label='${\\mathrm{ln} U}$', ms=3, fmt='o-', elinewidth=2, capsize=4, color='b', ecolor='r')

    plt.xlabel('${p}$\n', fontsize = 13)
    plt.ylabel('${\\mathrm{ln} U}$', fontsize = 13)

    plt.legend(fontsize = 12)
    plt.show()

    plt.figure()

    plt.plot(omega, n, ms=3, label="n(x)")
    plt.xlabel("${\\omega}$", fontsize = 13)
    plt.ylabel("${\\rho (\\omega)}$", fontsize = 13)

    # plt.plot(omega, data[80][1])

    plt.show()

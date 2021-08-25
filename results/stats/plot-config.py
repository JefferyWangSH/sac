import numpy as np
from matplotlib import pyplot as plt


def readData(filename):

    theta = []
    chi_2 = []
    err_chi_2 = []
    s = []
    err_s = []

    dirtData = {}
    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split()
                if float(strList[1]) > -15:
                    data = [float(strList[i+2]) for i in range(2)]
                    err = [float(strList[i+4]) for i in range(2)]
                    dirtData[float(strList[1])] = [data, err]
    file.close()
    
    theta = [key for key in dirtData.keys()]
    theta.sort(reverse = False)
    for key in dirtData.keys():
        chi_2.append(dirtData[key][0][0])
        s.append(dirtData[key][0][1])
        err_chi_2.append(dirtData[key][1][0])
        err_s.append(dirtData[key][1][1])

    return theta, chi_2, err_chi_2, s, err_s


def plotFigure(theta, obs, err, label):
    plt.figure()
    plt.grid(linestyle='-.')

    # plt.plot(omega, A, 's', ms=4, label="${A(\\omega)}$", ls = ":")
    plt.errorbar(theta, obs, err, label=label, ms=3, fmt='o', ecolor='r', color='b', elinewidth=1.5, capsize=4)

    plt.xlabel('${ln(1/\\Theta)}$\n')
    plt.ylabel(label)

    plt.legend()
    plt.show()


if __name__ == "__main__":

    theta, chi_2, err_chi_2, s, err_s = readData("output-cst2.txt")

    plotFigure(theta, chi_2, err_chi_2, "${\\chi^{2}}$")
    plotFigure(theta, s, err_s, "${S}$")

import numpy as np
from matplotlib import pyplot as plt


def readDate(filename):

    p = []
    lnU = []

    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split()
                p.append(float(strList[0]))
                lnU.append(float(strList[1]))
    file.close()

    return p, lnU


if __name__ == "__main__":


    p, lnU = readDate("test.txt")

    plt.figure()
    plt.grid(linestyle='-.')

    plt.plot(p, lnU, 's', ms=4, label='${\\mathrm{ln} U}$', ls=":")

    plt.xlabel('${p}$\n', fontsize = 13)
    plt.ylabel('${\\mathrm{ln} U}$', fontsize = 13)

    plt.legend(fontsize = 12)
    plt.show()

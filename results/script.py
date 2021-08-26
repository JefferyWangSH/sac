import numpy as np
from matplotlib import pyplot as plt


def read_data(filename):
    x = []
    y = []

    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                data = line.split()
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
    plt.savefig('fig.pdf')
    plt.show()


if __name__ == "__main__":
    freq, spec = read_data(filename="spec.dat")
    print(sum(spec))
    plot(freq, spec, label="${A(\omega)}$")
    
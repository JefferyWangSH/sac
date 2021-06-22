import numpy as np
from matplotlib import pyplot as plt

def read_configs_from_file(filename):

    x = []
    n = []

    with open(filename, mode='r') as file:
        for line in file:
            if len(line) != 0:
                strList = line.split()
                x.append(float(strList[0]))
                n.append(float(strList[1]))
    file.close()

    return x, n


if __name__ == "__main__":
    
    x, n = read_configs_from_file(filename = "./configs/alpha_-29.80.txt")

    omega_min = -8.0
    omega_max = +8.0

    omega = [ omega_min + (omega_max - omega_min) * xi / len(x) for xi in x]

    plt.figure()
    
    plt.title("${\\mathbf{k} = (\\pi/2, \\pi/2)}$", fontsize = 12)

    plt.plot(omega, n)
    plt.xlabel("${\\omega}$", fontsize = 13)
    plt.ylabel("${\\mathrm{A}(\\mathrm{k}, \\omega)}$", fontsize = 13)

    plt.legend()
    plt.show()

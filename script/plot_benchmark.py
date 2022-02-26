import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

if "__main__":
    
    ###############################################################
    ##          Read spectral functions from file
    ###############################################################
    exact_spec_path = "data/exact_spec.dat"
    sac_spec_path = "../output/benchmark/spec.dat"

    exact_data = np.loadtxt(fname=exact_spec_path, dtype=float, skiprows=(0), usecols=(1,2))
    sac_data = np.loadtxt(fname=sac_spec_path, dtype=float, skiprows=(0), usecols=(1,2))
    exact_freq, exact_spec = exact_data[:,0], exact_data[:,1]
    sac_freq, sac_spec = sac_data[:,0], sac_data[:,1]
    

    ###############################################################
    ##          Plot the benchmark results
    ###############################################################

    # font convention
    font = {'family': 'Times New Roman', # or 'Arial' 'serif'
            'weight': 'regular', # or 'normal' 'bold' 
            # 'style': 'italic',
            'size': 16, 
            }
    # mpl.rcParams['mathtext.default'] = 'regular'
    # mpl.rcParams['xtick.direction'] = 'in'
    # mpl.rcParams['ytick.direction'] = 'in'
    # plt.style.use('seaborn-paper')
    mpl.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Times"],
    })

    # create figure and axis
    fig, ax = plt.subplots(1, 1, figsize=(8,6))
    ax.tick_params(axis='both', which='major', direction='in', length=5, width=1.5, labelsize=16)

    frame_size = 1.5
    ax.spines['left'].set_linewidth(frame_size)
    ax.spines['right'].set_linewidth(frame_size)
    ax.spines['top'].set_linewidth(frame_size)
    ax.spines['bottom'].set_linewidth(frame_size)
    
    pic_exact, = ax.plot(exact_freq, exact_spec, linewidth=2.5, label=r"Synthetic")
    pic_sac, = ax.plot(sac_freq, sac_spec, linewidth=2.5, label=r"SAC")

    ax.set_xlabel(r"$\omega$", font, fontsize=23)
    ax.set_ylabel(r"$A(\omega)$", font, fontsize=23)
    ax.set_title(r"Benchmark", font, fontsize=18)
    ax.legend( loc='best', labelspacing=0.4, markerfirst=True,
               frameon=True, fancybox=False, edgecolor='black', handlelength=1.5, 
               ncol=1, columnspacing=0.8, prop=font )
    # ax.grid(linestyle='-.')
    fig.tight_layout()
    plt.savefig("figure/benchmark.pdf", format='pdf', dpi=1200)
    plt.show()

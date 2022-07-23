# Stochastic Analytic Continuation
![workflow](https://github.com/JefferyWangSH/sac/actions/workflows/main.yml/badge.svg?branch=master)
![commit](https://img.shields.io/github/last-commit/JefferyWangSH/sac/master?color=blue)

We present in this repository a general C++ implementation of 
the Stochastic Analytic Continuation `(SAC)` algorithm, proposed by Anders W. Sandvik.

This program aims to extract the spectral information from the imaginary-time correlation functions calculated by quantum Monte Carlo `(QMC)` simulations 
in strongly correlated quantum systems.
Details about the `SAC` algorithm can be found in the [References](#references) section. 

---

## Installation ##

### Prerequisite ###
* `gcc/g++` `( version >= 7.1, support C++17 standard )` and `cmake` `( version >= 3.21 )` installed.
* `Boost C++ libraries` `( version >= 1.71 )` installed.
* `Eigen library` `( version >= 3.3.9 )` providing a friendly API of matrix structure
* `Intel Math Kernel Library (MKL)` for high-accuracy linear algebra, working as the backend of the `Eigen` library.

### Usage ###
1. Download the source code from github:
    ``` shel
    $ # download the source code
    $ git clone -b master https://github.com/JefferyWangSH/sac.git {program_root}
    ```

2. Enter the [`build`](build/) directory and run [`runcmake.sh`](build/runcmake.sh) to config the program with cmake.
    ``` shell
    $ # config cmake
    $ cd {program_root}/build && sh runcmake.sh
    ```
    
3. Enter the [`run`](run/) directory and compile the codes. 
    ``` shell
    $ # build the program
    $ cd {program_root}/run && sh make.sh
    ```

4. 
   We use `slurm` for task managements, and some example scripts for the submissions are presented in [`run/submit.slurm`](run/submit.slurm) and [`run/run.sh`](run/run.sh). One can also run the code directly from the command line.
   ```shell
   $ # run the program with slurm
   $ cd {program_root}/run
   $ # edit simulating parameters before running the code 
   $ vim submit.slurm
   $ sh submit.slurm
   
   $ # run directly from the command line
   $ # using program option `--help` to see the optional params and the helping messages.
   $ {program_root}/build/sac  --help
   ``` 

5. 
   The parameters of the SAC algorithm are organized in the TOML (v1.0.0) format, 
   and one should provide a TOML configuration file to the program through the option `--config` ( See [`benchmark/config.toml`](benchmark/config.toml) for example ).
   
   And to make the input QMC data detectable for the program, one should pass the path of the input QMC corrlation functions and the corresponding imaginary-time grids through the command line options `--corr` and `--tgrids`.
    ```shell
    $ # specify paths of input and output files
    $ {program_root}/build/sac \
        --config=${config_file} --tgrids=${tgrids_file} --corr=${corr_file} \
        --log=${log_file} --spec=${spec_file} --report=${report_file}
    ```

   During the SAC simulation, the history of the simulated annealing will be recorded in the output log file, which can be assigned through the option `--log`.
   After the simulation is finished, results of the recovered spectrum ( `--spec` ), together with a quality report of the SAC recovery ( `--report` ), will be automatically generated.
   ```shell
   $ # by default, the output files would be put under the [benchmark] folder unless explicitly assigned
   $ cd {program_root}/benchmark
   
   $ # monitor the process of annealing in log.out
   $ vim log.out

   $ # if someone runs the program using slurm, the screen output will be stored in stdout.out
   $ vim stdout.out
   
   $ # recovered spectral functions
   $ vim spec.out

   $ # quality report of the recovered spectrum to quantify the goodness of the SAC recovery
   $ vim report.out
   ```


## Repository Structure ##
* [`include/`](include) & [`src/`](src) - header files and source codes 
  * [`sac_core.h`](include/sac_core.h)/[`sac_core.cpp`](src/sac_core.cpp) 
    \- the core class for the SAC simulations, implementing the MonteCarlo sampling, simulated annealing, and the recovery of spectral functions
  * [`qmc_reader.h`](include/qmc_reader.h)/[`qmc_reader.cpp`](src/qmc_reader.cpp)
    \- independent module to analyse and preprocess the input QMC correlation functions
  * [`freq_grids.h`](include/freq_grids.h)/[`freq_grids.cpp`](src/freq_grids.cpp)
    \- grids class for the discretization of ( both hyperfine and spectral ) grids in the frequency domain
  * [`sac_kernel.h`](include/sac_kernel.h)/[`sac_kernel.cpp`](src/sac_kernel.cpp)
    \- the integration kernel of SAC, which connects the spectral functions and the QMC correlations. The kernel may vary for different physical quantities, and basic kernels for fermionic and bosonic quantum systems are provided by default. One can also add customized types of kernels to this class.
  * [`sac_annealing.h`](include/sac_annealing.h)/[`sac_annealing.cpp`](src/sac_annealing.cpp)
    \- the chain structure for the storage of metadata during the simulated annealing process of SAC
  * [`sac_measure.h`](include/sac_measure.h)/[`sac_measure.cpp`](src/sac_measure.cpp)
    \- class for measuring the averaged fitting goodness and the accepting ratio
  * [`sac_writer.h`](include/sac_writer.h)/[`sac_writer.cpp`](src/sac_writer.cpp)
    \- writer module to output the simulation information, e.g. the log info, the recovered spectrum and the quality report, to files
  * [`random.h`](include/random.h)/[`random.cpp`](src/random.cpp)
    \- independent random module for the Monte Carlo process
  * [`utils/`](include/utils) - Utilities
    * [`linear_algebra.hpp`](include/utils/linear_algebra.hpp) 
      \- Eigen-style interfaces for the high-efficiency diagonalization of matrices, using MKL LAPACK as the backend
    * [`toml.hpp`](include/utils/toml.hpp)
      \- [toml++ library](https://github.com/marzer/tomlplusplus), a header-only TOML config parser for C++ 
  * [`sac_main.cpp`](src/sac_main.cpp) - the main program

* [`benchmark/`](benchmark) - benchmark example
  * [`config.toml`](benchmark/config.toml)
    \- toml configuration file including the SAC parameters
  * [`tgrids.in`](benchmark/tgrids.in) - input QMC imaginary-time grids
  * [`corr.in`](benchmark/corr.in) - input QMC correlation functions
  * [`log.out`](benchmark/log.out)
    \- log file generated during the simulated annealing process of SAC. The metadata, including the sampling temperature, fitting goodness and the averaged accepting ratio, are recorded for the tracking of the annealing and equilibrium.
  * [`spec.out`](benchmark/spec.out)
    \- data of the accumulated spectral functions recovered by SAC
  * [`report.out`](benchmark/report.out) 
    \- quality report of the recovered spectral functions, compared with the input QMC data
  * [`stdout.out`](benchmark/stdout.out) - log file of the standard screen output during the running of the program, which is generated by slurm
  * [`scripts/`](benchmark/scripts) - python scripts
    * [`generate_data.py`](benchmark/scripts/generate_data.py)
      \- generation of the synthetic spectral functions and the biased correlation functions for the usage of benchmark
    * [`plot_benchmark.py`](benchmark/scripts/plot_benchmark.py)
      \- visualization of benchmark results
    * [`data/exact_spec.dat`](benchmark/scripts/data/exact_spec.dat)
      \- synthetic spectral functions
    * [`figure/benchmark.pdf`](benchmark/scripts/figure/benchmark.pdf)
      \- benchmark figure
* [`build/`](build) - building directory
  * [`runcmake.sh`](build/runcmake.sh)- bash script for the configuration of cmake
* [`run/`](run) - running directory
  * [`make.sh`](run/make.sh) - bash script for the compilation of the project
  * [`run.sh`](run/run.sh)
    \- bash script for the running of project with customized program options
  * [`submit.slurm`](run/submit.slurm) 
    \- slurm script for the submission of program tasks
* [`cmake/`](cmake) - FindXXX.cmake files for cross-platform compilation


## References ##
<span id="reference"></span>

1. Anders W. Sandvik, 
   Stochastic method for analytic continuation of quantum Monte Carlo data, 
   *Phys. Rev. B 57, 10287*, 1998. [doi](https://doi.org/10.1103/PhysRevB.57.10287)
2. K. S. D. Beach, 
   Identifying the maximum entropy method as a special limit of stochastic analytic continuation,
   *arXiv:cond-mat/040305*, 2004. [doi](http://arxiv.org/abs/cond-mat/0403055)  
3. Hui Shao, Yan Qi Qin, Sylvain Capponi et al.,
   Nearly deconfined spinon excitations in the square-lattice spin-1/2 Heisenberg antiferromagnet,
   *Phys. Rev. X 7, 041072*, 2017. [doi](https://link.aps.org/doi/10.1103/PhysRevX.7.041072)


## License & Support ##

The code is open source under the GPL-3.0 License.

I would like to sincerely thank Prof. Hui Shao and Zheng Zhou for the generous sharing of their source codes in Fortran.

If any questions, feel free to contact me via email 17307110117@fudan.edu.cn.
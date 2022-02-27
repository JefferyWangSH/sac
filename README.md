# Stochastic Analytic Continuation
![workflow](https://github.com/JefferyWangSH/sac/actions/workflows/main.yml/badge.svg?branch=master)
![commit](https://img.shields.io/github/last-commit/JefferyWangSH/sac/master?color=blue)

A general C++ implementation of Stochastic Analytic Continuation method `(SAC)`,
proposed by Anders W. Sandvik, is presented in this repository. 

The program aims to extract information of fermion spectral functions from imaginary-time correlations calculated by quantum Monte Carlo `(QMC)` simulations. 
Details about the `SAC` algorithm can be found in [References](#references) section. 

---

## Installation ##

### Prerequisite ###
* `gcc/g++` `( version >= 7.1, support C++17 standard )` and `cmake` `( version >= 3.21 )` installed.
* `Boost C++ libraries` `( version >= 1.71 )` installed.
* `Eigen library` `( version >= 3.3.9 )` providing a user-friendly interface of matrices.
* `Intel Math Kernel Library (MKL)` for high-accuracy linear algebra, alse working as the backend of `Eigen`.

### Usage ###
1. Download source codes from github:
    ``` shell
    $ # download the source code
    $ git clone https://github.com/JefferyWangSH/sac.git {program_root}
    ```

2. Enter directory [`build`](build/) and run [`runcmake.sh`](build/runcmake.sh) for the analysis of program using cmake.
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
   We use `slurm` for task management, and some sample scripts for program running are presented in [`submit.slurm`](run/submit.slurm) and [`run.sh`](run/run.sh). One can also run the code directly from the command line.
   ```shell
   $ # run the program using slurm
   $ cd {program_root}/run
   $ # edit simulating parameters before running the code 
   $ vim submit.slurm
   $ sh submit.slurm
   
   $ # run directly from command line
   $ # using program option `--help` to see optional params and helping messages.
   $ {program_root}/build/sac  --help
   ``` 

5. 
   To make the input QMC data detectable for the program, one should pass the path of the input corrlation function and the corresponding imaginary-time grids through command line options. 
   
   It is recommanded that all these input files should be put in certain user customized input folder, say "benchmark", under the [`input/`](input/) directory. 
   An example of the input format can be found in [`input/benchmark`](input/benchmark/), and the corresponding format of options are shown below.
    ```shell
    $ # specify paths of input and output files
    $ {program_root}/build/sac --tau-file-path=${tau_file_path}   \
                               --corr-file-path=${corr_file_path} \
                               --log-file-path=${log_file_path}   \
                               --spec-file-path=${spec_file_path} \
                               --report-file-path=${report_file_path}
    ```

6. 
   During the simulation, the progress rate of annealing can be checked in the output log file, which is assigned through option `--log-file-path`.
   And after the simulation is finished, results of the recovered spectrum `(--spec-file-path)`, together with a quality report of the recovery`(--report-file-path)`, will be automatically generated.
   ```shell
   $ # by default, the output files would be under output/benchmark unless explicitly assigned
   $ cd {program_root}/output/benchmark
   
   $ # monitor the process of annealing in log.log
   $ vim log.log

   $ # if someone runs the program using slurm, the screen output will get stored in out.log
   $ vim out.log
   
   $ # recovered spectral function
   $ vim spec.dat

   $ # quality report of recovered spectrum to quantify the goodness of the SAC recovery
   $ vim report.dat
   ```


## Repository Structure ##
* [`include`](include) & [`src`](src) - header files and source codes 
  * [`sac.h`](include/sac.h)/[`sac.cpp`](src/sac.cpp) 
    \- main class for SAC calculations, implementing MC sampling, simulated annealing, and the recovery of spectrum
  * [`qmc_data_reader.h`](include/qmc_data_reader.h)/[`qmc_data_reader.cpp`](src/qmc_data_reader.cpp)
    \- interface module for analysing and preprocessing the input correlation functions from QMC
  * [`freq_grids.h`](include/freq_grids.h)/[`freq_grids.cpp`](src/freq_grids.cpp)
    \- class for the discretization of grids in frequency domain
  * [`kernel.h`](include/kernel.h)/[`kernel.cpp`](src/kernel.cpp)
    \- integral kernels which connects spectral function and QMC correlations
  * [`annealing_chain.h`](include/annealing_chain.h)/[`annealing_chain.cpp`](src/annealing_chain.cpp)
    \- chain structure for the storage of metadata during the SAC annealing process
  * [`measure.h`](include/measure.h)/[`measure.cpp`](src/measure.cpp)
    \- class for the measurements of averaged fitting goodness and accepting radio
  * [`random.h`](include/random.h)/[`random.cpp`](src/random.cpp)
    \- independent random module for the Monte Carlo process
  * [`matrix_util.hpp`](src/matrix_util.hpp) 
    \- Eigen-style C++ interface for high-efficiency matrix diagonalization, using MKL LAPACK as backend
  * [`main.cpp`](src/main.cpp) 
    \- the main program
* [`input`](input) - input of SAC simulations
  * [`benchmark`](input/benchmark) - example input for benchmark
      * [`tau.dat`](input/benchmark/tau.dat) 
        \- imaginary-time grids of input correlation functions
      * [`cor.dat`](input/benchmark/cor.dat) 
        \- correlation functions measured by QMC
* [`output`](output) - simulation outputs
  * [`benchmark`](output/benchmark) - results of benchmark calculations
    * [`out.log`](output/benchmark/out.log)
      \- file records of the standard screen output when running code with slurm 
    * [`log.log`](output/benchmark/log.log) 
      \- log file generated during the simulated annealing process of SAC. 
      Sampling temperature, fitting goodness and average accepting radio are recorded for the tracking of annealing and equilibrium.
    * [`spec.dat`](output/benchmark/spec.dat) 
      \- data of accumulated spectral functions recovered by SAC
    * [`report.dat`](output/benchmark/report.dat)
      \- quality report of recovered spectral functions, compared with input data
* [`build`](build) - building directory
  * [`runcmake.sh`](build/runcmake.sh)- bash script for the configuration of cmake
* [`run`](run) - running directory
  * [`make.sh`](run/make.sh) - bash script for the compilation of the project
  * [`run.sh`](run/run.sh) 
    \- bash script for the running of project with customized program options
  * [`submit.slurm`](run/submit.slurm) 
    \- sbatch script for the submittment of program tasks using slurm
* [`script`](script) - python scripts 
  * [`generate_data.py`](script/generate_data.py)
    \- generation of synthetic spectral functions and biased correlation functions for benchmark use
  * [`plot_benchmark.py`](script/plot_benchmark.py)
    \- visualization of benchmark results
  * [`data/exact_spec.dat`](script/data/exact_spec.dat)
    \- synthetic spectral functions
  * [`figure/benchmark.pdf`](script/figure/benchmark.pdf)
    \- benchmark figure
* [`cmake`](cmake) - FindXXX.cmake files for cross-platform compilation


## References ##
<span id="reference"></span>

1. Anders W. Sandvik, 
   Stochastic method for analytic continuation of quantum Monte, 
   *Phys. Rev. B 57, 10287*, 1998. [doi](https://doi.org/10.1103/PhysRevB.57.10287)
2. K. S. D. Beach, 
   Identifying the maximum entropy method as a special limit of stochastic analytic continuation,
   *arXiv:cond-mat/040305*, 2004. [doi](http://arxiv.org/abs/cond-mat/0403055)  
3. Hui Shao, Yan Qi Qin, Sylvain Capponi et al.,
   Nearly Deconfined Spinon Excitations in the Square-Lattice Spin-1/2 Heisenberg Antiferromagnet,
   *Phys. Rev. X 7, 041072*, 2017. [doi](https://link.aps.org/doi/10.1103/PhysRevX.7.041072)


## License & Support ##

The code is free software under GPL-3.0 License.

I would like to sincerely thank Prof. Hui Shao and Zheng Zhou for the generous sharing of their source codes in Fortran.

If any questions, please feel free to contact me via email 17307110117@fudan.edu.cn.
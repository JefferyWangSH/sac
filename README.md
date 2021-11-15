# Stochastic Analytic Continuation [![forthebadge](https://forthebadge.com/images/badges/works-on-my-machine.svg)](https://forthebadge.com)


A primary C++ implementation of Stochastic Analytic Continuation method `(SAC)`,
proposed by Anders W. Sandvik, is presented in this repository. 

The program aims to extract information of fermion spectrum functions from imaginary-time quantum Monte Carlo `(QMC)` data. 
More details about the `SAC` algorithm can be found in [References](#references) section. 

---

## Installation ##

### Prerequisite ###
* `gcc/g++` `( version > 4.8.5 )` and `cmake` `( version > 2.8.12 )` installed.
* `Boost C++ libraries` `( version > 1.71 )` installed.
* `Eigen library` `( version > 3.3.9 )` providing a user-friendly interface of matrices.
* `Intel Math Kernel Library (MKL)` for high-accuracy linear algebra, alse working as the backend of `Eigen`.
* `Intel implementation of Message Passing Interface (Intel MPI) (optional)` for large range of distributed parallel.

### Usage ###
1. Download source codes from github:
    ``` shell
    $ # download the source code
    $ git clone https://github.com/JefferyWangSH/SAC.git {PROGRAM_ROOT}
    ```

2. Enter directory [`build/`](build/) and and run [`runcmake.sh`](build/runcmake.sh) for the analysis of program using cmake.
    ``` shell
    $ # initialize cmake
    $ cd {PROGRAM_ROOT}/build && ./runcmake.sh
    ```
    
3. Enter the [`run/`](run/) directory and compile the codes. 
    ``` shell
    $ # build the program
    $ cd {PROGRAM_ROOT}/run && ./make.sh
    ```

4. We use `slurm` for task management, and some sample scripts for program running are presented in [`batch.sh`](run/batch.sh) and [`run.sh`](run/run.sh). Of course, one can also run the code directly from the command line.
   ```shell
   $ # run the program using slurm
   $ cd {PROGRAM_ROOT}/run
   $ ./batch.sh
   $ # edit simulating parameters in run.sh
   $ vim run.sh
   
   $ # run directly from command line
   $ # using program option `--help` to see optional params and helping messages.
   $ cd {PROGRAM_ROOT}/build
   $ ./sac  --help --input-folder-name {INPUT_FOLDER_NAME}
   ``` 

5. To organize the input QMC data and make them detectable for the program, one should create a customized input folder and put it under the [`input/`](input/) directory. 
   
   The input should include at least two files named `tau.dat` and `cor.dat`, which contains the time sequence and correlation functions respectively.
   Some examples can be found in [`input/benchmark/`](input/benchmark/) to show the format of input. 

6. During the simulating process, an output folder, with the same name as the original input, will be created under the [`results/`](results/) directory, and after the simulation is finished, results of recovered spectrum will be automatically generated there.
   ```shell
   $ cd {PROGRAM_ROOT}/results/{INPUT_FOLDER_NAME}
   
   $ # monitor the process of annealing in log.log
   $ vim log.log
   
   $ # if someone runs the program using slurm, the screen output will get stored in out.out
   $ vim out.out
   
   $ # recovered spectrum
   $ vim spec.dat
   ```


## Repository Structure ##
* [`include/`](include) & [`src/`](src) - header file declarations and source codes 
  * [`SAC.h`](include/SAC.h)/[`SAC.cpp`](src/SAC.cpp) 
    \- main class for SAC calculations, implementing MC sampling, simulated annealing, and spectrum recovery
  * [`ReadInModule.h`](include/ReadInModule.h)/[`ReadInModule.cpp`](src/ReadInModule.cpp)
    \- interface module for analysis and pretreatments of input QMC data
  * [`FrequencyGrid.h`](include/FrequencyGrid.h)/[`FrequencyGrid.cpp`](src/FrequencyGrid.cpp)
    \- helping class for the division of grids in frequency domain, both grids for MC sampling and spectrum collection are supported
  * [`Kernel.h`](include/Kernel.h)/[`Kernel.cpp`](src/Kernel.cpp)
    \- integral kernels connecting spectrum function to QMC correlations
  * [`AnnealChain.h`](include/AnnealChain.h)/[`AnnealChain.cpp`](src/AnnealChain.cpp)
    \- structure and class to store simulating information during the SAC annealing process
  * [`Measure.h`](include/Measure.h)/[`Measure.cpp`](src/Measure.cpp)
    \- class for the measurements of average fitting goodness and accepting radio
  * [`Random.h`](include/Random.h)/[`Random.cpp`](src/Random.cpp)
    \- random module for the Monte Carlo process
  * [`MatrixUtil.hpp`](src/MatrixUtil.hpp) 
    \- C++ and Eigen style interface file for high-efficiency matrix diagonalization, using MKL LAPACK as backend
  * [`main.cpp`](src/main.cpp) 
    \- the main program
* [`input/`](input) - input files for SAC simulations
  * [`benchmark/`](input/benchmark) - benchmark input for the correctness testing of codes
      * [`tau.dat`](input/benchmark/tau.dat) 
        \- sequence of imaginary time points
      * [`cor.dat`](input/benchmark/cor.dat) 
        \- synthetic correlation functions with Gaussian noise introduced
      * [`spec_base.dat`](input/benchmark/spec_base.dat)
        \- standard spectral function to be benchmarked
* [`results/`](results) - simulation results
  * [`benchmark/`](results/benchmark) - results of benchmark calculation
    * [`out.out`](results/benchmark/out.out)
      \- file records of the screen output when running the code using slurm 
    * [`log.log`](results/benchmark/log.log) 
      \- log file generated during the simulated annealing process of SAC. 
      Sampling temperature, fitting goodness and average accepting radio are recorded for tracking of annealing and equilibrium.
    * [`spec.dat`](results/benchmark/spec.dat) 
      \- spectrum data generated by SAC
    * [`benchmark.pdf`](results/benchmark/benchmark.pdf) - figure of test spectrum
  * [`script.py`](results/script.py)
      \- python script to plot the recovered spectrum
* [`build/`](build) - building directory
  * [`runcmake.sh`](build/runcmake.sh)- bash script for cmake building 
* [`run/`](run) - running directory
  * [`make.sh`](run/make.sh) - bash script for program building
  * [`run.sh`](run/run.sh) 
    \- bash script for program running with customized program options
  * [`batch.sh`](run/batch.sh) 
    \- sbatch script for upload program tasks using slurm
* [`references/`](references) - references and some reading materials about SAC algorithms


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
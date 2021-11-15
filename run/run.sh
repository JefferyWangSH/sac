#!/bin/bash

#SBATCH --partition=v6_384
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

exe="../build/sac"
input_name="benchmark"
lt=150
beta=6.0
nbin_qmc=1000
theta=1e8

mpirun ${exe} --lt=${lt} --beta=${beta} --nbin-qmc=${nbin_qmc} --theta=${theta} --input-folder-name=${input_name}

exit 0
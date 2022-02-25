#!/bin/bash

# submitting params
partition="v6_384"
nodes=1
ntasks_per_node=1
cpus_per_task=8
mem_per_cpu=4G

# program params
exe="../build/sac"
input_name="L8b8.00u-4.00"
lt=160
beta=8.0
nbin_qmc=1000
theta=1e6

# create output folder if not exist
output_folder="../output/"$input_name
if [ ! -d $output_folder ]; then
    mkdir -p $output_folder
fi

# set up jobname and log output name
jobname=$input_name
output=$output_folder"/out.log"
error=$out_folder"/err.log"

# submit mission
sbatch --job-name=$jobname --output=$output --partition=$partition \
       --nodes=$nodes --ntasks-per-node=$ntasks_per_node --cpus-per-task=$cpus_per_task --mem-per-cpu=$mem_per_cpu \
       --export=exe=$exe,input_name=$input_name,lt=$lt,beta=$beta,nbin_qmc=$nbin_qmc,theta=$theta ./run.sh

exit 0
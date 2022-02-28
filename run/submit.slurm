#!/bin/bash

########################################################### submitting params
partition="v6_384"
nodes=1
ntasks_per_node=1
cpus_per_task=8
mem_per_cpu=4G

########################################################### program params
exe="../build/sac"

# input info of qmc data
lt=160
beta=8.0
nbin_qmc=1000
rebin_pace=1
nbootstrap=5000

kernel_type="fermion"
update_type="single"

# grids params
freq_min=-10.0
freq_max=+10.0
freq_interval=0.00001
spec_interval=0.01000

# annealing params
ndelta=1000
theta=1000000
annealing_pace=0.90
stablization_pace=10

# sac measurements and colloection of spetral functions
bin_size=4000
bin_num=5
collecting_steps=100000

# file names
folder_name="benchmark"
output_folder="../output/"$folder_name
tau_file_path="../input/"$folder_name"/tau.dat"
corr_file_path="../input/"$folder_name"/cor.dat"
log_file_path="../output/"$folder_name"/log.log"
spec_file_path="../output/"$folder_name"/spec.dat"
report_file_path="../output/"$folder_name"/report.dat"


########################################################### preprocessing
# check input files
if [ ! -f $tau_file_path ]; then
    echo "intput file "$tau_file_path" not exist. "
    exit 1
fi
if [ ! -f $corr_file_path ]; then
    echo "intput file "$corr_file_path" not exist. "
    exit 1
fi
# create output folder if not exist
if [ ! -d $output_folder ]; then
    mkdir -p $output_folder
fi
# delete previous logging file if exist
if [ -f $log_file_path ]; then
    rm $log_file_path && touch $log_file_path
fi

# set up jobname and log output name
jobname=$folder_name
output=$output_folder"/out.log"
error=$out_folder"/err.log"


########################################################### submit mission
sbatch --job-name=$jobname --output=$output --partition=$partition \
--nodes=$nodes --ntasks-per-node=$ntasks_per_node --cpus-per-task=$cpus_per_task --mem-per-cpu=$mem_per_cpu \
--export=exe=$exe,kernel_type=$kernel_type,update_type=$update_type,\
tau_file_path=$tau_file_path,corr_file_path=$corr_file_path,\
log_file_path=$log_file_path,spec_file_path=$spec_file_path,report_file_path=$report_file_path,\
lt=$lt,beta=$beta,nbin_qmc=$nbin_qmc,rebin_pace=$rebin_pace,nbootstrap=$nbootstrap,\
freq_min=$freq_min,freq_max=$freq_max,freq_interval=$freq_interval,spec_interval=$spec_interval,\
ndelta=$ndelta,theta=$theta,annealing_pace=$annealing_pace,stablization_pace=$stablization_pace,\
bin_size=$bin_size,bin_num=$bin_num,collecting_steps=$collecting_steps ./run.sh

exit 0
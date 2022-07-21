#!/bin/bash


#######################################################################################
#                                Submitting params
#######################################################################################

partition="v6_384"
nodes=1
ntasks_per_node=1
cpus_per_task=8
mem_per_cpu=4G


#######################################################################################
#                                  Program params
#######################################################################################

exe="../build/sac"

# file names
folder_name="benchmark"
output_folder="../"$folder_name
config_file="../"$folder_name"/config.toml"
tgrids_file="../"$folder_name"/tgrids.in"
corr_file="../"$folder_name"/corr.in"
log_file="../"$folder_name"/log.out"
spec_file="../"$folder_name"/spec.out"
report_file="../"$folder_name"/report.out"


#######################################################################################
#                             Pre-check the input files
#######################################################################################

# check the QMC input files
if [ ! -f $tgrids_file ]; then
    echo "intput file "$tgrids_file" do not exist. "
    exit 1
fi
if [ ! -f $corr_file ]; then
    echo "intput file "$corr_file" do not exist. "
    exit 1
fi
# delete previous log file if exists
if [ -f $log_file ]; then
    rm $log_file && touch $log_file
fi

# set up jobname and log output name
jobname=$folder_name
output=$output_folder"/stdout.out"
error=$out_folder"/stderr.out"


#######################################################################################
#                                Submit the mission
#######################################################################################

sbatch --job-name=$jobname --output=$output --partition=$partition \
--nodes=$nodes --ntasks-per-node=$ntasks_per_node --cpus-per-task=$cpus_per_task --mem-per-cpu=$mem_per_cpu \
--export=exe=$exe,config_file=$config_file,tgrids_file=$tgrids_file,corr_file=$corr_file,\
log_file=$log_file,spec_file=$spec_file,report_file=$report_file ./run.sh

exit 0
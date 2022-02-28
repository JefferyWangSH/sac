#!/bin/bash

# execute program
$exe --kernel-type=$kernel_type --update-type=$update_type \
     --tau-file-path=$tau_file_path --corr-file-path=$corr_file_path \
     --log-file-path=$log_file_path --spec-file-path=$spec_file_path  --report-file-path=$report_file_path \
     --lt=$lt --beta=$beta --nbin-qmc=$nbin_qmc --rebin-pace=$rebin_pace --nbootstrap=$nbootstrap \
     --freq-min=$freq_min --freq-max=$freq_max --freq-interval=$freq_interval --spec-interval=$spec_interval --ndelta=$ndelta \
     --theta=$theta --annealing-pace=$annealing_pace --stablization-pace=$stablization_pace \
     --sbin-sac=$bin_size --nbin-sac=$bin_num --collecting-steps=$collecting_steps

exit 0
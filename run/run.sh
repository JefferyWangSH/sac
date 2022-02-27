#!/bin/bash

# execute program
$exe --kernel-type=$kernel_type --update-type=$update_type \
     --tau-file-path=$tau_file_path --corr-file-path=$corr_file_path \
     --log-file-path=$log_file_path --spec-file-path=$spec_file_path  --report-file-path=$report_file_path \
     --lt=$lt --beta=$beta --nbin-qmc=$nbin_qmc --theta=$theta

exit 0
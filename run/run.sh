#!/bin/bash

# execute the program
$exe --config=$config_file --tgrids=$tgrids_file --corr=$corr_file \
     --log=$log_file --spec=$spec_file  --report=$report_file

exit 0
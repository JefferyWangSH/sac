#!/bin/bash

# name of folder containing input data
input_name="benchmark"

# create output folder if not exist
output_folder="../results/"${input_name}
if [ ! -d ${output_folder} ]; then
  mkdir ${output_folder}
fi

# set up jobname and log output name
jobname=${input_name}
output=${output_folder}"/out.out"

sbatch --job-name ${jobname} --output ${output} ./run.sh

exit 0
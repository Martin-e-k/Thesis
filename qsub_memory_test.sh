#!/bin/sh
### Every line in this header section should start with ### for a comment
### or #PBS for an option for qsub
### Note: No unix commands may be executed until after the last #PBS line
###
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
###
### Send mail when job is aborted or terminates abnormally
#PBS -M s193125@student.dtu.dk
#PBS -m abe
###
### Compute resources, here x core on x node
#PBS -l nodes=1:ppn=40:fatnode
###
### Required RAM in GB
#PBS -l mem=1200GB
###
### How long (max) will the job take
#PBS -l walltime=00:05:00
###

### Output files - not required to be specified
### Comment out the next 2 lines to use the job id instead in the file names
#PBS -e /home/projects/ku_00041/people/markri/qsub_logs/memory_test2.stderr
#PBS -o /home/projects/ku_00041/people/markri/qsub_logs/memory_test2.stdout
###
### Job name - not required to be specified
### It is often easier just to use the job id instead for recognition
#PBS -N memory_test
###
### More qsub options can be added here

# This part is the real job script
# Here follows the user commands:

# Load all required modules for the job

# Go to the directory to work from
cat /proc/sys/vm/max_map_count


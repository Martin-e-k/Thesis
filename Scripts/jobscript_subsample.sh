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
### Compute resources, here 1 core on 1 node
#PBS -l nodes=1:ppn=40:thinnode
###
### Required RAM in GBls

#PBS -l mem=120GB
###
### How long (max) will the job take
#PBS -l walltime=4:00:00
###

### Output files - not required to be specified
### Comment out the next 2 lines to use the job id instead in the file names
#PBS -e /home/projects/ku_00041/people/markri/qsub_logs/subsampling.err
#PBS -o /home/projects/ku_00041/people/markri/qsub_logs/subsampling.log
###
### Job name - not required to be specified
### It is often easier just to use the job id instead for recognition
#PBS -N Martin_subsampling
###
### More qsub options can be added here

# This part is the real job script
# Here follows the user commands:

# Load all required modules for the job
module load tools seqtk/1.3 parallel/20210722 

cd /home/projects/ku_00041/people/markri/data/02_trimmed/

ls /home/projects/ku_00041/people/markri/data/02_trimmed | grep 'trim.fq.gz$' > preprocessed.list

parallel -j 38 "seqtk sample {} 10000 | gzip - -c > /home/projects/ku_00041/people/markri/data/subset/02_trimmed/{}" :::: preprocessed.list



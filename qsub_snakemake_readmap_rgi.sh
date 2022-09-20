#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N snakemake_readmap_tox2
### Output files (comment out the next 2 lines to get the job name used instead)
###PBS -e $PBS_JOBID.snakemake_readmap_tox2.err
###PBS -o $PBS_JOBID.snakemake_readmap_tox2.log
### Number of nodes
#PBS -l nodes=1:ppn=40:thinnode,walltime=10:00:00:00
### Add current shell environment to job (comment out if not needed)
###Only send mail when job is aborted or terminates abnormally
#PBS -m n

##PBS -t 1-7

#Load bash modules
source /home/projects/ku_00041/people/markri/scripts/bash_modules


snakemake --snakefile /home/projects/ku_00041/people/markri/scripts/readmap_vir.smk --resources mem_gb=5880 --cluster 'qsub -W group_list=ku_00041 -A ku_00041 -l nodes=1:ppn={threads},mem={resources.mem_gb}gb,walltime={resources.runtime} -V -d /home/projects/ku_00041/people/markri/scripts -e /home/projects/ku_00041/people/markri/qsub_logs/snakelogs/ -o /home/projects/ku_00041/people/markri/qsub_logs/snakelogs/' --jobs 100 --printshellcmds --latency-wait 60 

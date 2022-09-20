#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N qsub_diamond_bacmet
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.qsub_diamond_vfdb_bins.stderr
#PBS -o $PBS_JOBID.qsub_diamond_vfdb_bins.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes	
#PBS -l nodes=1:ppn=4:thinnode
### Memory
#PBS -l mem=16GB

### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=1:00:00:00

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/data/

module load tools diamond/2.0.13

#parallel qsub -F '{1}' /home/projects/ku_00041/people/markri/scripts/qsub_diamond_vir_bins.sh :::: /home/projects/ku_00041/people/markri/lists/proteins_bins.list

F1=$1

#Run commands
diamond blastp -k 1 -e 0.00001 -d /home/projects/dtu_00009/data/vfdb/VFDB_setA_pro.fas.gz -q /home/projects/ku_00041/people/markri/data/prodigal/bins/protein_files/"$F1" -o /home/projects/ku_00041/people/markri/data/vfdb/bins/"$F1"

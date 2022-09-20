#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N qsub_diamond_bacmet
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.qsub_diamond_vfdb_contig.stderr
#PBS -o $PBS_JOBID.qsub_diamond_vfdb_contig.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes	
#PBS -l nodes=1:ppn=40:thinnode
### Memory

### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=10:00:00:00

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/data/

module load tools diamond/2.0.13

#Run commands
diamond blastp -k 1 -e 0.00001 -d /home/projects/dtu_00009/data/vfdb/VFDB_setA_pro.fas.gz -q /home/projects/ku_00041/people/markri/data/mmseqs2/all_ORF_rep_aaseq.fasta -o /home/projects/ku_00041/people/markri/data/vfdb/contigs/results_tophits.tsv

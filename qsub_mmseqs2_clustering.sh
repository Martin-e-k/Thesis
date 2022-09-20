#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N qsub_mmseqs2_clustering
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.qsub_mmseqs2_clustering.stderr
#PBS -o $PBS_JOBID.qsub_mmseqs2_clustering.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40:thinnode
### Memory

### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=00:10:00:00

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/data/prodigal/contig/dna_files/new_headers/

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/projects/ku_00041/people/markri/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/projects/ku_00041/people/markri/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/projects/ku_00041/people/markri/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/projects/ku_00041/people/markri/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<


######COMMAND used on command line######
#qsub qsub_mmseqs2_clustering.sh

module load tools mmseqs2/release_12-113e3

#Creates temp file of all ORF
cat *.fna > temp_all_ORF.fasta

#Run commands
mmseqs easy-linclust --min-seq-id 0.95 -c 0.95 --threads 24 temp_all_ORF.fasta ../../../../mmseqs2/all_ORF ../../../../mmseqs2/tmp;

#Removes file of all ORF..
rm temp_all_ORF.fasta



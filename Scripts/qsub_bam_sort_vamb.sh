#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N qsub_bam_sort_vamb
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.qsub_bam_sort_vamb.stderr
#PBS -o $PBS_JOBID.qsub_bam_sort_vamb.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40:thinnode
### Memory
#PBS -l mem=80gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=10:00:00:00

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/data/vamb2/bam

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
#commandline runing qsub:
#parallel "qsub -F '{1}' /home/projects/ku_00041/people/markri/scripts/qsub_bam_sort_vamb.sh" :::: /home/projects/ku_00041/people/markri/lists/samplenames.list

#Run commands
sample=$1

conda activate vamb_env
module load tools
module load samtools/1.14

samtools sort --threads 8 /home/projects/ku_00041/people/markri/data/vamb2/bam/"$sample".bam -O BAM -o /home/projects/ku_00041/people/markri/data/vamb2/bam/sorted/"$sample".sorted.bam 

conda deactivate

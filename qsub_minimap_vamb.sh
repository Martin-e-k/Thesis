#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N qsub_minimap_vamb
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.qsub_minimap_vamb.stderr
#PBS -o $PBS_JOBID.qsub_minimap_vamb.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40:thinnode
### Memory
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=02:00:00:00

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/data/vamb2

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
#parallel "qsub -F '{1}' /home/projects/ku_00041/people/markri/scripts/qsub_minimap_vamb.sh" :::: /home/projects/ku_00041/people/markri/lists/samplenames.list

#Run commands

sample=$1

conda activate vamb_env
module load tools
module load samtools/1.14

minimap2 -I 128g -t 40 -N 50 -ax sr catalogue.mmi /home/projects/ku_00041/people/markri/data/02_trimmed/"$sample"R1_001.trim.fq.gz /home/projects/ku_00041/people/markri/data/02_trimmed/"$sample"R2_001.trim.fq.gz | samtools view -F 3584 -b --threads 8 > /home/projects/ku_00041/people/markri/data/vamb2/bam/"$sample".bam

conda deactivate

#minimap2 -I 64g -t 8 -N 50 -ax sr catalogue.mmi /home/projects/ku_00041/people/markri/data/02_trimmed/Nuuk_ID90_S1_StV23C_5_mid1_S0_L001_R1_001.trim.fq.gz /home/projects/ku_00041/people/markri/data/02_trimmed/Nuuk_ID90_S1_StV23C_5_mid1_S0_L001_R2_001.trim.fq.gz | samtools view -F 3584 -b --threads 8 > /home/projects/ku_00041/people/markri/data/vamb/bam/Nuuk_ID90_S1_StV23C_5_mid1_S0_L001_.bam

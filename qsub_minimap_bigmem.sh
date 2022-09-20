#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N qsub_minimap2
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.qsub_minimap2_bigmem.stderr
#PBS -o $PBS_JOBID.qsub_minimap2_bigmem.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40:fatnode
### Memory
#PBS -l mem=1200gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=10:00:00:00

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/data/mmseqs2/

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
#parallel "qsub -F '{1}' /home/projects/ku_00041/people/markri/scripts/qsub_minimap_contigs.sh" :::: /home/projects/ku_00041/people/markri/lists/samplenames.list

#Run commands

sample=$1

module load tools samtools/1.14 minimap2/2.24r1122

minimap2 -I 1000g -d all_ORF_rep_index.mmi all_ORF_rep_seq.fasta

#minimap2 -I 180g -t 40 -N 50 -ax sr all_ORF_rep_index.mmi /home/projects/ku_00041/people/markri/data/02_trimmed/"$sample"R1_001.trim.fq.gz /home/projects/ku_00041/people/markri/data/02_trimmed/"$sample"R2_001.trim.fq.gz | samtools view -F 3584 -b --threads 40 > /home/projects/ku_00041/people/markri/data/minimap/bam/"$sample".bam

minimap2 -I 1000g -t 40 -N 50 -ax sr all_ORF_rep_index.mmi /home/projects/ku_00041/people/markri/data/02_trimmed/Nuuk_ID94_S1_StV23C_35_mid5_S0_L001_R1_001.trim.fq.gz /home/projects/ku_00041/people/markri/data/02_trimmed/Nuuk_ID94_S1_StV23C_35_mid5_S0_L001_R2_001.trim.fq.gz | samtools view -F 3584 -b --threads 40 > /home/projects/ku_00041/people/markri/data/minimap/bam/Nuuk_ID94_S1_StV23C_35_mid5_S0_L001_.bam

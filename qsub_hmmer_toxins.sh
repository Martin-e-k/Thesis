#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N qsub_hmmer_toxinst40
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.hmmer_toxinst40.stderr
#PBS -o $PBS_JOBID.hmmer_toxinst40.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40:thinnode
### Memory

### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=10:00:00:00

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/data/

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

#Load modules
module load tools hmmer/3.2.1


#Running hmmer on non_redudant prodigal gene catalogue agianst custom mobile genetic element hmmer database 
hmmsearch -T 40 --cpu 40 --domtblout /home/projects/ku_00041/people/markri/data/toxins/contigs/hmmerT40.results /home/projects/ku_00041/people/markri/PathoFact/databases/toxins/combined_Toxin.hmm /home/projects/ku_00041/people/markri/data/mmseqs2/all_ORF_rep_aaseq.fasta

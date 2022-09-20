#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N qsub_rgi4
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.qsub_rgi4.stderr
#PBS -o $PBS_JOBID.qsub_rgi4.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40:thinnode
### Memory

### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=10:00:00:00

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/PathoFact/localDB/


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
#qsub /home/projects/ku_00041/people/markri/scripts/qsub_rgi.sh

conda activate rgi

rgi load --card_json /home/projects/ku_00041/people/markri/PathoFact/localDB/card.json --local

#Run commands

rgi main --input_sequence /home/projects/ku_00041/people/markri/data/mmseqs2/all_ORF_rep_aaseq.fasta --output_file /home/projects/ku_00041/people/markri/data/rgi/contigs/all_protein_rep_rgi --input_type protein --alignment_tool DIAMOND --num_threads 40 --local --clean


#Deactivate conda enviroment
conda deactivate


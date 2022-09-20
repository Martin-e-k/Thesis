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
### Required RAM in GB
#PBS -l mem=184GB
###
### How long (max) will the job take
#PBS -l walltime=20:00:00:00
###

### Output files - not required to be specified
### Comment out the next 2 lines to use the job id instead in the file names
#PBS -e /home/projects/ku_00041/people/markri/qsub_logs/patho_172.err
#PBS -o /home/projects/ku_00041/people/markri/qsub_logs/patho_172.log
###
### Job name - not required to be specified
### It is often easier just to use the job id instead for recognition
#PBS -N patho_172

###
### More qsub options can be added here

# This part is the real job script
# Here follows the user commands:

# Load all required modules for the job

module load tools signalp/5.0b

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/PathoFact

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

conda activate PathoFact


#Get input filenames from commandline

#Single node run
snakemake -s Snakefile --use-conda --reason --cores 30 -p --latency-wait 60



#Batch job
#snakemake -s Snakefile --configfile config.yaml \
#--use-conda --cores 20 --reason -p \
#--cluster "qsub -W group_list=ku_00041 -A ku_00041 -l walltime=48:00:00 -l nodes=1:ppn=40:thinnode -l mem=180gb" \
#--latency-wait 60 --rerun-incomplete



#Deactivate conda enviroment
conda deactivate



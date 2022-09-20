#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N qsub_drep_compare
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.qsub_drep_compare.stderr
#PBS -o $PBS_JOBID.qsub_drep_compare.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40:thinnode
### Memory
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=30:00:00:00
##contact incase of event

# Go to the directory from where the job was submitted (initial directory is $HOME)
#echo Working directory is "$PBS_O_WORKDIR"
#cd "$PBS_O_WORKDIR"

# or define specific path
cd /home/projects/ku_00041/people/markri/data/drep/


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

# commandline command for running this script
# qsub -F '/home/projects/ku_00041/people/markri/data/drep/compare/' /home/projects/ku_00041/people/markri/scripts/qsub_drep_compare.sh

# activate the correct conda environment
conda activate drep

#code for running dRep
outputDir=$1
#inputFiles=$2

dRep compare "$outputDir" -g /home/projects/ku_00041/people/markri/data/vamb/vamb_output/bins/* -d 

conda deactivate


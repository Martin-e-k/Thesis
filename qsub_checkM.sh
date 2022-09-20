#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N qsub_checkM
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.qsub_checkM.stderr
#PBS -o $PBS_JOBID.qsub_checkM.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40:thinnode
### Memory
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=10:00:00:00

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/data/checkM


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

#Run commands
# commandline for running this script: 
# qsub -F '/home/projects/ku_00041/people/markri/data/checkM/result_file /home/projects/ku_00041/people/markri/data/vamb2/vamb_output/bins /home/projects/ku_00041/people/markri/data/checkM' /home/projects/ku_00041/people/markri/scripts/qsub_checkM.sh

conda activate vamb3

results_file=$1
input_bins=$2
output_dir=$3

checkm lineage_wf -t 32 -f "$results_file" --tab_table "$input_bins" "$output_dir"  

conda deactivate

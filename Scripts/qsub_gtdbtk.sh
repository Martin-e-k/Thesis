#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N GTDBtk
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.GTDBtk_vamb.stderr
#PBS -o $PBS_JOBID.GTDBtk_vamb.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
### Number of nodes
#PBS -l nodes=1:ppn=20:fatnode
### Memory
#PBS -l mem=1000gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=30:00:00:00

# Go to the directory from where the job was submitted (initial directory is $HOME)
cd /home/projects/ku_00041/people/markri/qsub_logs/
 
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

#commandline code to runs this script: 
#qsub -F '/home/projects/ku_00041/people/markri/data/vamb/vamb_output/bins /home/projects/ku_00041/people/markri/data/gtdbtk_vamb' /home/projects/ku_00041/people/markri/scripts/qsub_gtdbtk.sh

### Here follows the user commands:
conda activate gtdbtk

INPUT=$1
OUT_DIR=$2

gtdbtk classify_wf --genome_dir "$INPUT" --out_dir "$OUT_DIR" --cpus 20 --force --debug

conda deactivate

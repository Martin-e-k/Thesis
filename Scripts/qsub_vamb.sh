#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N vamb_binning
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.vamb2_binning_module.stderr
#PBS -o $PBS_JOBID.vamb2_binning_module.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:gpus=1:ppn=40
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=30:00:00:00

# Go to the directory from where the job was submitted (initial directory is $HOME)
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

#commandline code to runs this script: 
#qsub -F '/home/projects/ku_00041/people/markri/data/vamb2/vamb_output /home/projects/ku_00041/people/markri/data/vamb2/catalogue.fna.gz /home/projects/ku_00041/people/markri/data/vamb2/jgi_depthMatrix' /home/projects/ku_00041/people/markri/scripts/qsub_vamb.sh

### Here follows the user commands:
conda activate vamb3
module load tools
module load vamb/3.0.3

JGI=$3
FASTA=$2
OUTDIR=$1

vamb -o C --cuda --outdir "$OUTDIR" --fasta "$FASTA" --jgi "$JGI"  --minfasta 200000

#conda deactivate

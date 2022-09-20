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
### Compute resources, here x core on x node
#PBS -l nodes=1:ppn=40:fatnode
###
### Required RAM in GB
#PBS -l mem=1200GB
###
### How long (max) will the job take
#PBS -l walltime=20:00:00:00
###

### Output files - not required to be specified
### Comment out the next 2 lines to use the job id instead in the file names
#PBS -e $PBS_JOBID.spades.stderr
#PBS -o $PBS_JOBID.spades.stdout
###
### Job name - not required to be specified
### It is often easier just to use the job id instead for recognition
#PBS -N spades_two
###
### More qsub options can be added here

# This part is the real job script
# Here follows the user commands:

# Load all required modules for the job

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/data/02_trimmed

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
#parallel --xapply "qsub -F '{1} {2} {3}' /home/projects/ku_00041/people/markri/scripts/jobscript_spades.sh " :::: /home/projects/ku_00041/people/markri/lists/r1.list :::: /home/projects/ku_00041/people/markri/lists/r2.list :::: /home/projects/ku_00041/people/markri/lists/singles.list

conda activate assembly

#Get input filenames from commandline
R1=$1
R2=$2
S1=$3

#Run commands 
output_dir="${R1}"

spades.py --meta -k 21,33,55,77 -t 8 -o /home/projects/ku_00041/people/markri/data/05_assembly3/"$output_dir" -1 "$R1" -2 "$R2" -s "$S1"
#Deactivate conda enviroment
conda deactivate



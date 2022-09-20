#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N qsub_hmm_tox_sort
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.hmm_tox_sort.stderr
#PBS -o $PBS_JOBID.hmm_tox_sort.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=1:thinnode
### Memory
#PBS -l mem=8GB

### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 30 days)
#PBS -l walltime=00:10:00

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/data/toxins/bins2

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

#parallel qsub -F '{1}' /home/projects/ku_00041/people/markri/scripts/sort_hmm_output_batch.sh :::: /home/projects/ku_00041/people/markri/lists/hmm_tox_files_test.list

F1=$1

#Command to sort out hits with e-value below 10e-5
head -n -10 "$F1" | tail -n+4 | sed 's/ \+ /\t/g' | awk -F"\t" '$6<10e-5' > /filtered/"$F1"
#command to retain top hits for ORF with multiple hits
sort -rnk6 /filtered/"$F1" | awk '!x[$1]++' > /filtered/sorted/"$F1"

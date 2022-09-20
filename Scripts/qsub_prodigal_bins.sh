#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00041 -A ku_00041
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N qsub_prodigal_bins
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e $PBS_JOBID.qsub_prodigal_bins.stderr
#PBS -o $PBS_JOBID.qsub_prodigal_bins.stdout
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=20:thinnode
### Memory

### Requesting time - format is <days>:<hours>:<minutes>:<seconds>
#PBS -l walltime=01:00:00:00

# Go to the directory to work from
cd /home/projects/ku_00041/people/markri/data/vamb/vamb_output/bins/


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
#parallel qsub -F '{1}' /home/projects/ku_00041/people/markri/scripts/qsub_prodigal_bins.sh :::: /home/projects/ku_00041/people/markri/lists/bins.list

conda activate prodigal


#Get input filenames from commandline
F1=$1

#create outfile name from input by removing everything after the underscores
outfile_gbk="${F1%%.*}"".gbk"

outfile_protein="${F1%%.*}"".faa"

outfile_dna="${F1%%.*}"".fna"


#Run commands
prodigal -p meta -i "$F1" -o ../../../prodigal/bins/gbk_files/"$outfile_gbk" -a ../../../prodigal/bins/protein_files/"$outfile_protein" -d ../../../prodigal/bins/dna_files/"$outfile_dna"

#Deactivate conda enviroment
conda deactivate

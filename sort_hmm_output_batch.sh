#!/bin/sh

#parallel -F '{1}' /home/projects/ku_00041/people/markri/scripts/sort_hmm_output_batch.sh :::: /home/projects/ku_00041/people/markri/lists/hmm_tox_files_test.list

for FILE in *.results; 
do touch filtered/$FILE;
touch filtered/sorted/$FILE;
head -n -10 $FILE | tail -n+4 | sed 's/ \+ /\t/g' | awk -F"\t" '$6<10e-5' > filtered/$FILE;
sort -rgk6 filtered/$FILE | awk '!x[$1]++' > filtered/sorted/$FILE;
echo $FILE; done

F1=$1

#Command to sort out hits with e-value below 10e-5
head -n -10 "$F1" | tail -n+4 | sed 's/ \+ /\t/g' | awk -F"\t" '$6<10e-5' > /filtered/"$F1"
#command to retain top hits for ORF with multiple hits
sort -rnk6 /filtered/"$F1" | awk '!x[$1]++' > /filtered/sorted/"$F1"


for FILE in *.results; 
do touch best_hits/$FILE;
head -n -10 $FILE | tail -n+4 | awk '!x[$1]++' | sed 's/ \+ /\t/g' > best_hits/$FILE;
echo $FILE; done



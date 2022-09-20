#!/usr/bin/python3

#import libraries
import sys
import subprocess
import os

#Give paths on commandline
#file_dir = sys.argv[1]
#outfilename = sys.argv[2]

#Hard coded paths
file_list = open("/home/projects/ku_00041/people/markri/lists/hmm_outputfiles.list", "r")
outfilename = "/home/projects/ku_00041/people/markri/data/toxins/bins_combined2.txt"

#files = os.listdir(file_dir)

#Remove newlines
converted_list = []
for element in file_list:
    converted_list.append(element.strip())

#Write outputfile
outfile = open(outfilename, 'w')

#write_header = True

for element in converted_list:
    print(element)
    infile = open(element, 'r')
    genome = element.split('.')[0]
    #first_line = True
    
    for line in infile:
    #    if first_line != True:
            outfile.write(genome + "\t" + line)
    #    else:
    #        first_line = False
    #        if write_header == True:
    #            outfile.write("Bin\t" + line)
    #            write_header = False
    
    infile.close()

outfile.close()

#clean up temporary directory
#subprocess.run(['rm', '-rf', file_dir])

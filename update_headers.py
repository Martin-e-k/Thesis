import os
os.chdir("/home/projects/ku_00041/people/markri/data/prodigal/contig/protein_files")

prodigal_files = os.listdir()


#Update all headers to include the sample ID

for file in prodigal_files:
    if os.path.isfile(file):
        print(file)
        infile = open(file, "r")
        outfile = open("/home/projects/ku_00041/people/markri/data/prodigal/contig/protein_files/new_headers/" + file, "w+")
        sample_ID = file.split("MG_")[1]
        sample_ID = sample_ID.split(".")[0]
        for line in infile:
            if line.startswith(">"):
                new_header = line.split(" #")[0] + "_" + sample_ID + line[line.find(" #"):]
                outfile.write(new_header)
                continue
            outfile.write(line)
        infile.close()
        outfile.close()



outfile = open("/home/projects/ku_00041/people/markri/data/quast_output/assembly_summary.tsv", "w+")
list_file = open("/home/projects/ku_00041/people/markri/lists/contigs.list", "r")
first_file = True
list_ID = list_file.readlines()
for ID in list_ID:
	if first_file:
		#Handles first ID by adding in the header information
		infile = open("/home/projects/ku_00041/people/markri/data/quast_output/" + ID.strip() + "/transposed_report.tsv", "r")
		first_line = infile.readline().replace("Assembly\t", "")

		values = infile.readline().replace("contigs\t", "")

		outfile.write("ID" + "\t" + first_line)
		outfile.write(ID.strip() + "\t" + values)

		first_file = False
		infile.close()
		continue

	#Adds the values to the outputfile
	infile = open("/home/projects/ku_00041/people/markri/data/quast_output/" + ID.strip() + "/transposed_report.tsv", "r")
	first_line = infile.readline()

	values = infile.readline().replace("contigs\t", "")

	outfile.write(ID.strip() + "\t" + values)

	infile.close()

outfile.close()


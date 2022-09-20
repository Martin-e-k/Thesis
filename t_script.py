import sys
import csv

#Function to sum to lists
def add_counts(first_list, second_list):
	first_list = [int(i) for i in first_list]
	second_list = [int(i) for i in second_list]
	summed_counts = list(map(sum, zip(first_list, second_list)))

	summed_counts = [str(i) for i in summed_counts]
	return summed_counts

def sum_list(input_list):
	input_list = [int(i) for i in input_list]
	return sum(input_list)

gene_count_matrix = csv.reader(open(sys.argv[1],"r"), delimiter="\t")
metadata = open(sys.argv[2], "r")
outfile_matrix = open(sys.argv[3], "w")

#dict_prod_count = {prod_ID : countlist}
dict_prod_count = {}

first_line = True
for row in gene_count_matrix:
	if first_line:
		outfile_matrix.write("\t".join(row) + "\n")
		first_line = False
		continue

	dict_prod_count[row[0]] = row[1:]


#creat dict of rgi_ID and prodigal ID, structure is: {rgi_7 : [node1, node4, node3], rgi_2 : [node2]}
dict_rgi_prod_ID = {}
first_line = True
for line in metadata:

	if first_line:
		first_line = False
		continue

	prodigal_ID = line.split(";")[0]
	rgi_ID = line.split(";")[1]

	if rgi_ID not in dict_rgi_prod_ID:
		dict_rgi_prod_ID[rgi_ID] = [prodigal_ID]
	else:
		dict_rgi_prod_ID[rgi_ID] = dict_rgi_prod_ID[rgi_ID] + [prodigal_ID]

for key in dict_rgi_prod_ID:
	#print(key)
	#print(dict_rgi_prod_ID[key])
	new_readcount_row = dict_rgi_prod_ID
	if len(dict_rgi_prod_ID[key])>1:
		first_prod_ID = True 
		count_sum_all_sampels = []
		for prod_id in dict_rgi_prod_ID[key]:
			#get toal readcounts across sampels
			count_sum_all_sampels.append(sum_list(dict_prod_count[prod_id]))

			if first_prod_ID:
				new_readcount_row = dict_prod_count[prod_id]
				first_prod_ID = False
			else:				
				new_readcount_row = add_counts(new_readcount_row, dict_prod_count[prod_id])

		#Get prodigal ID with highest count across sampels.
		winning_prod_ID = dict_rgi_prod_ID[key][count_sum_all_sampels.index(max(count_sum_all_sampels))]
		try:
			outfile_matrix.write(winning_prod_ID + "\t" + "\t".join(new_readcount_row) + "\n")
		except:
			print("failed at winning_prod_ID")
			print(type(winning_prod_ID))
			print(key)
			print()
			sys.exit()

	else:


		prodigal_ID_to_write = dict_rgi_prod_ID[key][0]
		readcount_to_write = "\t".join(dict_prod_count[dict_rgi_prod_ID[key][0]])
		try:
			outfile_matrix.write(prodigal_ID_to_write + "\t" + readcount_to_write + "\n")
		except:
			print(prodigal_ID_to_write)
			print(readcount_to_write)
			print(type(prodigal_ID_to_write))
			print(type(readcount_to_write))
			print(dict_rgi_prod_ID[key])
			print(key)
			sys.exit()

		#outfile_matrix.write(dict_rgi_prod_ID[key][0] + "\t" + "\t".join(dict_prod_count[dict_rgi_prod_ID[key][0]] + "\n"))

print("DONE!")
outfile_matrix.close()
#gene_count_matrix.close()
metadata.close()

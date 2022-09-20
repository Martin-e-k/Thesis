library(readr)

read_counts <- read_delim("/home/projects/ku_00041/people/markri/data/readmap/virulence/gene_cat/nonredundant_gene_clusters/gene_count_matrix.tsv", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

read_counts[is.na(read_counts)] <- 0

#write.table(read_counts, file='test.tsv', quote=FALSE, sep='\t')
write_delim(read_counts, file='test.tsv', quote=NULL, delim='\t')

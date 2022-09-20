library(DESeq2)
library(tidyverse)
library(readxl)
library(patchwork)
library(reshape2)

#Set directories
work_dir = "C:/Users/marti/Desktop/Speciale/scripts/"
data_dir = "C:/Users/marti/Desktop/Speciale/data/"
result_dir = "C:/Users/marti/Desktop/Speciale/results/"
setwd(work_dir)

set.seed(28031991)

###Chose input type###
gene_type <- "amr" #amr, vir or tox#

#Load data
read_counts <- read.table(paste(data_dir, "gene_count_matrix_rgi.tsv", sep=""), header = TRUE)
gene_meta <- read_delim("C:/Users/marti/Desktop/Speciale/data/rgi_gene_meta.tsv", delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
gene_meta$ORF_ID <-word(gene_meta$ORF_ID)
gene_meta <- gene_meta %>% 
  select(ORF_ID, `Drug Class`)

metadata <- read_excel(paste(data_dir, "metadata.xlsx", sep = ""))
samples <- read_excel(paste(data_dir, "projekt_samples.xlsx", sep = ""), col_names = FALSE)


#Fix sample names in count matrix
read_count_matrix <- as.matrix(read_counts[ , -1])
colnames(read_count_matrix) <- str_sub(colnames(read_count_matrix), end = -10) 
colnames(read_count_matrix) <- str_replace(colnames(read_count_matrix), pattern = "_\\d+\\_\\d\\d", replacement = "_")
row.names(read_count_matrix) <- read_counts[ , 1]

#Change 'Greenland' location in metadata
metadata$Location <- str_replace(metadata$Location, pattern = "Greenland", replacement = "Qajaa")


rm(read_counts)
meta_location = metadata[1:2]

#Check forkerte navne i count matrice!!!!
ncol(read_count_matrix)
nrow(metadata)
colnames(read_count_matrix)[!(colnames(read_count_matrix) %in% metadata$SampleName)]
setdiff(metadata$SampleName, colnames(read_count_matrix))

#Sort read counts as the metadata
read_count_matrix <- read_count_matrix[, meta_location$SampleName]

gene_melt <- melt(read_count_matrix)
colnames(gene_melt) <- c("ORF_ID", "SampleName", "Count")
as_tibble(gene_melt)
gene_melt <- inner_join(gene_melt, meta_location, by = "SampleName")

next_tibl <- gene_melt %>% 
  group_by(Location, ORF_ID) %>% 
  summarise_at(vars(Count), list(Abundance = sum)) %>% 
  pivot_wider(names_from = Location, values_from = Abundance)

nexter_tibl <- inner_join(next_tibl, gene_meta, by = "ORF_ID")

#Get all drug classes
drug_classes <- gene_meta$`Drug Class` %>% 
  str_split("; ") %>% 
  unlist() %>% unique()

for(drug in drug_classes){
  nexter_tibl[, paste(drug, "_Nuuk", sep = "")] <- grepl(drug, nexter_tibl$`Drug Class`, fixed = TRUE) * nexter_tibl$Nuuk
  nexter_tibl[, paste(drug, "_Qajaa", sep = "")] <- grepl(drug, nexter_tibl$`Drug Class`, fixed = TRUE) * nexter_tibl$Qajaa
  nexter_tibl[, paste(drug, "_Sermermiut", sep = "")] <- grepl(drug, nexter_tibl$`Drug Class`, fixed = TRUE) * nexter_tibl$Sermermiut
}

nextest_tbl <- nexter_tibl %>% 
  select(-ORF_ID, -Nuuk, -Qajaa, -Sermermiut, -`Drug Class`) %>% 
  colSums()

chunklength=3

# split the vector by length

det <- split(nextest_tbl,ceiling(seq_along(nextest_tbl) / chunklength))

Nuuk <- c()
Qajaa <- c()
Sermermiut <- c()
for(f in det){
  Nuuk <- c(Nuuk, f[[1]]/58)
  Qajaa <- c(Qajaa, f[[2]]/7)
  Sermermiut <- c(Sermermiut, f[[3]]/12)
}
poppoo <- tibble(drug_classes, Nuuk, Qajaa, Sermermiut)
poppoo <- as_tibble(melt(poppoo))
colnames(poppoo) <- c("Drug_class", "Location", "Counts")

ggplot(data = poppoo, aes(x = Location, y = Drug_class, fill = Counts)) +
  geom_tile()

ggplot(data = poppoo, aes(x = Counts, y = Drug_class, fill = Location)) +
  geom_col()

abundant_drug_classes <- poppoo %>% 
  filter(Drug_class %in% c("tetracycline antibiotic", "triclosan", "rifamycin antibiotic", "phenicol antibiotic", "peptide antibiotic", "penem", "penam", 
                           "macrolide antibiotic", "glycylcycline", "glycopeptide antibiotic", "fosfomycin", "fluoroquinolone antibiotic",
                           "elfamycin antibiotic", "diaminopyrimidine antibiotic", "cephamycin", "cephalosporin", "carbapenem", "aminoglycoside antibiotic"
                           , "aminocoumarin antibiotic", "acridine dye"))

ggplot(data = abundant_drug_classes, aes(x = Counts, y = Drug_class, fill = Location)) +
  geom_col() +
  scale_fill_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Geographic location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  labs(y = "Drug class",
       x = "Normalized counts")

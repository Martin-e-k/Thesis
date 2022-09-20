library(tidyverse)
library(patchwork)

work_dir = "C:/Users/marti/Desktop/Speciale/scripts/"
data_dir = "C:/Users/marti/Desktop/Speciale/data/"
result_dir = "C:/Users/marti/Desktop/Speciale/results/"
setwd(work_dir)

#Load data
checkM_result <- read_delim("C:/Users/marti/Desktop/Speciale/data/bins/checkM_result.tsv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

gtdbtk_ar122_summary <- read_delim("C:/Users/marti/Desktop/Speciale/data/bins/gtdbtk.ar122.summary.tsv", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

gtdbtk_bac120_summary <- read_delim("C:/Users/marti/Desktop/Speciale/data/bins/gtdbtk.bac120.summary.tsv", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)


#Load gene information
amr_genes <- read_delim("C:/Users/marti/Desktop/Speciale/data/bins/rgi_combined_bins.tsv", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)


tox_genes <- read_delim("C:/Users/marti/Desktop/Speciale/data/bins/tox_combined_bins.tsv", 
                                delim = "\t", escape_double = FALSE, 
                                col_names = FALSE, trim_ws = TRUE)

tox_genes <- subset(tox_genes, select = -c(X2, X10, X11, X23))


vir_genes <- read_delim("C:/Users/marti/Desktop/Speciale/data/bins/vir_combined_bins.tsv", 
                                delim = "\t", escape_double = FALSE, 
                                col_names = FALSE, trim_ws = TRUE)

#This is stupid
hmmsearch_header <- read_delim("C:/Users/marti/Desktop/Speciale/data/bins/hmmsearch_header.tsv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
header <- colnames(hmmsearch_header)

colnames(tox_genes) <- header

#How many genes are found in each group across all bins
length(amr_genes$Bin)
length(tox_genes$X1)
length(vir_genes$X1)


#Check if any bins are identified as both både archea og bacteria
sum(!(gtdbtk_bac120_summary$user_genome %in% gtdbtk_ar122_summary$user_genome)) 


#Create bin quality groups
checkM_result <- checkM_result %>% 
  mutate(Quality = case_when(Completeness >= 90 & Contamination < 5 ~ 'High',
                             Completeness >= 50 & Completeness < 90 & Contamination < 10 ~ 'Medium',
                             Completeness >= 90 & Contamination >= 5 & Contamination < 10 ~'Medium',
                             TRUE ~ 'Low'))

checkM_result$Quality <- factor(checkM_result$Quality, levels = c('Low', 'Medium', 'High'))

#Save to use in other script
write.csv(checkM_result, file = paste(data_dir, "checkM_tibble.csv", sep = ""))

p1 <- ggplot(data = checkM_result,
             aes(x = Completeness,
                 y = Contamination,
                 color = Quality)) +
  geom_point() +
  scale_color_manual(values = c("firebrick", "orange", "darkgreen"), 
                     name = "Quality", 
                     labels = c("Low", "Medium", "High")) +
  theme_minimal() +
  ggtitle("Completeness v Contamination")

p2 <- ggplot(data = checkM_result,
             aes(x = Quality,
                 fill = Quality)) +
  geom_bar() +
  scale_fill_manual(values = c("firebrick", "orange", "darkgreen"), 
                     name = "Quality", 
                     labels = c("Low", "Medium", "High")) +
  labs(y = "Bin count") +
  theme_minimal() +
  ggtitle("Quality distribution")

p1 + p2  + plot_annotation(tag_levels = "a")

#Combine some with something else


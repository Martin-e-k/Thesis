library(treeio)
library(tidyverse)
library(ggtree)
library(ggplot2)
library(reshape2)
library(ggstance)
library(ape)
library(ggnewscale)

#library(rphylopic) cool but not used

#Set directories
work_dir = "C:/Users/marti/Desktop/Speciale/scripts/"
data_dir = "C:/Users/marti/Desktop/Speciale/data/"
result_dir = "C:/Users/marti/Desktop/Speciale/results/"
setwd(work_dir)

#Load bin quality info
checkM_tibble <- read_csv("C:/Users/marti/Desktop/Speciale/data/checkM_tibble.csv")

checkM_tibble <- checkM_tibble %>% 
   select(user_genome = `Bin Id`, Quality) # %>% 
  # rename(user_genome = `Bin Id`)

#Load gene count data
amr_count <- read_table2("C:/Users/marti/Desktop/Speciale/data/amr_count_bin.txt")
amr_count <- amr_count %>% 
  rename(user_genome = Bin, amr = `1`)
tox_count <- read_table2("C:/Users/marti/Desktop/Speciale/data/tox_count_bin.txt", 
                             col_names = FALSE)
tox_count <- tox_count %>% 
  rename(user_genome = X2, tox = X1)
vir_count <- read_table2("C:/Users/marti/Desktop/Speciale/data/vf_counts.txt", 
                         col_names = FALSE)
vir_count <- vir_count %>% 
  rename(user_genome = X2, vir = X1)

vir_count$user_genome <- str_sub(vir_count$user_genome, start = 1, end = -5)


#Load Tree files
bac_tree <- read.tree(file = paste(data_dir, "bac120.user.tree", sep = ""))
arc_tree <- read.tree(file = paste(data_dir, "ar122.user.tree", sep = ""))

#Load summary files
bac_summary <- read_table2("C:/Users/marti/Desktop/Speciale/data/gtdbtk.bac120.summary.tsv")
bac_summary <- bac_summary %>% 
  select(user_genome, classification)


arc_summary <- read_table2("C:/Users/marti/Desktop/Speciale/data/gtdbtk.ar122.summary.tsv")
arc_summary <- arc_summary %>% 
  select(user_genome, classification)

#Sort according to tree file
bac_summary <- bac_summary[match(bac_tree$tip.label, bac_summary$user_genome), ]
arc_summary <- arc_summary[match(arc_tree$tip.label, arc_summary$user_genome), ]

#Split taxa
bac_summary <- bac_summary %>% 
  mutate(kingdom = str_sub(sapply(bac_summary$classification, function(x) strsplit(x, ".__")[[1]][2]), start = 1, end = -2),
               phylum = str_sub(sapply(bac_summary$classification, function(x) strsplit(x, ";.__")[[1]][2]), start = 1),
               class = str_sub(sapply(bac_summary$classification, function(x) strsplit(x, ";.__")[[1]][3]), start = 1),
               order = str_sub(sapply(bac_summary$classification, function(x) strsplit(x, ";.__")[[1]][4]), start = 1),
               family = str_sub(sapply(bac_summary$classification, function(x) strsplit(x, ";.__")[[1]][5]), start = 1),
               genus = str_sub(sapply(bac_summary$classification, function(x) strsplit(x, ";.__")[[1]][6]), start = 1),
               species = str_sub(sapply(bac_summary$classification, function(x) strsplit(x, ";.__")[[1]][7]), start = 1))

arc_summary <- arc_summary %>% 
  mutate(kingdom = str_sub(sapply(arc_summary$classification, function(x) strsplit(x, ".__")[[1]][2]), start = 1, end = -2),
         phylum = str_sub(sapply(arc_summary$classification, function(x) strsplit(x, ";.__")[[1]][2]), start = 1),
         class = str_sub(sapply(arc_summary$classification, function(x) strsplit(x, ";.__")[[1]][3]), start = 1),
         order = str_sub(sapply(arc_summary$classification, function(x) strsplit(x, ";.__")[[1]][4]), start = 1),
         family = str_sub(sapply(arc_summary$classification, function(x) strsplit(x, ";.__")[[1]][5]), start = 1),
         genus = str_sub(sapply(arc_summary$classification, function(x) strsplit(x, ";.__")[[1]][6]), start = 1),
         species = str_sub(sapply(arc_summary$classification, function(x) strsplit(x, ";.__")[[1]][7]), start = 1))



#Filter checkM output
bac_checkM <- checkM_tibble %>% 
  filter(user_genome %in% bac_summary$user_genome)

arc_checkM <- checkM_tibble %>% 
  filter(user_genome %in% arc_summary$user_genome)

bac_summary <- inner_join(bac_summary, bac_checkM, by = "user_genome")

arc_summary <- inner_join(arc_summary, arc_checkM, by = "user_genome")


#Group OTUs by quality
bac_high <- bac_summary %>% 
  filter(Quality == "High") %>% 
  select(user_genome)

bac_med <- bac_summary %>% 
  filter(Quality == "Medium") %>% 
  select(user_genome)

bac_low <- bac_summary %>% 
  filter(Quality == "Low") %>% 
  select(user_genome)

bac_qual_groups <- c(bac_low, bac_med, bac_high)
names(bac_qual_groups) <- c("Low", "Medium", "High")

arc_high <- arc_summary %>% 
  filter(Quality == "High") %>% 
  select(user_genome)

arc_med <- arc_summary %>% 
  filter(Quality == "Medium") %>% 
  select(user_genome)

arc_low <- arc_summary %>% 
  filter(Quality == "Low") %>% 
  select(user_genome)

arc_qual_groups <- c(arc_low, arc_med, arc_high)
names(arc_qual_groups) <- c("Low", "Medium", "High")



bac_tree <- groupOTU(bac_tree, bac_qual_groups, "Quality")

arc_tree <- groupOTU(arc_tree, arc_qual_groups, "Quality")


#Group otus by phylum
bac_phyla <- unique(bac_summary$phylum)
bac_groups  <- c()
for (phy in bac_phyla){
  temp <- filter(bac_summary, phylum == phy)[1]
  bac_groups <- c(bac_groups, temp)
}
names(bac_groups) <- bac_phyla


#Reduce the number of phyla
bac_groups$Firmicutes <- c(bac_groups$Firmicutes, bac_groups$Firmicutes_A, bac_groups$Firmicutes_B)
bac_groups$Firmicutes_A <- NULL
bac_groups$Firmicutes_B <- NULL
bac_groups$Desulfobacterota <- c(bac_groups$Desulfobacterota, bac_groups$Desulfobacterota_B, bac_groups$Desulfobacterota_E, bac_groups$Desulfobacterota_F)
bac_groups$Desulfobacterota_B <- NULL
bac_groups$Desulfobacterota_E <- NULL
bac_groups$Desulfobacterota_F <- NULL

bac_groups$Other <- c(bac_groups$Eisenbacteria, bac_groups$Planctomycetota, bac_groups$Dependentiae, bac_groups$Armatimonadota)
bac_groups$Eisenbacteria <- NULL
bac_groups$Planctomycetota <- NULL
bac_groups$Dependentiae <- NULL
bac_groups$Armatimonadota <- NULL


arc_phyla <- unique(arc_summary$phylum)
arc_groups  <- c()
for (phy in arc_phyla){
 temp <- filter(arc_summary, phylum == phy)[1]
 arc_groups <- c(arc_groups, temp)
}
names(arc_groups) <- arc_phyla

bac_tree <- groupOTU(bac_tree, bac_groups, "phylum")
arc_tree <- groupOTU(arc_tree, arc_groups, "phylum")


#Cleanup
#rm(temp, arc_checkM, bac_checkM, arc_low, arc_med, arc_high, bac_low, bac_med, bac_high, bac_groups, arc_groups, bac_qual_groups, arc_qual_groups)


#Plot quality trees

#bac_tree$tip.label <- bac_summary$species
# ggtree(bac_tree, layout = "circular") +
#   geom_tiplab(size = 1, aes(color=Quality), offset = 0.2) +
#   scale_color_manual(values = c("darkgreen", "firebrick", "orange")) +
#   new_scale_color() +
#   geom_tippoint(aes(color= phylum), size = 1)
# 



# ggtree(bac_tree, aes(color = phylum), layout = "circular") +
#   new_scale_color() +
#   geom_tippoint(aes(color = Quality), size = 1, alpha = 1) +
#   scale_color_manual(values = c("darkgreen", "firebrick", "orange"))

# 
# ggtree(bac_tree, layout = "circular") +
#   geom_tippoint(aes(color = Quality, shape = Quality), size = 2, alpha = .8) +
#   scale_color_manual(values = c("darkgreen", "firebrick", "orange"))

ggtree(bac_tree, layout = "circular", branch.length = "none") +
  geom_tippoint(aes(color = Quality, shape = Quality, alpha = Quality, size = Quality)) +
  scale_color_manual(values = c("darkgreen", "firebrick", "orange")) +
  scale_alpha_manual(values = c(0.9, 0.5, 0.5)) +
  scale_size_manual(values = c(2, 1.5, 1.5))


ggtree(arc_tree, layout = "circular", branch.length = "none") +
  geom_tippoint(aes(color = Quality, shape = Quality, alpha = Quality, size = Quality)) +
  scale_color_manual(values = c("darkgreen", "firebrick", "orange")) +
  scale_alpha_manual(values = c(0.9, 0.5, 0.5)) +
  scale_size_manual(values = c(3, 2, 2)) 


#Remove low quality bins
bac_subset <- drop.tip(bac_tree, tip = bac_low$user_genome)
arc_subset <- drop.tip(arc_tree, tip = arc_low$user_genome)

#Get gene counts for remaining bins
bac_tox_count <- tox_count %>% 
  filter(user_genome %in% bac_subset$tip.label)
bac_amr_count <- amr_count %>% 
  filter(user_genome %in% bac_subset$tip.label)
bac_vir_count <- vir_count %>% 
  filter(user_genome %in% bac_subset$tip.label)

#
bac_gene_counts <- inner_join(bac_tox_count, bac_vir_count, by = "user_genome")
bac_gene_counts <- left_join(bac_gene_counts, bac_amr_count, by = "user_genome")
bac_gene_counts$amr <- bac_gene_counts$amr %>%  replace_na(0)
bac_gene_counts <- bac_gene_counts %>% 
  relocate(user_genome)



arc_tox_count <- tox_count %>% 
  filter(user_genome %in% arc_subset$tip.label)
arc_amr_count <- amr_count %>% 
  filter(user_genome %in% arc_subset$tip.label)
arc_vir_count <- vir_count %>% 
  filter(user_genome %in% arc_subset$tip.label)



#
arc_gene_counts <- inner_join(arc_tox_count, arc_vir_count, by = "user_genome")
arc_gene_counts <- left_join(arc_gene_counts, arc_amr_count, by = "user_genome")
arc_gene_counts$amr <- arc_gene_counts$amr %>%  replace_na(0)
arc_gene_counts <- arc_gene_counts %>% 
  relocate(user_genome)

#Select pathogens
bac_pathogens <- bac_gene_counts %>% 
  filter(vir > quantile(vir)[4] & tox > quantile(tox)[4] )

bac_non_pathogens <- bac_gene_counts %>% 
  filter(!user_genome %in% bac_pathogens$user_genome)

patho_groups_bac <- c(list(bac_pathogens$user_genome), list(bac_non_pathogens$user_genome))
names(patho_groups_bac) <- c("pathogens", "non_pathogens")

arc_pathogens <- arc_gene_counts %>% 
  filter(vir > quantile(vir)[4] & tox > quantile(tox)[4] )

arc_non_pathogens <- arc_gene_counts %>% 
  filter(!user_genome %in% arc_pathogens$user_genome)

patho_groups_arc <- c(list(arc_pathogens$user_genome), list(arc_non_pathogens$user_genome))
names(patho_groups_arc) <- c("pathogens", "non_pathogens")

bac_subset <- groupOTU(bac_subset, patho_groups_bac, "Pathogens")
arc_subset <- groupOTU(arc_subset, patho_groups_arc, "Pathogens")


#Prepare variables for heatmaps
vir_heat_bac <- data.frame(row.names = bac_gene_counts$user_genome, "Vir_count" = bac_gene_counts$vir)
tox_heat_bac <- data.frame(row.names = bac_gene_counts$user_genome, "Tox_count" = bac_gene_counts$tox)
amr_heat_bac <- data.frame(row.names = bac_gene_counts$user_genome, "AMR_count" = bac_gene_counts$amr)

vir_heat_arc <- data.frame(row.names = arc_gene_counts$user_genome, "Vir_count" = arc_gene_counts$vir)
tox_heat_arc <- data.frame(row.names = arc_gene_counts$user_genome, "Tox_count" = arc_gene_counts$tox)
amr_heat_arc <- data.frame(row.names = arc_gene_counts$user_genome, "AMR_count" = arc_gene_counts$amr)


#Create treeplots with gene heat maps
p_bac <- ggtree(bac_subset, layout = "circular", branch.length = "none") %<+% bac_gene_counts +
  geom_tiplab(size = 1, aes(color = Pathogens))+
  scale_color_manual(values = c("black", "red"), guide = "none")

h1_bac <- p_bac + new_scale_fill()
h2_bac <- gheatmap(h1_bac, vir_heat_bac,
                   offset = 11,
                   width = 0.10,
                   colnames = FALSE) +
  scale_fill_continuous(name ="Virulence genes",
                        low = "yellow", high = "red")

h3_bac <- h2_bac + new_scale_fill()

h4_bac <- gheatmap(h3_bac, tox_heat_bac,
                   offset = 6,
                   width = 0.10,
                   colnames = FALSE) +
  scale_fill_continuous(name = "Toxin genes",
                        low = "turquoise1", high = "purple")

h5_bac <- h4_bac + new_scale_fill()
h6_bac <- gheatmap(h5_bac, amr_heat_bac,
                   offset = 1,
                   width = 0.10,
                   colnames = FALSE) +
  scale_fill_continuous(name = "AMR genes",
                        low = "seagreen1", high = "seagreen4")

h6_bac


p_arc <- ggtree(arc_subset, layout = "circular", branch.length = "none") %<+% arc_gene_counts +
  geom_tiplab(size = 3, aes(color = Pathogens))+
  scale_color_manual(values = c("black", "red"), guide = "none") 

h1_arc <- p_arc + new_scale_fill()
h2_arc <- gheatmap(h1_arc, vir_heat_arc,
                   offset = 6,
                   width = 0.10,
                   colnames = FALSE) +
  scale_fill_continuous(name ="Virulence genes",
                        low = "yellow", high = "red")

h3_arc <- h2_arc + new_scale_fill()

h4_arc <- gheatmap(h3_arc, tox_heat_arc,
                   offset = 5,
                   width = 0.10,
                   colnames = FALSE) +
  scale_fill_continuous(name = "Toxin genes",
                        low = "turquoise1", high = "purple")

h5_arc <- h4_arc + new_scale_fill()
h6_arc <- gheatmap(h5_arc, amr_heat_arc,
                   offset = 4,
                   width = 0.10,
                   colnames = FALSE) +
  scale_fill_continuous(name = "AMR genes",
                        low = "seagreen1", high = "seagreen4")

h6_arc











# IQR(bac_gene_counts$vir)
# median
# median(bac_gene_counts$vir)
# median(bac_gene_counts$amr)
# IQR(bac_gene_counts$amr)
# quantile(bac_gene_counts$amr)[4]
# 



#Trines kode
# e_tree_reroot <- root(curated_tree, outgroup="GB_GCA_009694025.1",resolve.root = FALSE)
# p_reroot <- ggtree(e_tree_reroot, layout="circular")#,branch.length="none")
# groups <- list(infant= e_tree_reroot$tip.label[(grepl("18097D", e_tree_reroot$tip.label))], adult= e_tree_reroot$tip.label[(grepl("plate", e_tree_reroot$tip.label))], outgroup=e_tree_reroot$tip.label[!(grepl("18097D", e_tree_reroot$tip.label) | grepl("plate", e_tree_reroot$tip.label))])
# groupOTU(p_reroot, groups, "Cohort") + geom_tippoint(aes(color=Cohort),size=1)
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

THE_tibl <- inner_join(checkM_tibble, vir_count, by = "user_genome")
THE_tibl <- left_join(THE_tibl, tox_count, by = "user_genome")
THE_tibl <- left_join(THE_tibl, amr_count, by = "user_genome")

THE_tibl <- THE_tibl %>% 
  mutate(classification = case_when(user_genome %in% arc_summary$user_genome ~ "Archaea",
                                    user_genome %in% bac_summary$user_genome ~ "Bacteria",
                                    TRUE ~ "Unclassified"))


THE_tibl$amr <- THE_tibl$amr %>% replace_na(0)
THE_tibl$tox <- THE_tibl$tox %>% replace_na(0)


THE_tibl %>% 
  filter(Quality != "Low", classification == "Unclassified") %>% 
  select(user_genome, vir) %>% 
  ggplot(aes(x = user_genome, y=vir)) +
  geom_col()

test <- THE_tibl %>% 
  filter(Quality != "Low")
  

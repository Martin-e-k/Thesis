library(tidyverse)
library(vegan)
library(ape)
library(phyloseq)
library(pairwiseAdonis)
library(microbiome)
library(stats)
#Set directories
work_dir = "C:/Users/marti/Desktop/Speciale/scripts/"
data_dir = "C:/Users/marti/Desktop/Speciale/data/"
result_dir = "C:/Users/marti/Desktop/Speciale/results/"
setwd(work_dir)

data = read.table(paste(data_dir, "unknown.virus.all_samples.tsv", sep = ""), sep = "\t", header = T, row.names = 1)
data$NCBI_tax_id = NULL
unknown <- data[1,]
data <- data[-1,]
colnames(unknown) <- str_sub(colnames(unknown), start = 15, end = -27)

# #Prepare some stuff to use later with phyloseq
taxa <- tibble(kingdom = str_sub(sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][1]), start = 4),
               phylum = str_sub(sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][2]), start = 4),
               class = str_sub(sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][3]), start = 4),
               order = str_sub(sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][4]), start = 4),
               family = str_sub(sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][5]), start = 4),
               genus = str_sub(sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][6]), start = 4),
               species = str_sub(sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][7]), start = 4))

# taxa <- tibble(kingdom = sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][1]),
#                phylum = sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][2]),
#                class = sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][3]),
#                order = sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][4]),
#                family = sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][5]),
#                genus = sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][6]),
#                species = sapply(rownames(data), function(x) strsplit(x, "\\|")[[1]][7]))
# 


colnames(data) <- str_sub(colnames(data), start = 15, end = -27)
otu_mat <-  as.matrix(data)
otu_mat = otu_table(otu_mat, taxa_are_rows = TRUE)

tax_mat <- as.matrix(taxa)
rownames(tax_mat) <- rownames(data)
tax_mat <- tax_table(tax_mat)


#Load metadata (fastqc output)
multiqc_general_stats <- read_table2("C:/Users/marti/Desktop/Speciale/data/multiqc_general_stats.txt")

multiqc_general_stats <- multiqc_general_stats[grepl("R1", multiqc_general_stats$Sample), ]

multiqc_general_stats$Sample <- str_sub(multiqc_general_stats$Sample, start = 23, end = -16)
colnames(multiqc_general_stats) <- c("Sample", "Dublicate_%", "GC_content_%", "Avg_seq_length", "Fails_%", "Total_reads")

total_reads <- multiqc_general_stats %>% 
  select(Sample, Total_reads)

data = data[grepl("g__", row.names(data)) & !grepl("s__", row.names(data)),]

datadist = vegdist(t(data))
datadistfit = pcoa(datadist)


#Extract principal components
PCs <-datadistfit$vectors
#Convert to dataframe for manipulation
PCs_df <- as.data.frame(PCs)

#Compute varriance explained by each PC to be used as axis labels in plot
vars_transformed <- apply(PCs_df, 2, var)
var_explained_PC <- vars_transformed/sum(vars_transformed)

datajoint = data.frame(datadistfit$vectors[,1:2], Label = sub("unknown\\.virus\\.", "", sapply(row.names(datadistfit$vectors), function(x) strsplit(x, "_")[[1]][1])))
 
datajoint[3][datajoint[3] == "Mid"] <- "Qajaa"
datajoint[3][datajoint[3] == "Ser"] <- "Sermermiut"
pcoa_plot <- ggplot(datajoint, aes(x = Axis.1, y = Axis.2, color = Label )) + 
  geom_point() + stat_ellipse() +
  theme_minimal() +
  xlab(paste("PC1 [", as.character(round(var_explained_PC[1]*100, 2)),"%]", sep = "")) +
  ylab(paste("PC2 [", as.character(round(var_explained_PC[2]*100, 2)),"%]", sep = "")) +
  scale_color_manual(values = c("mediumaquamarine",  "coral2", "darkorchid3"), 
                    name = "Geographic location", 
                    labels = c( "Nuuk", "Qajaa", "Sermermiut")) 



#Adonis test
adonis_test <- adonis(formula = datadist  ~ Label, data = datajoint, permutations = 10000)

#p-value only
adonis_test$aov.tab$`Pr(>F)`[1]

#Pairwise
pair_adonis <- pairwise.adonis(datadist, factors = datajoint$Label, p.adjust.m = "holm", perm = 10000)

#Add p-value to plot

pcoa_plot + annotate("text", x=0.45, y=0.65, label = paste("p-value: 9.999e-05"))


#Alpha diversity
shannon = diversity(data, index = "shannon")

shannon_tibl = tibble(Sample = rownames(datajoint), shannon, datajoint[3])

shannon_tibl <- left_join(shannon_tibl, total_reads, by = "Sample")

shannon_library_lm <- shannon_tibl %>% lm(formula=shannon ~ Total_reads)

shannon_tibl$norm_shannon <- shannon_library_lm$residuals + shannon_library_lm$coefficients[1]


ggplot(shannon_tibl, aes(y = norm_shannon, x = Label, fill = Label)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0) +
  geom_jitter() +
  theme_minimal() +
  scale_fill_manual(values = c( "mediumaquamarine", "coral2", "darkorchid3"), 
                    name = "Geographic location",
                    labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ylab("Shannon index") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggplot(shannon_tibl, aes(y = norm_shannon, x = Label, fill = Label)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0) +
  geom_jitter() +
  theme_minimal() +
  scale_fill_manual(values = c( "mediumaquamarine", "coral2", "darkorchid3"), 
                    name = "Geographic location",
                    labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ylab("Shannon index") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggplot(shannon_tibl, aes(y = norm_shannon, x = Label, fill = Label)) +
  geom_violin(alpha = 0.7) +
  geom_jitter() +
  theme_minimal() +
  scale_fill_manual(values = c( "mediumaquamarine", "coral2", "darkorchid3"), 
                    name = "Geographic location",
                    labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ylab("Shannon index") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#Wilcoxon test between diversity indexies
wilcox <- pairwise.wilcox.test(shannon_tibl$norm_shannon, shannon_tibl$Label, p.adj = "none")
wilcox_adj <- pairwise.wilcox.test(shannon_tibl$norm_shannon, shannon_tibl$Label, p.adj = "holm")

wilcox$p.value



#Unknown DNA in each sample

shannon_tibl$unknown_dna <-t(unknown)

ggplot(shannon_tibl, aes(x = Sample, y = unknown_dna, color = Label)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("mediumaquamarine",  "coral2", "darkorchid3"), 
                     name = "Geographic location", 
                     labels = c( "Nuuk", "Qajaa", "Sermermiut")) +
  labs(y = "Percentage unknown DNA") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggplot(shannon_tibl, aes(x = Label, y = unknown_dna, fill = Label)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0) +
  geom_jitter() +
  theme_minimal() +
  scale_fill_manual(values = c("mediumaquamarine",  "coral2", "darkorchid3"), 
                     name = "Geographic location", 
                     labels = c( "Nuuk", "Qajaa", "Sermermiut")) +
  labs(y = "Percentage unknown DNA") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


wilcox_unknown <- pairwise.wilcox.test(shannon_tibl$unknown_dna, shannon_tibl$Label, p.adj = "none")
wilcox_unknown_adj <- pairwise.wilcox.test(shannon_tibl$unknown_dna, shannon_tibl$Label, p.adj = "holm")

wilcox_unknown$p.value
wilcox_unknown_adj$p.value


###Stacked bar plots###

#Create phyloseq

meta_data <- shannon_tibl[,-1]
rownames(meta_data) <- shannon_tibl$Sample
phylo_meta <- sample_data(meta_data)

sample_names(phylo_meta) <- sample_names(otu_mat)

phylo_object <- phyloseq(otu_mat, tax_mat, phylo_meta)

#Relative abundance on kingdom level SORTED!

physeq_relativ = transform_sample_counts(phylo_object, function(x) x / sum(x) )
glom <- tax_glom(physeq_relativ, taxrank = 'kingdom')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$kingdom <- as.character(data_glom$kingdom) #convert to character


clostridia <- data_glom[data_glom$kingdom == "Bacteria",]
lvls <- clostridia[order(clostridia$Abundance, decreasing = TRUE), c(2)]

ggplot(data_glom, aes(factor(Sample, levels=lvls), y = Abundance, fill = fct_reorder(kingdom, Abundance))) +
  facet_grid(~Label, scales = "free") +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_brewer("Kingdom", palette = "Dark2") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Relative abundance")

#Relative abundance phylum
physeq_relativ = transform_sample_counts(phylo_object, function(x) x / sum(x) )
glom <- tax_glom(physeq_relativ, taxrank = 'phylum')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$phylum <- as.character(data_glom$phylum) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$phylum[data_glom$Abundance < 0.1] <- "< 10% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$phylum))
Count

unique(data_glom$phylum)
unique(data_glom$class)


#data_glom$class <- factor(data_glom$class, levels = c("Alphaproteobacteria", "Betaproteobacteria", "Deltaproteobacteria", "Gammaproteobacteria", "Flavobacteriia", "Verrucomicrobiae", "Planctomycetia", "Saprospirae", "Sphingobacteriia", "BD1-5", "Cytophagia", "Acidobacteriia", "OM190", "TA18", "Unknown", "< 1% abund."))
data_glom$phylum <- factor(data_glom$phylum, levels = c("Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Deinococcus_Thermus", "Euryarchaeota", "Firmicutes", "Negarnaviricota", "Nitrospirae", "Planctomycetes", "Proteobacteria", "Verrucomicrobia", "Viruses_unclassified", "< 10% abund."))

#plot with condensed phyla into "unknown" category
spatial_plot <- ggplot(data=data_glom, aes(x=Sample, y=Abundance, fill=phylum)) + facet_grid(~Label, scales = "free")
spatial_plot + geom_bar(aes(), stat="identity", position="stack") +
  theme_classic() +
  scale_fill_brewer("Phylum", palette = "Dark2") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Relative abundance")



#Relative abundance bacteria only
bacteria <- subset_taxa(phylo_object, kingdom == "Bacteria")
bacteria_relativ = transform_sample_counts(bacteria, function(x) x / sum(x) )
glom_bac <- tax_glom(bacteria_relativ, taxrank = 'phylum')
glom_bac # should list # taxa as # phyla
data_glom_bac<- psmelt(glom_bac) # create dataframe from phyloseq object
data_glom_bac$phylum <- as.character(data_glom_bac$phylum) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom_bac$phylum[data_glom_bac$Abundance < 0.1] <- "< 10% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom_bac$phylum))
Count

unique(data_glom_bac$phylum)


#data_glom$class <- factor(data_glom$class, levels = c("Alphaproteobacteria", "Betaproteobacteria", "Deltaproteobacteria", "Gammaproteobacteria", "Flavobacteriia", "Verrucomicrobiae", "Planctomycetia", "Saprospirae", "Sphingobacteriia", "BD1-5", "Cytophagia", "Acidobacteriia", "OM190", "TA18", "Unknown", "< 1% abund."))
data_glom_bac$phylum <- factor(data_glom_bac$phylum, levels = c("Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Deinococcus_Thermus", "Euryarchaeota", "Firmicutes", "Negarnaviricota", "Nitrospirae", "Planctomycetes", "Proteobacteria", "Verrucomicrobia", "< 10% abund."))


#plot with condensed phyla into "unknown" category
spatial_plot_bac <- ggplot(data=data_glom_bac, aes(x=Sample, y=Abundance, fill=phylum)) + facet_grid(~Label, scales = "free")
spatial_plot_bac + geom_bar(aes(), stat="identity", position="stack") +
  theme_classic() +
  scale_fill_brewer("Phylum", palette = "Dark2") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Relative abundance")



#All phyla sorted by proteiobacteria abundance
physeq_relativ = transform_sample_counts(phylo_object, function(x) x / sum(x) )
glom <- tax_glom(physeq_relativ, taxrank = 'phylum')
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$phylum <- as.character(data_glom$phylum) #convert to character
data_glom$phylum[data_glom$Abundance < 0.1] <- "< 10% abund."


clostridia <- data_glom[data_glom$phylum == "Proteobacteria",]
lvls <- clostridia[order(clostridia$Abundance, decreasing = TRUE), c(2)]
remaining_lvls <- setdiff(unique(data_glom$Sample), lvls)
lvls <- c(lvls, remaining_lvls)

ggplot(data_glom, aes(factor(Sample, levels=lvls), y = Abundance, fill = fct_reorder(phylum, Abundance))) +
  facet_grid(~Label, scales = "free") +
  geom_bar(stat = "identity") +
  theme_classic() +
  #scale_fill_manual(values = c("red", "blue", "green", "yellow", "purple", "darkgreen")) +
  scale_fill_brewer("Phylum", palette="Dark2") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Relative abundance")


  
#Bacteria phylum sorted by Firmicutes abundance
physeq_relativ = transform_sample_counts(bacteria, function(x) x / sum(x) )
glom <- tax_glom(physeq_relativ, taxrank = 'phylum')
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$phylum <- as.character(data_glom$phylum) #convert to character
data_glom$phylum[data_glom$Abundance < 0.1] <- "< 10% abund."


clostridia <- data_glom[data_glom$phylum == "Firmicutes",]
lvls <- clostridia[order(clostridia$Abundance, decreasing = TRUE), c(2)]
remaining_lvls <- setdiff(unique(data_glom$Sample), lvls)
lvls <- c(lvls, %`)remaining_lvls)

ggplot(data_glom, aes(factor(Sample, levels=lvls), y = Abundance, fill = fct_reorder(phylum, Abundance))) +
  facet_grid(~Label, scales = "free") +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_brewer("Phylum", palette="Dark2") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Relative abundance")


#Archaea

ps_archaea <- subset_taxa(phylo_object, kingdom == "Archaea")
archaea_relativ = transform_sample_counts(ps_archaea, function(x) x / sum(x) )



#General data stats for rapport
multiqc_general_stats$Sequence_depth <- multiqc_general_stats$Avg_seq_length * multiqc_general_stats$Total_reads

mean(multiqc_general_stats[70:77,]$`Dublicate_%`)
multiqc_general_stats[70:77,]

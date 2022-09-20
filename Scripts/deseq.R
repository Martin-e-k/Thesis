library(DESeq2)
library(tidyverse)
library(readxl)
library(stringr)
library(pheatmap)
library(vegan)
library(ape)
library(patchwork)
library(readr)
library(pairwiseAdonis)


#Set directories
work_dir = "C:/Users/marti/Desktop/Speciale/scripts/"
data_dir = "C:/Users/marti/Desktop/Speciale/data/"
result_dir = "C:/Users/marti/Desktop/Speciale/results/"
setwd(work_dir)

set.seed(28031991)

###Chose input type###
gene_type <- "amr" #amr, vir or tox#

#Load data

if(gene_type == "tox"){
  read_counts <- read.table(paste(data_dir, "new_gene_counts_tox.tsv", sep=""), header = TRUE)
  gene_meta <- read_delim("C:/Users/marti/Desktop/Speciale/data/tox_gene_meta.csv", 
                          delim = ";", escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE)
  gene_meta <- gene_meta %>% 
    select(-X4)
  colnames(gene_meta) <- c("Accession", "ORF_ID", "Gene symbol", "Database", "Gene name")
  
} else if(gene_type == "amr"){
  read_counts <- read.table(paste(data_dir, "gene_count_matrix_rgi.tsv", sep=""), header = TRUE)
  gene_meta <- read_delim("C:/Users/marti/Desktop/Speciale/data/rgi_gene_meta.tsv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
  gene_meta$ORF_ID <-word(gene_meta$ORF_ID)
  gene_meta <- gene_meta %>% 
    select(ARO, ORF_ID, Best_Hit_ARO, `Drug Class`, `Resistance Mechanism`, `AMR Gene Family`)
  
  
} else if (gene_type == "vir"){
  read_counts <- read.table(paste(data_dir, "vfdb_gene_counts.tsv", sep=""), header = TRUE)
  vfdb_meta <- read_excel("C:/Users/marti/Desktop/Speciale/data/vfdb_meta.xls", skip = 1)
  
  vfdb_prodID_acc <- read_delim("C:/Users/marti/Desktop/Speciale/data/vfdb_prodID_acc.tsv", 
                                delim = "\t", escape_double = FALSE, 
                                col_names = c("ORF_ID", "VFGID"), trim_ws = TRUE)
  
  vfdb_fasta_headers <- read_csv("C:/Users/marti/Desktop/Speciale/data/vfdb_fasta_headers.txt", 
                                 col_names = FALSE)
  
  VFGID <- word(vfdb_fasta_headers$X1, 1)
  VFGID <- substring(VFGID, 2)
  
  VFID <- str_extract(vfdb_fasta_headers$X1, pattern = "VF\\d+")
  
  vfg_vfid_fromheaders <- tibble(VFGID, VFID)
  
  joined_meta <- left_join(vfg_vfid_fromheaders, vfdb_meta, by = "VFID")
  
  gene_meta <- left_join(vfdb_prodID_acc, joined_meta, by = "VFGID")
  rm(joined_meta, vfdb_fasta_headers, vfdb_meta, vfdb_prodID_acc, VFGID, VFID,vfg_vfid_fromheaders)
} else{
  stop(paste("Learn to type fool"))
}

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

colnames(read_count_matrix) == meta_location$SampleName


#Create DESeq2 object
dds = DESeqDataSetFromMatrix(countData = read_count_matrix,
                             colData = meta_location,
                             design = ~ Location)

#rm(read_count_matrix, metadata, df)
#Run deseq
dds <- DESeq(dds)

###Differential abundance###
res <- results(dds)
resultsNames(dds)

summary(res)
plotMA(res)

#Vulcano plot
res1 = as.data.frame(res)

res1 = mutate(res1, Significant=ifelse(res1$padj<0.1, "FDR<0.1", "Not significant"))

res1[which(abs(res1$log2FoldChange)<1.0), "Significant"] = "Not significant"

ggplot(res1, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(col=Significant)) +
  scale_color_manual(values = c("darkgreen", "red"))

#Inspect genes with highest fold change (Both directions)
res = res[order(res$log2FoldChange, decreasing = TRUE), ]

upGenes = rownames(head(res, 5))

res = res[order(res$log2FoldChange, decreasing = FALSE), ]
downGenes = rownames(head(res, 5))

outlier_fold <- rownames(head(res, 1))

#Plot normalized count of specific gene (left out of presentation)
plotCounts(dds, gene=which.min(res$padj), intgroup="Location")

#Plot of only the 100 most significant genes
sig100 <- res[order(res$padj),][1:100,]

plotMA(sig100)

res05 <- results(dds, alpha=0.05)
summary(res05)
plotMA(res05)


#### FIND UD AF HVILKEN 
#Count transformation strategier
#vsd <- vst(dds) #virker ikke med små matricer
ntd <- normTransform(dds)
rld <- rlog(dds) #Virker ikke med store matricer


###Tutorial plots###

head( assay(rld) )

#par( mfrow = c( 1, 2 ) )
plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )


sampleDists_rld <- dist( t( assay(rld) ) )

sampleDistMatrix <- as.matrix( sampleDists_rld )
rownames(sampleDistMatrix) <- paste( rld$Location)
colnames(sampleDistMatrix) <- NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

#PCA
plotPCA( rld, intgroup = c( "Location"))

#heatmap
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )

heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))



###Difference between samples###
#Create dendrogram
norm.data = assay(rld)
sampleDists_ntd <- dist(t(norm.data),  method = "euclidean")
sampleDists_bray <- vegdist(t(norm.data),  method = "bray")

clusters=hclust(sampleDists_ntd)
clusters_bray=hclust(sampleDists_bray)

plot(clusters, labels = FALSE)
plot(clusters_bray, label = FALSE)
#PCA plot
plotPCA(ntd, intgroup=c("Location"))
plotPCA(rld, intgroup=c("Location"))

#plotSparsity(dds)
#plot(hclust(sampleDists_bray, method = "average"), label = FALSE)


#Compare PLC between normalized (Not used for presentation)
ntd <- normTransform(dds)
plotPCA(ntd, intgroup=c("Location"))
plotPCA(rld, intgroup=c("Location"))

#Heatmap (Not used in presentation)
select <- order(rowMeans(counts(dds, normalized=TRUE)),
               decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,"Location"])
pheatmap(assay(dds)[select,], cluster_rows = TRUE, show_rownames = FALSE,
        cluster_cols = TRUE, show_colnames = FALSE)


#Heatmap of sample distance (Not used)
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists_ntd)
rownames(sampleDistMatrix) <- dds$Location
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
        clustering_distance_rows=sampleDists_ntd,
        clustering_distance_cols=sampleDists_ntd,
        col=colors)



###PCoA####
#Compute bray curtis distance of ARG count matrix
distance3 <- vegdist(x=t(norm.data), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)

#plot(hclust(distance3, method = "average"), labels = colnames(countData_rgi_2),
#     main="TPM transformed RGI read counts\ndistance: Bray-curtis ; linkage: average")

#dendogram_rgi <- hclust(distance3, method = "average")
#Write dendrogram to file
#write(hc2Newick(dendogram_rgi),file='dendogram_rgi_TPM.newick')


#Perform PCoA using ape
PCoA_ARG <- pcoa(distance3)

#Extract principal components
PCs <- PCoA_ARG$vectors
#Convert to dataframe for manipulation
PCs_df <- as.data.frame(PCs)


#Extract principal component 1 and 2 and convert to tibble structure 
PC1_2_df <- rownames_to_column(.data = PCs_df[ , c("Axis.1", "Axis.2")], var = "SampleName")

#Join PC1, PC2 with sample metadata
PCoA_df <- inner_join(x=PC1_2_df, y=meta_location, by="SampleName")

#Convert to tibble object
PCoA_tibble <- as_tibble(PCoA_df)
#Rename PC columns

PCoA_tibble <- PCoA_tibble %>%
  rename("PC1" = "Axis.1",
         "PC2" = "Axis.2")
# PCoA_tibble <- PCoA_tibble %>% 
#   rename("Axis.1" = "PC1", 
#          "Axis.2" = "PC2")


#Compute varriance explained by each PC to be used as axis labels in plot
vars_transformed <- apply(PCs_df, 2, var)
var_explained_PC <- vars_transformed/sum(vars_transformed)

#Plot the samples according to PC1 and PC2
PCOA_ARG_plot <- ggplot(data = PCoA_tibble,aes(x=PC1,
                                               y=PC2,
                                               color=Location)) +
  geom_point(size=2) +
  stat_ellipse(size=0.8, geom="polygon", aes(fill = Location), alpha = 0.1) + 
  theme_minimal() +
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Geographic location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  scale_fill_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                    name = "Geographic location", 
                    labels = c("Nuuk", "Qajaa", "Sermermiut")) + 
  xlab(paste("PC1 [", as.character(round(var_explained_PC[1]*100, 2)),"%]", sep = "")) +
  ylab(paste("PC2 [", as.character(round(var_explained_PC[2]*100, 2)),"%]", sep = ""))

PCOA_ARG_plot



#Adonis test
adonis_test <- adonis(formula = distance3  ~ Location, data = meta_location, permutations = 10000)

#p-value only
p_value <- adonis_test$aov.tab$`Pr(>F)`[1]

#Pairwise
pair_adonis <- pairwise.adonis(distance3, factors = meta_location$Location, p.adjust.m = "holm", perm = 10000)

PCOA_ARG_plot + annotate("text", x=0.25, y=0.25, label = paste("p-value: 9.999e-05"))

###Boxplots of significant genes###
#List of significant genes
sig_genes <- subset(res1, Significant == "FDR<0.1")

sig_genes <- sig_genes[order(sig_genes$padj, decreasing = FALSE), ]

#top_sig_genes <- head(sig_genes, 5)

temp = assays(dds)
counts = temp$counts
rm(temp)

sig_genes_counts <-  counts[row.names(counts) %in% row.names(sig_genes),]
#top_sig_genes_counts <- counts[row.names(counts) %in% row.names(top_sig_genes),]

#Plot normalized count of specific gene (left out of presentation)
plotCounts(dds, gene=which.min(res$padj), intgroup="Location")

#There is a cleaner way to do this
norm_df <- as.data.frame(t(norm.data))
tibl <- row.names(norm_df)

tibl2 <- as_tibble(norm_df)

tibl3 <- tibl2 %>% 
  add_column(SampleName = tibl, .before = 1)

data_tibl<- left_join(meta_location, tibl3, by = "SampleName")

rm(norm_df, tibl, tibl2, tibl3)

#Plot
gene1 <- rownames(sig_genes)[1]
gene_symbol1 <- filter(gene_meta, ORF_ID == gene1)[3]
if(gene_type == "vir"){
  gene_symbol1 <-  filter(gene_meta, ORF_ID == gene1)[4]
}
p1 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene1))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene1),
                  color = Location)) +
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  labs(x = "", y = "Normalized counts") +
  ggtitle(gene_symbol1)

gene2 = rownames(sig_genes)[2]
gene_symbol2 <- filter(gene_meta, ORF_ID == gene2)[3]
if(gene_type == "vir"){
  gene_symbol2 <-  filter(gene_meta, ORF_ID == gene2)[4]
}

p2 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene2))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene2),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  labs(x = "", y = "Normalized counts") +
  ggtitle(gene_symbol2)

gene3 = rownames(sig_genes)[3]
gene_symbol3 <- filter(gene_meta, ORF_ID == gene3)[3]
if(gene_type == "vir"){
  gene_symbol3 <- filter(gene_meta, ORF_ID == gene3)[4]
}
p3 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene3))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene3),
                  color = Location)) +
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  labs(x = "", y = "Normalized counts") +
  ggtitle(gene_symbol3)

gene4 = rownames(sig_genes)[4]
gene_symbol4 <- filter(gene_meta, ORF_ID == gene4)[3]
if(gene_type == "vir"){
  gene_symbol4 <-  filter(gene_meta, ORF_ID == gene4)[4]
}
p4 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene4))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene4),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ggtitle(gene_symbol4) +
  labs(x = "", y = "Normalized counts")


gene5 = rownames(sig_genes)[5]
gene_symbol5 <- filter(gene_meta, ORF_ID == gene5)[3]
if(gene_type == "vir"){
  gene_symbol5 <-  filter(gene_meta, ORF_ID == gene5)[4]
}
p5 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene5))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene5),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ggtitle(gene_symbol5) +
  labs(x = "", y = "Normalized counts")

gene6 = rownames(sig_genes)[6]
gene_symbol6 <- filter(gene_meta, ORF_ID == gene6)[3]
if(gene_type == "vir"){
  gene_symbol6 <-  filter(gene_meta, ORF_ID == gene6)[4]
}
p6 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene6))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene6),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ggtitle(gene_symbol6) +
  labs(x = "", y = "Normalized counts")

gene7 = rownames(sig_genes)[7]
gene_symbol7 <- filter(gene_meta, ORF_ID == gene7)[3]
if(gene_type == "vir"){
  gene_symbol7 <-  filter(gene_meta, ORF_ID == gene7)[4]
}
p7 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene7))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene7),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ggtitle(gene_symbol7) +
  labs(x = "", y = "Normalized counts")

gene8 = rownames(sig_genes)[8]
gene_symbol8 <- filter(gene_meta, ORF_ID == gene8)[3]
if(gene_type == "vir"){
  gene_symbol8 <-  filter(gene_meta, ORF_ID == gene8)[4]
}
p8 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene8))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene8),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ggtitle(gene_symbol8) +
  labs(x = "", y = "Normalized counts")

gene9 = rownames(sig_genes)[9]
gene_symbol9 <- filter(gene_meta, ORF_ID == gene9)[3]
if(gene_type == "vir"){
  gene_symbol9 <-  filter(gene_meta, ORF_ID == gene9)[4]
}
p9 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene9))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene9),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ggtitle(gene_symbol9) +
  labs(x = "", y = "Normalized counts")




#Plot med patchwork
(p1 + p2) / (p3 + p4) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = "a")

(p1 + p2 + p3) / (p4 + p5 + p6) / (p7+ p8+ p9) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = "a")



#Misc used for thesis

genes_to_look_up <- c(gene1, gene2, gene3, gene4, gene5, gene6, gene7, gene8, gene9)

look_up_table <- filter(gene_meta, ORF_ID %in% genes_to_look_up) %>% 
  select(-ORF_ID)


sig_genes
sig_genes_tibl <- tibble(ORF_ID = rownames(sig_genes), sig_genes)
test <- rownames(sig_genes)

test2 <- filter(gene_meta, ORF_ID %in% test) %>% 
  filter(Best_Hit_ARO != "adeF") %>% 
  select(ORF_ID)

sig_genes_tibl <- filter(sig_genes_tibl, ORF_ID %in% test2$ORF_ID)



gene1 <- sig_genes_tibl$ORF_ID[1]
gene_symbol1 <- filter(gene_meta, ORF_ID == gene1)[3]
p1 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene1))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene1),
                  color = Location)) +
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  labs(x = "", y = "Normalized counts") +
  ggtitle(gene_symbol1)

gene2 = sig_genes_tibl$ORF_ID[2]
#gene_symbol2 <- filter(gene_meta, ORF_ID == gene2)[3]
p2 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene2))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene2),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  labs(x = "", y = "Normalized counts") +
  ggtitle(gene_symbol2)

gene3 = sig_genes_tibl$ORF_ID[3]
gene_symbol3 <- filter(gene_meta, ORF_ID == gene3)[3]
p3 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene3))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene3),
                  color = Location)) +
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  labs(x = "", y = "Normalized counts") +
  ggtitle(gene_symbol3)

gene4 = sig_genes_tibl$ORF_ID[4]
gene_symbol4 <- filter(gene_meta, ORF_ID == gene4)[3]
p4 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene4))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene4),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ggtitle(gene_symbol4) +
  labs(x = "", y = "Normalized counts")


gene5 = sig_genes_tibl$ORF_ID[5]
gene_symbol5 <- filter(gene_meta, ORF_ID == gene5)[3]
p5 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene5))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene5),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ggtitle(gene_symbol5) +
  labs(x = "", y = "Normalized counts")

gene6 = sig_genes_tibl$ORF_ID[6]
gene_symbol6 <- filter(gene_meta, ORF_ID == gene6)[3]
p6 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene6))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene6),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ggtitle(gene_symbol6) +
  labs(x = "", y = "Normalized counts")

gene7 = sig_genes_tibl$ORF_ID[7]
gene_symbol7 <- filter(gene_meta, ORF_ID == gene7)[3]
p7 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene7))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene7),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ggtitle(gene_symbol7) +
  labs(x = "", y = "Normalized counts")

gene8 = sig_genes_tibl$ORF_ID[8]
gene_symbol8 <- filter(gene_meta, ORF_ID == gene8)[3]
p8 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene8))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene8),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ggtitle(gene_symbol8) +
  labs(x = "", y = "Normalized counts")

gene9 = sig_genes_tibl$ORF_ID[9]
gene_symbol9 <- filter(gene_meta, ORF_ID == gene9)[3]
p9 <- ggplot(data = data_tibl,
             aes(x = Location,
                 y = !!ensym(gene9))) +
  theme_minimal() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.25,
              aes(x = Location,
                  y = !!ensym(gene9),
                  color = Location)) + 
  scale_color_manual(values = c("mediumaquamarine", "coral2", "darkorchid3"), 
                     name = "Location", 
                     labels = c("Nuuk", "Qajaa", "Sermermiut")) +
  ggtitle(gene_symbol9) +
  labs(x = "", y = "Normalized counts")

(p1 + p2 + p3) / (p4 + p5 + p6) / (p7+ p8+ p9) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = "a")
#Write file of genenames
#lapply(row.names(sig_genes), write, paste(data_dir, "sig_genes_vir.txt", sep = ""), append=TRUE)

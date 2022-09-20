##############################################################################################

#This Snakefile produces a gene count matrix of the nonredundant genes for each cluster of antibiotic resistance genes


# The output is:
# /home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/gene_count_matrix.tsv: Gene count matrix of the readcounts genes in each sample


##############################################################################################

#!python

#!/usr/bin/env python3
import os
import sys
import pandas
import re
from itertools import product


#Get wildcards for each sample ID

SAMPLES, = glob_wildcards("/home/projects/ku_00041/people/markri/data/02_trimmed/{sample}R1_001.trim.fq.gz")


#Helper rule
rule all:
        input: "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/gene_count_matrix.tsv"

# Extract the DNA sequences of the hmm virulence protein hits
rule get_rep_sequences:
    input:
        vir_results = "/home/projects/ku_00041/people/markri/data/vfdb/contigs/results_tophits.tsv",
        prodigal_clusters = "/home/projects/ku_00041/people/markri/data/mmseqs2/all_ORF_rep_seq.fasta"
    output:
        vir_ID_list = "/home/projects/ku_00041/people/markri/data/vfdb/contigs/vir_ID.list",
        gene_cat = "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/non_redundant_vir_gene_catalogue.fasta"
    threads: 1
    resources:
        mem_gb = 20,
        runtime = 1800 #30min
    shell:
        "cut -f1 {input.vir_results} > {output.vir_ID_list}; seqkit grep -f {output.vir_ID_list} {input.prodigal_clusters} -o {output.gene_cat}"

# The reads are mapped to the nonredundant gene catalogue with minimap2
rule minimap_index:
    input:
        "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/non_redundant_vir_gene_catalogue.fasta"
    output:
        index = "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/nonredundant.mmi",
        gene_lengths = "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/all_genes_nonredundant.fasta.fai"
    threads: 1
    resources:
        mem_gb = 50,
        runtime = 86400
    shell:
        "minimap2 -d {output.index} {input}; touch {output.gene_lengths}"


rule minimap_align:
    input:
        index = "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/nonredundant.mmi",
        gene_cat = "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/non_redundant_vir_gene_catalogue.fasta",
        fastq1 = "/home/projects/ku_00041/people/markri/data/02_trimmed/{sample}R1_001.trim.fq.gz",
       	fastq2 = "/home/projects/ku_00041/people/markri/data/02_trimmed/{sample}R2_001.trim.fq.gz"
    output:
        bam = "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/mapped_reads/gene_counts_{sample}.bam"
    threads: 40
    resources:
        mem_gb = 188,
        runtime = 86400 #24h
    shell:
        "minimap2 -v 1 -t 40 -N 50 -ax sr {input.index} {input.fastq1} {input.fastq2} | samtools view -T {input.gene_cat} -F 3584 -b --threads {threads} | samtools sort --threads {threads} > {output.bam}; "


# Use samtools to count the genes
rule count_genes:
    input:
        "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/mapped_reads/gene_counts_{sample}.bam"
    output:
        count = "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/counts/{sample}.counts"
    threads: 1
    resources:
        mem_gb = 20,
        runtime = 7200 #2h
    shell:
        "samtools idxstats {input} | cut -f3 > {output.count}"

# Create the gene index
rule gene_names:
    input:
        "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/mapped_reads/gene_counts_{sample}.bam"
    output:
        "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/mapped_reads/gene_names_{sample}"
    wildcard_constraints:
        sample=SAMPLES[0]
    threads: 1
    resources:
        mem_gb = 20,
        runtime = 3600 #1h
    shell:
        "samtools idxstats {input} | cut -f1 > {output}"


# Create the header for the matrix
rule create_header:
    input: 
        expand("/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/counts/{sample}.counts", sample=SAMPLES)
    output:
        header = "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/mapped_reads/header.txt"
    threads: 1
    resources:
        mem_gb = 188,
        runtime = 3600 #1h
    run:
        header = "Gene"
        for f in input:
            sample_name = f.split("/")[-1].split(".")[0]
            header = header + "\t" + sample_name
        header = header + "\n"
        with open(output.header, "w") as out:
            out.write(header)


#Combining the gene names with the counts of all samples
rule gene_count_matrix:
    input:
        counts = expand("/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/counts/{sample}.counts", sample=SAMPLES),
        gene_names = expand("/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/mapped_reads/gene_names_{sample}", sample=SAMPLES[0]),
        header = "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/mapped_reads/header.txt"
    output:
        "/home/projects/ku_00041/people/markri/data/readmap/vfdb/gene_cat/nonredundant_gene_clusters/gene_count_matrix.tsv"
    threads: 1
    resources:
        mem_gb = 188,
        runtime = 1800 #30min
    shell:
        "paste {input.gene_names} {input.counts} | cat {input.header} - > {output}; sed -i '$d' {output}"

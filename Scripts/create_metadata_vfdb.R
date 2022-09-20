library(readxl)
library(readr)
library(tidyverse)
vfdb_meta <- read_excel("C:/Users/marti/Desktop/Speciale/data/vfdb_meta.xls", skip = 1)

vfdb_prodID_acc <- read_delim("C:/Users/marti/Desktop/Speciale/data/vfdb_prodID_acc.tsv", 
                              delim = "\t", escape_double = FALSE, 
                              col_names = c("ProdID", "VFGID"), trim_ws = TRUE)

vfdb_fasta_headers <- read_csv("C:/Users/marti/Desktop/Speciale/data/vfdb_fasta_headers.txt", 
                               col_names = FALSE)

VFGID <- word(vfdb_fasta_headers$X1, 1)
VFGID <- substring(VFGID, 2)

VFID <- str_extract(vfdb_fasta_headers$X1, pattern = "VF\\d+")

vfg_vfid_fromheaders <- tibble(VFGID, VFID)

joined_meta <- left_join(vfg_vfid_fromheaders, vfdb_meta, by = "VFID")

final_table <- left_join(vfdb_prodID_acc, joined_meta, by = "VFGID")


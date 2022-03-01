## code to prepare `beau_annots` dataset goes here

library(tidyverse, warn.conflicts = F)

beau_names <-
  read.csv("./data-raw/beau_protein_gene_names_ncbi.csv",
           header=T, stringsAsFactors = F, na.strings = c(NA,""," ")) %>%
  as_tibble() %>%
  select(gene_name = locus_tag, gene_desc=protein_name)

beau_annots <-
  read.csv(paste0("./data-raw/beau_annots_robin_ncbi.csv"),
           header=T, stringsAsFactors = F, na.strings = c(NA,""," ")) %>%
  as_tibble() %>%
  left_join(beau_names, by="gene_name") %>%
  select(gene_name, gene_desc, everything())

write.csv(beau_annots,"data-raw/beau_annots.csv")
usethis::use_data(beau_annots, overwrite = T)

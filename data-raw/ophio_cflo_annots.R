## code to prepare `ophio_cflo_annots` dataset goes here

library(tidyverse)

ophio_cflo_annots <-
  read.csv(paste0("./data-raw/ophio_cflo_annots_robin_ncbi.csv"),
                              header=T, stringsAsFactors = F, na.strings = c(NA,""," ")) %>%
  as_tibble() %>%
  select(gene_name, gene_desc=blast_annot, everything())

usethis::use_data(ophio_cflo_annots, overwrite = T)

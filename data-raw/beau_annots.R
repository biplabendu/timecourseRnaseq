## code to prepare `beau_annots` dataset goes here

library(tidyverse)

beau_annots <-
  read.csv(paste0("./data-raw/beau_annots_robin_ncbi.csv"),
           header=T, stringsAsFactors = F, na.strings = c(NA,""," ")) %>%
  as_tibble() %>%
  mutate(gene_desc="not_available")


usethis::use_data(beau_annots)

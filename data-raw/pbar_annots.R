## code to prepare `pbar_annots` dataset goes here

library(tidyr)
library(dplyr)

pbar_annots <-
  read.csv(paste0("./data-raw/pbar_annots_2May23.csv"),
                 header=T, stringsAsFactors = F, na.strings = c(NA,""," ")) %>%
  as_tibble() %>%
  select(gene_name, gene_desc=old_annotation, everything())

usethis::use_data(pbar_annots)

## code to prepare `cflo_annots` dataset goes here

library(tidyr)
library(dplyr)

cflo_annots <-
  read.csv(paste0("./data-raw/cflo_annots.csv"),
                 header=T, stringsAsFactors = F, na.strings = c(NA,""," ")) %>%
  as_tibble() %>%
  select(gene_name, gene_desc=old_annotation, everything())

usethis::use_data(cflo_annots)

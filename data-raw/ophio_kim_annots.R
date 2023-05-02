## code to prepare `ophio_kim_annots` dataset goes here

library(tidyverse)

ophio_kim_annots <-
  read.csv(paste0("./data-raw/ophio_kim_annots_robin_ncbi.csv"),
           header=T, stringsAsFactors = F, na.strings = c(NA,""," ")) %>%
  as_tibble() %>%
  mutate(gene_desc="not_available")

usethis::use_data(ophio_kim_annots, overwrite = T)

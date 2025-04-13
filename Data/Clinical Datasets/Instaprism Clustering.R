#libraries####
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(purrr)
library(ggplot2)
library(ggpubr)
library(mclust)
library(NbClust)

# Load the data####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell Types (not batch corrected).txt")
combined_standard <- combined_standard %>% mutate_at(vars(25:33), as.numeric)
immune <- combined_standard %>%
  filter(!(platform == "microarray")) %>%
  filter(!(primary == "FALSE")) %>%
  select(c(patient_id, ER_status, HER2_status, "B-cells", "T-cells", Myeloid, DC))

#ER negative patients
er_neg <- immune %>%
  filter(ER_status == "neg") %>%
  filter(is.na(Endothelial) == FALSE) %>%
  select(-ER_status, -HER2_status) %>%
  column_to_rownames(var = "patient_id") %>%
  mutate_all(as.numeric)

mclust_er_neg <- Mclust(er_neg) #mclust clustering
plot(mclust_er_neg)

nbclust_er_neg <- NbClust(er_neg, distance = "euclidean", min.nc = 2, max.nc = 15, 
        method = "single") #nbclust clustering

#ER positive patients
er_pos <- immune %>%
  filter(ER_status == "pos") %>%
  filter(is.na(Endothelial) == FALSE) %>%
  select(-ER_status, -HER2_status) %>%
  column_to_rownames(var = "patient_id") %>%
  mutate_all(as.numeric)

mclust_er_pos <- Mclust(er_pos) #mclust clustering
plot(mclust_er_pos)

nbclust_er_pos <-NbClust(er_pos, distance = "euclidean", min.nc = 2, max.nc = 15, 
        method = "average") #nbclust clustering











  




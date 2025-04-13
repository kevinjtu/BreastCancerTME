#libraries####
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(xlsx)
library(reshape2)
library(ggplot2)
library(forcats)

#import data
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell Types (not batch corrected).txt")
combined_standard <- rename(combined_standard, c(B_cells = "B-cells", T_cells = "T-cells", Normal_Epithelial = "Normal Epithelial", Cancer_Epithelial = "Cancer Epithelial"))
combined_standard <- combined_standard %>% mutate_at(vars(25:33), as.numeric)

#make proportion graph####ß
prop <- combined_standard %>%
  filter(primary == "TRUE") %>%
  filter(platform == "rnaseq") %>%
  filter(!(study == "MBC")) %>%
  filter(!(study == "ispy2")) %>%
  select(patient_id, study, pam50, intclust, ER_status, Endothelial, CAFs, PVL, B_cells, T_cells, Myeloid, Normal_Epithelial, Cancer_Epithelial, DC)

prop_long <- melt(prop, id.vars = c("patient_id", "study", "pam50", "intclust", "ER_status"), variable.name = "cell_type", value.name = "proportion")

prop_long$patient_id <- factor(prop_long$patient_id, 
                               levels = sum_cancer$patient_id[order(sum_cancer$sum_Cancer_Epithelial)])

ggplot(prop_long, aes(fill = cell_type, y = proportion, x = patient_id)) +  
  geom_bar(stat="identity") +
  facet_wrap(~ intclust, scales = "free") +
  labs(x = "Patient ID", y = "Proportion", title = "Stacked Bar Chart of Cell Type Proportions by intclust") +
  theme(axis.text.x = element_blank())  # To hide patient IDs on x-axis


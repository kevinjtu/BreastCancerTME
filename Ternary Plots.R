#Libraries####
library(data.table)
library(dplyr)
library(tidyr)
library(writexl)
library(readxl)
library(ggplot2)
library(ggtern)
library(ggpubr)
require(ggtern_themes)

###data prep####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
imc <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 1)
instaprism <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 17)
digpath <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 13)

setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
clinical <- fread('metabric_clinical.txt')
clinical <- clinical %>%
  rename_with(~"metabric_id", .cols = "METABRIC.ID") %>%
  select(metabric_id, Cellularity)

imc_instaprism <- imc %>%
  merge(instaprism, by = "metabric_id") %>%
  merge(clinical, by = "metabric_id") %>%
  filter(!is.na(Cellularity))

imc_digpath <- imc %>%
  merge(digpath, by = "metabric_id") %>%
  merge(clinical, by = "metabric_id") %>%
  filter(!is.na(Cellularity))

all <- merge(imc, instaprism, by = "metabric_id")
all <- merge(all, digpath, by = "metabric_id")
all <- merge(all, clinical, by = "metabric_id")

imc <- all[, c("metabric_id", "IMC_immune", "IMC_stroma", "IMC_epithelial", "Cellularity")]
imc <- imc %>%
  filter(!is.na(Cellularity))
instaprism <- all[, c("metabric_id", "instaprism_immune", "instaprism_stroma", "instaprism_epithelial", "Cellularity")]
digpath <- all[, c("metabric_id", "digpath_immune", "digpath_stroma", "digpath_epithelial", "Cellularity")]

#Ternary plots####
imc$Cellularity <- factor(imc$Cellularity, levels = c("low", "moderate", "high"))
instaprism$Cellularity <- factor(instaprism$Cellularity, levels = c("low", "moderate", "high"))
digpath$Cellularity <- factor(digpath$Cellularity, levels = c("low", "moderate", "high"))

ggtern(data = imc, aes(x = IMC_immune, y = IMC_stroma, z = IMC_epithelial)) +
  geom_point(aes(fill = Cellularity), shape = 21, size = 2.5) +  
  scale_fill_manual(values = c("#A0C4FF", "#CAFFBF", "#FFADAD")) +
  theme_linedraw() +
  ggtitle("Image Mass Cytometry") +
  xlab("Imm") +
  ylab("Str") +
  zlab("Epi")
 
ggtern(data = instaprism, aes(x = instaprism_immune, y = instaprism_stroma, z = instaprism_epithelial)) +
  geom_point(aes(fill = Cellularity), shape = 21, size = 2.5) +  
  scale_fill_manual(values = c("#A0C4FF", "#CAFFBF", "#FFADAD")) +
  theme_linedraw() +
  ggtitle("InstaPrism") +
  xlab("Imm") +
  ylab("Str") +
  zlab("Epi")

ggtern(data = digpath, aes(x = digpath_immune, y = digpath_stroma, z = digpath_epithelial)) +
  geom_point(aes(fill = Cellularity), shape = 21, size = 2.5) +  
  scale_fill_manual(values = c("#A0C4FF", "#CAFFBF", "#FFADAD")) +
  theme_linedraw() +
  ggtitle("Digital Pathology") +
  xlab("Imm") +
  ylab("Str") +
  zlab("Epi")

#IMC-Instaprism Data modality plots####
long_form <- imc_instaprism %>%
  gather(key = "method_tissue", value = "proportion", -metabric_id) %>%
  left_join(imc_instaprism[, c("metabric_id", "Cellularity")], by = "metabric_id") %>%
  na.omit() %>%
  filter(method_tissue != "Cellularity")

# Split the method_tissue column into two columns: method and tissue_type
long_form <- long_form %>%
  separate(method_tissue, into = c("method", "tissue_type"), sep = "_")

# Spread the method column into separate columns
long_form <- long_form %>%
  spread(key = method, value = proportion)

# Rename the columns
colnames(long_form) <- c("metabric_id", "tissue_type", "Cellularity", "IMC", "instaprism")

ggplot(long_form, aes(x = IMC_value, y = instaprism_value, color = tissue_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  ggtitle("RNA") +
  xlab("IMC") +
  ylab("InstaPrism") +
  stat_cor(method = "pearson", show.legend = TRUE)

#digpath-Instaprism Data modality plots####
imc_digpath <- merge(imc, digpath, by = "metabric_id")
long_form <- imc_digpath %>%
  gather(key = "method_tissue", value = "proportion", -metabric_id)

# Split the method_tissue column into two columns: method and tissue_type
long_form <- long_form %>%
  separate(method_tissue, into = c("method", "tissue_type"), sep = "_")

# Spread the method column into separate columns
long_form <- long_form %>%
  spread(key = method, value = proportion)

# Rename the columns
colnames(long_form) <- c("metabric_id", "tissue_type", "IMC_value", "instaprism_value")

ggplot(long_form, aes(x = IMC_value, y = instaprism_value, color = tissue_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  ggtitle("Pathology") +
  xlab("IMC") +
  ylab("Digital Pathology") +
  stat_cor(method = "pearson", show.legend = TRUE)



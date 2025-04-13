# libraries ####
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(stringr)
library(readxl)
library(DescTools)
library(Metrics)
library(ComplexHeatmap)
library(circlize)
library(hydroGOF)
library(LaplacesDemon)
library(pheatmap)

#import data ####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Bulk Expression Cell Type Quantification")
imagedata <- fread("SingleCells.csv")

##find number of Treg_Fox3P and Endothelial_CXCL12
#identify cancer cycling cells
imagedata[is_tumour == 1 & cellPhenotype == "Ep Ki67^{+}", cellPhenotype := "Cancer_Cycling"]

#treg fox3p
subset_data <- subset(imagedata, cellPhenotype == "T_{Reg} & T_{Ex}" & !is.na(FOXP3))
foxp3_values <- as.numeric(subset_data$FOXP3)
median_fox3p <- median(foxp3_values, na.rm = TRUE)

imagedata$cellPhenotype <- ifelse(
  imagedata$cellPhenotype == "T_{Reg} & T_{Ex}" & imagedata$FOXP3 > median_fox3p,
  "T_{Reg} & T_{Ex}_fox3p_hi",
  ifelse(
    imagedata$cellPhenotype == "T_{Reg} & T_{Ex}",
    "T_{Reg} & T_{Ex}_fox3p_low",
    imagedata$cellPhenotype
  )
)

#Endothelial_CXCL12
subset_data <- subset(imagedata, cellPhenotype == "Endothelial" & !is.na(CXCL12))
CXCL12_values <- as.numeric(subset_data$CXCL12)
median_CXCL12 <- median(CXCL12_values, na.rm = TRUE)

imagedata$cellPhenotype <- ifelse(
  imagedata$cellPhenotype == "Endothelial" & imagedata$CXCL12 > median_CXCL12,
  "Endothelial_CXCL12_hi",
  ifelse(
    imagedata$cellPhenotype == "Endothelial",
    "Endothelial_CXCL12_low",
    imagedata$cellPhenotype
  )
)

##find proportion of cell types
Proportionresult <- imagedata %>%
  dplyr::group_by(metabric_id, cellPhenotype) %>%
  dplyr::summarize(Total_Cells = n()) %>%
  dplyr::group_by(metabric_id) %>%
  dplyr::mutate(Proportion = Total_Cells/sum(Total_Cells))

#filter for HER2^{+}, Basal, Ep Ki67^{+}, Myofibroblasts, Endothelial_CXCL12_hi, T_{Reg} & T_{Ex}_fox3p_hi
Proportionresult <- Proportionresult %>%
  dplyr::filter(
    cellPhenotype %in% c(
      "HER2^{+}",
      "Basal",
      "Cancer_Cycling",
      "Endothelial_CXCL12_hi",
      "T_{Reg} & T_{Ex}_fox3p_hi"))

Proportionresult_wide <- Proportionresult %>%
  dplyr::select(metabric_id, cellPhenotype, Proportion) %>%
  tidyr::pivot_wider(names_from = cellPhenotype, values_from = Proportion) %>%
  dplyr::rename(
    imc_Cancer_Her2 = "HER2^{+}",
    imc_Cancer_Basal = "Basal",
    imc_Cancer_Cycling = "Cancer_Cycling",
    imc_Endothelial_CXCL12 = "Endothelial_CXCL12_hi",
    imc_Treg_FOX3P = "T_{Reg} & T_{Ex}_fox3p_hi",
    patient_id = "metabric_id"
  )

#benchmark against metabric instaprism results
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
instaprism <- read_xlsx("METABRIC InstaPrism results.xlsx", sheet = 3)
instaprism$instaprism_Cancer_Cycling <- rowSums(instaprism[, c("Cancer Cycling_1", "Cancer Cycling_2", "Cancer Cycling_3", "Cancer Cycling_4", "Cancer Cycling_5", "Cancer Cycling_6", "Cancer Cycling_7", "Cancer Cycling_8")], na.rm = TRUE)
instaprism$instaprism_Cancer_Her2 <- rowSums(instaprism[, c("Cancer Her2 SC_1", "Cancer Her2 SC_2", "Cancer Her2 SC_3", "Cancer Her2 SC_4", "Cancer Her2 SC_5", "Cancer Her2 SC_6")], na.rm = TRUE)
instaprism$instaprism_Cancer_Basal <- rowSums(instaprism[, c("Cancer Basal SC_1", "Cancer Basal SC_2", "Cancer Basal SC_3", "Cancer Basal SC_4", "Cancer Basal SC_5", "Cancer Basal SC_6", "Cancer Basal SC_7", "Cancer Basal SC_8")], na.rm = TRUE)
instaprism <- instaprism %>% 
  dplyr::select(-starts_with("Cancer Cycling_"), -starts_with("Cancer Her2 SC_"), -starts_with("Cancer LumB SC_"), -starts_with("Cancer LumA SC_"), -starts_with("Cancer Basal SC_"))
instaprism <- instaprism %>%
  dplyr::rename(
    instaprism_Treg_FOX3P = `T_cells_c2_CD4+_T-regs_FOXP3`,
    instaprism_Endothelial_CXCL12 = `Endothelial CXCL12`) %>%
  dplyr::select(patient_id, `instaprism_Cancer_Her2`, `instaprism_Cancer_Basal`, `instaprism_Cancer_Cycling`, `instaprism_Treg_FOX3P`, `instaprism_Endothelial_CXCL12`)

setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
clinical <- fread('metabric_clinical.txt')
clinical <- clinical %>%
  dplyr::rename(patient_id = `METABRIC.ID`) %>%
  select(patient_id, Cellularity)
instaprism <- merge(instaprism, clinical, by = "patient_id", all = TRUE)

allimmune <- merge(Proportionresult_wide, instaprism, by = "patient_id", all = TRUE)

###benchmarking results heatmap####
#mdae
mdae_values <- c(
  mdae(
    na.omit(allimmune[c("imc_Endothelial_CXCL12", "instaprism_Endothelial_CXCL12")])$imc_Endothelial_CXCL12, 
    na.omit(allimmune[c("imc_Endothelial_CXCL12", "instaprism_Endothelial_CXCL12")])$instaprism_Endothelial_CXCL12),
  mdae(
    na.omit(allimmune[c("imc_Cancer_Cycling", "instaprism_Cancer_Cycling")])$imc_Cancer_Cycling, 
    na.omit(allimmune[c("imc_Cancer_Cycling", "instaprism_Cancer_Cycling")])$instaprism_Cancer_Cycling),
  mdae(
    na.omit(allimmune[c("imc_Cancer_Basal", "instaprism_Cancer_Basal")])$imc_Cancer_Basal, 
    na.omit(allimmune[c("imc_Cancer_Basal", "instaprism_Cancer_Basal")])$instaprism_Cancer_Basal),
  mdae(
    na.omit(allimmune[c("imc_Cancer_Her2", "instaprism_Cancer_Her2")])$imc_Cancer_Her2, 
    na.omit(allimmune[c("imc_Cancer_Her2", "instaprism_Cancer_Her2")])$instaprism_Cancer_Her2),
  mdae(
    na.omit(allimmune[c("imc_Treg_FOX3P", "instaprism_Treg_FOX3P")])$imc_Treg_FOX3P, 
    na.omit(allimmune[c("imc_Treg_FOX3P", "instaprism_Treg_FOX3P")])$instaprism_Treg_FOX3P))

ordered_values <- sort(mdae_values)
mdae_matrix <- matrix(ordered_values, nrow = 1, ncol = 5)

colnames(mdae_matrix) <- c("Endothelial [CXCL12]", "Cancer Cyling", "Cancer Basal", "Cancer HER2", "Treg [FOX3P]")

mdae_matrix_t <- t(mdae_matrix)
mdae_matrix_t <- as.data.frame(mdae_matrix_t)
colnames(mdae_matrix_t)[colnames(mdae_matrix_t) == "V1"] <- "MAE"

col_fun <- colorRamp2(c(0, 0.1), c("#FFFFFF", "#FFADAD"))
pheatmap(mdae_matrix_t,
         cellwidth = 30,
         display_numbers = TRUE,    
         fontsize_number = 12,
         cluster_rows = FALSE,          # Do not cluster rows
         cluster_cols = FALSE,          # Do not cluster columns
         color = col_fun(seq(0, 0.1, length.out=100)),
         legend = TRUE,
         ) 

#RMSE
rmse_values <- c(
  hydroGOF::rmse(allimmune$imc_Endothelial_CXCL12, allimmune$instaprism_Endothelial_CXCL12, na.rm = TRUE),
  hydroGOF::rmse(allimmune$imc_Cancer_Cycling, allimmune$instaprism_Cancer_Cycling, na.rm = TRUE),
  hydroGOF::rmse(allimmune$imc_Cancer_Her2, allimmune$instaprism_Cancer_Her2, na.rm = TRUE),
  hydroGOF::rmse(allimmune$imc_Cancer_Basal, allimmune$instaprism_Cancer_Basal, na.rm = TRUE),
  hydroGOF::rmse(allimmune$imc_Treg_FOX3P, allimmune$instaprism_Treg_FOX3P, na.rm = TRUE))
ordered_values <- sort(rmse_values)
rmse_matrix <- matrix(ordered_values, nrow = 1, ncol = 5)

colnames(rmse_matrix) <- c("Endothelial [CXCL12]", "Cancer Cyling", "Cancer Basal", "Cancer HER2", "Treg [FOX3P]")

rmse_matrix_t <- t(rmse_matrix)
rmse_matrix_t <- as.data.frame(rmse_matrix_t)
colnames(rmse_matrix_t)[colnames(rmse_matrix_t) == "V1"] <- "RMSE"

col_fun <- colorRamp2(c(0, 0.1), c("#FFFFFF", "#A0C4FF"))
pheatmap(rmse_matrix_t,
         cellwidth = 30,
         display_numbers = TRUE, 
         fontsize_number = 12,
         cluster_rows = FALSE,          # Do not cluster rows
         cluster_cols = FALSE,          # Do not cluster columns
         color = col_fun(seq(0, 0.1, length.out=100)),
         legend = TRUE,
) 

medae %v% rmse

#KDL
KLD_values <- c(
  KLD(
    na.omit(allimmune[c("imc_Endothelial_CXCL12", "instaprism_Endothelial_CXCL12")])$imc_Endothelial_CXCL12, 
    na.omit(allimmune[c("imc_Endothelial_CXCL12", "instaprism_Endothelial_CXCL12")])$instaprism_Endothelial_CXCL12)[["mean.sum.KLD"]],
  KLD(
    na.omit(allimmune[c("imc_Cancer_Cycling", "instaprism_Cancer_Cycling")])$imc_Cancer_Cycling, 
    na.omit(allimmune[c("imc_Cancer_Cycling", "instaprism_Cancer_Cycling")])$instaprism_Cancer_Cycling)[["mean.sum.KLD"]],
  KLD(
    na.omit(allimmune[c("imc_Cancer_Basal", "instaprism_Cancer_Basal")])$imc_Cancer_Basal, 
    na.omit(allimmune[c("imc_Cancer_Basal", "instaprism_Cancer_Basal")])$instaprism_Cancer_Basal)[["mean.sum.KLD"]],
  KLD(
    na.omit(allimmune[c("imc_Cancer_Her2", "instaprism_Cancer_Her2")])$imc_Cancer_Her2, 
    na.omit(allimmune[c("imc_Cancer_Her2", "instaprism_Cancer_Her2")])$instaprism_Cancer_Her2)[["mean.sum.KLD"]],
  KLD(
    na.omit(allimmune[c("imc_Treg_FOX3P", "instaprism_Treg_FOX3P")])$imc_Treg_FOX3P, 
    na.omit(allimmune[c("imc_Treg_FOX3P", "instaprism_Treg_FOX3P")])$instaprism_Treg_FOX3P)[["mean.sum.KLD"]])

ordered_values <- sort(KLD_values)
KLD_matrix <- matrix(ordered_values, nrow = 1, ncol = 5)

colnames(KLD_matrix) <- c("Endothelial [CXCL12]", "Cancer Cyling", "Cancer Basal", "Cancer [Her2]", "Treg [FOX3P]")
col_fun = colorRamp2(c(0, 0.20, 0.40, 0.60, 0.80), c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"))

kld <- Heatmap(
  KLD_matrix,
  col = col_fun,
  na_col = "grey",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", KLD_matrix[i, j]), x, y, gp = gpar(fontsize = 16))},
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  rect_gp = gpar(col = "white", lwd = 2),
  row_title = "KLD",
  name = "Metric",
  width = unit(10, "cm"),
  height = unit(2, "cm"),
  column_names_side = "bottom",
)

medae %v% rmse

#spearmanR
spearman_values <- c(
  hydroGOF::rSpearman(allimmune$imc_Endothelial_CXCL12, allimmune$instaprism_Endothelial_CXCL12, na.rm = TRUE),
  hydroGOF::rSpearman(allimmune$imc_Cancer_Cycling, allimmune$instaprism_Cancer_Cycling, na.rm = TRUE),
  hydroGOF::rSpearman(allimmune$imc_Cancer_Her2, allimmune$instaprism_Cancer_Her2, na.rm = TRUE),
  hydroGOF::rSpearman(allimmune$imc_Cancer_Basal, allimmune$instaprism_Cancer_Basal, na.rm = TRUE),
  hydroGOF::rSpearman(allimmune$imc_Treg_FOX3P, allimmune$instaprism_Treg_FOX3P, na.rm = TRUE))
ordered_values <- sort(spearman_values)
spearman_matrix <- matrix(ordered_values, nrow = 1, ncol = 5)

colnames(spearman_matrix) <- c("Endothelial [CXCL12]", "Cancer Cyling", "Cancer Basal", "Cancer HER2", "Treg [FOX3P]")

spearman_matrix_t <- t(spearman_matrix)
spearman_matrix_t <- as.data.frame(spearman_matrix_t)
colnames(spearman_matrix_t)[colnames(spearman_matrix_t) == "V1"] <- "SpearmanR"

col_fun <- colorRamp2(c(0, 0.1), c("#FFFFFF", "#80EF80"))
pheatmap(spearman_matrix_t,
         cellwidth = 30,
         display_numbers = TRUE, 
         fontsize_number = 12,
         cluster_rows = FALSE,          # Do not cluster rows
         cluster_cols = FALSE,          # Do not cluster columns
         color = col_fun(seq(0, 0.1, length.out=100)),
         legend = TRUE,
) 

medae %v% rmse
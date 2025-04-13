#libraries####
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(PCAtools)
library(ggplot2)
library(compositions)
library(patchwork)

# Load the data####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
combined_standard <- combined_standard %>%
  select(patient_id, study, platform, everything()[25:73])
combined_standard_trans <- as.data.frame(clr(acomp(combined_standard[,4:52])))
combined_standard_trans <- cbind(combined_standard[,1:3], combined_standard_trans)
combined_standard_trans <- distinct(combined_standard_trans, patient_id, .keep_all = TRUE)

#pca analysis of data####
data_subset <- as.data.frame(t(combined_standard_trans[, 4:52, with = FALSE]))
colnames(data_subset) <- combined_standard_trans$patient_id
metadata_df <- combined_standard_trans[, c("patient_id", "study", "platform")]
metadata_df <- column_to_rownames(metadata_df, var = "patient_id")
p <- pca(as.matrix(data_subset), 
         metadata = metadata_df)
           
biplot(p, colby = "study", lab = NULL, 
       legendPosition = 'right')
biplot(p, colby = "platform", lab = NULL, 
       legendPosition = 'right')


#pca analysis plot each study indivdually, microarray####
microarray_studies <- c("HatzisValidation", "HatzisDiscovery", "METABRIC", "ispy2")
plots <- list() 
for (currentstudy in microarray_studies) {
  combined_standard_trans_microarray <- combined_standard_trans %>%
    filter(platform == "microarray")
  combined_standard_trans_ordered <- combined_standard_trans_microarray[order(combined_standard_trans_microarray$study != currentstudy, decreasing = TRUE), ]
  data_subset <- as.data.frame(t(combined_standard_trans_ordered[, 4:52, with = FALSE]))
  colnames(data_subset) <- combined_standard_trans_ordered$patient_id
  metadata_df <- combined_standard_trans_ordered[, c("patient_id", "study", "platform")]
  metadata_df <- column_to_rownames(metadata_df, var = "patient_id")
  
  metadata_df$color <- ifelse(metadata_df$study == currentstudy, "Study", "AllOthers")
  p <- pca(as.matrix(data_subset), 
           metadata = metadata_df)
  
  plot <- biplot(p, 
                 colby = "color", # Don't use the default colby argument
                 x = "PC1",
                 y = "PC2",
                 lab = NULL, 
                 legendPosition = "none",
                 pointSize = 0.5, 
                 shape = 'color', shapekey = c('Study' = 16, 'AllOthers' = 1),
                 ellipse = TRUE,
                 ellipseLevel = 0.95,
                 ellipseFill = TRUE,
                 ellipseLineSize = 0.5,
                 ellipseLineCol = "black"
                 ) + 
    ggtitle(paste(currentstudy))+
    labs(x = "", y = "")+
    theme(plot.title = element_text(margin = margin(b = -20)))
    plot
    plots[[currentstudy]] <- plot
    setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Breast cancer microenvironment paper/Supplementary")
    ggsave(paste0("PCA_Plot_", currentstudy, ".png"), plot = plot, width = 3, height = 3, dpi = 300)
}

setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")

#pca analysis plot each study indivdually, rnaseq####
rnaseq_studies <- c("SCANB", "SMC", "TCGA", "MBC", "Derouane", "CALGB", "Matador", "aurora", "transneo", "newton")
plots <- list() 
for (currentstudy in rnaseq_studies) {
  combined_standard_trans_rnaseq <- combined_standard_trans %>%
    filter(platform == "rnaseq")
  combined_standard_trans_ordered <- combined_standard_trans_rnaseq[order(combined_standard_trans_rnaseq$study != currentstudy, decreasing = TRUE), ]
  data_subset <- as.data.frame(t(combined_standard_trans_ordered[, 4:52, with = FALSE]))
  colnames(data_subset) <- combined_standard_trans_ordered$patient_id
  metadata_df <- combined_standard_trans_ordered[, c("patient_id", "study", "platform")]
  metadata_df <- column_to_rownames(metadata_df, var = "patient_id")
  
  metadata_df$color <- ifelse(metadata_df$study == currentstudy, "Study", "AllOthers")
  p <- pca(as.matrix(data_subset), 
           metadata = metadata_df)
  
  # Invert PC1 and PC2 for SCANB to fix axis inversion
  if (currentstudy == "SCANB") {
    p$rotated$PC1 <- -p$rotated$PC1
    p$rotated$PC2 <- -p$rotated$PC2
  }
  
  plot <- biplot(p, 
                 colby = "color", # Don't use the default colby argument
                 x = "PC1",
                 y = "PC2",
                 lab = NULL, 
                 legendPosition = "none",
                 pointSize = 0.5, 
                 shape = 'color', shapekey = c('Study' = 16, 'AllOthers' = 1),
                 ellipse = TRUE,
                 ellipseLevel = 0.95,
                 ellipseFill = TRUE,
                 ellipseLineSize = 0.5,
                 ellipseLineCol = "black"
  ) + 
    ggtitle(paste(currentstudy))+
    labs(x = "", y = "")+
    theme(plot.title = element_text(margin = margin(b = -20)))
  plot
  plots[[currentstudy]] <- plot
  setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Breast cancer microenvironment paper/Supplementary")
  ggsave(paste0("PCA_Plot_", currentstudy, ".png"), plot = plot, width = 3, height = 3, dpi = 300)
}

setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")

#Bargraphs of batch effect####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
combined_standard <- combined_standard %>%
  select(patient_id, study, platform, everything()[25:73])

combined_summed <- combined_standard %>%
  mutate(Cancer_Epithelial = rowSums(select(., matches("Cancer"))),
         Normal_Epithelial = rowSums(select(., matches("Luminal|epithelial"))),
         T_cells = rowSums(select(., matches("T_cells"))),
         B_cells = rowSums(select(., matches("B cells"))),
         Myeloid = rowSums(select(., matches("Myeloid"))),
         Endothelial = rowSums(select(., matches("Endothelial"))),
         PVL = rowSums(select(., matches("PVL"))),
         CAFs = rowSums(select(., matches("CAF"))),
         Plasmablasts = rowSums(select(., matches("Plasmablasts"))))%>%
  dplyr::select(patient_id, study, platform, Cancer_Epithelial, Normal_Epithelial, T_cells,
         B_cells, Myeloid, Endothelial, PVL, CAFs, Plasmablasts)

combined_summed <- melt(combined_summed)

ggplot(combined_summed, aes(x = study, y = value, fill = platform)) +
  geom_boxplot() +
  stat_summary(fun.y=median, geom="point", shape=20, size=3, color="red", fill="red") +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  labs(x = "Study", y = "Proportion", title = "Batch Effects") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






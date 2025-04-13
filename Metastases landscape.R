#library####
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(ggplot2)
library(purrr)
library(rstatix)
library(ggpubr)
library(ggridges)

# Load the data####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
combined_standard$pcr <- gsub("NA", NA, combined_standard$pcr)

new_names <- names(combined_standard)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "T-regs", "Tregs")  # Replace T-regs with Tregs
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_standard <- setNames(combined_standard, new_names)

#rename cell types
combined_standard_new <- setDT(combined_standard)
column_mapping <- data.table(
  old_name = c(
    "Endothelial_ACKR1", "Endothelial_RGS5", "Endothelial_CXCL12", 
    "CAFs_MSC_iCAFlike_s1", "CAFs_MSC_iCAFlike_s2", "CAFs_Transitioning_s3",
    "CAFs_myCAF_like_s4", "CAFs_myCAF_like_s5", "PVL_Differentiated_s3", 
    "PVL_Immature_s2", "PVL_Immature_s1", "Endothelial_Lymphatic_LYVE1", 
    "B_cells_Memory", "B_cells_Naive", "T_cells_c4_CD8_ZFP36", 
    "T_cells_c6_IFIT1", "T_cells_c5_CD8_GZMK", "T_cells_c7_CD8_IFNG", 
    "T_cells_c8_CD8_LAG3", "T_cells_c0_CD4_CCR7", "T_cells_c1_CD4_IL7R", 
    "T_cells_c2_CD4_Tregs_FOXP3", "T_cells_c3_CD4_Tfh_CXCL13", "T_cells_c9_NK_cells_AREG", 
    "T_cells_c11_MKI67", "T_cells_c10_NKT_cells_FCGR3A", "Myeloid_c10_Macrophage_1_EGR1", 
    "Myeloid_c12_Monocyte_1_IL1B", "Myeloid_c2_LAM2_APOE", "Myeloid_c1_LAM1_FABP5", 
    "Myeloid_c8_Monocyte_2_S100A9", "Myeloid_c7_Monocyte_3_FCGR3A", "Myeloid_c9_Macrophage_2_CXCL10", 
    "Cycling_Myeloid", "Myeloid_c11_cDC2_CD1C", "Myeloid_c4_DCs_pDC_IRF7", 
    "Myeloid_c3_cDC1_CLEC9A", "Myeloid_c0_DC_LAMP3", "Myoepithelial", 
    "Luminal_Progenitors", "Mature_Luminal", "Plasmablasts", "Cycling_PVL", 
    "Myeloid_c5_Macrophage_3_SIGLEC1"
  ),
  new_name = c(
    "Venous Vessel", "Tip-like Vessel", "Tip-like Vessel", 
    "iCAF", "iCAF", "Transitioning CAF", "myCAF", "myCAF", 
    "Differentiated PVL", "Immature PVL", "Immature PVL", "Lymphatic Vessel", 
    "Memory B Cells", "Naive B Cells", "Chemokine Expressing CD8 T Cells", 
    "Type 1 Interferon T Cells", "Effector Cytotoxic CD8 T Cells", 
    "Activated Cytotoxic CD8 T Cells", "Inhibitory CD8 T Cells", 
    "Central Memory CD4 T Cells", "Th1 Cells", "Regulatory T Cells", 
    "T Follicular Helper Cells", "Natural Killer Cells", "Proliferating T Cells", 
    "Natural killer T Cells", "Macrophage M2", "Monocyte [IL1B]", 
    "Lipid-Associated Macrophage 2", "Lipid-Associated Macrophage 1", 
    "Monocyte [S100A9]", "Monocyte [FCGR3A]", "Macrophage M1", 
    "Proliferating Myeloid Cells", "Conventional DC", "Plasmacytoid DC", 
    "Conventional DC", "Plasmacytoid DC", "Myoepithelial Cells", 
    "Luminal Progenitor Cells", "Mature Luminal Cells", "Plasmablast Cells", 
    "Proliferating PVL", "Macrophage M2"
  )
)
for (new_col in unique(column_mapping$new_name)) {
  cols_to_sum <- column_mapping[old_name %in% names(combined_standard_new) & new_name == new_col, old_name]
  if (length(cols_to_sum) > 0) {
    combined_standard_new[, (new_col) := rowSums(.SD, na.rm = TRUE), .SDcols = cols_to_sum]
  }
}
combined_standard_new[, (column_mapping$old_name) := NULL]
combined_standard_new <- as.data.frame(combined_standard_new)
combined_standard_new[, c(25:29, 32:68)] <- lapply(combined_standard_new[, c(25:29, 32:68)], function(x) as.numeric(as.character(x)))
rowSums(combined_standard_new[, c(25:29, 32:68)])
combined_standard <- combined_standard_new

#make metastases landscape plots (DEPRACATED)####
metastases <- combined_standard %>%
  select(patient_id, sample_id, study, platform, metastases_location, primary, 25:73) %>%
  filter(study %in% c("aurora",))

metastases_long <- melt(metastases, id.vars = c("patient_id", "sample_id", "study", "platform", "primary", "metastases_location"))
metastases_long$value <- as.numeric(metastases_long$value)
metastases_long <- metastases_long %>% 
  #filter(metastases_location %in% c("Bone", "Brain", "Breast", "Liver", "Lung")) %>%
  group_by(study, metastases_location) %>%
  mutate(proportion = as.numeric(value) / sum(as.numeric(value)))

unique_patient_counts <- metastases_long %>%
  group_by(metastases_location, variable, platform) %>%
  summarise(unique_patient_count = n_distinct(patient_id))

ggplot(metastases_long, aes(x = metastases_location, y = proportion, fill = variable)) +
  geom_bar(stat = "identity") +
  facet_wrap(~platform) +
  labs(x = "Metastases Location", y = "Value")


#Compare values between pooled metastases samples (unpaired t test) (DEPRACATED)####
metastases <- combined_standard %>%
  filter(study %in% c("MBC", "aurora"))%>%
  select(patient_id, metastases_location, Endothelial, CAFs, PVL, B_cells, T_cells, Myeloid, Normal_Epithelial, Cancer_Epithelial, DC) %>%
  na.omit() %>%
  filter(metastases_location %in% c("Bone", "Brain", "Breast", "Chest", "Liver", "Lung", "LymphNode", "SoftTissue")) #only looking at metastases_location with over 5 samples

metastases_long <- melt(metastases, id.vars = c("patient_id", "metastases_location"))

cell_type <- colnames(combined_standard)[25:33]
results <- data.frame()

for (cell in cell_type) {
metastases_cell <- metastases_long %>%
  filter(variable == cell)

pwc <- metastases_cell %>% 
  pairwise_t_test(value ~ metastases_location, pool.sd = FALSE, p.adjust.method = "BH")

result <- pwc %>% 
  mutate(cell_type = cell)

results <- rbind(results, result)
}

ggplot(metastases_long, aes(x = metastases_location, y = value, fill = metastases_location)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_y") +
  labs(title = "Metastases Location vs. Cell Type",
       x = "Metastases Location",
       y = "Value")
  
#Compare values between metastases samples (paired t test)####
metastases <- combined_standard %>%
  filter(study %in% c("aurora"))%>%
  select(patient_id, metastases_location, c(25:29, 32:68)) %>%
  na.omit() %>%
  filter(metastases_location %in% c("Brain", "Breast", "Chest", "Liver", "Lung", "LymphNode", "Skin", "SoftTissue")) #only looking at metastases_location with a least 3 paired samples

breast_patients <- metastases %>%
  filter(grepl("Breast", metastases_location)) %>%
  select(patient_id)

breast_patients <- unique(breast_patients$patient_id)
metastases_filtered <- metastases %>% #makes metastases_filtered only have paired samples
  filter(patient_id %in% breast_patients)%>%
  group_by(patient_id) %>%
  filter(n() >= 2) %>%
  ungroup()

metastases_long <- melt(metastases_filtered, id.vars = c("patient_id", "metastases_location"))

met_loc <- list("Liver", "Brain", "Skin", "Lung", "SoftTissue", "LymphNode", "Chest")
cell_type <- colnames(metastases)[3:44]
results <- data.frame()

for (loc in met_loc){
  patients <- metastases_long %>%
    filter(grepl(loc, metastases_location)) %>%
    select(patient_id)
  patients <- unique(patients$patient_id)
  metastases_loc <- metastases_long %>% #makes dataset with only paired samples at location
    filter(patient_id %in% patients)%>%
    filter(metastases_location == "Breast" | metastases_location == loc)
  
  for (cell in cell_type) {
  metastases_loc <- metastases_loc %>%
    distinct(patient_id, metastases_location, variable, .keep_all = TRUE)
  
  metastases_cell <- metastases_loc %>%
    filter(variable == cell) %>%
    select(metastases_location, value)

  pwc <- metastases_cell %>% 
    pairwise_t_test(value ~ metastases_location, pool.sd = FALSE, paired = TRUE, p.adjust.method = "BH")
  pwc <- pwc %>% add_xy_position(x = "metastases_location")
  
  result <- pwc %>% 
    mutate(cell_type = cell)
  
  results <- rbind(results, result)
  }}

#ggplot(metastases_long, aes(x = metastases_location, y = value, fill = metastases_location)) +
  #geom_boxplot() +
  #facet_wrap(~variable, scales = "free_y") +
  #labs(title = "Metastases Location vs. Cell Type",
       #x = "Metastases Location",
       #y = "Value")

#plot the metastases location with breast as baseline of 0####
non_breast_values <- metastases_filtered %>%
  filter(metastases_location != "Breast")
breast_values <- metastases_filtered %>%
  filter(metastases_location == "Breast")
merged_values <- merge(non_breast_values, breast_values, by = "patient_id")

# Subtract values from breast_values columns 3:78 from non_breast_values columns 3:78
for (i in 3:44) {
  merged_values[,i] <- log2((merged_values[,i]/merged_values[,i+43]))
}

merged_values <- merged_values[, c(1:44)]
colnames(merged_values) <- colnames(non_breast_values)
normalized <- melt(merged_values, id.vars = c("patient_id", "metastases_location"))
normalized <- normalized %>%
  dplyr::rename("Cell_State" = "variable")

normalized$Cell_Type <- NA
normalized$Cell_Type[grep("Cancer", normalized$Cell_State)] <- "Cancer Epithelium"
normalized$Cell_Type[grep("Luminal", normalized$Cell_State)] <- "Normal Epithelium"
normalized$Cell_Type[grep("epithelial", normalized$Cell_State)] <- "Epithelium"
normalized$Cell_Type[grep("Plasmablast", normalized$Cell_State)] <- "Lymphoid"
normalized$Cell_Type[grep("Th1", normalized$Cell_State)] <- "Lymphoid"
normalized$Cell_Type[grep("Follicular", normalized$Cell_State)] <- "Lymphoid"
normalized$Cell_Type[grep("Natural Killer", normalized$Cell_State)] <- "Lymphoid"
normalized$Cell_Type[grep("PVL", normalized$Cell_State)] <- "Vasculature"
normalized$Cell_Type[grep("T Cells", normalized$Cell_State)] <- "Lymphoid"
normalized$Cell_Type[grep("B Cell", normalized$Cell_State)] <- "Lymphoid"
normalized$Cell_Type[grep("Monocyte", normalized$Cell_State)] <- "Myeloid"
normalized$Cell_Type[grep("Myeloid", normalized$Cell_State)] <- "Myeloid"
normalized$Cell_Type[grep("Macrophage", normalized$Cell_State)] <- "Myeloid"
normalized$Cell_Type[grep("DC", normalized$Cell_State)] <- "Myeloid"
normalized$Cell_Type[grep("Vessel", normalized$Cell_State)] <- "Vasculature"
normalized$Cell_Type[grep("CAF", normalized$Cell_State)] <- "Fibroblasts"

#add paired test p-values to the normalized dataset
results <- subset(results, grepl("Breast", group1) | grepl("Breast", group2))
results$metloc <- ifelse(grepl("Breast", results$group1), results$group2, results$group1)
normalized <- normalized %>%
  merge(results, by.x = c("Cell_State", "metastases_location"), by.y = c("cell_type", "metloc"), all.x = TRUE) %>%
  select(c("patient_id", "metastases_location", "Cell_State", "value", "Cell_Type", "p.adj"))

mlocations <- unique(normalized$metastases_location)
for(location in mlocations){
normalized_subset <- dplyr::filter(normalized, metastases_location == location)
normalized_subset <- dplyr::filter(normalized_subset,(p.adj < 0.05))

plot <- ggplot(normalized_subset, aes(x = reorder(Cell_State, value), y = value, fill = Cell_Type)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red", size = 1.5) +  # Add dotted red line at y = 0
  labs(title = location,
       x = "Cell State",
       y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1))

print(plot)
}

normalized_subset <- dplyr::filter(normalized,(p.adj < 0.05))
normalized_subset <- normalized_subset %>%
  group_by(metastases_location, Cell_State) %>%
  mutate(average_value = mean(value)) %>%
  ungroup()
  #filter(abs(average_value) >= 0.00)
normalized_subset <- filter(normalized_subset, metastases_location %in% c("Liver"))

#BOXPLOT - DEPRECATED
#plot <- ggplot(normalized_subset, aes(x = reorder(Cell_State, value), y = value, fill = Cell_Type)) +
  #geom_boxplot() +
  #geom_hline(yintercept = 0, linetype = "dotted", color = "red", size = 1.5) +  # Add dotted red line at y = 0
  #labs(title = "Microenvironment changes in metastases",
   #    x = "Cell State",
    #   y = "Percentage Change") +
  #facet_wrap(~metastases_location, scales = "free_y", ncol = 2) +
  #theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1))+
  #coord_cartesian(ylim = c(-100, 100))

x_position <- max(normalized_subset$value)+1

plot <- ggplot(normalized_subset, aes(x = value, y = reorder(Cell_State, -value), fill = Cell_Type)) +
  geom_density_ridges(alpha = 0.8, scale = 1.25) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.5) +  # Add dotted red line at x = 0
  labs(title = "Liver Metastases",
       x = "Log(FC)",
       y = "'") +
  #facet_wrap(~metastases_location, scales = "free_x", ncol = 2) +
  theme(axis.text.y = element_text(size = 8))+
  scale_fill_manual(values = c("Cancer Epithelium" = "#d4a373",
                               "Normal Epithelium" = "#BDB2FF",
                               "Vasculature" = "#A0C4FF",
                               "Lymphoid" = "#FFADAD",
                               "Myeloid" = "#FDFFB6",
                               "Fibroblasts" = "#CAFFBF")) +
  theme_minimal()+
  geom_text(aes(x = x_position, 
                label = round(average_value, 2)), 
            hjust = 0, size = 3)+
  coord_cartesian(clip = "off")+
  theme(legend.position = "bottom", 
        plot.margin = margin(5.5, 15, 5.5, 5.5))  

print(plot)
#Compare values between metastases samples (paired t test, samples that did not subtype switch)####
#list of patient_ids who did not switch
metastases <- combined_standard %>%
  select(patient_id, sample_id, study, metastases_location, primary, pam50) %>%
  filter(study %in% c("aurora")) %>%
  filter(metastases_location %in% c("Brain", "Chest", "Liver", "Lung", "Skin", "SoftTissue", "LymphNode"))

primaries <- combined_standard %>%
  select(patient_id, sample_id, study, metastases_location, primary, pam50) %>%
  filter(study %in% c("aurora")) %>%
  filter(metastases_location %in% c("Breast"))

merged_df <- merge(primaries, metastases, by = "patient_id", suffixes = c("_primary", "_met"))
merged_df <- merged_df[, c("patient_id", "pam50_primary", "pam50_met", "metastases_location_met")]

matching_patients <- merged_df %>%
  filter(pam50_primary == pam50_met) %>%
  pull(patient_id)

#data analysis
metastases <- combined_standard %>%
  filter(study %in% c("aurora"))%>%
  select(patient_id, metastases_location, 25:73) %>%
  na.omit() %>%
  filter(metastases_location %in% c("Brain", "Breast", "Chest", "Liver", "Lung", "LymphNode", "Skin", "SoftTissue")) #only looking at metastases_location with a least 3 paired samples

breast_patients <- metastases %>%
  filter(grepl("Breast", metastases_location)) %>%
  select(patient_id)

breast_patients <- unique(breast_patients$patient_id)
breast_patients <- breast_patients[breast_patients %in% matching_patients]

metastases_filtered <- metastases %>% #makes metastases_filtered only have paired samples
  filter(patient_id %in% breast_patients)%>%
  group_by(patient_id) %>%
  filter(n() >= 2) %>%
  ungroup()

metastases_long <- melt(metastases_filtered, id.vars = c("patient_id", "metastases_location"))

met_loc <- list("Liver", "Brain", "Skin", "Lung", "SoftTissue", "LymphNode", "Chest")
cell_type <- colnames(combined_standard)[25:73]
results <- data.frame()

for (loc in met_loc){
  patients <- metastases_long %>%
    filter(grepl(loc, metastases_location)) %>%
    select(patient_id)
  patients <- unique(patients$patient_id)
  metastases_loc <- metastases_long %>% #makes dataset with only paired samples at location
    filter(patient_id %in% patients)%>%
    filter(metastases_location == "Breast" | metastases_location == loc)
  
  for (cell in cell_type) {
    metastases_loc <- metastases_loc %>%
      distinct(patient_id, metastases_location, variable, .keep_all = TRUE)
    
    metastases_cell <- metastases_loc %>%
      filter(variable == cell) %>%
      select(metastases_location, value)
    
    pwc <- metastases_cell %>% 
      pairwise_t_test(value ~ metastases_location, pool.sd = FALSE, paired = TRUE, p.adjust.method = "BH")
    pwc <- pwc %>% add_xy_position(x = "metastases_location")
    
    result <- pwc %>% 
      mutate(cell_type = cell)
    
    results <- rbind(results, result)
  }}

#plot the metastases location with breast as baseline of 0 - only samples that did not switch subtypes####
#list of patient_ids who did not switch
metastases <- combined_standard %>%
  select(patient_id, sample_id, study, metastases_location, primary, pam50) %>%
  filter(study %in% c("aurora")) %>%
  filter(metastases_location %in% c("Brain", "Chest", "Liver", "Lung", "Skin", "SoftTissue", "LymphNode"))

primaries <- combined_standard %>%
  select(patient_id, sample_id, study, metastases_location, primary, pam50) %>%
  filter(study %in% c("aurora")) %>%
  filter(metastases_location %in% c("Breast"))

merged_df <- merge(primaries, metastases, by = "patient_id", suffixes = c("_primary", "_met"))
merged_df <- merged_df[, c("patient_id", "pam50_primary", "pam50_met", "metastases_location_met")]

matching_patients <- merged_df %>%
  filter(pam50_primary == pam50_met) %>%
  pull(patient_id)

#get metastases values
non_breast_values <- metastases_filtered %>%
  filter(metastases_location != "Breast")
breast_values <- metastases_filtered %>%
  filter(metastases_location == "Breast")
merged_values <- merge(non_breast_values, breast_values, by = "patient_id")
merged_values <- merged_values %>%
  filter(patient_id %in% matching_patients)

# Subtract values from breast_values columns 3:78 from non_breast_values columns 3:78
for (i in 3:51) {
  merged_values[,i] <- log2((merged_values[,i]/merged_values[,i+50]))
}

merged_values <- merged_values[, c(1:51)]
colnames(merged_values) <- colnames(non_breast_values)
normalized <- melt(merged_values, id.vars = c("patient_id", "metastases_location"))
normalized <- normalized %>%
  dplyr::rename("Cell_State" = "variable")

normalized$Cell_Type[grep("Cancer", normalized$Cell_State)] <- "Cancer Epithelium"
normalized$Cell_Type[grep("Luminal", normalized$Cell_State)] <- "Normal Epithelium"
normalized$Cell_Type[grep("epithelial", normalized$Cell_State)] <- "Normal Epithelium"
normalized$Cell_Type[grep("Plasmablasts", normalized$Cell_State)] <- "Lymphoid"
normalized$Cell_Type[grep("PVL", normalized$Cell_State)] <- "Vasculature"
normalized$Cell_Type[grep("T_cells", normalized$Cell_State)] <- "Lymphoid"
normalized$Cell_Type[grep("B_cells", normalized$Cell_State)] <- "Lymphoid"
normalized$Cell_Type[grep("Myeloid", normalized$Cell_State)] <- "Myeloid"
normalized$Cell_Type[grep("Endothelial", normalized$Cell_State)] <- "Vasculature"
normalized$Cell_Type[grep("CAFs", normalized$Cell_State)] <- "Fibroblasts"

#add paired test p-values to the normalized dataset
results <- subset(results, grepl("Breast", group1) | grepl("Breast", group2))
results$metloc <- ifelse(grepl("Breast", results$group1), results$group2, results$group1)
normalized <- normalized %>%
  merge(results, by.x = c("Cell_State", "metastases_location"), by.y = c("cell_type", "metloc"), all.x = TRUE) %>%
  select(c("patient_id", "metastases_location", "Cell_State", "value", "Cell_Type", "p.adj"))

mlocations <- unique(normalized$metastases_location)
for(location in mlocations){
  normalized_subset <- dplyr::filter(normalized, metastases_location == location)
  normalized_subset <- dplyr::filter(normalized_subset,(p.adj < 0.05))
  
  plot <- ggplot(normalized_subset, aes(x = reorder(Cell_State, value), y = value, fill = Cell_Type)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dotted", color = "red", size = 1.5) +  # Add dotted red line at y = 0
    labs(title = location,
         x = "Cell State",
         y = "Proportion") +
    theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1))
  
  #print(plot)
}

normalized_subset <- dplyr::filter(normalized,(p.adj < 0.05))
normalized_subset <- normalized_subset %>%
  group_by(metastases_location, Cell_State) %>%
  mutate(average_value = mean(value)) %>%
  ungroup()
#filter(abs(average_value) >= 0.00)
normalized_subset <- filter(normalized_subset, metastases_location %in% c("LymphNode"))

x_position <- max(normalized_subset$value)+0.25

plot <- ggplot(normalized_subset, aes(x = value, y = reorder(Cell_State, -value), fill = Cell_Type)) +
  geom_density_ridges(alpha = 0.8, scale = 1.25) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.5) +  # Add dotted red line at x = 0
  labs(title = "Lymph Node Metastases",
       x = "Log(FC)",
       y = "'") +
  #facet_wrap(~metastases_location, scales = "free_x", ncol = 2) +
  theme(axis.text.y = element_text(size = 8))+
  scale_fill_manual(values = c("Cancer Epithelium" = "#d4a373",
                               "Normal Epithelium" = "#BDB2FF",
                               "Vasculature" = "#A0C4FF",
                               "Lymphoid" = "#FFADAD",
                               "Myeloid" = "#FDFFB6",
                               "Fibroblasts" = "#CAFFBF")) +
  theme_minimal()+
  geom_text(aes(x = x_position, 
                label = round(average_value, 2)), 
            hjust = 0, size = 3)+
  coord_cartesian(clip = "off")+
  theme(legend.position = "bottom", 
        plot.margin = margin(5.5, 15, 5.5, 5.5))  

print(plot)
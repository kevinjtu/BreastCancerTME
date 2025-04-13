#libraries####
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggsignif)
library(stringr)
library(ggpmisc)

#data input####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
cell_types <- fread("Clinical Associations with Cell Types.txt")
cell_types <- filter(cell_types, primary == TRUE)
cell_types <- filter(cell_types, platform == "rnaseq")
cell_types <- filter(cell_types, study != "MBC")
cell_types <- filter(cell_types, !is.na(intclust) & !is.na(pam50) & !is.na(ER_status) & !is.na(HER2_status))
cell_types <- select(cell_types, c(1,4,5,10,11,25:33))
cell_types$intclust <- gsub("4ER\\+|4ER-", "4", cell_types$intclust)

cell_types$tme_total <- 1 - cell_types$"Cancer Epithelial"
columns_to_update <- c("Endothelial", "CAFs", "PVL", "B-cells", "T-cells", "Myeloid", "Normal Epithelial", "DC")
for (col in columns_to_update) {
  cell_types[[col]] <- cell_types[[col]] / cell_types$tme_total
}
cell_types$tme_total <- NULL

cell_types_long <- reshape2::melt(cell_types, id.vars = c("patient_id", "ER_status", "HER2_status", "intclust", "pam50", "Cancer Epithelial"))
cell_types_long <- rename(cell_types_long, Cancer_Epithelial = "Cancer Epithelial")

#make stacked barchart plot separated by intclust####
cell_types_long <- cell_types_long %>%
  arrange(intclust)

# Extract the ordered patient IDs
ordered_data <- cell_types_long %>%
  arrange(as.numeric(intclust), Cancer_Epithelial)
ordered_patient_ids <- unique(ordered_data$patient_id)

#color pallette
intclust_colors <- c(
  "1" = "#EB411E",
  "2" = "#59B05E",
  "3" = "#BA1D65",
  "4" = "#2BAFAA",
  "5" = "#770714",
  "6" = "#F4E83B",
  "7" = "#25348A",
  "8" = "#F49D19",
  "9" = "#C570AD",
  "10" = "#8c62bd"
)

pam50_colors <- c(
  "Basal" = "#E41A1C",  
  "Her2" = "#FB9A99",       
  "LumA" = "#1F78B4",         
  "LumB" = "#A6CEE3",         
  "Normal" = "#66A61E")

boolean <- c(
  "pos" = "#6d98ba", 
  "neg" = "#e98472"
)

# intclust tile plot
intclust <- ggplot(ordered_data, aes(x = factor(patient_id, levels = ordered_patient_ids), y = factor(1), fill = intclust)) +
  geom_tile() +
  labs(y = "IntClust", fill = "IntClust") +
  scale_fill_manual(values = intclust_colors, breaks = 1:10) +  # Specify breaks from 1 to 10
  theme_void()

# pam50 tile plot
pam <- ggplot(ordered_data, aes(x = factor(patient_id, levels = ordered_patient_ids), y = factor(1), fill = pam50)) +
  geom_tile() +
  labs(y = "Pam50", fill = "Pam50") +
  scale_fill_manual(values = pam50_colors) +
  theme_void()

# ER status tile plot
er <- ggplot(ordered_data, aes(x = factor(patient_id, levels = ordered_patient_ids), y = factor(1), fill = ER_status)) +
  geom_tile() +
  labs(y = "ER status") +
  scale_fill_manual(values = boolean) +
  theme_void()+
  labs(fill = "ER/HER2 status")
  
# her tile plot
her <- ggplot(ordered_data, aes(x = factor(patient_id, levels = ordered_patient_ids), y = factor(1), fill = HER2_status)) +
  geom_tile() +
  labs(y = "HER2 status") +
  scale_fill_manual(values = boolean) +
  theme_void()+
  guides(fill = "none")

# Create stacked cell_type plot
desired_order <- c("DC", "B-cells", "T-cells", "Myeloid", "PVL", "Endothelial", "CAFs", "Normal Epithelial", "Cancer Epithelial")
set1_palette <- brewer.pal(9, "Set1")
bar_plot <- ggplot(ordered_data, aes(x = factor(patient_id, levels = ordered_patient_ids), y = value, fill = factor(variable, levels = desired_order))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Proportion") +
  scale_fill_manual(values = set1_palette, breaks = desired_order) +  # Set the desired order of variables in legend
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 16)  
  ) +
  labs(fill = "Cell Type")

# cellularity tile plot
cancer_epithelial_gradient <- ggplot(ordered_data, aes(x = factor(patient_id, levels = ordered_patient_ids), y = factor(1), fill = Cancer_Epithelial)) +
  geom_tile() +
  labs(y = "Cancer Epithelial") +
  scale_fill_gradient(low = "white", high = "black") +  # Black to white gradient
  theme_void() +
  guides(fill = guide_colorbar(title = "Cellularity"))

# Arrange plots side by side
combined_plot <- intclust + pam + er + her + bar_plot + cancer_epithelial_gradient +
  plot_layout(ncol = 1, heights = c(0.02, 0.02, 0.02, 0.02, 0.9, 0.02), guides = "collect")
combined_plot


#make boxplot####

cell_types_long$intclust <- factor(cell_types_long$intclust, levels = sort(as.numeric(unique(cell_types_long$intclust))))
gg <- ggplot(cell_types_long, aes(x = intclust, y = value, fill = intclust)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_y") +
  scale_fill_manual(values = intclust_colors) +
  labs(x = "Intclust", y = "Proportion", title = "Cell Types in IntClusters") +
  theme_minimal()

#t test
perform_t_tests <- function(data) {
  intclust_combinations <- combn(unique(data$intclust), 2)
  results <- matrix(NA, nrow = ncol(intclust_combinations), ncol = 4,
                    dimnames = list(NULL, c("Intclust1", "Intclust2", "Variable", "p_value")))
  for (i in seq_len(ncol(intclust_combinations))) {
    intclust1 <- intclust_combinations[1, i]
    intclust2 <- intclust_combinations[2, i]
    for (variable in unique(data$variable)) {
      subset_data <- data %>%
        filter(intclust %in% c(intclust1, intclust2), variable == variable)
      t_test <- t.test(value ~ intclust, data = subset_data)
      results[i, ] <- c(intclust1, intclust2, variable, t_test$p.value)
    }
  }
  return(results)
}

# Perform t-tests for each variable
t_test_results <- cell_types_long %>%
  group_by(variable) %>%
  do(data.frame(perform_t_tests(.))) %>%
  ungroup()

t_test_results$p_value <- as.numeric(t_test_results$p_value)

#bonferroni adjustment
t_test_results$adjusted_p <- p.adjust(t_test_results$p_value, method = "bonferroni")


gg <- ggplot(cell_types_long, aes(x = intclust, y = value)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "Intclust", y = "Proportion", title = "Cell Types in IntClusters")

# Filter the t-test results to only include significant comparisons after Bonferroni adjustment
significant_results <- t_test_results %>%
  filter(adjusted_p < 0.05)

#cohen's d
calculate_cohens_d <- function(x, y) {
  mean_diff <- mean(x) - mean(y)
  pooled_sd <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / (length(x) + length(y) - 2))
  cohens_d <- mean_diff / pooled_sd
  return(cohens_d)
}

# Function to perform pairwise comparisons using Cohen's d
perform_cohens_d <- function(data) {
  intclust_combinations <- combn(unique(data$intclust), 2)
  results <- matrix(NA, nrow = ncol(intclust_combinations), ncol = 4,
                    dimnames = list(NULL, c("Intclust1", "Intclust2", "Variable", "Cohens_d")))
  for (i in seq_len(ncol(intclust_combinations))) {
    intclust1 <- intclust_combinations[1, i]
    intclust2 <- intclust_combinations[2, i]
    for (variable in unique(data$variable)) {
      subset_data <- data %>%
        filter(intclust %in% c(intclust1, intclust2), variable == variable)
      group1 <- subset_data$value[subset_data$intclust == intclust1]
      group2 <- subset_data$value[subset_data$intclust == intclust2]
      cohens_d <- calculate_cohens_d(group1, group2)
      results[i, ] <- c(intclust1, intclust2, variable, cohens_d)
    }
  }
  return(results)
}

# Perform pairwise comparisons using Cohen's d for each variable
cohens_d_results <- cell_types_long %>%
  group_by(variable) %>%
  do(data.frame(perform_cohens_d(.))) %>%
  ungroup()

# Convert Cohen's d to numeric type for proper handling
cohens_d_results$Cohens_d <- as.numeric(cohens_d_results$Cohens_d)

# Display the first few rows of the results
head(cohens_d_results)


significant_results <- cohens_d_results %>%
  filter(abs(Cohens_d) > 0.80)

# Create the ggplot object with the specified colors
filtered_data <- cell_types_long %>%
  filter(variable %in% c("Myeloid", "CAFs", "Normal Epithelial"))

# Create the ggplot object for the filtered data
gg <- ggplot(filtered_data, aes(x = intclust, y = value, fill = intclust)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "fixed") +
  scale_fill_manual(values = intclust_colors) +
  labs(x = "Intclust", y = "Proportion", title = "Cell Types in IntClusters (Filtered)") +
  theme_minimal() +
  theme(
    panel.spacing = unit(2, "lines"), # Space between facets
    plot.title = element_text(size = 18), # Title text size
    axis.title.x = element_text(size = 14), # X-axis title text size
    axis.title.y = element_text(size = 14), # Y-axis title text size
    axis.text = element_text(size = 12), # Axis text size
    strip.text = element_text(size = 14) # Facet label text size
  )

# Display the plot
print(gg)


#test myeloid effect size when stratifying ER- and ER+, use this block directly after "#make stacked barchart plot" block####
# Update `intclust` to include "4ER-" and "4ER+"
cell_types_long$intclust <- ifelse(
  cell_types_long$intclust == 4 & cell_types_long$ER_status == "neg", "4ER-",
  ifelse(cell_types_long$intclust == 4 & cell_types_long$ER_status == "pos", "4ER+",
         cell_types_long$intclust))

cell_types_long$intclust <- factor(
  cell_types_long$intclust,
  levels = c("1", "2", "3", "4ER+", "4ER-", "5", "6", "7", "8", "9", "10")
)

ic4_erneg_myeloid <- cell_types_long %>%
  filter(intclust == "4ER-", variable == "Myeloid")

ic4_erpos_myeloid <- cell_types_long %>%
  filter(intclust == "4ER+", variable == "Myeloid")

ic10_myeloid <- cell_types_long %>%
  filter(intclust == "10", variable == "Myeloid")

# Compute Cohen's d for the `value` column
calculate_cohens_d <- function(x, y) {
  mean_diff <- mean(x) - mean(y)
  pooled_sd <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / 
                      (length(x) + length(y) - 2))
  cohens_d <- mean_diff / pooled_sd
  return(cohens_d)
}

# Compute Cohen's d for the `value` column
ic4neg_ic10 <- calculate_cohens_d(
  ic4_erneg_myeloid$value,
  ic10_myeloid$value
)
ic4pos_ic10 <- calculate_cohens_d(
  ic4_erpos_myeloid$value,
  ic10_myeloid$value
)
ic4neg_ic10
ic4pos_ic10

#plot with ER-
gg <- ggplot(filtered_data, aes(x = intclust, y = value, fill = intclust)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "fixed") +
  scale_fill_manual(values = intclust_colors) +
  labs(x = "Intclust", y = "Proportion", title = "Cell Types in IntClusters (Filtered)") +
  theme_minimal() +
  theme(
    panel.spacing = unit(2, "lines"), # Space between facets
    plot.title = element_text(size = 18), # Title text size
    axis.title.x = element_text(size = 14), # X-axis title text size
    axis.title.y = element_text(size = 14), # Y-axis title text size
    axis.text = element_text(size = 12), # Axis text size
    strip.text = element_text(size = 14) # Facet label text size
  )

#make boxplot separated by ER/HER2 subtypes####
cell_types_long <- cell_types_long %>%
  mutate(ER_HER2_Subtype = case_when(
    ER_status == "pos" & HER2_status == "neg" ~ "ER+/HER2-",
    ER_status == "pos" & HER2_status == "pos" ~ "ER+/HER2+",
    ER_status == "neg" & HER2_status == "pos" ~ "ER-/HER2+",
    ER_status == "neg" & HER2_status == "neg" ~ "ER-/HER2-",
    TRUE ~ NA_character_
  ))

cell_types_long <- cell_types_long %>%
  arrange(ER_HER2_Subtype)

# Create stacked cell_type plot, stratefiied by tumor cellularity
desired_order <- c("DC", "B-cells", "T-cells", "Myeloid", "PVL", "Endothelial", "CAFs", "Normal Epithelial")
set1_palette <- brewer.pal(9, "Set1")

ordered_data <- cell_types_long %>%
  arrange(ER_HER2_Subtype, Cancer_Epithelial)
ordered_data$variable <- factor(ordered_data$variable, levels = desired_order)

ordered_data$Cellularity <- cut(
  ordered_data$Cancer_Epithelial * 100,
  breaks = c(0, 25, 50, 75, 100),
  labels = c("0-25%", "25-50%", "50-75%", "75-100%"),
  include.lowest = TRUE
)

ordered_data$ER_HER2_Subtype <- factor(ordered_data$ER_HER2_Subtype)

ggplot(ordered_data, aes(x = Cellularity, y = value, color = ER_HER2_Subtype, group = ER_HER2_Subtype)) +
  stat_summary(
    fun = mean, 
    geom = "point", 
    size = 3, 
    shape = 21,
    fill = "white",
  ) +
  stat_summary(
    fun.data = ~ data.frame(
      y = mean(.),
      ymin = mean(.) - sd(.) / sqrt(length(.)),
      ymax = mean(.) + sd(.) / sqrt(length(.))
    ),
    geom = "errorbar",
    width = 0.2,
  ) +
  stat_summary(
    aes(group = ER_HER2_Subtype), 
    fun = mean, 
    geom = "line"
  ) +
  facet_wrap(~variable, nrow = 2, scales = "free_y") +
  scale_color_manual(values = c("#E41A1C", "#FB9A99", "#1F78B4", "#66A61E")) +
  theme_bw() +
  labs(x = "Tumor Cellularity", y = "Proportion of TME", color = "Subtype") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

###t tests with sig labels on plot DEPRECATED
#unique_variables <- unique(ordered_data$variable)
#plot_list <- list()
#for (var in unique_variables) {
#  p <- ggplot(ordered_data[ordered_data$variable == var, ], aes(x = ER_HER2_Subtype, y = value, fill = cellularity)) +
#    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, notch = TRUE) +
 #   scale_fill_manual(values = c("grey90", "grey70", "grey50", "grey30")) +
#    scale_color_manual(values = c("grey90", "grey70", "grey50", "grey30")) +
#    theme_bw() +
#    theme(
#      axis.text.x = element_text(angle = 45, hjust = 1),
#      axis.title = element_text(size = 14),
#      axis.text = element_text(size = 12),
##      legend.text = element_text(size = 12),
#      legend.title = element_text(size = 12),
##      plot.title = element_text(size = 16)
#   ) +
  #  geom_signif(
#      comparisons = list(
#        c("ER", "ER_HER"),
#        c("ER", "HER"),
#        c("ER", "TNBC"),
#        c("ER_HER", "HER"),
#        c("HER", "TNBC"),
#        c("ER_HER", "TNBC")
#      ),
#      map_signif_level = TRUE,
#      step_increase = 0.1
#    ) +
#    ggtitle(paste(var))+
#    labs(x = "", y = "Proportion")
  
  # Store the plot in the list
#  plot_list[[var]] <- p
#}


#(grobs = plot_list, nrow = 2)
#make boxplot separated by pam50 subtypes####
desired_order <- c("DC", "B-cells", "T-cells", "Myeloid", "PVL", "Endothelial", "CAFs", "Normal Epithelial")

ordered_data <- cell_types_long %>%
  arrange(ER_HER2_Subtype, Cancer_Epithelial) %>%
  filter(pam50 != "Claudin")
ordered_data$variable <- factor(ordered_data$variable, levels = desired_order)

ordered_data$Cellularity <- cut(
  ordered_data$Cancer_Epithelial * 100,
  breaks = c(0, 25, 50, 75, 100),
  labels = c("0-25%", "25-50%", "50-75%", "75-100%"),
  include.lowest = TRUE
)

ggplot(ordered_data, aes(x = Cellularity, y = value, color = pam50, group = pam50)) +
  stat_summary(
    fun = mean, 
    geom = "point", 
    size = 3, 
    shape = 21,
    fill = "white",
  ) +
  stat_summary(
    fun.data = ~ data.frame(
      y = mean(.),
      ymin = mean(.) - sd(.) / sqrt(length(.)),
      ymax = mean(.) + sd(.) / sqrt(length(.))
    ),
    geom = "errorbar",
    width = 0.2,
  ) +
  stat_summary(
    aes(group = pam50), 
    fun = mean, 
    geom = "line"
  ) +
  facet_wrap(~variable, nrow = 2, scales = "free_y") +
  scale_color_manual(values = c("#E41A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#66A61E")) +
  theme_bw() +
  labs(x = "Tumor Cellularity", y = "Proportion of TME", color = "Pam50 Subtype") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

#T cell state increases in ER+/HER2- tumors####
# Load the cell type data
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
print(names(combined_standard))
new_names <- names(combined_standard)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "T-regs", "Tregs")  # Replace T-regs with Tregs
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_standard <- setNames(combined_standard, new_names)
print(names(combined_standard))

combined_standard <- combined_standard %>% mutate_at(vars(12:19), as.numeric) #set survival data to numeric
cox_results <- data.frame(Cell_Type = character(),
                          Variable = character(),
                          Hazard_Ratio = numeric(),
                          CI_Lower = numeric(),
                          CI_Upper = numeric(),
                          P_Value = numeric(),
                          stringsAsFactors = FALSE)

cancer_columns <- combined_standard[, grep("Cancer", names(combined_standard)), with = FALSE]
cancer_sum <- rowSums(cancer_columns, na.rm = TRUE)
combined_standard$Cancer_Epithelial <- cancer_sum
cols_to_normalize <- names(combined_standard)[25:68]
combined_standard[, cols_to_normalize] <- lapply(combined_standard[, ..cols_to_normalize, drop = FALSE], function(x) {
  x / (1 - cancer_sum)
})
rowSums(combined_standard[, 25:68], na.rm = TRUE)
combined_standard <- combined_standard[,-c(69:73)]

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
combined_standard_new[, 27:63] <- lapply(combined_standard_new[, 27:63], function(x) as.numeric(as.character(x)))
rowSums(combined_standard_new[, 27:63])
combined_standard <- combined_standard_new

combined_standard <- combined_standard %>%
  mutate(ER_HER2_Subtype = case_when(
    ER_status == "pos" & HER2_status == "neg" ~ "ER+/HER2-",
    ER_status == "pos" & HER2_status == "pos" ~ "ER+/HER2+",
    ER_status == "neg" & HER2_status == "pos" ~ "ER-/HER2+",
    ER_status == "neg" & HER2_status == "neg" ~ "ER-/HER2-",
    TRUE ~ NA_character_
  ))

#select just the t cell columns
combined_standard_tcell <- combined_standard %>%
  select(patient_id, pam50, ER_HER2_Subtype, Cancer_Epithelial, `Chemokine Expressing CD8 T Cells`, 
         `Type 1 Interferon T Cells`, `Effector Cytotoxic CD8 T Cells`, 
         `Activated Cytotoxic CD8 T Cells`, `Inhibitory CD8 T Cells`, 
         `Central Memory CD4 T Cells`, `Th1 Cells`, `Regulatory T Cells`, 
         `T Follicular Helper Cells`, `Natural killer T Cells`, 
         `Proliferating T Cells`)

combined_standard_tcell$Cellularity <- cut(
  combined_standard_tcell$Cancer_Epithelial * 100,
  breaks = c(0, 25, 50, 75, 100),
  labels = c("0-25%", "25-50%", "50-75%", "75-100%"),
  include.lowest = TRUE
)

combined_standard_tcell_long <- melt(combined_standard_tcell, id.vars = c("patient_id", "pam50", "Cellularity", "ER_HER2_Subtype", "Cancer_Epithelial"))

#scatterplot with lines
#combined_standard_tcell_long <- combined_standard_tcell_long %>%
# filter(!is.na(pam50) & pam50 != "Claudin" & pam50 != "")

combined_standard_tcell_long <- combined_standard_tcell_long %>%
  filter(!is.na(ER_HER2_Subtype))

ggplot(combined_standard_tcell_long, aes(x = Cellularity, y = value, color = ER_HER2_Subtype, group = ER_HER2_Subtype)) +
  stat_summary(
    fun = mean, 
    geom = "point", 
    size = 3, 
    shape = 21,
    fill = "white",
  ) +
  stat_summary(
    fun.data = ~ data.frame(
      y = mean(.),
      ymin = mean(.) - sd(.) / sqrt(length(.)),
      ymax = mean(.) + sd(.) / sqrt(length(.))
    ),
    geom = "errorbar",
    width = 0.2,
  ) +
  stat_summary(
    aes(group = ER_HER2_Subtype), 
    fun = mean, 
    geom = "line"
  ) +
  facet_wrap(~variable, nrow = 2, scales = "free_y") +
  scale_color_manual(values = c("#E41A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#66A61E")) +
  theme_bw() +
  labs(x = "Tumor Cellularity", y = "Proportion of TME", color = "Pam50 Subtype") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

#T cell states associated with celarultin in LumA
combined_standard_tcell_long_subset <- combined_standard_tcell_long %>%
  filter(variable %in% c("Type 1 Interferon T Cells", "Regulatory T Cells"))

ggplot(combined_standard_tcell_long_subset, aes(x = Cellularity, y = value, color = ER_HER2_Subtype, group = ER_HER2_Subtype)) +
  stat_summary(
    fun = mean, 
    geom = "point", 
    size = 3, 
    shape = 21,
    fill = "white",
  ) +
  stat_summary(
    fun.data = ~ data.frame(
      y = mean(.),
      ymin = mean(.) - sd(.) / sqrt(length(.)),
      ymax = mean(.) + sd(.) / sqrt(length(.))
    ),
    geom = "errorbar",
    width = 0.2,
  ) +
  stat_summary(
    aes(group = ER_HER2_Subtype), 
    fun = mean, 
    geom = "line"
  ) +
  facet_wrap(~variable, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("#E41A1C", "#FB9A99", "#1F78B4", "#66A61E")) +
  theme_bw() +
  labs(x = "Tumor Cellularity", y = "Proportion of TME", color = "Subtype") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )
#Endothelial state increases in ER-/HER2- tumors####
# Load the cell type data
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
print(names(combined_standard))
new_names <- names(combined_standard)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "T-regs", "Tregs")  # Replace T-regs with Tregs
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_standard <- setNames(combined_standard, new_names)
print(names(combined_standard))

combined_standard <- combined_standard %>% mutate_at(vars(12:19), as.numeric) #set survival data to numeric
cox_results <- data.frame(Cell_Type = character(),
                          Variable = character(),
                          Hazard_Ratio = numeric(),
                          CI_Lower = numeric(),
                          CI_Upper = numeric(),
                          P_Value = numeric(),
                          stringsAsFactors = FALSE)

cancer_columns <- combined_standard[, grep("Cancer", names(combined_standard)), with = FALSE]
cancer_sum <- rowSums(cancer_columns, na.rm = TRUE)
combined_standard$Cancer_Epithelial <- cancer_sum
cols_to_normalize <- names(combined_standard)[25:68]
combined_standard[, cols_to_normalize] <- lapply(combined_standard[, ..cols_to_normalize, drop = FALSE], function(x) {
  x / (1 - cancer_sum)
})
rowSums(combined_standard[, 25:68], na.rm = TRUE)
combined_standard <- combined_standard[,-c(69:73)]

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
combined_standard_new[, 27:63] <- lapply(combined_standard_new[, 27:63], function(x) as.numeric(as.character(x)))
rowSums(combined_standard_new[, 27:63])
combined_standard <- combined_standard_new

combined_standard <- combined_standard %>%
  mutate(ER_HER2_Subtype = case_when(
    ER_status == "pos" & HER2_status == "neg" ~ "ER+/HER2-",
    ER_status == "pos" & HER2_status == "pos" ~ "ER+/HER2+",
    ER_status == "neg" & HER2_status == "pos" ~ "ER-/HER2+",
    ER_status == "neg" & HER2_status == "neg" ~ "ER-/HER2-",
    TRUE ~ NA_character_
  ))

#select just the vessel columns
combined_standard_vessel <- combined_standard %>%
  select(patient_id, pam50, ER_HER2_Subtype, Cancer_Epithelial,
         `Venous Vessel`, `Tip-like Vessel`, `Lymphatic Vessel`)

combined_standard_vessel$Cellularity <- cut(
  combined_standard_vessel$Cancer_Epithelial * 100,
  breaks = c(0, 25, 50, 75, 100),
  labels = c("0-25%", "25-50%", "50-75%", "75-100%"),
  include.lowest = TRUE
)

combined_standard_vessel_long <- melt(combined_standard_vessel, 
                                     id.vars = c("patient_id", "pam50", "Cellularity", 
                                                 "ER_HER2_Subtype", "Cancer_Epithelial"))

#scatterplot with lines
#combined_standard_tcell_long <- combined_standard_tcell_long %>%
# filter(!is.na(pam50) & pam50 != "Claudin" & pam50 != "")

combined_standard_vessel_long <- combined_standard_vessel_long %>%
  filter(!is.na(ER_HER2_Subtype))

ggplot(combined_standard_vessel_long, aes(x = Cellularity, y = value, color = ER_HER2_Subtype, group = ER_HER2_Subtype)) +
  stat_summary(
    fun = mean, 
    geom = "point", 
    size = 3, 
    shape = 21,
    fill = "white",
  ) +
  stat_summary(
    fun.data = ~ data.frame(
      y = mean(.),
      ymin = mean(.) - sd(.) / sqrt(length(.)),
      ymax = mean(.) + sd(.) / sqrt(length(.))
    ),
    geom = "errorbar",
    width = 0.2,
  ) +
  stat_summary(
    aes(group = ER_HER2_Subtype), 
    fun = mean, 
    geom = "line"
  ) +
  facet_wrap(~variable, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("#E41A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#66A61E")) +
  theme_bw() +
  labs(x = "Tumor Cellularity", y = "Proportion of TME", color = "Pam50 Subtype") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

#Vessel states associated with celarultin in LumA
combined_standard_vessel_long_subset <- combined_standard_vessel_long %>%
  filter(variable %in% c("Tip-like Vessel", "Lymphatic Vessel"))

ggplot(combined_standard_vessel_long_subset, aes(x = Cellularity, y = value, color = ER_HER2_Subtype, group = ER_HER2_Subtype)) +
  stat_summary(
    fun = mean, 
    geom = "point", 
    size = 3, 
    shape = 21,
    fill = "white",
  ) +
  stat_summary(
    fun.data = ~ data.frame(
      y = mean(.),
      ymin = mean(.) - sd(.) / sqrt(length(.)),
      ymax = mean(.) + sd(.) / sqrt(length(.))
    ),
    geom = "errorbar",
    width = 0.2,
  ) +
  stat_summary(
    aes(group = ER_HER2_Subtype), 
    fun = mean, 
    geom = "line"
  ) +
  facet_wrap(~variable, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("#E41A1C", "#FB9A99", "#1F78B4", "#66A61E")) +
  theme_bw() +
  labs(x = "Tumor Cellularity", y = "Proportion of TME", color = "Subtype") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )
#B Cells state increases in ER- tumors####
# Load the cell type data
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
print(names(combined_standard))
new_names <- names(combined_standard)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "T-regs", "Tregs")  # Replace T-regs with Tregs
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_standard <- setNames(combined_standard, new_names)
print(names(combined_standard))

combined_standard <- combined_standard %>% mutate_at(vars(12:19), as.numeric) #set survival data to numeric
cox_results <- data.frame(Cell_Type = character(),
                          Variable = character(),
                          Hazard_Ratio = numeric(),
                          CI_Lower = numeric(),
                          CI_Upper = numeric(),
                          P_Value = numeric(),
                          stringsAsFactors = FALSE)

cancer_columns <- combined_standard[, grep("Cancer", names(combined_standard)), with = FALSE]
cancer_sum <- rowSums(cancer_columns, na.rm = TRUE)
combined_standard$Cancer_Epithelial <- cancer_sum
cols_to_normalize <- names(combined_standard)[25:68]
combined_standard[, cols_to_normalize] <- lapply(combined_standard[, ..cols_to_normalize, drop = FALSE], function(x) {
  x / (1 - cancer_sum)
})
rowSums(combined_standard[, 25:68], na.rm = TRUE)
combined_standard <- combined_standard[,-c(69:73)]

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
combined_standard_new[, 27:63] <- lapply(combined_standard_new[, 27:63], function(x) as.numeric(as.character(x)))
rowSums(combined_standard_new[, 27:63])
combined_standard <- combined_standard_new

combined_standard <- combined_standard %>%
  mutate(ER_HER2_Subtype = case_when(
    ER_status == "pos" & HER2_status == "neg" ~ "ER+/HER2-",
    ER_status == "pos" & HER2_status == "pos" ~ "ER+/HER2+",
    ER_status == "neg" & HER2_status == "pos" ~ "ER-/HER2+",
    ER_status == "neg" & HER2_status == "neg" ~ "ER-/HER2-",
    TRUE ~ NA_character_
  ))

#select just the b cell columns
combined_standard_bcell <- combined_standard %>%
  select(patient_id, pam50, ER_HER2_Subtype, Cancer_Epithelial,
         `Memory B Cells`, `Naive B Cells`, `Plasmablast Cells`)

combined_standard_bcell$Cellularity <- cut(
  combined_standard_bcell$Cancer_Epithelial * 100,
  breaks = c(0, 25, 50, 75, 100),
  labels = c("0-25%", "25-50%", "50-75%", "75-100%"),
  include.lowest = TRUE
)

combined_standard_bcell_long <- melt(combined_standard_bcell, 
                                      id.vars = c("patient_id", "pam50", "Cellularity", 
                                                  "ER_HER2_Subtype", "Cancer_Epithelial"))

#scatterplot with lines
#combined_standard_tcell_long <- combined_standard_tcell_long %>%
# filter(!is.na(pam50) & pam50 != "Claudin" & pam50 != "")

combined_standard_bcell_long <- combined_standard_bcell_long %>%
  filter(!is.na(ER_HER2_Subtype))

ggplot(combined_standard_bcell_long, aes(x = Cellularity, y = value, color = ER_HER2_Subtype, group = ER_HER2_Subtype)) +
  stat_summary(
    fun = mean, 
    geom = "point", 
    size = 3, 
    shape = 21,
    fill = "white",
  ) +
  stat_summary(
    fun.data = ~ data.frame(
      y = mean(.),
      ymin = mean(.) - sd(.) / sqrt(length(.)),
      ymax = mean(.) + sd(.) / sqrt(length(.))
    ),
    geom = "errorbar",
    width = 0.2,
  ) +
  stat_summary(
    aes(group = ER_HER2_Subtype), 
    fun = mean, 
    geom = "line"
  ) +
  facet_wrap(~variable, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("#E41A1C", "#FB9A99", "#1F78B4", "#66A61E")) +
  theme_bw() +
  labs(x = "Tumor Cellularity", y = "Proportion of TME", color = "Pam50 Subtype") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

#b cell states associated with celarultin in er-
combined_standard_bcell_long_subset <- combined_standard_bcell_long %>%
  filter(variable %in% c("Plasmablast Cells"))

ggplot(combined_standard_vessel_long_subset, aes(x = Cellularity, y = value, color = ER_HER2_Subtype, group = ER_HER2_Subtype)) +
  stat_summary(
    fun = mean, 
    geom = "point", 
    size = 3, 
    shape = 21,
    fill = "white",
  ) +
  stat_summary(
    fun.data = ~ data.frame(
      y = mean(.),
      ymin = mean(.) - sd(.) / sqrt(length(.)),
      ymax = mean(.) + sd(.) / sqrt(length(.))
    ),
    geom = "errorbar",
    width = 0.2,
  ) +
  stat_summary(
    aes(group = ER_HER2_Subtype), 
    fun = mean, 
    geom = "line"
  ) +
  facet_wrap(~variable, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("#E41A1C", "#FB9A99", "#1F78B4", "#66A61E")) +
  theme_bw() +
  labs(x = "Tumor Cellularity", y = "Proportion of TME", color = "Subtype") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )
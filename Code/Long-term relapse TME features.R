#libraries####
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(survival)
library(gridExtra)
library(ggplot2)
library(rms)
library(forestmodel)
library(stringr)
library(survminer)
library(ggsurvfit)

# Load the metabric cell state data####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
combined_standard <- filter(combined_standard, study == "METABRIC", primary == TRUE)

print(names(combined_standard))
new_names <- names(combined_standard)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "T-regs", "Tregs")  # Replace T-regs with Tregs
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_standard <- setNames(combined_standard, new_names)
print(names(combined_standard))

combined_standard <- combined_standard %>% mutate_at(vars(12:19), as.numeric) #set survival data to numeric

cancer_columns <- combined_standard[, grep("Cancer", names(combined_standard)), with = FALSE]
cancer_sum <- rowSums(cancer_columns, na.rm = TRUE)
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

# Assess long-term relapse candidates in Kaplan Meier Curves####
# Candidates (selected from cox regression model on just long-term relapse samples): 
# ER+: Conventional DC, Myoepithelial Cells, Lipid−Associated Macrophage 2, Monocyte [IL1B], Monocyte [FCGR3A], Monocyte [S100A9], Immature PVL
# ER-: Naive B Cells
cell_types <- names(combined_standard)[28:63]
combined_standard$RFS_time <- combined_standard$RFS_time / 365.25

# Split data by ER status
ERpos_data <- combined_standard %>% filter(ER_status == "pos")
ERneg_data <- combined_standard %>% filter(ER_status == "neg")

for (ct in cell_types) {
  tertiles <- quantile(ERpos_data[[ct]], 
                       probs = c(1/3, 2/3), na.rm = TRUE)
  ERpos_data <- ERpos_data %>%
    mutate(!!paste0(ct, "_tertile") := cut(
      .data[[ct]],
      breaks = c(-Inf, tertiles, Inf),
      labels = c("Low", "Medium", "High"),
      include.lowest = TRUE
    ))
}

for (ct in cell_types) {
  tertiles <- quantile(ERneg_data[[ct]], 
                       probs = c(1/3, 2/3), na.rm = TRUE)
  ERneg_data <- ERneg_data %>%
    mutate(!!paste0(ct, "_tertile") := cut(
      .data[[ct]],
      breaks = c(-Inf, tertiles, Inf),
      labels = c("Low", "Medium", "High"),
      include.lowest = TRUE
    ))
}

# Define cell types for each group
ERpos_cells <- colnames(combined_standard[28:63])
  c("Conventional DC", "Myoepithelial Cells",
                 "Lipid-Associated Macrophage 2", "Monocyte [IL1B]",
                 "Monocyte [FCGR3A]", "Monocyte [S100A9]", "Immature PVL")
ERneg_cells <- "Naive B Cells"

# Function to generate KM plots
generate_km_plots <- function(data, cell_types, er_status) {
  for (cell in cell_types) {
    tert_col <- paste0(cell, "_tertile")
    
    surv_obj <- surv_fit(Surv(RFS_time, RFS_event) ~ get(tert_col),
                         data = data)
    
    plot <- ggsurvplot(
      surv_obj,
      data = data,
      title = paste(er_status, ":", cell),
      pval = TRUE,
      risk.table = TRUE,
      legend.title = "Tertile",
      legend.labs = levels(data[[tert_col]]),
      palette = "nejm",
      ggplot2.component = list(geom_vline(xintercept = 10, 
                                          linetype = "dashed",
                                          color = "black"))
    )
    print(plot)
  }
}

#generate KM plots with p values for long term relapse
generate_km_plots <- function(data, cell_types, er_status) {
  for (cell in cell_types) {
    tert_col <- paste0(cell, "_tertile")
    
    # Original survival object for plotting
    surv_obj <- surv_fit(Surv(RFS_time, RFS_event) ~ get(tert_col), data = data)
    
    # Calculate p-value for time <10
    data_before <- data %>%
      mutate(
        time_before = pmin(RFS_time, 10),
        event_before = ifelse(RFS_time < 10 & RFS_event == 1, 1, 0)
      )
    
    surv_before <- survdiff(Surv(time_before, event_before) ~ get(tert_col), data = data_before)
    p_before <- 1 - pchisq(surv_before$chisq, df = length(surv_before$n) - 1)
    
    # Calculate p-value for time >10 (left-truncated)
    data_after <- data %>%
      dplyr::filter(RFS_time >= 10) %>%
      mutate(
        time_after = RFS_time - 10,
        event_after = RFS_event
      )
    
    if(nrow(data_after) > 0) {
      surv_after <- survdiff(Surv(time_after, event_after) ~ get(tert_col), data = data_after)
      p_after <- 1 - pchisq(surv_after$chisq, df = length(surv_after$n) - 1)
    } else {
      p_after <- NA
    }
    
    # Convert p-values to scientific notation
    p_before_sci <- formatC(p_before, format = "e", digits = 2)
    p_after_sci <- ifelse(!is.na(p_after), formatC(p_after, format = "e", digits = 2), "NA")
    
    # Create plot with both p-values
    plot <- ggsurvplot(
      surv_obj,
      data = ERpos_data,
      title = paste(er_status, ":", cell),
      pval = FALSE,  # Remove original p-value
      risk.table = TRUE,
      legend.title = "Tertile",
      legend.labs = levels(data[[tert_col]]),
      palette = "nejm"
      )
    plot$plot <- plot$plot +
      geom_vline(xintercept = 10, linetype = "dashed", color = "black") +
      annotate("text", x = 3, y = 0.2, label = paste("p =", p_before_sci)) +
      annotate("text", x = 16, y = 0.2, label = paste("p =", p_after_sci))
    print(plot)
  }
}


# Generate plots for ER- cases
generate_km_plots(ERpos_data, ERpos_cells, "ER+")
generate_km_plots(ERneg_data, ERpos_cells, "ER-")


# Depracated Compare IntClust 4ER- vs 10 ####
intclustfour <- combined_standard %>%
  dplyr::filter(ER_status == "neg") %>%
  dplyr::filter(intclust %in% c("4ER-", "10"))

variables <- colnames(intclustfour[, 27:63])
t_test_results <- data.frame(
  Column = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
for (variable in variables) {
  t_test <- t.test(intclustfour[[variable]] ~ intclustfour$intclust) # Perform t-test
  t_test_results <- rbind(
    t_test_results,
    data.frame(Column = variable, p_value = t_test$p.value)
  )
}
t_test_results$p_value_bh <- p.adjust(t_test_results$p_value, method = "BH")

significant_variables <- t_test_results %>%
  filter(p_value_bh < 0.05) %>%
  pull(Column)

intclustfour_long <- melt(intclustfour, id.vars = c("patient_id", "intclust"), measure.vars = significant_variables)
summary_data <- intclustfour_long %>%
  group_by(intclust, variable) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    se_value = sd(value, na.rm = TRUE) / sqrt(n())
  )

# Create the plot
ggplot(summary_data, aes(x = intclust, y = mean_value, group = variable)) +
  geom_point(aes(color = intclust), size = 3) +  # Add dots for means
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.2) +  # Add error bars
  geom_line(aes(group = variable), color = "black") +  # Draw line between points
  facet_wrap(~variable, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Intclust", y = "Mean Value ± SE")

plasmablast_data <- summary_data %>%
  filter(variable == "Naive B Cells")
t_test_result <- t.test(value ~ intclust, data = plasmablast_data)
p_value <- t_test_result$p.value

ggplot(plasmablast_data, aes(x = intclust, y = mean_value, group = variable)) +
  geom_point(aes(color = intclust), size = 3) +  # Add dots for means
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.2) +  # Add error bars
  geom_line(aes(group = variable), color = "black") +  # Draw line between points
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    plot.title = element_text(size = 18, hjust = 0.5),  # Increase plot title size and center it
    legend.text = element_text(size = 14),  # Increase legend text size (if applicable)
    legend.title = element_text(size = 16)  # Increase legend title size (if applicable)
  ) +
  labs(x = "Intclust", y = "Mean Proportion ± SE") +
  ggtitle("Naive B Cells")


# DEPRACATED Identify which cell states are associated with long-term relapse (10+ years)####
longtermrelapse_candidates <- intclustfour %>%
  dplyr::select(patient_id, ER_status, intclust, RFS_event, RFS_time, `Memory B Cells`,
                `Effector Cytotoxic CD8 T Cells`, `Natural killer T Cells`, `Monocyte [IL1B]`, `Lipid-Associated Macrophage 2`)

longtermrelapse_candidates <- longtermrelapse_candidates %>% #turn RFS_time from days to years
  mutate(RFS_time = RFS_time / 365.25)

longtermrelapse_candidates <- longtermrelapse_candidates %>%
  mutate(
    Memory_tertile = ntile(`Memory B Cells`, 2),
    Effector_tertile = ntile(`Effector Cytotoxic CD8 T Cells`, 2),
    NKT_tertile = ntile(`Natural killer T Cells`, 2),
    Monocyte_tertile = ntile(`Monocyte [IL1B]`, 2),
    Lam_tertile = ntile(`Lipid-Associated Macrophage 2`, 2),
  )

#kaplan meier curves
longtermrelapse_candidates_er <- longtermrelapse_candidates %>%
  filter(intclust == "4ER-")
km_fit <- survfit(Surv(RFS_time, RFS_event) ~ Lam_tertile, data = longtermrelapse_candidates_er)
logrank <- survdiff(Surv(RFS_time, RFS_event) ~ Lam_tertile, data = longtermrelapse_candidates_er)
survfit2(Surv(RFS_time, RFS_event) ~ Lam_tertile, 
                 data = longtermrelapse_candidates_er) %>% 
  ggsurvfit() +
  labs(x = "Years", y = "RFS probability", color = "Tertile") +
  geom_text(x = Inf, y = Inf, 
            label = paste("Logrank P =",  format(logrank$pvalue, digits = 2, scientific = TRUE)),
            hjust = 1, vjust = 1, 
            size = 4) +
  ggtitle("Lam_tertile, IntClust 4ER-")+
  theme_bw()+
  theme(
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 16)  
  ) 


#kaplan meier curve for the entire METABRIC dataset, ER-
longtermrelapse_metabric <- combined_standard %>%
  dplyr::select(patient_id, ER_status, intclust, RFS_event, RFS_time, 
                `Monocyte [S100A9]`)

longtermrelapse_metabric <- longtermrelapse_metabric %>% #turn RFS_time from days to years
  mutate(RFS_time = RFS_time / 365.25)

longtermrelapse_metabric <- longtermrelapse_metabric %>%
  mutate(
    mono_tertile = ntile(`Monocyte [S100A9]`, 2),) %>%
  filter(ER_status == "pos")

km_fit <- survfit(Surv(RFS_time, RFS_event) ~ mono_tertile, data = longtermrelapse_metabric)
logrank <- survdiff(Surv(RFS_time, RFS_event) ~ mono_tertile, data = longtermrelapse_metabric)
survfit2(Surv(RFS_time, RFS_event) ~ mono_tertile, 
         data = longtermrelapse_metabric) %>% 
  ggsurvfit() +
  labs(x = "Years", y = "RFS probability", color = "Tertile") +
  geom_text(x = Inf, y = Inf, 
            label = paste("Logrank P =",  format(logrank$pvalue, digits = 2, scientific = TRUE)),
            hjust = 1, vjust = 1, 
            size = 4) +
  ggtitle("mono_tertile, ER+")+
  theme_bw()+
  theme(
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 16)  
  ) 



factors <- c("Memory_tertile", "Effector_tertile", "NKT_tertile", "Monocyte_tertile", "Lam_tertile")
for (factor in factors) {
    # Subset the data based on ER_status
    data_subset <- subset(longtermrelapse_candidates, ER_status == er_status)
    
    # Fit survival model
    km_fit <- survfit(Surv(RFS_time, RFS_event) ~ get(factor), data = data_subset)
    logrank <- survdiff(Surv(RFS_time, RFS_event) ~ get(factor), data = data_subset)
    
    # Create the plot
    plot <- survfit2(Surv(RFS_time, RFS_event) ~ get(factor), 
                     data = data_subset) %>% 
      ggsurvfit() +
      labs(x = "Years", y = "RFS probability", color = factor) +
      geom_text(x = Inf, y = Inf, 
                label = paste("Logrank P =",  format(logrank$pvalue, digits = 2, scientific = TRUE)),
                hjust = 1, vjust = 1, 
                size = 4) +
      ggtitle(paste(factor, er_status)) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 16)  
      )
    print(plot)
  }
# DEPRACATED plasmablasts associations in all intclusts####
longtermrelapse_plasmablast <- combined_standard %>%
  dplyr::select(patient_id, ER_status, intclust, RFS_event, 
                RFS_time, `Plasmablast Cells`) %>%
  #filter(intclust %in% c("1", "6", "9", "2"))

longtermrelapse_plasmablast <- longtermrelapse_plasmablast %>% #turn RFS_time from days to years
  mutate(RFS_time = RFS_time / 365.25)

longtermrelapse_plasmablast <- longtermrelapse_plasmablast %>%
  mutate(
    Plasmablast_tertile = ntile(`Plasmablast Cells`, 3)
  )

# Loop through each level of ER_status 
er_statuses <- c("pos", "neg")
for (er_status in er_statuses) {
    data_subset <- subset(longtermrelapse_plasmablast, ER_status == er_status)
    
    # Fit survival model
    km_fit <- survfit(Surv(RFS_time, RFS_event) ~ Plasmablast_tertile, data = data_subset)
    logrank <- survdiff(Surv(RFS_time, RFS_event) ~ Plasmablast_tertile, data = data_subset)
    
    # Create the plot
    plot <- survfit2(Surv(RFS_time, RFS_event) ~ Plasmablast_tertile, 
                     data = data_subset) %>% 
      ggsurvfit() +
      labs(x = "Years", y = "RFS probability", color = "Plasmablast_tertile") +
      geom_text(x = Inf, y = Inf, 
                label = paste("Logrank P =",  format(logrank$pvalue, digits = 2, scientific = TRUE)),
                hjust = 1, vjust = 1, 
                size = 4) +
      ggtitle("Plasmablast", paste0(er_status)) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 16)  
      )
    print(plot)
}

# DEPRACATED amount of plasmablasts in long vs short term relapse groups####
longtermrelapse_plasmablast <- combined_standard %>%
  dplyr::select(patient_id, ER_status, intclust, RFS_event, 
                RFS_time, `Plasmablast Cells`)

longtermrelapse_plasmablast$relapse <- ifelse(longtermrelapse_plasmablast$intclust %in% c(1, 6, 9, 2), 
                                              "Long-term relapse", 
                                              "Short-term relapse")

longtermrelapse_plasmablast$"Plasmablast Cells" <- as.numeric(longtermrelapse_plasmablast$"Plasmablast Cells")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggsignif)

# Calculate summary statistics: mean and SEM
summary_data <- longtermrelapse_plasmablast %>%
  group_by(relapse) %>%
  summarise(
    mean_value = mean(`Plasmablast Cells`, na.rm = TRUE),
    sem = sd(`Plasmablast Cells`, na.rm = TRUE) / sqrt(n())
  )

# Create the dot plot with error bars
ggplot(longtermrelapse_plasmablast, aes(x = relapse, y = `Plasmablast Cells`, color = relapse)) +
  geom_point(data = summary_data, aes(y = mean_value), size = 4) +  # Add points for the means
  geom_errorbar(data = summary_data, aes(y = mean_value, ymin = mean_value - sem, ymax = mean_value + sem), 
                width = 0.2, size = 0.8) +  # Add error bars
  labs(
    title = "Mean Plasmablast Cells by Relapse Type",
    x = "Relapse Type",
    y = "Plasmablast Cells (± SEM)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  ggsignif::geom_signif(
    comparisons = list(c("Short-term relapse", "Long-term relapse")),  # Specify groups to compare
    test = "t.test"
  ) +
  scale_y_continuous(
    limits = c(0, max(summary_data$mean_value + summary_data$sem) * 1.1),  # Adjust the upper limit
    expand = c(0, 0)  # Remove additional padding
  )


# DEPRACATED memory bcell associations in all intclusts####
longtermrelapse_plasmablast <- combined_standard %>%
  dplyr::select(patient_id, ER_status, intclust, RFS_event, 
                RFS_time, `Memory B Cells`) %>%
  filter(intclust %in% c("1", "6", "9", "2"))

  longtermrelapse_plasmablast <- longtermrelapse_plasmablast %>% #turn RFS_time from days to years
  mutate(RFS_time = RFS_time / 365.25)

longtermrelapse_plasmablast <- longtermrelapse_plasmablast %>%
  mutate(
    Plasmablast_tertile = ntile(`Memory B Cells`, 3)
  )

# Loop through each level of ER_status 
er_statuses <- c("pos", "neg")
for (er_status in er_statuses) {
  data_subset <- subset(longtermrelapse_plasmablast, ER_status == er_status)
  
  # Fit survival model
  km_fit <- survfit(Surv(RFS_time, RFS_event) ~ Plasmablast_tertile, data = data_subset)
  logrank <- survdiff(Surv(RFS_time, RFS_event) ~ Plasmablast_tertile, data = data_subset)
  
  # Create the plot
  plot <- survfit2(Surv(RFS_time, RFS_event) ~ Plasmablast_tertile, 
                   data = data_subset) %>% 
    ggsurvfit() +
    labs(x = "Years", y = "RFS probability", color = "Plasmablast_tertile") +
    geom_text(x = Inf, y = Inf, 
              label = paste("Logrank P =",  format(logrank$pvalue, digits = 2, scientific = TRUE)),
              hjust = 1, vjust = 1, 
              size = 4) +
    ggtitle("Memory B Cell", paste0(er_status)) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14), 
      axis.text = element_text(size = 12), 
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      plot.title = element_text(size = 16)  
    )
  print(plot)
}

# DEPRACATED amount of memory bcell in long vs short term relapse groups####
longtermrelapse_plasmablast <- combined_standard %>%
  dplyr::select(patient_id, ER_status, intclust, RFS_event, 
                RFS_time, `Memory B Cells`) %>%
  filter(ER_status == "pos")

longtermrelapse_plasmablast$relapse <- ifelse(longtermrelapse_plasmablast$intclust %in% c(1, 6, 9, 2), 
                                              "Long-term relapse", 
                                              "Short-term relapse")

longtermrelapse_plasmablast$"Memory B Cells" <- as.numeric(longtermrelapse_plasmablast$"Memory B Cells")

# Calculate summary statistics: mean and SEM
summary_data <- longtermrelapse_plasmablast %>%
  group_by(relapse) %>%
  summarise(
    mean_value = mean(`Memory B Cells`, na.rm = TRUE),
    sem = sd(`Memory B Cells`, na.rm = TRUE) / sqrt(n())
  )

# Create the dot plot with error bars
ggplot(longtermrelapse_plasmablast, aes(x = relapse, y = `Memory B Cells`, color = relapse)) +
  geom_point(data = summary_data, aes(y = mean_value), size = 4) +  # Add points for the means
  geom_errorbar(data = summary_data, aes(y = mean_value, ymin = mean_value - sem, ymax = mean_value + sem), 
                width = 0.2, size = 0.8) +  # Add error bars
  labs(
    title = "Mean Memory B Cells Cells by Relapse Type",
    x = "Relapse Type",
    y = "Memory B Cells Cells (± SEM)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  ggsignif::geom_signif(
    comparisons = list(c("Short-term relapse", "Long-term relapse")),  # Specify groups to compare
    test = "t.test"
  ) +
  scale_y_continuous(
    limits = c(0.05, 0.0525),  # Adjust the upper limit
    expand = c(0, 0)  # Remove additional padding
  )


# DEPRACATED compare long-term vs short-term relapse groups####
combined_standard$RFS_time <- combined_standard$RFS_time / 365.25

# Create long-term and short-term relapse groups
long_term <- combined_standard[!is.na(combined_standard$RFS_time) & combined_standard$RFS_time > 10, ]
short_term <- combined_standard[!is.na(combined_standard$RFS_time) & combined_standard$RFS_time <= 10, ]

# Initialize results dataframe
results <- data.frame(
  Cell_Type = character(),
  P_Value = numeric(),
  Cohen_D = numeric(),
  stringsAsFactors = FALSE
)

# Loop through cell type columns (27 to 63)
for (i in 27:63) {
  col_name <- names(combined_standard)[i]
  
  # Extract data for both groups
  long_data <- long_term[[i]]
  short_data <- short_term[[i]]
  
  # Skip if not enough observations
  if (sum(!is.na(long_data)) < 2 | sum(!is.na(short_data)) < 2) next
  
  # Perform t-test
  t_test <- t.test(long_data, short_data, na.action = na.omit)
  
  # Calculate Cohen's D
  mean_long <- mean(long_data, na.rm = TRUE)
  mean_short <- mean(short_data, na.rm = TRUE)
  sd_long <- sd(long_data, na.rm = TRUE)
  sd_short <- sd(short_data, na.rm = TRUE)
  n_long <- sum(!is.na(long_data))
  n_short <- sum(!is.na(short_data))
  
  pooled_sd <- sqrt(((n_long - 1)*sd_long^2 + (n_short - 1)*sd_short^2) / (n_long + n_short - 2))
  cohen_d <- (mean_long - mean_short) / pooled_sd
  
  # Add to results
  results <- rbind(results, data.frame(
    Cell_Type = col_name,
    P_Value = t_test$p.value,
    Cohen_D = cohen_d
  ))
}

# Filter significant results (p < 0.05) and sort by p-value
significant_results <- results[results$P_Value < 0.05, ]
significant_results <- significant_results[order(significant_results$P_Value), ]

# Optional: Add adjusted p-values (Benjamini-Hochberg)
significant_results$Adjusted_P <- p.adjust(significant_results$P_Value, method = "BH")

# View results
print(significant_results)

# DEPRACATED stratify by pam50####
# Get unique PAM50 subtypes (excluding NA)
pam50_subtypes <- na.omit(unique(combined_standard$ER_status))

# Initialize results dataframe
stratified_results <- data.frame(
  PAM50_Subtype = character(),
  Cell_Type = character(),
  P_Value = numeric(),
  Cohen_D = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each PAM50 subtype
for (subtype in pam50_subtypes) {
  # Subset data for current PAM50 subtype
  subtype_data <- combined_standard[combined_standard$pam50 == subtype & !is.na(combined_standard$pam50), ]
  
  # Split into long/short-term groups
  long_term <- subtype_data[subtype_data$RFS_time > 10 & !is.na(subtype_data$RFS_time), ]
  short_term <- subtype_data[subtype_data$RFS_time <= 10 & !is.na(subtype_data$RFS_time), ]
  
  # Skip subtype if either group has <2 samples
  if (nrow(long_term) < 2 | nrow(short_term) < 2) next
  
  # Loop through cell types (columns 27-63)
  for (i in 27:63) {
    col_name <- names(subtype_data)[i]
    
    # Extract data
    long_data <- long_term[[i]]
    short_data <- short_term[[i]]
    
    # Skip if too many NAs
    if (sum(!is.na(long_data)) < 2 | sum(!is.na(short_data)) < 2) next
    
    # Perform t-test
    t_test <- t.test(long_data, short_data, na.action = na.omit)
    
    # Calculate Cohen's D
    mean_diff <- mean(long_data, na.rm = TRUE) - mean(short_data, na.rm = TRUE)
    pooled_sd <- sqrt(((sum(!is.na(long_data))-1)*var(long_data, na.rm = TRUE) + 
                         (sum(!is.na(short_data))-1)*var(short_data, na.rm = TRUE)) / 
                        (sum(!is.na(c(long_data, short_data))) - 2))
    
    cohen_d <- mean_diff / pooled_sd
    
    # Store results
    stratified_results <- rbind(stratified_results, data.frame(
      PAM50_Subtype = subtype,
      Cell_Type = col_name,
      P_Value = t_test$p.value,
      Cohen_D = cohen_d
    ))
  }
}

# Add adjusted p-values within each PAM50 subtype
stratified_results <- stratified_results %>%
  group_by(PAM50_Subtype) %>%
  mutate(Adjusted_P = p.adjust(P_Value, method = "BH")) %>%
  ungroup()

# Filter significant results (raw p < 0.05)
significant_stratified <- stratified_results %>%
  filter(P_Value < 0.05) %>%
  arrange(PAM50_Subtype, P_Value)

# View results
print(significant_stratified)
# DEPRACATED stratfiy by er/her2####
# Create ER/HER2 status combinations
er_her_combinations <- expand.grid(
  ER = c("pos", "neg"),
  HER2 = c("pos", "neg"),
  stringsAsFactors = FALSE
)

# Initialize results dataframe
er_her_results <- data.frame(
  ER_HER2_Group = character(),
  Cell_Type = character(),
  P_Value = numeric(),
  Cohen_D = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each ER/HER2 combination
for (j in 1:nrow(er_her_combinations)) {
  current_er <- er_her_combinations$ER[j]
  current_her2 <- er_her_combinations$HER2[j]
  
  # Subset data for current ER/HER2 combination
  subgroup_data <- combined_standard[
    combined_standard$ER_status == current_er &
      combined_standard$HER2_status == current_her2 &
      !is.na(combined_standard$ER_status) &
      !is.na(combined_standard$HER2_status),
  ]
  
  # Skip if subgroup has too few samples
  if (nrow(subgroup_data) < 5) next  # Minimum sample size threshold
  
  # Create long/short-term groups
  long_term <- subgroup_data[subgroup_data$RFS_time > 10 & !is.na(subgroup_data$RFS_time), ]
  short_term <- subgroup_data[subgroup_data$RFS_time <= 10 & !is.na(subgroup_data$RFS_time), ]
  
  # Skip if either group has <2 samples
  if (nrow(long_term) < 2 | nrow(short_term) < 2) next
  
  # Loop through cell types (columns 27-63)
  for (i in 27:63) {
    col_name <- names(subgroup_data)[i]
    
    # Extract data
    long_data <- long_term[[i]]
    short_data <- short_term[[i]]
    
    # Skip if insufficient data
    if (sum(!is.na(long_data)) < 2 | sum(!is.na(short_data)) < 2) next
    
    # Perform t-test
    t_test <- t.test(long_data, short_data, na.action = na.omit)
    
    # Calculate Cohen's D
    mean_diff <- mean(long_data, na.rm = TRUE) - mean(short_data, na.rm = TRUE)
    pooled_sd <- sqrt(
      ((length(long_data[!is.na(long_data)]) - 1) * var(long_data, na.rm = TRUE) +
         (length(short_data[!is.na(short_data)]) - 1) * var(short_data, na.rm = TRUE)) /
        (length(c(long_data, short_data)[!is.na(c(long_data, short_data))]) - 2)
    )
    
    cohen_d <- mean_diff / pooled_sd
    
    # Store results
    er_her_results <- rbind(er_her_results, data.frame(
      ER_HER2_Group = paste0("ER", current_er, "_HER2", current_her2),
      Cell_Type = col_name,
      P_Value = t_test$p.value,
      Cohen_D = cohen_d
    ))
  }
}

# Add adjusted p-values within each ER/HER2 group
er_her_results <- er_her_results %>%
  group_by(ER_HER2_Group) %>%
  mutate(Adjusted_P = p.adjust(P_Value, method = "BH")) %>%
  ungroup()

# Filter significant results (raw p < 0.05) and sort
significant_er_her <- er_her_results %>%
  filter(P_Value < 0.05) %>%
  arrange(ER_HER2_Group, P_Value)

# View results
print(significant_er_her)
# DEPRACATED Get unique ER statuses (excluding NA) ####
er_statuses <- na.omit(unique(combined_standard$ER_status))

# Initialize results dataframe
er_results <- data.frame(
  ER_Group = character(),
  Cell_Type = character(),
  P_Value = numeric(),
  Cohen_D = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each ER status
for (er_status in er_statuses) {
  # Subset data for current ER group
  er_data <- combined_standard[
    combined_standard$ER_status == er_status & 
      !is.na(combined_standard$ER_status),
  ]
  
  # Skip if subgroup has too few samples (optional threshold)
  if (nrow(er_data) < 5) next
  
  # Create long/short-term groups
  long_term <- er_data[er_data$RFS_time > 10 & !is.na(er_data$RFS_time), ]
  short_term <- er_data[er_data$RFS_time <= 10 & !is.na(er_data$RFS_time), ]
  
  # Skip if either group has <2 samples
  if (nrow(long_term) < 2 | nrow(short_term) < 2) next
  
  # Loop through cell types (columns 27-63)
  for (i in 27:63) {
    col_name <- names(er_data)[i]
    
    # Extract data
    long_data <- long_term[[i]]
    short_data <- short_term[[i]]
    
    # Skip if insufficient data
    if (sum(!is.na(long_data)) < 2 | sum(!is.na(short_data)) < 2) next
    
    # Perform t-test
    t_test <- t.test(long_data, short_data, na.action = na.omit)
    
    # Calculate Cohen's D
    mean_diff <- mean(long_data, na.rm = TRUE) - mean(short_data, na.rm = TRUE)
    pooled_sd <- sqrt(
      ((sum(!is.na(long_data))-1)*var(long_data, na.rm = TRUE) + 
         ((sum(!is.na(short_data))-1)*var(short_data, na.rm = TRUE)
         ) / (sum(!is.na(c(long_data, short_data))) - 2)))
       
       cohen_d <- mean_diff / pooled_sd
       
       # Store results
       er_results <- rbind(er_results, data.frame(
         ER_Group = er_status,
         Cell_Type = col_name,
         P_Value = t_test$p.value,
         Cohen_D = cohen_d
       ))
  }
}

# Add adjusted p-values within each ER group
er_results <- er_results %>%
  group_by(ER_Group) %>%
  mutate(Adjusted_P = p.adjust(P_Value, method = "BH")) %>%
  ungroup()

# Filter significant results (raw p < 0.05) and sort
significant_er <- er_results %>%
  filter(P_Value < 0.05) %>%
  arrange(ER_Group, P_Value)

# View results
print(significant_er)
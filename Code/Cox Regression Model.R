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
 
surv_metric <- "RFS" # options are RFS, BCSS, MFS, OS

# Load the cell type data####
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

# Cox Regression Models for cell types DEPRACATED####
cell_types <- unique(colnames(combined_standard[,25:68]))
ER_status <- c("pos", "neg")

for (ER_sub in ER_status) {
  print(ER_sub)
  combined_standard_er <- combined_standard[combined_standard$ER_status == ER_sub, ]
  
for (cell in cell_types) {
  print(cell)
  formula <- as.formula(paste("Surv(", paste0(surv_metric, "_time"), ", ", paste0(surv_metric, "_event)"), " ~ ", 
                              cell," + grade + stage + node_status + strata(study)"))
  cox_model <- coxph(formula, data = combined_standard_er, na.action = na.exclude)
  summary_cox <- summary(cox_model)
  hr <- coef(cox_model) # Hazard Ratio
  ci <- confint(cox_model, level = .95) # Confidence Interval
  p_value <- summary_cox$coefficients[, "Pr(>|z|)"] # P-Value
  
  # Store results in the data frame
  variable_names <- names(coef(cox_model)) # Exclude Intercept
  cox_results_cell <- data.frame(ER_status = ER_sub,
                                 Cell_Type = cell,
                                 Variable = variable_names,
                                 Hazard_Ratio = hr,
                                 CI_Lower = ci[, 1],
                                 CI_Upper = ci[, 2],
                                 P_Value = p_value)
  
  cox_results <- rbind(cox_results, cox_results_cell)}
}
rownames(cox_results) <- NULL

#plot cox regression model results for cell types DEPRACATED####
filtered_cox_results <- cox_results %>%
  filter(Cell_Type == Variable)

create_forest_plot <- function(data, er_status) {
  er_data <- data[data$ER_status == er_status, ]
  
  # Reorder levels of Variable based on Hazard_Ratio
  er_data$Variable <- factor(er_data$Variable, levels = er_data$Variable[order(er_data$Hazard_Ratio)])
  
  # Calculate size of points based on p-value
  point_size <- -log10(er_data$P_Value)
  
  # Create forest plot
  p <- ggplot(er_data, aes(x = Hazard_Ratio, y = Variable, fill = Cell_Type)) +
    geom_point(aes(size = point_size), shape = 22) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0) +
    scale_size_continuous(range = c(2, 10)) +
    labs(title = paste("Forest Plot for", er_status, "ER status,", surv_metric)) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  return(p)
}

# Create and display forest plot for ER status "pos"
pos <- create_forest_plot(filtered_cox_results, "pos")
neg <- create_forest_plot(filtered_cox_results, "neg")

# Arrange plots side by side
grid.arrange(pos, neg, ncol = 2)





# cox regression model test at cell state with 4-quantile ####
cell_types <- unique(colnames(combined_standard[,27:63]))
ER_status <- c("pos", "neg")
#combined_standard <- combined_standard %>% #have samples where RFS_time > 10 for long term relapse patient
#  filter(RFS_time > 3652.5)

for (ER_sub in ER_status) {
  print(ER_sub)
  combined_standard_er <- combined_standard[combined_standard$ER_status == ER_sub, ]
  
  for (cell in cell_types) {
    print(cell)
    
    # Create quartiles for the 'cell' variable
    combined_standard_er$cell_quartiles <- cut(combined_standard_er[[cell]], 
                                               breaks = quantile(combined_standard_er[[cell]], probs = seq(0, 1, 0.25), na.rm = TRUE), 
                                               labels = FALSE)
    
    # Convert 'cell_quartiles' to a numeric
    combined_standard_er$cell_quartiles <- as.integer(combined_standard_er$cell_quartiles)
    
    # Define the formula with 'cell_quartiles'
    formula <- as.formula(paste("Surv(", paste0(surv_metric, "_time"), ", ", paste0(surv_metric, "_event)"), 
                                " ~ cell_quartiles + stage + node_status + grade + HER2_status + strata(study)"))
    
    
    # Fit the Cox proportional hazards model
    cox_model <- coxph(formula, data = combined_standard_er)
    summary_cox <- summary(cox_model)
    hr <- coef(cox_model)["cell_quartiles"] # Hazard Ratio
    ci <- confint(cox_model, level = 0.95)["cell_quartiles", ] # Confidence Interval
    p_value <- summary_cox$coefficients["cell_quartiles", "Pr(>|z|)"] # P-Value
    
    # Store results in the data frame
    cox_results_cell <- data.frame(ER_status = ER_sub,
                                   Cell_Type = cell,
                                   Hazard_Ratio = hr,
                                   CI_Lower = ci[1],
                                   CI_Upper = ci[2],
                                   P_Value = p_value)
    
    cox_results <- rbind(cox_results, cox_results_cell)
  }
}

cox_results$P_Value_BH <- p.adjust(cox_results$P_Value, method = "BH")

#plot cox regression model results####
#plot er statuses together
cox_results <- data.table(cox_results)
cox_results <- dplyr::rename(cox_results, c("Cell_State" = "Cell_Type"))
cox_results$Cell_Type[grep("Cancer", cox_results$Cell_State)] <- "Cancer Epithelial"
cox_results$Cell_Type[grep("Luminal", cox_results$Cell_State)] <- "Epithelium"
cox_results$Cell_Type[grep("epithelial", cox_results$Cell_State)] <- "Epithelium"
cox_results$Cell_Type[grep("Plasmablast", cox_results$Cell_State)] <- "Lymphoid"
cox_results$Cell_Type[grep("Th1", cox_results$Cell_State)] <- "Lymphoid"
cox_results$Cell_Type[grep("Follicular", cox_results$Cell_State)] <- "Lymphoid"
cox_results$Cell_Type[grep("Natural Killer", cox_results$Cell_State)] <- "Lymphoid"
cox_results$Cell_Type[grep("PVL", cox_results$Cell_State)] <- "Vascular"
cox_results$Cell_Type[grep("T Cells", cox_results$Cell_State)] <- "Lymphoid"
cox_results$Cell_Type[grep("B Cell", cox_results$Cell_State)] <- "Lymphoid"
cox_results$Cell_Type[grep("Monocyte", cox_results$Cell_State)] <- "Myeloid"
cox_results$Cell_Type[grep("Myeloid", cox_results$Cell_State)] <- "Myeloid"
cox_results$Cell_Type[grep("Macrophage", cox_results$Cell_State)] <- "Myeloid"
cox_results$Cell_Type[grep("DC", cox_results$Cell_State)] <- "Myeloid"
cox_results$Cell_Type[grep("Vessel", cox_results$Cell_State)] <- "Vascular"
cox_results$Cell_Type[grep("CAF", cox_results$Cell_State)] <- "Fibroblast"

#plot er statuses together, ordered
er_pos_data <- cox_results[ER_status == "pos"]
er_neg_data <- cox_results[ER_status == "neg"]
er_pos_data$Cell_State <- factor(er_pos_data$Cell_State, levels = er_pos_data$Cell_State[order(er_pos_data$Hazard_Ratio)])
combined_data <- rbind(er_pos_data, er_neg_data)

# Create forest plot
ggplot(combined_data, aes(x = Hazard_Ratio, y = Cell_State, fill = ER_status)) +
  geom_point(size = 2, shape = 22, position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, position = position_dodge(width = 0.5)) +
  labs(title = paste("Forest Plot for", surv_metric)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  facet_wrap(~Cell_Type, scales = "free_y")+
  guides(fill = guide_legend(override.aes = list(size = c(5, 5))))

#plot univariate model data that is er-specific####
#er+
  combined_data_sig <- combined_data#[combined_data$P_Value < 0.05, ]
  combined_data_er <- filter(combined_data_sig, ER_status == "pos")
  combined_data_er$Hazard_Ratio <- exp(combined_data_er$Hazard_Ratio)
  combined_data_er$CI_Lower <- exp(combined_data_er$CI_Lower)
  combined_data_er$CI_Upper <- exp(combined_data_er$CI_Upper)
  
  ggplot(combined_data_er, aes(x = Hazard_Ratio, y = reorder(Cell_State, Hazard_Ratio), fill = Cell_Type)) +
    geom_point(size = 5, shape = 22, position = position_dodge(width = 0.5)) + 
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, position = position_dodge(width = 0.5)) +
    geom_text(aes(x = max(Hazard_Ratio) + 0.05, 
                  label = paste(round(Hazard_Ratio, 2), " (", round(CI_Lower, 2), ",", round(CI_Upper, 2), ")", sep = "")),
              position = position_dodge(width = 0.5), hjust = 0, size = 3.5,
              color = ifelse(combined_data_er$P_Value_BH > 0.05, "grey", "black")) + 
    labs(title = paste("ER+,", surv_metric), y = "Cell Subtype", fill = "Cell Lineage") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          legend.position = "bottom",  
          plot.margin = unit(c(0.5, 3, 0.5, 0.5), "cm"),
          panel.grid = element_blank()) +  # Remove gridlines
    geom_vline(xintercept = 1, linetype = "solid", color = "black") +
    coord_cartesian(clip = "off") +
    scale_color_manual(values = c("black", "grey"))

#er-
combined_data_er <- filter(combined_data_sig, ER_status == "neg")
combined_data_er <- combined_data_er[order(combined_data_er$Hazard_Ratio), ]
combined_data_er$Hazard_Ratio <- exp(combined_data_er$Hazard_Ratio)
combined_data_er$CI_Lower <- exp(combined_data_er$CI_Lower)
combined_data_er$CI_Upper <- exp(combined_data_er$CI_Upper)

ggplot(combined_data_er, aes(x = Hazard_Ratio, y = reorder(Cell_State, Hazard_Ratio), fill = Cell_Type)) +
  geom_point(size = 5, shape = 22, position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, position = position_dodge(width = 0.5)) +
  geom_text(aes(x = max(Hazard_Ratio) + 0.05, 
                label = paste(round(Hazard_Ratio, 2), " (", round(CI_Lower, 2), ",", round(CI_Upper, 2), ")", sep = "")),
            position = position_dodge(width = 0.5), hjust = 0, size = 3.5,
            color = ifelse(combined_data_er$P_Value_BH > 0.05, "grey", "black")) + 
  labs(title = paste("ER-,", surv_metric), y = "Cell Subtype", fill = "Cell Lineage") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        legend.position = "bottom",  
        plot.margin = unit(c(0.5, 3, 0.5, 0.5), "cm"),
        panel.grid = element_blank()) +  # Remove gridlines
  geom_vline(xintercept = 1, linetype = "solid", color = "black") +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = c("black", "grey"))

#plot cox regression model results for myeloid cells####
#plot er statuses together
cox_results <- data.table(cox_results)

cox_results$Cell_Type[grep("Macrophage", cox_results$Cell_State)] <- "Macrophage"
cox_results$Cell_Type[grep("LAM", cox_results$Cell_State)] <- "Macrophage"
cox_results$Cell_Type[grep("Monocyte", cox_results$Cell_State)] <- "Monocyte"
cox_results$Cell_Type[grep("DC", cox_results$Cell_State)] <- "Dendritic Cell"
cox_results$Cell_Type[grep("T_cell", cox_results$Cell_State)] <- "T Cell"
cox_results$Cell_Type[grep("B_cell", cox_results$Cell_State)] <- "B Cell"
cox_results$Cell_Type[grep("Cycling", cox_results$Cell_State)] <- "Stem Cell"
cox_results$Cell_Type[grep("Endothel", cox_results$Cell_State)] <- "Endothelial"
cox_results$Cell_Type[grep("PVL", cox_results$Cell_State)] <- "Perivascular"

#plot er statuses together, ordered
er_pos_data <- cox_results[ER_status == "pos"]
er_neg_data <- cox_results[ER_status == "neg"]
er_pos_data$Cell_State <- factor(er_pos_data$Cell_State, levels = er_pos_data$Cell_State[order(er_pos_data$Hazard_Ratio)])
combined_data <- rbind(er_pos_data, er_neg_data)
combined_data <- combined_data %>%
  filter(grepl("Normal Epithelium", Cell_Type) | grepl("Dendritic Cell", Cell_Type) | grepl("T Cell", Cell_Type))

# Create forest plot
ggplot(combined_data, aes(x = Hazard_Ratio, y = Cell_State, fill = ER_status)) +
  geom_point(size = 3, shape = 22, position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, position = position_dodge(width = 0.5)) +
  labs(title = paste("Forest Plot for", surv_metric, ", stratified by cell type")) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black")+
  facet_wrap(~Cell_Type, scales = "free_y")+
  guides(fill = guide_legend(override.aes = list(size = c(5, 5))))


#plot adjusted multivariable model (ER-)####
combined_standard_er <- filter(combined_standard, ER_status == "neg")

combined_data_sig <- combined_data[combined_data$P_Value_BH < 0.05, ]
combined_data_er <- filter(combined_data_sig, ER_status == "neg")
multi_variables <- unique(combined_data_er$Cell_State)

cell_types <- colnames(combined_standard[,27:63])

# Function to calculate quantiles
calculate_quantiles <- function(df, cols, quantiles) {
  df %>%
    mutate(across(all_of(cols), ~ ntile(.x, quantiles), .names = "quantile_{col}"))
}

# Apply the quantile calculation
combined_standard_er_quantiles <- calculate_quantiles(combined_standard_er, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  dplyr::select(-all_of(cell_types))

colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
quantile_vars_escaped <- paste0("`", multi_variables, "`")

# Prepare the formula for the Cox model
cox_formula <- as.formula(paste(
  "Surv(", paste0(surv_metric, "_time"), ", ", paste0(surv_metric, "_event)"), " ~",
  paste(quantile_vars_escaped, collapse = " + "),
  "+ stage + node_status + grade + HER2_status + strata(study)"
))

# Fit the Cox proportional hazards model
cox_model <- coxph(cox_formula, data = combined_standard_er_quantiles)

# Extract coefficients and hazard ratios
hr <- exp(cox_model$coefficients)
hr_df <- data.frame(
  variable = names(hr),
  hr = hr
)

hr_df <- hr_df[order(hr_df$hr), ]
hr_df <- hr_df[!(hr_df$variable %in% c("node_status", "stage", "grade", "HER2_statuspos")), ]

hr_df_variables <- hr_df$variable
hr_df_variables <- gsub("`", "", hr_df_variables)
hr_df_variables <- paste0("`", hr_df_variables, "`")

# Prepare the formula for the Cox model
cox_formula <- as.formula(paste(
  "Surv(", paste0(surv_metric, "_time"), ", ", paste0(surv_metric, "_event)"), " ~",
  paste(hr_df_variables, collapse = " + "),
  "+ stage + node_status + grade + HER2_status + strata(study)"
))

# Fit the Cox proportional hazards model
cox_model_reorder_erneg <- coxph(cox_formula, data = combined_standard_er_quantiles)

# Summarize the model
summary <- summary(cox_model_reorder_erneg)

#forest plot
color_map <- function(variable) {
  if (grepl("Myeloid", variable)) {
    return("#FDFFB6")
  } else if (grepl("T_cells", variable) | grepl("B_cells", variable)) {
    return("#FFADAD")
  } else if (grepl("Endothelial", variable) | grepl("PVL", variable)) {
    return("#A0C4FF")
  } else if (grepl("CAFs", variable)) {
    return("#CAFFBF")
  } else if (grepl("Luminal", variable)) {
    return("#BDB2FF")
  } else {
    return("black")  # default color
  }
}

panels <- list(
  list(width = 0.01, display = ~variable, heading = paste("ER-,", surv_metric), hjust = 1),
  list(width = 0.01, item = "vline", hjust = 0.5),
  list(
    width = 0.3, item = "forest", hjust = 0.5, heading = "HR (95% CI)", linetype = "dashed", 
    line_x = 0
  ),
  list(width = 0.03),
  list(width = 0.12, display = ~ ifelse(reference, "Reference", sprintf(
    "%0.2f (%0.2f, %0.2f)",
    trans(estimate), trans(conf.low), trans(conf.high)
  )), display_na = NA),
  list(width = 0.01, item = "vline", hjust = 0.5),
  list(
    width = 0.01,
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 2, eps = 0.001)),
    display_na = NA, hjust = 1, heading = "p"
  ),
  list(width = 0.03)
)

forest_model_order <- forest_model(cox_model_reorder_erneg, covariates = multi_variables, panels)
print(forest_model_order)

#plot adjusted multivariable model (ER+)####
combined_standard_er <- filter(combined_standard, ER_status == "pos")

combined_data_sig <- combined_data[combined_data$P_Value_BH < 0.05, ]
combined_data_er <- filter(combined_data_sig, ER_status == "pos")
multi_variables <- unique(combined_data_er$Cell_State)

cell_types <- colnames(combined_standard[,27:63])

# Function to calculate quantiles
calculate_quantiles <- function(df, cols, quantiles) {
  df %>%
    mutate(across(all_of(cols), ~ ntile(.x, quantiles), .names = "quantile_{col}"))
}

# Apply the quantile calculation
combined_standard_er_quantiles <- calculate_quantiles(combined_standard_er, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  dplyr::select(-all_of(cell_types))

colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
quantile_vars_escaped <- paste0("`", multi_variables, "`")

# Prepare the formula for the Cox model
cox_formula <- as.formula(paste(
  "Surv(", paste0(surv_metric, "_time"), ", ", paste0(surv_metric, "_event)"), " ~",
  paste(quantile_vars_escaped, collapse = " + "),
  "+ stage + node_status + grade + HER2_status + strata(study)"
))

# Fit the Cox proportional hazards model
cox_model <- coxph(cox_formula, data = combined_standard_er_quantiles)

# Extract coefficients and hazard ratios
hr <- exp(cox_model$coefficients)
hr_df <- data.frame(
  variable = names(hr),
  hr = hr
)

hr_df <- hr_df[order(hr_df$hr), ]
hr_df <- hr_df[!(hr_df$variable %in% c("node_status", "stage", "grade", "HER2_statuspos")), ]

hr_df_variables <- hr_df$variable
hr_df_variables <- gsub("`", "", hr_df_variables)
hr_df_variables <- paste0("`", hr_df_variables, "`")

# Prepare the formula for the Cox model
cox_formula <- as.formula(paste(
  "Surv(", paste0(surv_metric, "_time"), ", ", paste0(surv_metric, "_event)"), " ~",
  paste(hr_df_variables, collapse = " + "),
  "+ stage + node_status + grade + HER2_status + strata(study)"
))

# Fit the Cox proportional hazards model
cox_model_reorder_erpos <- coxph(cox_formula, data = combined_standard_er_quantiles)

# Summarize the model
summary <- summary(cox_model_reorder_erpos)

#forest plot
panels <- list(
  list(width = 0.02),
  list(width = 0.01, display = ~variable, heading = paste("ER+,", surv_metric), hjust = 1),
  list(width = 0.01, item = "vline", hjust = 0.5),
  list(
    width = 0.2, item = "forest", hjust = 0.5, heading = "HR (95% CI)", linetype = "dashed", 
    line_x = 0
  ),
  list(width = 0.03),
  list(width = 0.12, display = ~ ifelse(reference, "Reference", sprintf(
    "%0.2f (%0.2f, %0.2f)",
    trans(estimate), trans(conf.low), trans(conf.high)
  )), display_na = NA),
  list(width = 0.01, item = "vline", hjust = 0.5),
  list(
    width = 0.01,
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 3, eps = 0.001)),
    display_na = NA, hjust = 1, heading = "p"
  ),
  list(width = 0.03)
)

forest_model_order <- forest_model(cox_model_reorder_erpos, covariates = multi_variables, panels)
print(forest_model_order)




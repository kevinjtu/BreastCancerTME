#LIBRARY####
library(data.table)
library(Biobase)
library(dplyr)
library(tidyr)
library(textshape)
library(tibble)
library(readxl)
library(ggplot2)
library(rms)
library(DescTools)
library(forestmodel)
library(stringr)
library(vcd)
library(cowplot)

# Load the cell type data####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
combined_standard$pcr <- gsub("NA", NA, combined_standard$pcr)

new_names <- names(combined_standard)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "T-regs", "Tregs")  # Replace T-regs with Tregs
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_standard <- setNames(combined_standard, new_names)
combined_standard <- filter(combined_standard, !is.na(pcr))

or_results <- data.frame(Cell_State = character(),
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

# Function to calculate quantiles
calculate_quantiles <- function(df, cols, quantiles) {
  df %>%
    dplyr::mutate(across(all_of(cols), ~ ntile(.x, quantiles), .names = "quantile_{col}"))
}

#combined_standard <- as.data.frame(combined_standard)
#combined_standard$Endothelial <- rowSums(combined_standard[, grep("Endothelial", names(combined_standard))], na.rm = TRUE)
#combined_standard$CAFs <- rowSums(combined_standard[, grep("CAFs", names(combined_standard))], na.rm = TRUE)
#combined_standard$PVL <- rowSums(combined_standard[, grep("PVL_I|PVL_D", names(combined_standard))], na.rm = TRUE)
#combined_standard$B_cells <- rowSums(combined_standard[, grep("B_cells", names(combined_standard))], na.rm = TRUE)
#combined_standard$CD4_T_cells <- rowSums(combined_standard[, grep("IL7R|CCR7", names(combined_standard))], na.rm = TRUE)
#combined_standard$CD8_T_cells <- rowSums(combined_standard[, grep("CD8", names(combined_standard))], na.rm = TRUE)
#combined_standard$Monocyte <- rowSums(combined_standard[, grep("Monocyte", names(combined_standard))], na.rm = TRUE)
#combined_standard$DC <- rowSums(combined_standard[, grep("DC", names(combined_standard))], na.rm = TRUE)
#combined_standard$Normal_Epithelial <- rowSums(combined_standard[, grep("Luminal|Myoepithelial", names(combined_standard))], na.rm = TRUE)

#rename columns T_cells_c6_IFIT1 -> IFIT_T_Cells
#combined_standard <- combined_standard %>% 
#  rename(T_cells_IFIT = T_cells_c6_IFIT1) %>%
#  rename(Tregs = T_cells_c2_CD4_Tregs_FOXP3) %>%
#  rename(Tfh = T_cells_c3_CD4_Tfh_CXCL13) %>%
#  rename(NK_cells = T_cells_c9_NK_cells_AREG) %>%
#  rename(Cycling_T_cells = T_cells_c11_MKI67) %>%
#  rename(NKT_cells = T_cells_c10_NKT_cells_FCGR3A) %>%
#  rename(M1_Macrophage = Myeloid_c10_Macrophage_1_EGR1) %>%
#  rename(LAM2 = Myeloid_c2_LAM2_APOE) %>%
#  rename(LAM1 = Myeloid_c1_LAM1_FABP5) %>%
#  rename(M2_Macrophage = Myeloid_c9_Macrophage_2_CXCL10) %>%
#  rename(M3_Macrophage = Myeloid_c5_Macrophage_3_SIGLEC1)

#remove columns 25-39, 41-45, 52, 55-56, 59-65 from combined_standard, as these are the old columns beore combining
#combined_standard <- combined_standard[, -c(25:39, 41:45, 52, 55:56, 59:65)]
#print(rowSums(combined_standard[, c(25:38, 41:49)], na.rm = TRUE))
#cols_to_move <- c("treatment_catagory", "platform")
#combined_standard <- combined_standard[, c(setdiff(names(combined_standard), cols_to_move), cols_to_move)]


# or regression model test at cell state with quartiles (CHEMOTHERAPY) ####
cell_types <- unique(colnames(combined_standard[,27:63]))
combined_standard <- combined_standard %>% 
  filter(treatment_catagory == "chemo") 

ER_status <- c("pos", "neg")

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
    
    
    # Fit the linear regression model
    model <- glm(as.integer(pcr) ~ cell_quartiles + age + HER2_status+ study, data = combined_standard_er)
    summary <- summary(model)

    odds_ratio <- exp(coef(model)["cell_quartiles"])
    conf_interval <- exp(confint(model)["cell_quartiles", ])
    p_value <- coef(summary)["cell_quartiles", "Pr(>|t|)"]
    
    # Store results in the data frame
    or_results_cell <- data.frame(ER_status = ER_sub,
                                   Cell_Type = cell,
                                   Odds_Ratio = odds_ratio,
                                   CI_Lower = conf_interval[1],
                                   CI_Upper = conf_interval[2],
                                   P_Value = p_value)
    
    or_results <- rbind(or_results, or_results_cell)
  }
}
or_results$P_Value_BH <- p.adjust(or_results$P_Value, method = "BH")

#plot regression model results####
#assign cell types
cox_results <- data.table(or_results)
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
cox_results_2 <- cox_results

#plot er statuses together, ordered
er_pos_data <- cox_results_2[ER_status == "pos"]
er_neg_data <- cox_results_2[ER_status == "neg"]
er_pos_data$Cell_State <- factor(er_pos_data$Cell_State, levels = er_pos_data$Cell_State[order(er_pos_data$Odds_Ratio)])
combined_data <- rbind(er_pos_data, er_neg_data)

ggplot(combined_data, aes(x = Odds_Ratio, y = Cell_State, color = ER_status)) +
  geom_point(shape = 16) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0) +
  facet_wrap(~Cell_Type) +
  scale_size_continuous(range = c(2, 10)) +
  labs(title = "Forest Plot by ER status") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

#plot er+, sig only
combined_sig <- combined_data[combined_data$P_Value < 1, ]
combined_data_er <- combined_sig[ER_status == "pos", ]
combined_data_er <- combined_data_er %>% arrange(Odds_Ratio)
  
ggplot(combined_data_er, aes(x = Odds_Ratio, y = Cell_State, fill = Cell_Type)) +
  geom_point(size = 5, shape = 22, color = "black") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, col = "black") +
  scale_size_continuous(range = c(2, 10)) +
  geom_text(aes(x = max(Odds_Ratio) + 0.075, 
                label = paste0(round(Odds_Ratio, 2), " (", round(CI_Lower, 2), ", ", round(CI_Upper, 2), ")"),
                col = ifelse(P_Value > 0.05, "grey", "black")), 
            position = position_dodge(width = 0.5), hjust = 0, size = 3.5) +
  scale_fill_manual(values = c("Epithelium" = "#BDB2FF",
                               "Vascular" = "#A0C4FF",
                               "Lymphoid" = "#FFADAD",
                               "Myeloid" = "#FDFFB6",
                               "Fibroblast" = "#CAFFBF")) +
  labs(title = "pCR, ER+", y = NULL, fill = "Cell Lineage") +
  theme_bw() +
  geom_vline(xintercept = 1, linetype = "solid", color = "black") +
  theme(axis.text.y = element_text(size = 8, color = ifelse(combined_data_er$P_Value > 0.05, "grey", "black")),
        legend.position = "bottom",
        plot.margin = margin(5.5, 90, 5.5, 5.5),
        panel.grid = element_blank())+
  coord_cartesian(clip = "off")+
  scale_color_identity()

#plot er-, sig only
combined_data_er <- combined_sig %>%
  filter(ER_status == "neg") %>%
  arrange(Odds_Ratio)
combined_data_er$Cell_State <- factor(combined_data_er$Cell_State, levels = combined_data_er$Cell_State)

ggplot(combined_data_er, aes(x = Odds_Ratio, y = Cell_State, fill = Cell_Type)) +
  geom_point(size = 5, shape = 22, color = "black") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, col = "black") +
  scale_size_continuous(range = c(2, 10)) +
  geom_text(aes(x = max(Odds_Ratio) + 0.1, 
                label = paste0(round(Odds_Ratio, 2), " (", round(CI_Lower, 2), ", ", round(CI_Upper, 2), ")"),
                col = ifelse(P_Value > 0.05, "grey", "black")), 
            position = position_dodge(width = 0.5), hjust = 0, size = 3.5) +
  scale_fill_manual(values = c("Epithelium" = "#BDB2FF",
                               "Vascular" = "#A0C4FF",
                               "Lymphoid" = "#FFADAD",
                               "Myeloid" = "#FDFFB6",
                               "Fibroblast" = "#CAFFBF")) +
  labs(title = "pCR, ER-", y = NULL, fill = "Cell Lineage") +
  theme_bw() +
  geom_vline(xintercept = 1, linetype = "solid", color = "black") +
  theme(axis.text.y = element_text(size = 8, color = ifelse(combined_data_er$P_Value > 0.05, "grey", "black")),
        legend.position = "bottom",
        plot.margin = margin(5.5, 90, 5.5, 5.5),
        panel.grid = element_blank())+
  coord_cartesian(clip = "off")+
  scale_color_identity()

#plot adjusted multivariable model (ER+ speicifc)####
combined_standard_er <- filter(combined_standard, ER_status == "pos")

or_results_sig <- or_results[or_results$P_Value < 0.05, ]
or_results_er <- filter(or_results_sig, ER_status == "pos")
multi_variables <- unique(or_results_er$Cell_Type)

cell_types <- colnames(combined_standard[,25:47])

# Function to calculate quantiles
calculate_quantiles <- function(df, cols, quantiles) {
  df %>%
    dplyr::mutate(across(all_of(cols), ~ ntile(.x, quantiles), .names = "quantile_{col}"))
}

# Apply the quantile calculation
combined_standard_er_quantiles <- calculate_quantiles(combined_standard_er, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  dplyr::select(-all_of(cell_types))

colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
quantile_vars_escaped <- paste0("`", multi_variables, "`")

# Prepare the formula for the Cox model
formula_str <- paste("as.integer(pcr) ~", paste(quantile_vars_escaped, collapse = " + "), " + age + study + HER2_status")
model_formula <- as.formula(formula_str)
model <- glm(model_formula, data = combined_standard_er_quantiles)

# Extract coefficients and hazard ratios
hr <- exp(model$coefficients)
hr_df <- data.frame(
  variable = names(hr),
  hr = hr
)

hr_df <- hr_df[order(hr_df$hr), ]
patterns_to_exclude <- c("HER2", "study", "intercept")
keep_rows <- !grepl(paste(patterns_to_exclude, collapse = "|"), hr_df$variable, ignore.case = TRUE)
hr_df <- hr_df[keep_rows, ]

keep_rows <- !str_detect(hr_df$variable, "\\b(HER2|study|age|intercept)\\b")
keep_rows <- !str_detect(hr_df$variable, "^age\\d{2}-\\d{2}$")
hr_df <- hr_df[keep_rows, ]

hr_df_variables <- hr_df$variable
hr_df_variables <- gsub("`", "", hr_df_variables)
hr_df_variables <- paste0("`", hr_df_variables, "`")

# Prepare the formula for the Cox model
formula_str <- paste("as.integer(pcr) ~", paste(hr_df_variables, collapse = " + "), "+ age + study + HER2_status")
model_formula <- as.formula(formula_str)
model <- glm(model_formula, data = combined_standard_er_quantiles)
summary <- summary(model)

#forest plot
panels <- list(
  list(width = 0.02),
  list(width = 0.01, display = ~variable, heading = "ER+, pCR", hjust = 1),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(
    width = 0.2, item = "forest", hjust = 0.5, heading = "OR (95% CI)", linetype = "dashed", 
    line_x = 0),
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

forest_model_order <- forest_model(model, exponentiate = TRUE, covariates = multi_variables, panels)
print(forest_model_order)



#plot adjusted multivariable model (ER- speicifc)####
combined_standard_er <- filter(combined_standard, ER_status == "neg")

or_results_sig <- or_results[or_results$P_Value < 0.05, ]
or_results_er <- filter(or_results_sig, ER_status == "neg")
multi_variables <- unique(or_results_er$Cell_Type)

cell_types <- colnames(combined_standard[,25:47])

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
formula_str <- paste("as.integer(pcr) ~", paste(quantile_vars_escaped, collapse = " + "), " + age + study + HER2_status")
model_formula <- as.formula(formula_str)
model <- glm(model_formula, data = combined_standard_er_quantiles)

# Extract coefficients and hazard ratios
hr <- exp(model$coefficients)
hr_df <- data.frame(
  variable = names(hr),
  hr = hr
)

hr_df <- hr_df[order(hr_df$hr), ]
patterns_to_exclude <- c("HER2", "study", "intercept")
keep_rows <- !grepl(paste(patterns_to_exclude, collapse = "|"), hr_df$variable, ignore.case = TRUE)
hr_df <- hr_df[keep_rows, ]

keep_rows <- !str_detect(hr_df$variable, "\\b(HER2|study|age|intercept)\\b")
keep_rows <- !str_detect(hr_df$variable, "^age\\d{2}-\\d{2}$")
hr_df <- hr_df[keep_rows, ]

hr_df_variables <- hr_df$variable
hr_df_variables <- gsub("`", "", hr_df_variables)
hr_df_variables <- paste0("`", hr_df_variables, "`")

# Prepare the formula for the Cox model
formula_str <- paste("as.integer(pcr) ~", paste(hr_df_variables, collapse = " + "), "+ age + study + HER2_status")
model_formula <- as.formula(formula_str)
model <- glm(model_formula, data = combined_standard_er_quantiles)
summary <- summary(model)

#forest plot
panels <- list(
  list(width = 0.02),
  list(width = 0.01, display = ~variable, heading = "ER-, pCR", hjust = 0.5),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(
    width = 0.2, item = "forest", hjust = 0.5, heading = "OR (95% CI)", linetype = "dashed", 
    line_x = 0),
  list(width = 0.03),
  list(width = 0.12, display = ~ ifelse(reference, "Reference", sprintf(
    "%0.2f (%0.2f, %0.2f)",
    trans(estimate), trans(conf.low), trans(conf.high)
  )), display_na = NA),
  list(width = 0.01, item = "vline", hjust = 0.5),
  list(
    width = 0.01,
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 2, eps = 0.001)),
    display_na = NA, hjust = 0.5, heading = "p"
  ),
  list(width = 0.03)
)

forest_model_order <- forest_model(model, exponentiate = TRUE, covariates = multi_variables, panels)
print(forest_model_order)




#spineplots of RCB for all cell subtypes, ER-####
combined_standard_spine <- combined_standard
combined_standard_spine$rcb[combined_standard_spine$pcr == 1] <- "pcr"
calculate_quantiles_rcb <- combined_standard_spine %>% filter(!is.na(rcb))
calculate_quantiles_rcb <- calculate_quantiles_rcb %>% filter(study != 'ispy2')
combined_standard_er_quantiles <- calculate_quantiles(calculate_quantiles_rcb, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  dplyr::select(-all_of(cell_types)) %>%
  filter(ER_status == 'neg')
colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
combined_standard_er_quantiles$rcb[combined_standard_er_quantiles$rcb == "RCB-0/I"] <- "RCB-I"

combined_standard_er_quantiles$Cycling_T_cells <- as.factor(combined_standard_er_quantiles$Cycling_T_cells)
combined_standard_er_quantiles$rcb <- as.factor(combined_standard_er_quantiles$rcb)

spine_plots <- lapply(cell_types, function(cell_type) {
  combined_standard_er_quantiles[[cell_type]] <- as.factor(combined_standard_er_quantiles[[cell_type]])
  ggplot(combined_standard_er_quantiles, aes_string(x = cell_type, fill = "rcb")) +
    geom_bar(position = "fill", color = "black") +
    labs(x = "quantile", y = "Proportion") +
    ggtitle(paste(cell_type)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")  # Remove the legend
})

# Arrange plots in a grid
library(gridExtra)
do.call("grid.arrange", c(spine_plots, ncol = 5))  # Adjust the number of columns as needed

#spineplots of RCB for all cell subtypes, ER+####
combined_standard_spine <- combined_standard
combined_standard_spine$rcb[combined_standard_spine$pcr == 1] <- "pcr"
calculate_quantiles_rcb <- combined_standard_spine %>% filter(!is.na(rcb))
calculate_quantiles_rcb <- calculate_quantiles_rcb %>% filter(study != 'ispy2')
combined_standard_er_quantiles <- calculate_quantiles(calculate_quantiles_rcb, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  dplyr::select(-all_of(cell_types)) %>%
  filter(ER_status == 'pos')
colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
combined_standard_er_quantiles$rcb[combined_standard_er_quantiles$rcb == "RCB-0/I"] <- "RCB-I"

combined_standard_er_quantiles$Cycling_T_cells <- as.factor(combined_standard_er_quantiles$Cycling_T_cells)
combined_standard_er_quantiles$rcb <- as.factor(combined_standard_er_quantiles$rcb)

spine_plots <- lapply(cell_types, function(cell_type) {
  combined_standard_er_quantiles[[cell_type]] <- as.factor(combined_standard_er_quantiles[[cell_type]])
  ggplot(combined_standard_er_quantiles, aes_string(x = cell_type, fill = "rcb")) +
    geom_bar(position = "fill", color = "black") +
    labs(x = cell_type, y = "RCB Fraction") +
    ggtitle(paste(cell_type)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")  # Remove the legend
})

# Arrange plots in a grid
library(gridExtra)
do.call("grid.arrange", c(spine_plots, ncol = 5))  # Adjust the number of columns as needed

#spineplots of pCR for all cell subtypes, ER-####
combined_standard_spine <- combined_standard
calculate_quantiles_pcr <- combined_standard_spine %>% filter(!is.na(pcr))
combined_standard_er_quantiles <- calculate_quantiles(calculate_quantiles_pcr, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  select(-all_of(cell_types)) %>%
  filter(ER_status == 'neg')
colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
combined_standard_er_quantiles$pcr <- ifelse(combined_standard_er_quantiles$pcr == 1, "PCR", "No PCR")

combined_standard_er_quantiles$Cycling_Myeloid <- as.factor(combined_standard_er_quantiles$Cycling_Myeloid)
combined_standard_er_quantiles$pcr <- as.factor(combined_standard_er_quantiles$pcr)

spine_plots <- lapply(cell_types, function(cell_type) {
  combined_standard_er_quantiles[[cell_type]] <- as.factor(combined_standard_er_quantiles[[cell_type]])
  ggplot(combined_standard_er_quantiles, aes_string(x = cell_type, fill = "pcr")) +
    geom_bar(position = "fill", color = "black") +
    labs(x = "quantile", y = "Proportion") +
    ggtitle(paste(cell_type)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

# Arrange plots in a grid
library(gridExtra)
do.call("grid.arrange", c(spine_plots, ncol = 5))  # Adjust the number of columns as needed

#as a line graph
melted_data <- melt(combined_standard_er_quantiles, 
                    id.vars = c("pcr"), 
                    measure.vars = cell_types, 
                    variable.name = "cell_type", 
                    value.name = "quartile")
pcr_summary <- melted_data %>%
  group_by(cell_type, quartile) %>%
  summarise(pct_pcr = mean(as.numeric(pcr == 1)) * 100, .groups = "drop")

ggplot(pcr_summary, aes(x = quartile, y = pct_pcr, color = cell_type, group = cell_type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "Quartile", y = "% pCR", title = "Percentage of pCR by Quartile for Each Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

#spineplots of PCR####
##Luminal Progenitor, ER-
combined_standard_spine <- combined_standard
calculate_quantiles_pcr <- combined_standard_spine %>% filter(!is.na(pcr))
combined_standard_er_quantiles <- calculate_quantiles(calculate_quantiles_pcr, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  select(-all_of(cell_types)) %>%
  filter(ER_status == 'neg')
colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
combined_standard_er_quantiles$pcr <- ifelse(combined_standard_er_quantiles$pcr == 1, "PCR", "No PCR")

combined_standard_er_quantiles$Luminal_Progenitors <- as.factor(combined_standard_er_quantiles$Luminal_Progenitors)
combined_standard_er_quantiles$pcr <- as.factor(combined_standard_er_quantiles$pcr)

bar_colors <- ifelse(levels(combined_standard_er_quantiles$pcr) == "PCR", "#491b59", "#BDB2FF")

spineplot(as.factor(pcr) ~ Luminal_Progenitors, data = combined_standard_er_quantiles, 
          main = "Luminal Progenitors, ER-", 
          xlab = "Quartile", 
          ylab = "",
          col = c("PCR" = "#491b59", "no PCR" = "#BDB2FF"))

##T cells MKI67, ER-
combined_standard_spine <- combined_standard
calculate_quantiles_pcr <- combined_standard_spine %>% filter(!is.na(pcr))
combined_standard_er_quantiles <- calculate_quantiles(calculate_quantiles_pcr, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  select(-all_of(cell_types)) %>%
  filter(ER_status == 'neg')
colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
combined_standard_er_quantiles$pcr <- ifelse(combined_standard_er_quantiles$pcr == 1, "PCR", "No PCR")

combined_standard_er_quantiles$T_cells_c11_MKI67 <- as.factor(combined_standard_er_quantiles$T_cells_c11_MKI67)
combined_standard_er_quantiles$pcr <- as.factor(combined_standard_er_quantiles$pcr)

spineplot(as.factor(pcr) ~ T_cells_c11_MKI67, data = combined_standard_er_quantiles, 
          main = "Progenitor T Cells, ER-", 
          xlab = "Quartile", 
          ylab = "",
          col = c("PCR" = "#FC5F5F", "no PCR" = "#FFADAD"))

##CAF s3, ER-
combined_standard_spine <- combined_standard
calculate_quantiles_pcr <- combined_standard_spine %>% filter(!is.na(pcr))
combined_standard_er_quantiles <- calculate_quantiles(calculate_quantiles_pcr, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  select(-all_of(cell_types)) %>%
  filter(ER_status == 'neg')
colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
combined_standard_er_quantiles$pcr <- ifelse(combined_standard_er_quantiles$pcr == 1, "PCR", "No PCR")

combined_standard_er_quantiles$CAFs_Transitioning_s3 <- as.factor(combined_standard_er_quantiles$CAFs_Transitioning_s3)
combined_standard_er_quantiles$pcr <- as.factor(combined_standard_er_quantiles$pcr)

spineplot(as.factor(pcr) ~ CAFs_Transitioning_s3, data = combined_standard_er_quantiles, 
          main = "Transitioning s3 CAFs, ER-", 
          xlab = "Quartile", 
          ylab = "",
          col = c("PCR" = "#2AFF6A", "no PCR" = "#CAFFBF"))

##Luminal progenitors, ER+
combined_standard_spine <- combined_standard
calculate_quantiles_pcr <- combined_standard_spine %>% filter(!is.na(pcr))
combined_standard_er_quantiles <- calculate_quantiles(calculate_quantiles_pcr, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  select(-all_of(cell_types)) %>%
  filter(ER_status == 'pos')
colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
combined_standard_er_quantiles$pcr <- ifelse(combined_standard_er_quantiles$pcr == 1, "PCR", "No PCR")

combined_standard_er_quantiles$Luminal_Progenitors <- as.factor(combined_standard_er_quantiles$Luminal_Progenitors)
combined_standard_er_quantiles$pcr <- as.factor(combined_standard_er_quantiles$pcr)

spineplot(as.factor(pcr) ~ Luminal_Progenitors, data = combined_standard_er_quantiles, 
          main = "Luminal Progenitors, ER+", 
          xlab = "Quartile", 
          ylab = "",
          col = c("PCR" = "#491b59", "no PCR" = "#BDB2FF"))

##T cells MKI67, ER+
combined_standard_spine <- combined_standard
calculate_quantiles_pcr <- combined_standard_spine %>% filter(!is.na(pcr))
combined_standard_er_quantiles <- calculate_quantiles(calculate_quantiles_pcr, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  select(-all_of(cell_types)) %>%
  filter(ER_status == 'pos')
colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
combined_standard_er_quantiles$pcr <- ifelse(combined_standard_er_quantiles$pcr == 1, "PCR", "No PCR")

combined_standard_er_quantiles$T_cells_c11_MKI67 <- as.factor(combined_standard_er_quantiles$T_cells_c11_MKI67)
combined_standard_er_quantiles$pcr <- as.factor(combined_standard_er_quantiles$pcr)

spineplot(as.factor(pcr) ~ T_cells_c11_MKI67, data = combined_standard_er_quantiles, 
          main = "Progenitor T Cells, ER+", 
          xlab = "Quartile", 
          ylab = "",
          col = c("PCR" = "#CBCE6D", "no PCR" = "#FDFFB6"))

#spineplots of RCB for specific cell subtypes, ER-####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
combined_standard$pcr <- gsub("NA", NA, combined_standard$pcr)

new_names <- names(combined_standard)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "T-regs", "Tregs")  # Replace T-regs with Tregs
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_standard <- setNames(combined_standard, new_names)
combined_standard <- filter(combined_standard, !is.na(pcr))

or_results <- data.frame(Cell_State = character(),
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

combined_standard <- combined_standard %>%
  filter(treatment_catagory == "chemo")
combined_standard <- as.data.frame(combined_standard)
combined_standard <- combined_standard %>% 
  dplyr::select(patient_id, ER_status, pcr, rcb, study, `Lipid-Associated Macrophage 2`, myCAF, `Transitioning CAF`)

cell_types <- colnames(combined_standard[,6:8])

combined_standard_spine <- combined_standard
combined_standard_spine$rcb[combined_standard_spine$pcr == 1] <- "pCR"
calculate_quantiles_rcb <- combined_standard_spine %>% filter(!is.na(rcb))
calculate_quantiles_rcb <- calculate_quantiles_rcb %>% filter(study != 'ispy2')
combined_standard_er_quantiles <- calculate_quantiles(calculate_quantiles_rcb, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  dplyr::select(-all_of(cell_types)) %>%
  filter(ER_status == 'neg')
colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
combined_standard_er_quantiles$rcb[combined_standard_er_quantiles$rcb == "RCB-0/I"] <- "RCB-I"

combined_standard_er_quantiles$rcb <- as.factor(combined_standard_er_quantiles$rcb)

# Specify the desired order of cell types
cell_types <- c("myCAF", "Transitioning CAF", "Lipid-Associated Macrophage 2")

# Generate spine plots in the specified order
spine_plots <- lapply(cell_types, function(cell_type) {
  combined_standard_er_quantiles[[cell_type]] <- as.factor(combined_standard_er_quantiles[[cell_type]])
  ggplot(combined_standard_er_quantiles, aes(x = .data[[cell_type]], fill = rcb)) +
    geom_bar(position = "fill", color = "black") +
    labs(x = paste0(cell_type), y = "Proportion", fill = "RCB Category") +  # Add legend title
    ggtitle("ER-") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank()
    )
})

# Combine plots and legend
combined_plot <- plot_grid(
  plot_grid(plotlist = spine_plots, nrow = 1),  # Adjust columns as needed
  ncol = 1,
  rel_heights = c(4, 0.5)  # Adjust relative heights
)

# Display the combined plot
print(combined_plot)


#spineplots of RCB for specific cell subtypes, ER+####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
combined_standard$pcr <- gsub("NA", NA, combined_standard$pcr)

new_names <- names(combined_standard)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "T-regs", "Tregs")  # Replace T-regs with Tregs
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_standard <- setNames(combined_standard, new_names)
combined_standard <- filter(combined_standard, !is.na(pcr))

or_results <- data.frame(Cell_State = character(),
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

combined_standard <- combined_standard %>%
  filter(treatment_catagory == "chemo")
combined_standard <- as.data.frame(combined_standard)
combined_standard <- combined_standard %>% 
  dplyr::select(patient_id, ER_status, pcr, rcb, study, 'Mature Luminal Cells', 'Monocyte [FCGR3A]', 'Proliferating PVL')

cell_types <- colnames(combined_standard[,6:8])

combined_standard_spine <- combined_standard
combined_standard_spine$rcb[combined_standard_spine$pcr == 1] <- "pCR"
calculate_quantiles_rcb <- combined_standard_spine %>% filter(!is.na(rcb))
calculate_quantiles_rcb <- calculate_quantiles_rcb %>% filter(study != 'ispy2')
combined_standard_er_quantiles <- calculate_quantiles(calculate_quantiles_rcb, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  dplyr::select(-all_of(cell_types)) %>%
  filter(ER_status == 'pos')
colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
combined_standard_er_quantiles$rcb[combined_standard_er_quantiles$rcb == "RCB-0/I"] <- "RCB-I"

combined_standard_er_quantiles$rcb <- as.factor(combined_standard_er_quantiles$rcb)

# Specify the desired order of cell types
cell_types <- c("Proliferating PVL", "Monocyte [FCGR3A]", "Mature Luminal Cells")

# Generate spine plots in the specified order
spine_plots <- lapply(cell_types, function(cell_type) {
  combined_standard_er_quantiles[[cell_type]] <- as.factor(combined_standard_er_quantiles[[cell_type]])
  ggplot(combined_standard_er_quantiles, aes(x = .data[[cell_type]], fill = rcb)) +
    geom_bar(position = "fill", color = "black") +
    labs(x = paste0(cell_type), y = "Proportion", fill = "RCB Category") +  # Add legend title
    ggtitle("ER+") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank()
    )
})

# Combine plots and legend
combined_plot <- plot_grid(
  plot_grid(plotlist = spine_plots, nrow = 1),  # Adjust columns as needed
  ncol = 1,
  rel_heights = c(4, 0.5)  # Adjust relative heights
)

# Display the combined plot
print(combined_plot)


#Chemo-immunotherapy, cell types transneo only####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell Types.txt")
combined_standard$pcr <- gsub("NA", NA, combined_standard$pcr)
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
transneo_clinical <- read_excel("TransNeo_Clinical.xlsx")
transneo_age <- transneo_clinical %>%
  select(Donor.ID, Age) %>%
  rename(patient_id = Donor.ID)

or_results <- data.frame(Cell_State = character(),
                         Variable = character(),
                         Hazard_Ratio = numeric(),
                         CI_Lower = numeric(),
                         CI_Upper = numeric(),
                         P_Value = numeric(),
                         stringsAsFactors = FALSE)

numeric_cols <- c(25:31, 33)
combined_standard[, (numeric_cols) := lapply(.SD, function(x) x / (1 - `Cancer Epithelial`)), .SDcols = numeric_cols]
rowSums(combined_standard[, ..numeric_cols], na.rm = TRUE)
combined_standard <- filter(combined_standard, treatment_catagory %in% c("chemoimmuno", "immuno"))
combined_standard <- filter(combined_standard, study %in% c("transneo"))
combined_standard <- merge(combined_standard, transneo_age, by = "patient_id")
combined_standard <- filter(combined_standard, !is.na(pcr))


cell_types <- unique(colnames(combined_standard[,..numeric_cols]))
ER_status <- c("pos", "neg")

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
    
    
    # Fit the linear regression model
    model <- glm(as.integer(pcr) ~ cell_quartiles + Age, data = combined_standard_er)
    summary <- summary(model)
    
    odds_ratio <- exp(coef(model)["cell_quartiles"])
    conf_interval <- exp(confint(model)["cell_quartiles", ])
    p_value <- coef(summary)["cell_quartiles", "Pr(>|t|)"]
    
    # Store results in the data frame
    or_results_cell <- data.frame(ER_status = ER_sub,
                                  Cell_Type = cell,
                                  Odds_Ratio = odds_ratio,
                                  CI_Lower = conf_interval[1],
                                  CI_Upper = conf_interval[2],
                                  P_Value = p_value)
    
    or_results <- rbind(or_results, or_results_cell)
  }
}
or_results$P_Value_BH <- p.adjust(or_results$P_Value, method = "BH")

#assign cell types
cox_results_2 <- data.table(or_results)
cox_results_2 <- dplyr::rename(cox_results_2, c("Cell_State" = "Cell_Type"))
cox_results_2$Cell_Type <- NA
cox_results_2$Cell_Type[grep("PVL", cox_results_2$Cell_State)] <- "Vascular"
cox_results_2$Cell_Type[grep("Endothelial", cox_results_2$Cell_State)] <- "Vascular"
cox_results_2$Cell_Type[grep("CAFs", cox_results_2$Cell_State)] <- "Fibroblast"
cox_results_2$Cell_Type[grep("T-cells", cox_results_2$Cell_State)] <- "Lymphoid"
cox_results_2$Cell_Type[grep("B-cells", cox_results_2$Cell_State)] <- "Lymphoid"
cox_results_2$Cell_Type[grep("Myeloid", cox_results_2$Cell_State)] <- "Myeloid"
cox_results_2$Cell_Type[grep("Normal", cox_results_2$Cell_State)] <- "Epithelium"
cox_results_2$Cell_Type[grep("DC", cox_results_2$Cell_State)] <- "Myeloid"

#plot er statuses together, ordered
er_pos_data <- cox_results_2[ER_status == "pos"]
er_neg_data <- cox_results_2[ER_status == "neg"]
er_pos_data$Cell_State <- factor(er_pos_data$Cell_State, levels = er_pos_data$Cell_State[order(er_pos_data$Odds_Ratio)])
combined_data <- rbind(er_pos_data, er_neg_data)

#plot er+, sig only
combined_sig <- combined_data[combined_data$P_Value < 1, ]
combined_data_er <- combined_sig[ER_status == "pos", ]
combined_data_er <- combined_data_er %>% arrange(Odds_Ratio)

ggplot(combined_data_er, aes(x = Odds_Ratio, y = Cell_State, fill = Cell_Type)) +
  geom_point(size = 5, shape = 22, color = "black") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, col = "black") +
  scale_size_continuous(range = c(2, 10)) +
  geom_text(aes(x = max(Odds_Ratio) + 0.075, 
                label = paste0(round(Odds_Ratio, 2), " (", round(CI_Lower, 2), ", ", round(CI_Upper, 2), ")"),
                col = ifelse(P_Value > 0.05, "grey", "black")), 
            position = position_dodge(width = 0.5), hjust = 0, size = 3.5) +
  scale_fill_manual(values = c("Epithelium" = "#BDB2FF",
                               "Vascular" = "#A0C4FF",
                               "Lymphoid" = "#FFADAD",
                               "Myeloid" = "#FDFFB6",
                               "Fibroblast" = "#CAFFBF")) +
  labs(title = "pCR, ER+", y = NULL, fill = "Cell Lineage") +
  theme_bw() +
  geom_vline(xintercept = 1, linetype = "solid", color = "black") +
  theme(axis.text.y = element_text(size = 8, color = ifelse(combined_data_er$P_Value > 0.05, "grey", "black")),
        legend.position = "bottom",
        plot.margin = margin(5.5, 90, 5.5, 5.5),
        panel.grid = element_blank())+
  coord_cartesian(clip = "off")+
  scale_color_identity()


#plot er-, sig only
combined_data_er <- combined_sig %>%
  filter(ER_status == "neg") %>%
  arrange(Odds_Ratio)
combined_data_er$Cell_State <- factor(combined_data_er$Cell_State, levels = combined_data_er$Cell_State)

ggplot(combined_data_er, aes(x = Odds_Ratio, y = Cell_State, fill = Cell_Type)) +
  geom_point(size = 5, shape = 22, color = "black") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, col = "black") +
  scale_size_continuous(range = c(2, 10)) +
  geom_text(aes(x = max(Odds_Ratio) + 0.075, 
                label = paste0(round(Odds_Ratio, 2), " (", round(CI_Lower, 2), ", ", round(CI_Upper, 2), ")"),
                col = ifelse(P_Value > 0.05, "grey", "black")), 
            position = position_dodge(width = 0.5), hjust = 0, size = 3.5) +
  scale_fill_manual(values = c("Epithelium" = "#BDB2FF",
                               "Vascular" = "#A0C4FF",
                               "Lymphoid" = "#FFADAD",
                               "Myeloid" = "#FDFFB6",
                               "Fibroblast" = "#CAFFBF")) +
  labs(title = "pCR, ER-", y = NULL, fill = "Cell Lineage") +
  theme_bw() +
  geom_vline(xintercept = 1, linetype = "solid", color = "black") +
  theme(axis.text.y = element_text(size = 8, color = ifelse(combined_data_er$P_Value > 0.05, "grey", "black")),
        legend.position = "bottom",
        plot.margin = margin(5.5, 90, 5.5, 5.5),
        panel.grid = element_blank())+
  coord_cartesian(clip = "off")+
  scale_color_identity()


#Chemo-immunotherapy, cell states spineplots transneo only in ER+####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
combined_standard$pcr <- gsub("NA", NA, combined_standard$pcr)

new_names <- names(combined_standard)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "T-regs", "Tregs")  # Replace T-regs with Tregs
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_standard <- setNames(combined_standard, new_names)
combined_standard <- filter(combined_standard, !is.na(pcr))

or_results <- data.frame(Cell_State = character(),
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

combined_standard <- combined_standard %>%
  filter(treatment_catagory %in% c("chemo", "chemoimmuno")) %>%
  filter(study %in% c("transneo"))
combined_standard <- as.data.frame(combined_standard)
combined_standard <- combined_standard %>% 
  dplyr::select(patient_id, ER_status, pcr, rcb, study, myCAF, iCAF, `Transitioning CAF`, `Plasmacytoid DC`)

cell_types <- colnames(combined_standard[,6:9])

combined_standard_spine <- combined_standard
combined_standard_spine$rcb[combined_standard_spine$pcr == 1] <- "pCR"
calculate_quantiles_rcb <- combined_standard_spine %>% filter(!is.na(rcb))
calculate_quantiles_rcb <- calculate_quantiles_rcb %>% filter(study != 'ispy2')
combined_standard_er_quantiles <- calculate_quantiles(calculate_quantiles_rcb, cell_types, 4)
combined_standard_er_quantiles <- combined_standard_er_quantiles %>%
  dplyr::select(-all_of(cell_types)) %>%
  filter(ER_status == 'pos')
colnames(combined_standard_er_quantiles) <- gsub("^quantile_", "", colnames(combined_standard_er_quantiles))
combined_standard_er_quantiles$rcb[combined_standard_er_quantiles$rcb == "RCB-0/I"] <- "RCB-I"

combined_standard_er_quantiles$rcb <- as.factor(combined_standard_er_quantiles$rcb)

# Specify the desired order of cell types
cell_types <- c("myCAF", "iCAF", "Transitioning CAF", "Plasmacytoid DC")

# Generate spine plots in the specified order
spine_plots <- lapply(cell_types, function(cell_type) {
  combined_standard_er_quantiles[[cell_type]] <- as.factor(combined_standard_er_quantiles[[cell_type]])
  ggplot(combined_standard_er_quantiles, aes(x = .data[[cell_type]], fill = rcb)) +
    geom_bar(position = "fill", color = "black") +
    labs(x = paste0(cell_type), y = "Proportion", fill = "RCB Category") +  # Add legend title
    ggtitle("ER+") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank()
    )
})

# Combine plots and legend
combined_plot <- plot_grid(
  plot_grid(plotlist = spine_plots, nrow = 1),  # Adjust columns as needed
  nrow = 1,
  rel_heights = c(4, 0.5)  # Adjust relative heights
)

# Display the combined plot
print(combined_plot)



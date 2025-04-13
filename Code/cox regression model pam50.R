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

surv_metric <- "RFS" # options are RFS, BCSS, MFS, OS

# Load the cell type data####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
new_names <- gsub("-", "_", names(combined_standard))
new_names <- gsub(" ", "_", names(combined_standard))
combined_standard <- setNames(combined_standard, new_names)
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

# cox regression model test at cell state with quartiles ####
cell_types <- unique(colnames(combined_standard[,25:68]))
pam50_status <- c("Normal", "Basal", "Her2", "LumA", "LumB")
combined_standard_filtered <- combined_standard[combined_standard$pam50 %in% pam50_status, ]

for (pam50_sub in pam50_status) {
  print(pam50_sub)
  combined_standard_er <- combined_standard[combined_standard$pam50 == pam50_status, ]
  
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
                                " ~ cell_quartiles + stage + node_status + strata(study)"))
    
    
    # Fit the Cox proportional hazards model
    cox_model <- coxph(formula, data = combined_standard_er)
    summary_cox <- summary(cox_model)
    hr <- coef(cox_model)["cell_quartiles"] # Hazard Ratio
    ci <- confint(cox_model, level = 0.95)["cell_quartiles", ] # Confidence Interval
    p_value <- summary_cox$coefficients["cell_quartiles", "Pr(>|z|)"] # P-Value
    
    # Store results in the data frame
    cox_results_cell <- data.frame(pam50_status = pam50_sub,
                                   Cell_Type = cell,
                                   Hazard_Ratio = hr,
                                   CI_Lower = ci[1],
                                   CI_Upper = ci[2],
                                   P_Value = p_value)
    
    cox_results <- rbind(cox_results, cox_results_cell)
  }
}

#plot cox regression model results####
#plot er statuses together
cox_results <- data.table(cox_results)
cox_results <- dplyr::rename(cox_results, c("Cell_State" = "Cell_Type"))
cox_results$Cell_Type[grep("Cancer", cox_results$Cell_State)] <- "Cancer Epithelial"
cox_results$Cell_Type[grep("Luminal", cox_results$Cell_State)] <- "Normal Epithelial"
cox_results$Cell_Type[grep("epithelial", cox_results$Cell_State)] <- "Normal Epithelial"
cox_results$Cell_Type[grep("Plasmablasts", cox_results$Cell_State)] <- "Plasmablasts"
cox_results$Cell_Type[grep("PVL", cox_results$Cell_State)] <- "PVL"
cox_results$Cell_Type[grep("T_cells", cox_results$Cell_State)] <- "T Cells"
cox_results$Cell_Type[grep("B_cells", cox_results$Cell_State)] <- "B Cells"
cox_results$Cell_Type[grep("Myeloid", cox_results$Cell_State)] <- "Myeloid"
cox_results$Cell_Type[grep("Endothelial", cox_results$Cell_State)] <- "Endothelial"
cox_results$Cell_Type[grep("CAFs", cox_results$Cell_State)] <- "CAFs"

#plot data that is pam50-specific and only significant####
#er+
combined_data_sig <- combined_data[combined_data$P_Value < 0.05, ]
combined_data_er <- filter(combined_data_sig, ER_status == "pos")

ggplot(combined_data_er, aes(x = Hazard_Ratio, y = Cell_State, fill = Cell_Type)) +
  geom_point(size = 5, shape = 22, position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, position = position_dodge(width = 0.5)) +
  geom_text(aes(x = max(Hazard_Ratio) + 0.1, label = paste(round(Hazard_Ratio, 2), "(", round(CI_Lower, 2), ",", round(CI_Upper, 2), ")", sep = "")),
            position = position_dodge(width = 0.5), hjust = 0, size = 3.5) + 
  labs(title = paste("ER+, Forest Plot for", surv_metric), y = "Cell Subtype", fill = "Cell Lineage") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  coord_cartesian(clip = "off")

#er-
combined_data_er <- filter(combined_data_sig, ER_status == "neg")
combined_data_er <- combined_data_er[order(combined_data_er$Hazard_Ratio), ]

ggplot(combined_data_er, aes(x = Hazard_Ratio, y = reorder(Cell_State, Hazard_Ratio), fill = Cell_Type)) +
  geom_point(size = 5, shape = 22, position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, position = position_dodge(width = 0.5)) +
  geom_text(aes(x = max(Hazard_Ratio) + 0.1, label = paste(round(Hazard_Ratio, 2), "(", round(CI_Lower, 2), ",", round(CI_Upper, 2), ")", sep = "")),
            position = position_dodge(width = 0.5), hjust = 0, size = 3.5) + 
  labs(title = paste("ER-, Forest Plot for", surv_metric), y = "Cell Subtype", fill = "Cell Lineage") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  coord_cartesian(clip = "off")


#plot cox regression model results for myloid cells####
#plot er statuses together
cox_results <- data.table(cox_results)
cox_results <- filter(cox_results, grepl("Myeloid", Cell_State))

cox_results$Cell_Type[grep("Macrophage", cox_results$Cell_State)] <- "Macrophage"
cox_results$Cell_Type[grep("LAM", cox_results$Cell_State)] <- "Macrophage"
cox_results$Cell_Type[grep("Monocyte", cox_results$Cell_State)] <- "Monocyte"
cox_results$Cell_Type[grep("DC", cox_results$Cell_State)] <- "Dendritic Cell"
cox_results$Cell_Type[grep("Cycling", cox_results$Cell_State)] <- "Stem Cell"

#plot er statuses together, ordered
er_pos_data <- cox_results[ER_status == "pos"]
er_neg_data <- cox_results[ER_status == "neg"]
er_pos_data$Cell_State <- factor(er_pos_data$Cell_State, levels = er_pos_data$Cell_State[order(er_pos_data$Hazard_Ratio)])
combined_data <- rbind(er_pos_data, er_neg_data)
combined_data <- combined_data %>%
  filter(grepl("Macrophage", Cell_Type) | grepl("Dendritic Cell", Cell_Type))

# Create forest plot
ggplot(combined_data, aes(x = Hazard_Ratio, y = Cell_State, fill = ER_status)) +
  geom_point(size = 6, shape = 22, position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, position = position_dodge(width = 0.5)) +
  labs(title = paste("Forest Plot for", surv_metric)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 14)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  facet_wrap(~Cell_Type, scales = "free_y")+
  guides(fill = guide_legend(override.aes = list(size = c(5, 5))))


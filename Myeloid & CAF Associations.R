#libraries####
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(broom)
library(ggplot2)

# Load the cell type data, set the filtered group of intclust####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
combined_standard <- fread("Clinical Associations with Cell States.txt")
names(combined_standard) <- gsub("[- ]", "_", names(combined_standard))
names(combined_standard)

filtered_combined_standard <- combined_standard %>%
  select(matches("^(patient_id|intclust)$|Myeloid|CAFs")) %>%
  filter(intclust %in% c("10")) %>%
  mutate(intclust = as.numeric(intclust))

myeloid_cols <- grep("^Myeloid", colnames(filtered_combined_standard), value = TRUE)
caf_cols <- grep("^CAFs", colnames(filtered_combined_standard), value = TRUE)

#assess spearman regression between CAFs and myloid populations####
# Initialize an empty list to store regression results
regression_results <- list()

# Loop over Myeloid and CAF combinations
for (myeloid in myeloid_cols) {
  for (caf in caf_cols) {
    # Perform Spearman correlation for each combination
    spearman_result <- cor.test(filtered_combined_standard[[myeloid]], 
                                filtered_combined_standard[[caf]], 
                                method = "spearman")
    
    # Extract relevant information (correlation coefficient, p-value)
    tidy_result <- tibble(
      Myeloid = myeloid,
      CAF = caf,
      estimate = spearman_result$estimate,  # Spearman correlation coefficient
      p.value = spearman_result$p.value     # p-value
    )
    
    # Append to the results list
    regression_results <- append(regression_results, list(tidy_result))
  }
}

# Combine all results into a single dataframe
regression_results_df <- bind_rows(regression_results)

# Rename columns for clarity
colnames(regression_results_df) <- c("Myeloid", "CAF", "Coefficient", "P_Value")

# Apply Bonferroni correction
num_tests <- nrow(regression_results_df)
regression_results_df <- regression_results_df %>%
  mutate(Adjusted_P_Value = pmin(P_Value * num_tests, 1)) # Apply Bonferroni and ensure p-value doesn't exceed 1

# View the resulting dataframe
print(regression_results_df)

#Lineplot for CAF and myeloid for intclust 10####
#Myeloid_c7_Monocyte_3_FCGR3A and CAFs_myCAF_like_s5
filtered_combined_standard <- filtered_combined_standard %>%
  mutate(
    CAFs_myCAF_like_s5_quartile = cut(
      CAFs_myCAF_like_s5,
      breaks = quantile(CAFs_myCAF_like_s5, probs = seq(0, 1, 0.25), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )
  )

summary_stats <- filtered_combined_standard %>%
  group_by(CAFs_myCAF_like_s5_quartile) %>%
  summarise(
    mean_value = mean(Myeloid_c7_Monocyte_3_FCGR3A, na.rm = TRUE),
    sem_value = sd(Myeloid_c7_Monocyte_3_FCGR3A, na.rm = TRUE) / sqrt(n())
  )

gg <- ggplot(summary_stats, aes(x = CAFs_myCAF_like_s5_quartile, y = mean_value, group = 1)) +
  geom_line(size = 1, color = "black") + # Line connecting means
  geom_point(size = 3, color = "red") + # Points at means
  geom_errorbar(aes(ymin = mean_value - sem_value, ymax = mean_value + sem_value), width = 0.2, color = "black") +
  labs(
    x = "myCAF-like (s5) Quartile",
    y = "Monocyte [FCGR3A] Proportion",
    title = "IntClust 10"
  ) +
  annotate(
    "text",
    x = 0.8, 
    y = max(summary_stats$mean_value + summary_stats$sem_value) * 1.05, # Y-coordinate above the error bars
    label = "ρ = -0.40****",
    size = 6, 
    hjust = 0.4 
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  )

gg


#Myeloid_c4_DCs_pDC_IRF7 and CAFs_myCAF_like_s5
filtered_combined_standard <- filtered_combined_standard %>%
  mutate(
    CAFs_myCAF_like_s5_quartile = cut(
      CAFs_myCAF_like_s5,
      breaks = quantile(CAFs_myCAF_like_s5, probs = seq(0, 1, 0.25), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )
  )

summary_stats <- filtered_combined_standard %>%
  group_by(CAFs_myCAF_like_s5_quartile) %>%
  summarise(
    mean_value = mean(Myeloid_c4_DCs_pDC_IRF7, na.rm = TRUE),
    sem_value = sd(Myeloid_c4_DCs_pDC_IRF7, na.rm = TRUE) / sqrt(n())
  )

gg <- ggplot(summary_stats, aes(x = CAFs_myCAF_like_s5_quartile, y = mean_value, group = 1)) +
  geom_line(size = 1, color = "black") + # Line connecting means
  geom_point(size = 3, color = "red") + # Points at means
  geom_errorbar(aes(ymin = mean_value - sem_value, ymax = mean_value + sem_value), width = 0.2, color = "black") +
  labs(
    x = "myCAF-like (s5) Quartile",
    y = "pDC [IRF7] Proportion",
    title = "IntClust 10"
  ) +
  annotate(
    "text",
    x = 0.8, 
    y = max(summary_stats$mean_value + summary_stats$sem_value) * 1.05, # Y-coordinate above the error bars
    label = "ρ = -0.48****",
    size = 6, 
    hjust = 0.4 
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  )

gg
#Lineplot for CAF and myeloid for intclust 3####
#Myeloid_c3_cDC1_CLEC9A and CAFs_myCAF_like_s4
filtered_combined_standard <- filtered_combined_standard %>%
  mutate(
    CAFs_myCAF_like_s4_quartile = cut(
      CAFs_myCAF_like_s4,
      breaks = quantile(CAFs_myCAF_like_s4, probs = seq(0, 1, 0.25), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )
  )

summary_stats <- filtered_combined_standard %>%
  group_by(CAFs_myCAF_like_s4_quartile) %>%
  summarise(
    mean_value = mean(Myeloid_c3_cDC1_CLEC9A, na.rm = TRUE),
    sem_value = sd(Myeloid_c3_cDC1_CLEC9A, na.rm = TRUE) / sqrt(n())
  )

gg <- ggplot(summary_stats, aes(x = CAFs_myCAF_like_s4_quartile, y = mean_value, group = 1)) +
  geom_line(size = 1, color = "black") + # Line connecting means
  geom_point(size = 3, color = "red") + # Points at means
  geom_errorbar(aes(ymin = mean_value - sem_value, ymax = mean_value + sem_value), width = 0.2, color = "black") +
  labs(
    x = "myCAF-like (s4) Quartile",
    y = "cDC1 [CLEC9A] Proportion",
    title = "IntClust 3"
  ) +
  annotate(
    "text",
    x = 0.8, 
    y = max(summary_stats$mean_value + summary_stats$sem_value) * 1.05, # Y-coordinate above the error bars
    label = "ρ = -0.32****",
    size = 6, 
    hjust = 0.4 
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  )

gg


#Myeloid_c7_Monocyte_3_FCGR3A and CAFs_myCAF_like_s5
filtered_combined_standard <- filtered_combined_standard %>%
  mutate(
    CAFs_myCAF_like_s5_quartile = cut(
      CAFs_myCAF_like_s5,
      breaks = quantile(CAFs_myCAF_like_s5, probs = seq(0, 1, 0.25), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )
  )

summary_stats <- filtered_combined_standard %>%
  group_by(CAFs_myCAF_like_s5_quartile) %>%
  summarise(
    mean_value = mean(Myeloid_c7_Monocyte_3_FCGR3A, na.rm = TRUE),
    sem_value = sd(Myeloid_c7_Monocyte_3_FCGR3A, na.rm = TRUE) / sqrt(n())
  )

gg <- ggplot(summary_stats, aes(x = CAFs_myCAF_like_s5_quartile, y = mean_value, group = 1)) +
  geom_line(size = 1, color = "black") + # Line connecting means
  geom_point(size = 3, color = "red") + # Points at means
  geom_errorbar(aes(ymin = mean_value - sem_value, ymax = mean_value + sem_value), width = 0.2, color = "black") +
  labs(
    x = "myCAF-like (s5) Quartile",
    y = "Monocyte [FCGR3A] Proportion",
    title = "IntClust 3"
  ) +
  annotate(
    "text",
    x = 0.8, 
    y = max(summary_stats$mean_value + summary_stats$sem_value) * 1.05, # Y-coordinate above the error bars
    label = "ρ = -0.33****",
    size = 6, 
    hjust = 0.4 
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  )

gg
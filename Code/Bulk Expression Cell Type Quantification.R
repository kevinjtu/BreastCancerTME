# Path: Bulk Expression Cell Type Quantification.R

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggplot2)
library(stringr)
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Bulk Expression Cell Type Quantification")
imagedata <- fread("SingleCells.csv")

#This groups the imagedata into groups that have the same metabric_id and cellPhenotype
#Then it counts the number of rows in each group
result <- imagedata %>%
  group_by(metabric_id, cellPhenotype) %>%
  summarize(Total_Cells = n())

#Now find the proportion of each cellPhenotype in each metabric_id
Proportionresult <- result %>%
  group_by(metabric_id) %>%
  mutate(Proportion = Total_Cells/sum(Total_Cells))

#Now we will represent the Proporitonresult in a heatmap
#First we need to arrange Proportionresult into alphabetical order
#Then make a matrix with the metabric_id as the rows and the cellPhenotype as the columns
#We will then fill the matrix with the Proportion values

Proportionresult <- Proportionresult %>%
  arrange(metabric_id)

ggplot(Proportionresult, aes(x = reorder(metabric_id, -Proportion), y = cellPhenotype, fill = Proportion)) +
  geom_tile() +
  scale_fill_viridis_c() +  
  theme_minimal() +
  labs(x = "Metabric_ID (numerical order from top to bottom)", y = "Cell Phenotype", fill = "Proportion") +
  coord_flip() + # To rotate x-axis labels for better readability
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # To rotate x-axis labels for better readability

#configure the image data into wide format
imageres <- Proportionresult %>%
  pivot_wider(
    id_cols = metabric_id,          # Column to use as the identifier
    names_from = cellPhenotype, # Column to spread into new columns
    values_from = Proportion    # Values to populate the new columns
  )

imageres <- imageres %>%
  select("metabric_id", "B cells", "CD8^{+} T cells", "CD4^{+} T cells", "Macrophages", "T_{Reg} & T_{Ex}")

#add together the different cell types to get a total for each cell type
imageres_comp <- imageres %>%
  mutate(
    Bcells = `B cells`,
    CD4Tcells = `CD4^{+} T cells`,
    macrophages = `Macrophages`,
    CD8Tcells = `CD8^{+} T cells`,
    Tregs = `T_{Reg} & T_{Ex}`,
  ) %>%
  select(metabric_id, Bcells, CD4Tcells, macrophages, `CD8Tcells`, `Tregs`)

#order the samples by Bcells
imageres_Bcells_ranked <- imageres_comp %>%
  as.data.frame() %>%
  arrange(Bcells)





###Clean up data for cibersortx analysis
# Remove the entrezid from the bulk expression data to format it correctly for cibersortx, export to txt
expressiondata <- fread("data_mrna_illumina_microarray.txt")
cibersortbulk <- expressiondata[, -2]

# Cibersort says there are strings in the dataframe, Create a custom function to check for non-empty strings
cibersortbulk[is.na(cibersortbulk)] <- 0

cibersortbulk <- cibersortbulk %>%
  distinct(Hugo_Symbol, .keep_all = TRUE)

cibersortbulk <- as.data.frame(cibersortbulk)

write.table(cibersortbulk, "METABRIC_bulkexpression.tsv", sep="\t", row.names = FALSE, quote = FALSE)

is.recursive(cibersortbulk)
is.atomic(cibersortbulk)

#Create a subset of the data (first 1000 columns) to have less data to work with, fit into Cibersort
cibersortbulk_subset1 <- cibersortbulk[, 1:500]
cibersortbulk_subset1 <- as.data.frame(cibersortbulk_subset1)

write.table(cibersortbulk_subset1, "METABRIC_bulkexpression_Subset1.tsv", sep="\t", row.names = FALSE, quote = FALSE)

is.recursive(cibersortbulk_subset1)
is.atomic(cibersortbulk_subset1)

#Create a subset of the data (last 981 columns) to have less data to work with, fit into Cibersort
cibersortbulk_subset2 <- cibersortbulk[, c(1, 499:1000)]
cibersortbulk_subset2 <- as.data.frame(cibersortbulk_subset2)

write.table(cibersortbulk_subset2, "METABRIC_bulkexpression_Subset2.tsv", sep="\t", row.names = FALSE, quote = FALSE)

is.recursive(cibersortbulk_subset2)
is.atomic(cibersortbulk_subset2)

#Create a subset of the data (last 981 columns) to have less data to work with, fit into Cibersort
cibersortbulk_subset3 <- cibersortbulk[, c(1, 999:1500)]
cibersortbulk_subset3 <- as.data.frame(cibersortbulk_subset3)

write.table(cibersortbulk_subset3, "METABRIC_bulkexpression_Subset3.tsv", sep="\t", row.names = FALSE, quote = FALSE)

is.recursive(cibersortbulk_subset3)
is.atomic(cibersortbulk_subset3)

#Create a subset of the data (last 981 columns) to have less data to work with, fit into Cibersort
cibersortbulk_subset4 <- cibersortbulk[, c(1, 1499:1981)]
cibersortbulk_subset4 <- as.data.frame(cibersortbulk_subset4)

write.table(cibersortbulk_subset4, "METABRIC_bulkexpression_Subset4.tsv", sep="\t", row.names = FALSE, quote = FALSE)

is.recursive(cibersortbulk_subset4)
is.atomic(cibersortbulk_subset4)


###Ran cibersortx @ https://cibersortx.stanford.edu/
#import the cibersort results

ciberres1 <- fread("CIBERSORTx_METABRIC Subset 1_Results.csv")
ciberres2 <- fread("CIBERSORTx_METABRIC Subset 2_Results.csv")
ciberres3 <- fread("CIBERSORTx_METABRIC Subset 3_Results.csv")
ciberres4 <- fread("CIBERSORTx_METABRIC Subset 4_Results.csv")

#combine subset results into one dataframe
ciberres_all <- bind_rows(ciberres1, ciberres2, ciberres3, ciberres4)

# Remove duplicate rows based on sample names
ciberres <- ciberres_all %>%
  distinct(Mixture, .keep_all = TRUE)

write.csv(ciberres, "Cibersort Results.csv", quote = FALSE)

###Run analysis between IMC and Cibersort data
#select for the most important immune cells
ciberres_comp <- ciberres %>%
  select("Mixture", "B cells naive", "B cells memory", "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated", "T cells regulatory (Tregs)", "Macrophages M0", "Macrophages M1", "Macrophages M2", "Dendritic cells resting", "Dendritic cells activated")

#add together the different cell types to get a total for each cell type
ciberres_comp <- ciberres_comp %>%
  mutate(
    Bcells = `B cells naive` + `B cells memory`,
    CD4Tcells = `T cells CD4 naive` + `T cells CD4 memory resting` + `T cells CD4 memory activated`,
    macrophages = `Macrophages M0` + `Macrophages M1` + `Macrophages M2`,
    CD8Tcells = `T cells CD8`,
    Tregs = `T cells regulatory (Tregs)`,
    metabric_id = `Mixture`
  ) %>%
  select(metabric_id, Bcells, CD4Tcells, macrophages, `CD8Tcells`, `Tregs`)

#order the samples by Bcells
ciberres_Bcells_ranked <- ciberres_comp %>%
  as.data.frame() %>%
  arrange(Bcells)



###Plan: rank the metabric ids by their ranking for each immune cell from the datasets, then regress them.
#start with bcells

ciberres_Bcells_ranked <- inner_join(ciberres_Bcells_ranked, imageres_Bcells_ranked, by = "metabric_id")

ciber_bcellrank <- data.frame(
  metabric_id = ciberres_Bcells_ranked$metabric_id,
  BcellRanking = 1:nrow(ciberres_Bcells_ranked)
)

imc_bcellrank <- data.frame(
  metabric_id = imageres_Bcells_ranked$metabric_id,
  BcellRanking = 1:nrow(imageres_Bcells_ranked)
)

merged_bcellrank <- merge(ciber_bcellrank, imc_bcellrank, by.x = "metabric_id", by.y = "metabric_id", all = FALSE)
new_column_names <- c("metabric_id", "cibersort_rank", "imc_rank")
colnames(merged_bcellrank) <- new_column_names

regression_model <- lm(imc_rank ~ cibersort_rank, data = merged_bcellrank) #Run a linear regression
summary(regression_model) #Print regression summary

ggplot(merged_bcellrank, aes(x = cibersort_rank, y = imc_rank)) +
  geom_point() +
  labs(x = "cibersort_rank", y = "imc_rank") +
  ggtitle("Bcell proportion ranking adjusted")

#now with CD8 T cells

ciberres_CD8Tcells_ranked <- ciberres_comp %>%
  as.data.frame() %>%
  arrange(CD8Tcells)

imageres_CD8Tcells_ranked <- imageres_comp %>%
  as.data.frame() %>%
  arrange(CD8Tcells)

ciberres_CD8Tcells_ranked <- inner_join(ciberres_CD8Tcells_ranked, imageres_CD8Tcells_ranked, by = "metabric_id")

ciber_CD8Tcellrank <- data.frame(
  metabric_id = ciberres_CD8Tcells_ranked$metabric_id,
  CD8TcellRanking = 1:nrow(ciberres_CD8Tcells_ranked)
)

imc_CD8Tcellrank <- data.frame(
  metabric_id = imageres_CD8Tcells_ranked$metabric_id,
  CD8TcellRanking = 1:nrow(imageres_CD8Tcells_ranked)
)

merged_CD8Tcellrank <- merge(ciber_CD8Tcellrank, imc_CD8Tcellrank, by.x = "metabric_id", by.y = "metabric_id", all = FALSE)
new_column_names <- c("metabric_id", "cibersort_rank", "imc_rank")
colnames(merged_CD8Tcellrank) <- new_column_names

regression_model <- lm(imc_rank ~ cibersort_rank, data = merged_CD8Tcellrank) #Run a linear regression
summary(regression_model) #Print regression summary

ggplot(merged_CD8Tcellrank, aes(x = cibersort_rank, y = imc_rank)) +
  geom_point() +
  labs(x = "cibersort_rank", y = "imc_rank") +
  ggtitle("CD8Tcell proportion ranking adjusted")

#now with CD4 T cells
ciberres_CD4Tcells_ranked <- ciberres_comp %>%
  as.data.frame() %>%
  arrange(CD4Tcells)

imageres_CD4Tcells_ranked <- imageres_comp %>%
  as.data.frame() %>%
  arrange(CD4Tcells)

ciberres_CD4Tcells_ranked <- inner_join(ciberres_CD4Tcells_ranked, imageres_CD4Tcells_ranked, by = "metabric_id")

ciber_CD4Tcellrank <- data.frame(
  metabric_id = ciberres_CD4Tcells_ranked$metabric_id,
  CD4TcellRanking = 1:nrow(ciberres_CD4Tcells_ranked)
)

imc_CD4Tcellrank <- data.frame(
  metabric_id = imageres_CD4Tcells_ranked$metabric_id,
  CD4TcellRanking = 1:nrow(imageres_CD4Tcells_ranked)
)

merged_CD4Tcellrank <- merge(ciber_CD4Tcellrank, imc_CD4Tcellrank, by.x = "metabric_id", by.y = "metabric_id", all = FALSE)
new_column_names <- c("metabric_id", "cibersort_rank", "imc_rank")
colnames(merged_CD4Tcellrank) <- new_column_names

regression_model <- lm(imc_rank ~ cibersort_rank, data = merged_CD4Tcellrank) #Run a linear regression
summary(regression_model) #Print regression summary

ggplot(merged_CD4Tcellrank, aes(x = cibersort_rank, y = imc_rank)) +
  geom_point() +
  labs(x = "cibersort_rank", y = "imc_rank") +
  ggtitle("CD4Tcell proportion ranking adjusted")

#now with macrophages
ciberres_macrophages_ranked <- ciberres_comp %>%
  as.data.frame() %>%
  arrange(macrophages)

imageres_macrophages_ranked <- imageres_comp %>%
  as.data.frame() %>%
  arrange(macrophages)

ciberres_macrophages_ranked <- inner_join(ciberres_macrophages_ranked, imageres_macrophages_ranked, by = "metabric_id")

ciber_macrophagesrank <- data.frame(
  metabric_id = ciberres_macrophages_ranked$metabric_id,
  macrophagesRanking = 1:nrow(ciberres_macrophages_ranked)
)

imc_macrophagesrank <- data.frame(
  metabric_id = imageres_macrophages_ranked$metabric_id,
  macrophagesRanking = 1:nrow(imageres_macrophages_ranked)
)

merged_macrophagesrank <- merge(ciber_macrophagesrank, imc_macrophagesrank, by.x = "metabric_id", by.y = "metabric_id", all = FALSE)
new_column_names <- c("metabric_id", "cibersort_rank", "imc_rank")
colnames(merged_macrophagesrank) <- new_column_names

regression_model <- lm(imc_rank ~ cibersort_rank, data = merged_macrophagesrank) #Run a linear regression
summary(regression_model) #Print regression summary

ggplot(merged_macrophagesrank, aes(x = cibersort_rank, y = imc_rank)) +
  geom_point() +
  labs(x = "cibersort_rank", y = "imc_rank") +
  ggtitle("macrophages proportion ranking adjusted")

#now with tregs
ciberres_tregs_ranked <- ciberres_comp %>%
  as.data.frame() %>%
  arrange(Tregs)

imageres_tregs_ranked <- imageres_comp %>%
  as.data.frame() %>%
  arrange(Tregs)

ciberres_tregs_ranked <- inner_join(ciberres_tregs_ranked, imageres_tregs_ranked, by = "metabric_id")

ciber_tregsrank <- data.frame(
  metabric_id = ciberres_tregs_ranked$metabric_id,
  tregsRanking = 1:nrow(ciberres_tregs_ranked)
)

imc_tregsrank <- data.frame(
  metabric_id = imageres_tregs_ranked$metabric_id,
  tregsRanking = 1:nrow(imageres_tregs_ranked)
)

merged_tregsrank <- merge(ciber_tregsrank, imc_tregsrank, by.x = "metabric_id", by.y = "metabric_id", all = FALSE)
new_column_names <- c("metabric_id", "cibersort_rank", "imc_rank")
colnames(merged_tregsrank) <- new_column_names

regression_model <- lm(imc_rank ~ cibersort_rank, data = merged_tregsrank) #Run a linear regression
summary(regression_model) #Print regression summary

ggplot(merged_tregsrank, aes(x = cibersort_rank, y = imc_rank)) +
  geom_point() +
  labs(x = "cibersort_rank", y = "imc_rank") +
  ggtitle("Tregs proportion ranking adjusted")



###did not output good data. Will now try just regressing tumor/nontumor data from the transNEO trial
transNEO_MCPcounter <- fread("transneo-diagnosis-immune-MCPcounter.tsv")
transNEO_digpath <- fread("transneo-diagnosis-DigPathology.tsv")
colnames(transNEO_MCPcounter)[1] <- "Trial.ID"

#first, digpath's tumor fraction vs mcpcounter's endothelial cells
TransNEO_merged <- inner_join(transNEO_digpath, transNEO_MCPcounter, by = "Trial.ID")

TransNEO_merged <- TransNEO_merged %>%
  arrange(TransNEO_merged$"Endothelial cells") %>%
  mutate(mcpRank = row_number())

TransNEO_merged <- TransNEO_merged %>%
  arrange(TransNEO_merged$fraction_cancer) %>%
  mutate(digpathRank = row_number())

regression_model <- lm(mcpRank ~ digpathRank, data = TransNEO_merged)
coef <- coef(summary(regression_model))
intercept <- coef["(Intercept)", "Estimate"]
slope <- coef["digpathRank", "Estimate"]
r_squared <- summary(regression_model)$r.squared

tumorplot <- ggplot(TransNEO_merged, aes(x = mcpRank, y = digpathRank)) +
  geom_point() +
  labs(x = "mcpRank", y = "digpath Rank") +
  ggtitle("Tumor cells proportion ranking (endothelial only)")

tumorplot <- tumorplot +
  geom_abline(intercept = intercept, slope = slope, color = "red", linetype = "dashed") +
  annotate("text", x = max(TransNEO_merged$mcpRank), y = max(TransNEO_merged$digpathRank) - 5, 
           label = paste("R-squared =", round(r_squared, 4)), hjust = 1)

print(tumorplot)

#next, digpath's tumor fraction vs mcpcounter's endothelial + fibroblast cells
TransNEO_merged <- inner_join(transNEO_digpath, transNEO_MCPcounter, by = "Trial.ID")

TransNEO_merged <- TransNEO_merged %>%
  mutate(endofibro = TransNEO_merged$"Endothelial cells" + Fibroblasts)

TransNEO_merged <- TransNEO_merged %>%         
  arrange(TransNEO_merged$endofibro) %>%
  mutate(mcpRank = row_number())

TransNEO_merged <- TransNEO_merged %>%
  arrange(TransNEO_merged$"fraction_cancer") %>%
  mutate(digpathRank = row_number())

regression_model <- lm(mcpRank ~ digpathRank, data = TransNEO_merged)
coef <- coef(summary(regression_model))
intercept <- coef["(Intercept)", "Estimate"]
slope <- coef["digpathRank", "Estimate"]
r_squared <- summary(regression_model)$r.squared

tumorplot <- ggplot(TransNEO_merged, aes(x = mcpRank, y = digpathRank)) +
  geom_point() +
  labs(x = "mcpRank", y = "digpath Rank") +
  ggtitle("Tumor cells proportion ranking (fibroblasts + endothelial)")

tumorplot <- tumorplot +
  geom_abline(intercept = intercept, slope = slope, color = "red", linetype = "dashed") +
  annotate("text", x = max(TransNEO_merged$mcpRank), y = max(TransNEO_merged$digpathRank) - 5, 
           label = paste("R-squared =", round(r_squared, 4)), hjust = 1)

print(tumorplot)

#next, digpath's tumor fraction vs mcpcounter's fibroblast cells
TransNEO_merged <- inner_join(transNEO_digpath, transNEO_MCPcounter, by = "Trial.ID")

TransNEO_merged <- TransNEO_merged %>%         
  arrange(TransNEO_merged$Fibroblasts) %>%
  mutate(mcpRank = row_number())

TransNEO_merged <- TransNEO_merged %>%
  arrange(TransNEO_merged$"fraction_cancer") %>%
  mutate(digpathRank = row_number())

regression_model <- lm(mcpRank ~ digpathRank, data = TransNEO_merged)
coef <- coef(summary(regression_model))
intercept <- coef["(Intercept)", "Estimate"]
slope <- coef["digpathRank", "Estimate"]
r_squared <- summary(regression_model)$r.squared

tumorplot <- ggplot(TransNEO_merged, aes(x = mcpRank, y = digpathRank)) +
  geom_point() +
  labs(x = "mcpRank", y = "digpath Rank") +
  ggtitle("Tumor cells proportion ranking (fibroblasts only)")

tumorplot <- tumorplot +
  geom_abline(intercept = intercept, slope = slope, color = "red", linetype = "dashed") +
  annotate("text", x = max(TransNEO_merged$mcpRank), y = max(TransNEO_merged$digpathRank) - 5, 
           label = paste("R-squared =", round(r_squared, 4)), hjust = 1)

print(tumorplot)

#try ranking enothelial and fibroblasts as part of stroma
TransNEO_merged <- inner_join(transNEO_digpath, transNEO_MCPcounter, by = "Trial.ID")

TransNEO_merged <- TransNEO_merged %>%
  mutate(endofibro = TransNEO_merged$"Endothelial cells" + Fibroblasts)

TransNEO_merged <- TransNEO_merged %>%         
  arrange(TransNEO_merged$endofibro) %>%
  mutate(mcpRank = row_number())

TransNEO_merged <- TransNEO_merged %>%
  arrange(TransNEO_merged$"fraction_stroma") %>%
  mutate(digpathRank = row_number())

regression_model <- lm(mcpRank ~ digpathRank, data = TransNEO_merged)
coef <- coef(summary(regression_model))
intercept <- coef["(Intercept)", "Estimate"]
slope <- coef["digpathRank", "Estimate"]
r_squared <- summary(regression_model)$r.squared

tumorplot <- ggplot(TransNEO_merged, aes(x = mcpRank, y = digpathRank)) +
  geom_point() +
  labs(x = "mcpRank", y = "digpath Rank") +
  ggtitle("Stroma cells proportion ranking (fibroblasts + endothelial)")

tumorplot <- tumorplot +
  geom_abline(intercept = intercept, slope = slope, color = "red", linetype = "dashed") +
  annotate("text", x = max(TransNEO_merged$mcpRank), y = max(TransNEO_merged$digpathRank) - 5, 
           label = paste("R-squared =", round(r_squared, 4)), hjust = 1)

print(tumorplot)

#try ranking just enothelial as part of stroma
TransNEO_merged <- inner_join(transNEO_digpath, transNEO_MCPcounter, by = "Trial.ID")

TransNEO_merged <- TransNEO_merged %>%         
  arrange(TransNEO_merged$"Endothelial cells") %>%
  mutate(mcpRank = row_number())

TransNEO_merged <- TransNEO_merged %>%
  arrange(TransNEO_merged$"fraction_stroma") %>%
  mutate(digpathRank = row_number())

regression_model <- lm(mcpRank ~ digpathRank, data = TransNEO_merged)
coef <- coef(summary(regression_model))
intercept <- coef["(Intercept)", "Estimate"]
slope <- coef["digpathRank", "Estimate"]
r_squared <- summary(regression_model)$r.squared

tumorplot <- ggplot(TransNEO_merged, aes(x = mcpRank, y = digpathRank)) +
  geom_point() +
  labs(x = "mcpRank", y = "digpath Rank") +
  ggtitle("Stroma cells proportion ranking (endothelial)")

tumorplot <- tumorplot +
  geom_abline(intercept = intercept, slope = slope, color = "red", linetype = "dashed") +
  annotate("text", x = max(TransNEO_merged$mcpRank), y = max(TransNEO_merged$digpathRank) - 5, 
           label = paste("R-squared =", round(r_squared, 4)), hjust = 1)

print(tumorplot)

#try ranking all immune cells as part of lymph
TransNEO_merged <- inner_join(transNEO_digpath, transNEO_MCPcounter, by = "Trial.ID")

TransNEO_merged <- TransNEO_merged %>%
  mutate(lymph = rowSums(TransNEO_merged[, 9:16]))

TransNEO_merged <- TransNEO_merged %>%         
  arrange(TransNEO_merged$lymph) %>%
  mutate(mcpRank = row_number())

TransNEO_merged <- TransNEO_merged %>%
  arrange(TransNEO_merged$"fraction_lymph") %>%
  mutate(digpathRank = row_number())

regression_model <- lm(mcpRank ~ digpathRank, data = TransNEO_merged)
coef <- coef(summary(regression_model))
intercept <- coef["(Intercept)", "Estimate"]
slope <- coef["digpathRank", "Estimate"]
r_squared <- summary(regression_model)$r.squared

tumorplot <- ggplot(TransNEO_merged, aes(x = mcpRank, y = digpathRank)) +
  geom_point() +
  labs(x = "mcpRank", y = "digpath Rank") +
  ggtitle("Lymph cells proportion ranking")

tumorplot <- tumorplot +
  geom_abline(intercept = intercept, slope = slope, color = "red", linetype = "dashed") +
  annotate("text", x = max(TransNEO_merged$mcpRank), y = max(TransNEO_merged$digpathRank) - 5, 
           label = paste("R-squared =", round(r_squared, 4)), hjust = 1)

print(tumorplot)




###Try regressions without ranking
#try all immune cells as part of lymph (unranked)
TransNEO_merged <- inner_join(transNEO_digpath, transNEO_MCPcounter, by = "Trial.ID")

TransNEO_merged <- TransNEO_merged %>%
  mutate(lymph = rowSums(TransNEO_merged[, 9:16]))

regression_model <- lm(lymph ~ fraction_lymph, data = TransNEO_merged)
coef <- coef(summary(regression_model))
intercept <- coef["(Intercept)", "Estimate"]
slope <- coef["fraction_lymph", "Estimate"]
r_squared <- summary(regression_model)$r.squared

tumorplot <- ggplot(TransNEO_merged, aes(x = fraction_lymph, y = lymph)) +
  geom_point() +
  labs(x = "fraction_lymph", y = "mcp_lymph score") +
  ggtitle("Lymph cells proportion regression (no rank)")+
  scale_x_continuous(limits = c(0, 0.6))

tumorplot <- tumorplot +
  geom_abline(intercept = intercept, slope = slope, color = "red", linetype = "dashed") +
  annotate("text", x = 0.5, y = 1, 
           label = paste("R-squared =", round(r_squared, 4)), hjust = 1)

print(tumorplot)

#try endothelial cells as part of cancer (unranked)
TransNEO_merged <- inner_join(transNEO_digpath, transNEO_MCPcounter, by = "Trial.ID")

regression_model <- lm(TransNEO_merged$`Endothelial cells` ~ fraction_cancer, data = TransNEO_merged)
coef <- coef(summary(regression_model))
intercept <- coef["(Intercept)", "Estimate"]
slope <- coef["fraction_cancer", "Estimate"]
r_squared <- summary(regression_model)$r.squared

tumorplot <- ggplot(TransNEO_merged, aes(x = fraction_cancer, y = TransNEO_merged$`Endothelial cells`)) +
  geom_point() +
  labs(x = "fraction_cancer", y = "mcp_endothelial score") +
  ggtitle("Cancer cells proportion regression (no rank)")+
  scale_x_continuous(limits = c(0, 0.65))

tumorplot <- tumorplot +
  geom_abline(intercept = intercept, slope = slope, color = "red", linetype = "dashed") +
  annotate("text", x = 0.5, y = 1, 
           label = paste("R-squared =", round(r_squared, 4)), hjust = 1)

print(tumorplot)



###Will now try the immunedeconv package, which contains a bunch of deconvultion methods in one package
###Have to use anaconda, paste this information into anaconda
install.packages("remotes")
remotes::install_github("omnideconv/immunedeconv")

setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Bulk Expression Cell Type Quantification")
expressiondata <- fread("data_mrna_illumina_microarray.txt")
metabricexp <- expressiondata[, -2]
metabricexp[is.na(metabricexp)] <- 0
metabricexp <- metabricexp %>%
  distinct(Hugo_Symbol, .keep_all = TRUE)
rownames(metabricexp) <- metabricexp$Hugo_Symbol
metabricexp <- data.matrix(metabricexp)
metabricexp <- metabricexp[, -1]

quantiseq_res <- immunedeconv::deconvolute(metabricexp, "quantiseq", arrays = TRUE)
write.csv(quantiseq_res, "quantiseq Results.csv", quote = FALSE)

timer_res <- immunedeconv::deconvolute(metabricexp, "timer", indications = rep( "brca" , 1980 ))
write.csv(timer_res, "TIMER Results.csv", quote = FALSE)

mcpcounter_res <- immunedeconv::deconvolute(metabricexp, "mcp_counter", feature_types = "HUGO_symbols")
write.csv(mcpcounter_res, "MCPCounter Results.csv", quote = FALSE)

xCell.data <- xCell::xCell.data
xcell_res <- immunedeconv::deconvolute(metabricexp, "xcell", arrays = TRUE, expected_cell_types = c('B cell','T cell CD4+','T cell CD8+','Myeloid dendritic cell','Eosinophil','Macrophage','Monocyte','Mast cell','Neutrophil','NK cell','Endothelial cell','Cancer associated fibroblast'))
write.csv(xcell_res, "XCell Results.csv", quote = FALSE)

abis_res <- immunedeconv::deconvolute(metabricexp, "abis", arrays = TRUE)
write.csv(abis_res, "Abis Results.csv", quote = FALSE)

estimate_res <- immunedeconv::deconvolute(metabricexp, "estimate")
write.csv(estimate_res, "Estimate Results.csv", quote = FALSE)

#set_cibersort_binary("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Bulk Expression Cell Type Quantification/CIBERSORT.R")
#set_cibersort_mat("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Bulk Expression Cell Type Quantification/LM22.txt")

#cibersort_res <- immunedeconv::deconvolute(metabricexp, "cibersort", array = TRUE)
#write.csv(cibersort_res, "Immunedeconv Cibersort Results.csv", quote = FALSE)

#cibersortabs_res <- immunedeconv::deconvolute(metabricexp, "cibersort_abs", array = TRUE)
#write.csv(cibersortabs_res, "CibersortAbsolute Results.csv", quote = FALSE)

epic_res <- immunedeconv::deconvolute(metabricexp, "epic", tumor = TRUE)
write.csv(epic_res, "EPIC Results.csv", quote = FALSE)

consensustme_res <- immunedeconv::deconvolute(metabricexp, "consensus_tme", indications = rep( "brca" , 1980 ))
write.csv(consensustme_res, "ConsensusTME Results.csv", quote = FALSE)



###analysis of deconvolution analyses
#quantiseq
quantiseq <- fread("quantiseq Results.csv")
quantiseq <- t(quantiseq)
quantiseq <- quantiseq[-1, ] # Remove the first row
colnames(quantiseq) <- quantiseq[1, ] # Make the first row the column names
quantiseq <- quantiseq[-1, ] # Remove the first row
quantiseq <- as.data.frame(quantiseq) 
quantiseq$metabric_id <- rownames(quantiseq)
rownames(quantiseq) <- NULL
quantiseq <- quantiseq %>%
  pivot_longer(cols = -metabric_id,  # Columns to be stacked into rows (excluding metabric_id)
               names_to = "cell_type",  # New column name for the cell type
               values_to = "quantiseq_score")  # New column name for the cell proportion

imageres_quantiseq_match <- Proportionresult %>%
  mutate(cellPhenotype = case_when(
    cellPhenotype == "B cells" ~ "B cell",
    cellPhenotype == "CD8^{+} T cells" ~ "T cell CD8+",
    cellPhenotype == "CD4^{+} T cells" ~ "T cell CD4+ (non-regulatory)",
    cellPhenotype == "T_{Reg} & T_{Ex}" ~ "T cell regulatory (Tregs)",
    cellPhenotype == "Macrophages" ~ "Macrophage M1",
    TRUE ~ cellPhenotype  # Keep other cell types unchanged
  ))

imageres_quantiseq_match <- imageres_quantiseq_match %>%
  filter(cellPhenotype %in% c("B cell", "T cell CD8+", "T cell CD4+ (non-regulatory)", "T cell regulatory (Tregs)", "Macrophage M1")) %>%
  rename(cell_type = cellPhenotype) %>%
  select(-Total_Cells) %>%
  rename(imc_proportion = Proportion)

quantiseq_imc <- merge(quantiseq, imageres_quantiseq_match, by = c("metabric_id", "cell_type"))

quantiseq_imc$quantiseq_score <- as.numeric(quantiseq_imc$quantiseq_score)
quantiseq_imc$imc_proportion <- as.numeric(quantiseq_imc$imc_proportion)

quantiseq_plot <- ggplot(quantiseq_imc, aes(x = quantiseq_score, y = imc_proportion)) +
  geom_point() +  # Add points for the scatterplot
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression lines
  facet_wrap(~cell_type, scales = "free") +  # Facet by cell type
  labs(title = "Quantiseq vs IMC proportions",
       x = "quantiseq_score", y = "imc_proportion") +  # Customize labels
  theme_minimal() +  # Customize the theme (optional)
  scale_x_continuous(limits = c(0, max(quantiseq_imc$quantiseq_score)))  # Set x-axis limits

quantiseq_rsquared <- quantiseq_imc %>%
  group_by(cell_type) %>%
  summarize(r_squared = cor(quantiseq_score, imc_proportion)^2)

quantiseq_plot <- quantiseq_plot +
  geom_text(data = quantiseq_rsquared, aes(label = paste("R-squared =", round(r_squared, 2))),
            x = Inf, y = -Inf, hjust = 1, vjust = 0, size = 4)

print(quantiseq_plot)

#mcpCounter
MCPCounter <- fread("MCPCounter Results.csv")
MCPCounter <- t(MCPCounter)
MCPCounter <- MCPCounter[-1, ] # Remove the first row
colnames(MCPCounter) <- MCPCounter[1, ] # Make the first row the column names
MCPCounter <- MCPCounter[-1, ] # Remove the first row
MCPCounter <- as.data.frame(MCPCounter) 
MCPCounter$metabric_id <- rownames(MCPCounter)
rownames(MCPCounter) <- NULL
MCPCounter <- MCPCounter %>%
  pivot_longer(cols = -metabric_id,  # Columns to be stacked into rows (excluding metabric_id)
               names_to = "cell_type",  # New column name for the cell type
               values_to = "MCPcounter_score")  # New column name for the cell proportion

imageres_MCPcounter_match <- Proportionresult %>%
  mutate(cellPhenotype = case_when(
    cellPhenotype == "B cells" ~ "B cell",
    cellPhenotype == "CD8^{+} T cells" ~ "T cell CD8+",
    cellPhenotype == "Macrophages" ~ "Macrophage/Monocyte",
    cellPhenotype == "Fibroblasts" ~ "Cancer associated fibroblast",
    TRUE ~ cellPhenotype  # Keep other cell types unchanged
  ))

imageres_MCPcounter_match <- imageres_MCPcounter_match %>%
  filter(cellPhenotype %in% c("B cell", "T cell CD8+", "Macrophage/Monocyte", "Cancer associated fibroblast")) %>%
  rename(cell_type = cellPhenotype) %>%
  select(-Total_Cells) %>%
  rename(imc_proportion = Proportion)

MCPcounter_imc <- merge(MCPCounter, imageres_MCPcounter_match, by = c("metabric_id", "cell_type"))

MCPcounter_imc$MCPcounter_score <- as.numeric(MCPcounter_imc$MCPcounter_score)
MCPcounter_imc$imc_proportion <- as.numeric(MCPcounter_imc$imc_proportion)

MCPcounter_plot <- ggplot(MCPcounter_imc, aes(x = MCPcounter_score, y = imc_proportion)) +
  geom_point() +  # Add points for the scatterplot
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression lines
  facet_wrap(~cell_type, scales = "free") +  # Facet by cell type
  labs(title = "MCPcounter vs IMC proportions",
       x = "MCPcounter_score", y = "imc_proportion") +  # Customize labels
  theme_minimal() +  # Customize the theme (optional)
  scale_x_continuous(limits = c(min(MCPcounter_imc$MCPcounter_score), max(MCPcounter_imc$MCPcounter_score)))  # Set x-axis limits

MCPcounter_rsquared <- MCPcounter_imc %>%
  group_by(cell_type) %>%
  summarize(r_squared = cor(MCPcounter_score, imc_proportion)^2)

MCPcounter_plot <- MCPcounter_plot +
  geom_text(data = MCPcounter_rsquared, aes(label = paste("R-squared =", round(r_squared, 2))),
            x = Inf, y = -Inf, hjust = 1, vjust = 0, size = 4)

print(MCPcounter_plot)

#xCell
XCell <- fread("XCell Results.csv")
XCell <- t(XCell)
XCell <- XCell[-1, ] # Remove the first row
colnames(XCell) <- XCell[1, ] # Make the first row the column names
XCell <- XCell[-1, ] # Remove the first row
XCell <- as.data.frame(XCell) 
XCell$metabric_id <- rownames(XCell)
rownames(XCell) <- NULL
XCell <- XCell %>%
  pivot_longer(cols = -metabric_id,  # Columns to be stacked into rows (excluding metabric_id)
               names_to = "cell_type",  # New column name for the cell type
               values_to = "XCell_score")  # New column name for the cell proportion

imageres_XCell_match <- Proportionresult %>%
  mutate(cellPhenotype = case_when(
    cellPhenotype == "B cells" ~ "B cell",
    cellPhenotype == "CD8^{+} T cells" ~ "T cell CD8+",
    cellPhenotype == "CD4^{+} T cells" ~ "T cell CD4+ (non-regulatory)",
    cellPhenotype == "Macrophages" ~ "Macrophage",
    cellPhenotype == "T_{Reg} & T_{Ex}" ~ "T cell regulatory (Tregs)",
    cellPhenotype == "Fibroblasts" ~ "Cancer associated fibroblast",
    TRUE ~ cellPhenotype  # Keep other cell types unchanged
  ))

imageres_XCell_match <- imageres_XCell_match %>%
  filter(cellPhenotype %in% c("B cell", "T cell CD8+", "Macrophage", "Cancer associated fibroblast", "T cell regulatory (Tregs)", "T cell CD4+ (non-regulatory)")) %>%
  rename(cell_type = cellPhenotype) %>%
  select(-Total_Cells) %>%
  rename(imc_proportion = Proportion)

XCell_imc <- merge(XCell, imageres_XCell_match, by = c("metabric_id", "cell_type"))

XCell_imc$XCell_score <- as.numeric(XCell_imc$XCell_score)
XCell_imc$imc_proportion <- as.numeric(XCell_imc$imc_proportion)

XCell_plot <- ggplot(XCell_imc, aes(x = XCell_score, y = imc_proportion)) +
  geom_point() +  # Add points for the scatterplot
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression lines
  facet_wrap(~cell_type, scales = "free") +  # Facet by cell type
  labs(title = "XCell vs IMC proportions",
       x = "XCell_score", y = "imc_proportion") +  # Customize labels
  theme_minimal() +  # Customize the theme (optional)
  scale_x_continuous(limits = c(min(XCell_imc$XCell_score), max(XCell_imc$XCell_score)))  # Set x-axis limits

XCell_rsquared <- XCell_imc %>%
  group_by(cell_type) %>%
  summarize(r_squared = cor(XCell_score, imc_proportion)^2)

XCell_plot <- XCell_plot +
  geom_text(data = XCell_rsquared, aes(label = paste("R-squared =", round(r_squared, 2))),
            x = Inf, y = -Inf, hjust = 1, vjust = 0, size = 4)

print(XCell_plot)

#abis
Abis <- fread("Abis Results.csv")
Abis <- t(Abis)
Abis <- Abis[-1, ] # Remove the first row
colnames(Abis) <- Abis[1, ] # Make the first row the column names
Abis <- Abis[-1, ] # Remove the first row
Abis <- as.data.frame(Abis) 
Abis$metabric_id <- rownames(Abis)
rownames(Abis) <- NULL

Abis$"B cell naive" <- as.numeric(Abis$"B cell naive")
Abis$"B cell memory" <- as.numeric(Abis$"B cell memory")
Abis$"B cell plasma immature" <- as.numeric(Abis$"B cell plasma immature")
Abis <- Abis %>%
  mutate(abis_bcell_sum = `B cell memory` + `B cell plasma immature`)

imc_bcell <- Proportionresult %>% filter(cellPhenotype == "B cells")
Abis_imc <- merge(Abis, imc_bcell, by = c("metabric_id"))

Abis_imc$abis_bcell_sum <- as.numeric(Abis_imc$abis_bcell_sum)
Abis_imc$Proportion <- as.numeric(Abis_imc$Proportion)

Abis_plot <- ggplot(Abis_imc, aes(x = abis_bcell_sum, y = Proportion)) +
  geom_point() +  # Add points for the scatterplot
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression lines
  labs(title = "Abis vs IMC Bcell proportions (B cell memory + B cell plasma immature)",
       x = "Abis_bcell_sum", y = "imc_proportion") +  # Customize labels
  theme_minimal() +  # Customize the theme (optional)
  scale_x_continuous(limits = c(min(Abis_imc$abis_bcell_sum), max(Abis_imc$abis_bcell_su)))  # Set x-axis limits

Abis_rsquared <- Abis_imc %>%
  summarize(r_squared = cor(abis_bcell_sum, Proportion)^2)

Abis_plot <- Abis_plot +
  geom_text(data = Abis_rsquared, aes(label = paste("R-squared =", round(r_squared, 2))),
            x = Inf, y = -Inf, hjust = 1, vjust = 0, size = 4)

print(Abis_plot)


#Estimate
Estimate <- fread("Estimate Results.csv")
Estimate <- t(Estimate)
Estimate <- Estimate[-1, ] # Remove the first row
colnames(Estimate) <- Estimate[1, ] # Make the first row the column names
Estimate <- Estimate[-1, ] # Remove the first row
Estimate <- as.data.frame(Estimate) 
Estimate$metabric_id <- rownames(Estimate)
rownames(Estimate) <- NULL

selected_cell_phenotypes <- c("B cells", "CD8^{+} T cells", "CD4^{+} T cells", "Macrophages")
filtered_df <- Proportionresult %>%
  filter(cellPhenotype %in% selected_cell_phenotypes)
result_df <- filtered_df %>%
  group_by(metabric_id) %>%
  summarize("imc" = sum(Proportion))

Estimate_imc <- merge(Estimate, result_df, by = c("metabric_id"))
Estimate_imc <- Estimate_imc %>%
  rename(immune_score = "immune score")

Estimate_imc$imc <- as.numeric(Estimate_imc$imc)
Estimate_imc$immune_score <- as.numeric(Estimate_imc$immune_score)

Estimate_plot <- ggplot(Estimate_imc, aes(x = immune_score, y = imc)) +
  geom_point() +  # Add points for the scatterplot
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression lines
  labs(title = "Estimate vs IMC immune proportions",
       x = "Estimate_score", y = "imc_proportion") +  # Customize labels
  theme_minimal() +  # Customize the theme (optional)
  scale_x_continuous(limits = c(min(Estimate_imc$immune_score), max(Estimate_imc$immune_score)))  # Set x-axis limits

Estimate_rsquared <- Estimate_imc %>%
  summarize(r_squared = cor(immune_score, imc)^2)

Estimate_plot <- Estimate_plot +
  geom_text(data = Estimate_rsquared, aes(label = paste("R-squared =", round(r_squared, 2))),
            x = Inf, y = -Inf, hjust = 1, vjust = 0, size = 4)

print(Estimate_plot)

#TIMER
TIMER <- fread("TIMER Results.csv")
TIMER <- t(TIMER)
TIMER <- TIMER[-1, ] # Remove the first row
colnames(TIMER) <- TIMER[1, ] # Make the first row the column names
TIMER <- TIMER[-1, ] # Remove the first row
TIMER <- as.data.frame(TIMER) 
TIMER$metabric_id <- rownames(TIMER)
rownames(TIMER) <- NULL
TIMER <- TIMER %>%
  pivot_longer(cols = -metabric_id,  # Columns to be stacked into rows (excluding metabric_id)
               names_to = "cell_type",  # New column name for the cell type
               values_to = "TIMER_score")  # New column name for the cell proportion

imageres_TIMER_match <- Proportionresult %>%
  mutate(cellPhenotype = case_when(
    cellPhenotype == "B cells" ~ "B cell",
    cellPhenotype == "CD8^{+} T cells" ~ "T cell CD8+",
    cellPhenotype == "CD4^{+} T cells" ~ "T cell CD4+",
    cellPhenotype == "Macrophages" ~ "Macrophage",
    TRUE ~ cellPhenotype  # Keep other cell types unchanged
  ))

imageres_TIMER_match <- imageres_TIMER_match %>%
  filter(cellPhenotype %in% c("B cell", "T cell CD8+", "Macrophage", "T cell CD4+")) %>%
  rename(cell_type = cellPhenotype) %>%
  select(-Total_Cells) %>%
  rename(imc_proportion = Proportion)

TIMER_imc <- merge(TIMER, imageres_TIMER_match, by = c("metabric_id", "cell_type"))

TIMER_imc$TIMER_score <- as.numeric(TIMER_imc$TIMER_score)
TIMER_imc$imc_proportion <- as.numeric(TIMER_imc$imc_proportion)

TIMER_plot <- ggplot(TIMER_imc, aes(x = TIMER_score, y = imc_proportion)) +
  geom_point() +  # Add points for the scatterplot
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression lines
  facet_wrap(~cell_type, scales = "free") +  # Facet by cell type
  labs(title = "TIMER vs IMC proportions",
       x = "TIMER_score", y = "imc_proportion") +  # Customize labels
  theme_minimal() +  # Customize the theme (optional)
  scale_x_continuous(limits = c(min(TIMER_imc$TIMER_score), max(TIMER_imc$TIMER_score)))  # Set x-axis limits

TIMER_rsquared <- TIMER_imc %>%
  group_by(cell_type) %>%
  summarize(r_squared = cor(TIMER_score, imc_proportion)^2)

TIMER_plot <- TIMER_plot +
  geom_text(data = TIMER_rsquared, aes(label = paste("R-squared =", round(r_squared, 2))),
            x = Inf, y = -Inf, hjust = 1, vjust = 0, size = 4)

print(TIMER_plot)


#EPIC
EPIC <- fread("EPIC Results.csv")
EPIC <- t(EPIC)
EPIC <- EPIC[-1, ] # Remove the first row
colnames(EPIC) <- EPIC[1, ] # Make the first row the column names
EPIC <- EPIC[-1, ] # Remove the first row
EPIC <- as.data.frame(EPIC) 
EPIC$metabric_id <- rownames(EPIC)
rownames(EPIC) <- NULL
EPIC <- EPIC %>%
  pivot_longer(cols = -metabric_id,  # Columns to be stacked into rows (excluding metabric_id)
               names_to = "cell_type",  # New column name for the cell type
               values_to = "EPIC_score")  # New column name for the cell proportion

imageres_EPIC_match <- Proportionresult %>%
  mutate(cellPhenotype = case_when(
    cellPhenotype == "B cells" ~ "B cell",
    cellPhenotype == "CD8^{+} T cells" ~ "T cell CD8+",
    cellPhenotype == "CD4^{+} T cells" ~ "T cell CD4+",
    cellPhenotype == "Macrophages" ~ "Macrophage",
    cellPhenotype == "Fibroblasts" ~ "Cancer associated fibroblast",
    TRUE ~ cellPhenotype  # Keep other cell types unchanged
  ))

imageres_EPIC_match <- imageres_EPIC_match %>%
  filter(cellPhenotype %in% c("B cell", "T cell CD8+", "Macrophage", "Cancer associated fibroblast", "T cell CD4+")) %>%
  rename(cell_type = cellPhenotype) %>%
  select(-Total_Cells) %>%
  rename(imc_proportion = Proportion)

EPIC_imc <- merge(EPIC, imageres_EPIC_match, by = c("metabric_id", "cell_type"))

EPIC_imc$EPIC_score <- as.numeric(EPIC_imc$EPIC_score)
EPIC_imc$imc_proportion <- as.numeric(EPIC_imc$imc_proportion)

EPIC_plot <- ggplot(EPIC_imc, aes(x = EPIC_score, y = imc_proportion)) +
  geom_point() +  # Add points for the scatterplot
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression lines
  facet_wrap(~cell_type, scales = "free") +  # Facet by cell type
  labs(title = "EPIC vs IMC proportions",
       x = "EPIC_score", y = "imc_proportion") +  # Customize labels
  theme_minimal() +  # Customize the theme (optional)
  scale_x_continuous(limits = c(min(EPIC_imc$EPIC_score), max(EPIC_imc$EPIC_score)))  # Set x-axis limits

EPIC_rsquared <- EPIC_imc %>%
  group_by(cell_type) %>%
  summarize(r_squared = cor(EPIC_score, imc_proportion)^2)

EPIC_plot <- EPIC_plot +
  geom_text(data = EPIC_rsquared, aes(label = paste("R-squared =", round(r_squared, 2))),
            x = Inf, y = -Inf, hjust = 1, vjust = 0, size = 4)

print(EPIC_plot)

#Consensus
Consensus <- fread("ConsensusTME Results.csv")
Consensus <- t(Consensus)
Consensus <- Consensus[-1, ] # Remove the first row
colnames(Consensus) <- Consensus[1, ] # Make the first row the column names
Consensus <- Consensus[-1, ] # Remove the first row
Consensus <- as.data.frame(Consensus) 
Consensus$metabric_id <- rownames(Consensus)
rownames(Consensus) <- NULL
Consensus <- Consensus %>%
  pivot_longer(cols = -metabric_id,  # Columns to be stacked into rows (excluding metabric_id)
               names_to = "cell_type",  # New column name for the cell type
               values_to = "Consensus_score")  # New column name for the cell proportion

imageres_Consensus_match <- Proportionresult %>%
  mutate(cellPhenotype = case_when(
    cellPhenotype == "B cells" ~ "B cell",
    cellPhenotype == "CD8^{+} T cells" ~ "T cell CD8+",
    cellPhenotype == "CD4^{+} T cells" ~ "T cell CD4+ (non-regulatory)",
    cellPhenotype == "Macrophages" ~ "Macrophage",
    cellPhenotype == "T_{Reg} & T_{Ex}" ~ "T cell regulatory (Tregs)",
    cellPhenotype == "Fibroblasts" ~ "Cancer associated fibroblast",
    TRUE ~ cellPhenotype  # Keep other cell types unchanged
  ))

imageres_Consensus_match <- imageres_Consensus_match %>%
  filter(cellPhenotype %in% c("B cell", "T cell CD8+", "Macrophage", "Cancer associated fibroblast", "T cell regulatory (Tregs)", "T cell CD4+ (non-regulatory)")) %>%
  rename(cell_type = cellPhenotype) %>%
  select(-Total_Cells) %>%
  rename(imc_proportion = Proportion)

Consensus_imc <- merge(Consensus, imageres_Consensus_match, by = c("metabric_id", "cell_type"))

Consensus_imc$Consensus_score <- as.numeric(Consensus_imc$Consensus_score)
Consensus_imc$imc_proportion <- as.numeric(Consensus_imc$imc_proportion)

Consensus_plot <- ggplot(Consensus_imc, aes(x = Consensus_score, y = imc_proportion)) +
  geom_point() +  # Add points for the scatterplot
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression lines
  facet_wrap(~cell_type, scales = "free") +  # Facet by cell type
  labs(title = "Consensus vs IMC proportions",
       x = "Consensus_score", y = "imc_proportion") +  # Customize labels
  theme_minimal() +  # Customize the theme (optional)
  scale_x_continuous(limits = c(min(Consensus_imc$Consensus_score), max(Consensus_imc$Consensus_score)))  # Set x-axis limits

Consensus_rsquared <- Consensus_imc %>%
  group_by(cell_type) %>%
  summarize(r_squared = cor(Consensus_score, imc_proportion)^2)

Consensus_plot <- Consensus_plot +
  geom_text(data = Consensus_rsquared, aes(label = paste("R-squared =", round(r_squared, 2))),
            x = Inf, y = -Inf, hjust = 1, vjust = 0, size = 4)

print(Consensus_plot)






#Libraries####
library(data.table)
library(dplyr)
library(tidyr)
library(writexl)
library(readxl)
library(ggplot2)
library(wCorr)
library(Metrics)
library(semEff)
library(DescTools)
library(pheatmap)
library(textshape)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(LaplacesDemon)

###Input data####
#deconvolution results
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Deconvolution Results")
abis <- read_xlsx("ABIS Results.xlsx", sheet = 3)
cibersort <- read_xlsx("Cibersort Results.xlsx", sheet = 3)
consensus <- read_xlsx("ConsensusTME Results.xlsx", sheet = 3)
epic <- read_xlsx("EPIC Results.xlsx", sheet = 3)
estimate <- read_xlsx("ESTIMATE Results.xlsx", sheet = 3) #tissue level deconvolution
instaprism <- read_xlsx("metabric InstaPrism results tpm.xlsx", sheet = 5)
kassandra <- read_xlsx("Kassandra Results.xlsx", sheet = 3)
mcp <- read_xlsx("MCPCounter Results.xlsx", sheet = 3)
quantiseq <- read_xlsx("quantiseq Results.xlsx", sheet = 3)
scaden <- read_xlsx("Scaden Results.xlsx", sheet = 3)
timer <- read_xlsx("TIMER Results.xlsx", sheet = 3)
xcell <- read_xlsx("XCell Results.xlsx", sheet = 3)

#convert scores to fraction (FOR RMSE and THEILS U)
df_list <- list(abis, mcp, xcell, consensus, timer, kassandra)
normalize_df <- function(df) {
  df_to_normalize <- df[, -1]
  normalized_df <- apply(df_to_normalize, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  normalized_df <- t(t(normalized_df) / rowSums(normalized_df))
  row_sums <- rowSums(normalized_df)
  normalized_df <- sweep(normalized_df, 1, row_sums, FUN = '/')
  normalized_df <- cbind(df[, 1, drop = FALSE], normalized_df)
  return(normalized_df)
}
normalized_dataframes <- lapply(df_list, normalize_df)

abis <- normalized_dataframes[[1]]
mcp <- normalized_dataframes[[2]]
xcell <- normalized_dataframes[[3]]
consensus <- normalized_dataframes[[4]]
timer <- normalized_dataframes[[5]]
kassandra <- normalized_dataframes[[6]]


###IMC results
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
imc <- read_xlsx("IMC Cell Count and Fraction Data.xlsx")
#turn to wide-format
imc <- imc %>%
  filter(cellPhenotype != "CD38^{+} lymphocytes") %>%
  filter(cellPhenotype != "Ki67^{+}") %>%
  select(metabric_id, cellPhenotype, FractionOfTotalImmuneAndStromaCells, total_objects) %>%
  pivot_wider(names_from = cellPhenotype, values_from = FractionOfTotalImmuneAndStromaCells) %>%
  dplyr::rename(imc_CD4_T_Cells = "CD4^{+} T cells", 
         imc_CD8_T_Cells = "CD8^{+} T cells", 
         imc_B_cells = "B cells", 
         imc_Granulocytes = "Granulocytes", 
         imc_Macrophages = "Macrophages", 
         imc_Tregs = "T_{Reg} & T_{Ex}",
         imc_Fibroblasts = "Fibroblasts",
         imc_Endothelial = "Endothelial",
         imc_Epithelial = "Epithelial",
         imc_NK_Cells = "CD57^{+}") %>%
  replace(is.na(.), 0)

###Clinical data
setwd("/Users/kevintu/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
clinical <- fread('metabric_clinical.txt')
clinical <- clinical %>%
  dplyr::rename(metabric_id = `METABRIC.ID`) %>%
  select(metabric_id, Cellularity, Pam50Subtype)

# Function to capitalize the first letter
capitalize_first <- function(x) {
  sapply(x, function(y) paste0(toupper(substring(y, 1, 1)), tolower(substring(y, 2))))
}

# Apply the function to the Cellularity column
clinical[, Cellularity := capitalize_first(Cellularity)]

allimmune <- merge(imc, abis, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, instaprism, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, cibersort, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, epic, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, mcp, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, quantiseq, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, scaden, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, timer, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, xcell, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, consensus, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, clinical, by = "metabric_id", all = TRUE)
allimmune <- na.omit(allimmune)

setwd("/Users/kevintu/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")

###benchmarking results####
methods <- c("abis", "instaprism", "cibersort", "epic", "mcp", "quantiseq", "scaden", "timer", "xcell", "consensus")
results <- data.frame(Method = character(), CellType = character(), Pearson_R = numeric(), Spearman_R = numeric(), rmse = numeric(), adjpR2 = numeric(), adjsR2 = numeric(), theilU = numeric(), meanae = numeric(), medianae = numeric(), bias = numeric(), kld = numeric())

#remove imc_epithelial, only use if using FractionOfTotalImmuneAndStromaCells in imc data
allimmune <- allimmune %>%
  select(-imc_Epithelial)

# Nested loops to iterate over methods and immune cell types
for (method in methods) {
  print(method)
  imc_cells <- grep("^imc", names(allimmune), value = TRUE)
  method_cells <- grep(method, names(allimmune), value = TRUE)
  shared_cells <- intersect(sub("imc", "", imc_cells), sub(method, "", method_cells)) # Extracting immune cell names shared between 'imc' and 'abis'
  print(shared_cells)
  
  for (cell_type in shared_cells) {
    # Extract columns for 'imc' and current method and filter only rows with non-missing values
    imc_column <- paste0("imc", cell_type)
    method_column <- paste0(method, cell_type)
    imc_data <- allimmune[, imc_column]
    method_data <- allimmune[, method_column]
    
    # Calculate Pearson correlation using weightedCorr
    Pcorr <- weightedCorr(imc_data, method_data, method = "Pearson")
    Scorr <- weightedCorr(imc_data, method_data, method = "Spearman")
    wPcorr <- weightedCorr(imc_data, method_data, method = "Pearson", weights = allimmune$total_objects)
    wScorr <- weightedCorr(imc_data, method_data, method = "Spearman", weights = allimmune$total_objects)
    rmse <- rmse(imc_data, method_data)
    data <- data.frame(imc_data, method_data)
    model <- lm(imc_data ~ method_data, data = data)
    pR2 <- R2(model, type = "pearson")
    adjpR2 <- pR2["R.squared.adj"]
    sR2 <- R2(model, type = "spearman")
    adjsR2 <- sR2["R.squared.adj"]
    theilu <- TheilU(imc_data, method_data, type = 2, na.rm = FALSE)
    meanae <- mae(imc_data, method_data)
    medianae <- mdae(imc_data, method_data)
    bias <- bias(imc_data, method_data)
    kld <- KLD(imc_data, method_data)
    
    # Store results in the dataframe
    results <- rbind(results, c(method, cell_type, Pcorr, Scorr, wPcorr, wScorr, rmse, adjpR2, adjsR2, theilu, meanae, medianae, bias, kld[["mean.sum.KLD"]]))
  }
}
colnames(results) <- c("Method", "Immune_Cell", "Unweighted Pearson", "Unweighted Spearman", "Weighted Pearson", "Weighted Spearman", "RMSE", "adjusted pearson R2", "adjusted spearman R2", "Theil's U Index", "Mean Absolute Error", "Median Absolute Error", "Bias", "KLD")

###Kass results
imc_kassandra <- merge(imc, kassandra, by = "metabric_id") #kass only has 440 samples, put back in

print(method)
imc_cells <- grep("^imc", names(imc_kassandra), value = TRUE)
method_cells <- grep("^kassandra", names(imc_kassandra), value = TRUE)
shared_cells <- intersect(sub("imc", "", imc_cells), sub("kassandra", "", method_cells)) # Extracting immune cell names shared between 'imc' and 'abis'
print(shared_cells)
  
for (cell_type in shared_cells) {
    # Extract columns for 'imc' and current method and filter only rows with non-missing values
    imc_column <- paste0("imc", cell_type)
    method_column <- paste0("kassandra", cell_type)
    imc_data <- imc_kassandra[, imc_column]
    method_data <- imc_kassandra[, method_column]
    
    # Calculate Pearson correlation using weightedCorr
    Pcorr <- weightedCorr(imc_data, method_data, method = "Pearson")
    Scorr <- weightedCorr(imc_data, method_data, method = "Spearman")
    wPcorr <- weightedCorr(imc_data, method_data, method = "Pearson", weights = imc_kassandra$total_objects)
    wScorr <- weightedCorr(imc_data, method_data, method = "Spearman", weights = imc_kassandra$total_objects)
    rmse <- rmse(imc_data, method_data)
    data <- data.frame(imc_data, method_data)
    model <- lm(imc_data ~ method_data, data = data)
    pR2 <- R2(model, type = "pearson")
    adjpR2 <- pR2["R.squared.adj"]
    sR2 <- R2(model, type = "spearman")
    adjsR2 <- sR2["R.squared.adj"]
    theilu <- TheilU(imc_data, method_data, type = 2, na.rm = FALSE)
    meanae <- mae(imc_data, method_data)
    medianae <- mdae(imc_data, method_data)
    bias <- bias(imc_data, method_data)
    kld <- KLD(imc_data, method_data)
    
    # Store results in the dataframe
    method <- "kassandra"
    results <- rbind(results, c(method, cell_type, Pcorr, Scorr, wPcorr, wScorr, rmse, adjpR2, adjsR2, theilu, meanae, medianae, bias, kld[["mean.sum.KLD"]]))
}
#make columns 3-13 of results numeric
results[,3:13] <- sapply(results[,3:13], as.numeric)

###ploting results in pheatmap####
#remove epithelial and abis for now
  #results <- results[results$Immune_Cell != "Epithelial",]
results <- results[results$Method != "abis",]

#list of R correlations methods

methods_list <- c("Unweighted Pearson", "Unweighted Spearman", "Weighted Pearson", 
                  "Weighted Spearman", "adjusted pearson R2", "adjusted spearman R2")

# Loop over each R correlations method
for (method in methods_list) {
  #make heatmap matrix
  mat <- results %>%
    select(Method, Immune_Cell, !!method) %>%
    pivot_wider(names_from = Immune_Cell, values_from = !!method)
  
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$Method
  mat <- mat[,-1]
  mat <- as.matrix(mat)
  
  mat <- t(mat)
  
  # Sort rows and columns by average
  row_order <- order(rowMedians(mat, na.rm=TRUE), decreasing = TRUE)
  col_order <- order(colMedians(mat, na.rm=TRUE), decreasing = TRUE)
  
  # Create a new heatmap with transposed and sorted data
  pheatmap(
    mat[row_order, col_order],
    color = colorRampPalette(c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"))(100),
    na_col = "grey",
    name = "Value",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    display_numbers = T,
    fontsize_number = 12,
    main = paste(method),
    fontsize = 14,
    angle_col = 45
  )
}

# list of "lower is better" methods

methods_list <- c("RMSE", "Theil's U Index", "Mean Absolute Error", 
                  "Median Absolute Error", "Bias")

# Loop over each method
for (method in methods_list) {
  #make heatmap matrix
  mat <- results %>%
    select(Method, Immune_Cell, !!method) %>%
    pivot_wider(names_from = Immune_Cell, values_from = !!method)
  
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$Method
  mat <- mat[,-1]
  mat <- as.matrix(mat)
  
  mat <- t(mat)
  
  # Sort rows and columns by average
  row_order <- order(rowMedians(mat, na.rm=TRUE), decreasing = TRUE)
  col_order <- order(colMedians(mat, na.rm=TRUE), decreasing = TRUE)
  
  # Create a new heatmap with transposed and sorted data
  pheatmap(
    mat[row_order, col_order],
    na_col = "grey",
    color = colorRampPalette(c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"))(100),
    name = "Value",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    display_numbers = T,
    fontsize_number = 12,
    main = paste(method),
    fontsize = 14,
    angle_col = 45
  )
}



###benchmarking results, stratified by cellularity####
##low cellularity
methods <- c("abis", "instaprism", "cibersort", "epic", "mcp", "quantiseq", "scaden", "timer", "xcell", "consensus")
results_cellularity <- data.frame(Method = character(), CellType = character(), cellularity = character(), Pearson_R = numeric(), Spearman_R = numeric(), rmse = numeric(), adjpR2 = numeric(), adjsR2 = numeric(), theilU = numeric(), meanae = numeric(), medianae = numeric(), bias = numeric(), kld = numeric())

# Nested loops to iterate over methods and immune cell types in low cellularity
allimmune_cellular <- allimmune %>%
    filter(Cellularity == "Low")
  for (method in methods) {
    print(method)
    imc_cells <- grep("^imc", names(allimmune_cellular), value = TRUE)
    method_cells <- grep(method, names(allimmune_cellular), value = TRUE)
    shared_cells <- intersect(sub("imc", "", imc_cells), sub(method, "", method_cells)) # Extracting immune cell names shared between 'imc' and 'abis'
    print(shared_cells)
    
    for (cell_type in shared_cells) {
      # Extract columns for 'imc' and current method and filter only rows with non-missing values
      imc_column <- paste0("imc", cell_type)
      method_column <- paste0(method, cell_type)
      imc_data <- allimmune_cellular[, imc_column]
      method_data <- allimmune_cellular[, method_column]
      
      # Calculate Pearson correlation using weightedCorr
      Pcorr <- weightedCorr(imc_data, method_data, method = "Pearson")
      Scorr <- weightedCorr(imc_data, method_data, method = "Spearman")
      wPcorr <- weightedCorr(imc_data, method_data, method = "Pearson", weights = allimmune_cellular$total_objects)
      wScorr <- weightedCorr(imc_data, method_data, method = "Spearman", weights = allimmune_cellular$total_objects)
      rmse <- rmse(imc_data, method_data)
      data <- data.frame(imc_data, method_data)
      model <- lm(imc_data ~ method_data, data = data)
      pR2 <- R2(model, type = "pearson")
      adjpR2 <- pR2["R.squared.adj"]
      sR2 <- R2(model, type = "spearman")
      adjsR2 <- sR2["R.squared.adj"]
      theilu <- TheilU(imc_data, method_data, type = 2, na.rm = FALSE)
      meanae <- mae(imc_data, method_data)
      medianae <- mdae(imc_data, method_data)
      bias <- bias(imc_data, method_data)
      kld <- KLD(imc_data, method_data)
      
      # Store results in the dataframe
      results_cellularity <- rbind(results_cellularity, c(method, cell_type, "Low", Pcorr, Scorr, wPcorr, wScorr, rmse, adjpR2, adjsR2, theilu, meanae, medianae, bias, kld[["mean.sum.KLD"]]))
    }}

colnames(results_cellularity) <- c("Method", "Immune_Cell", "Cellularity", "Unweighted Pearson", "Unweighted Spearman", "Weighted Pearson", "Weighted Spearman", "RMSE", "adjusted pearson R2", "adjusted spearman R2", "Theil's U Index", "Mean Absolute Error", "Median Absolute Error", "Bias", "KLD")

# Nested loops to iterate over methods and immune cell types in moderate cellularity
allimmune_cellular <- allimmune %>%
  filter(Cellularity == "Moderate")
for (method in methods) {
  print(method)
  imc_cells <- grep("^imc", names(allimmune_cellular), value = TRUE)
  method_cells <- grep(method, names(allimmune_cellular), value = TRUE)
  shared_cells <- intersect(sub("imc", "", imc_cells), sub(method, "", method_cells)) # Extracting immune cell names shared between 'imc' and 'abis'
  print(shared_cells)
  
  for (cell_type in shared_cells) {
    # Extract columns for 'imc' and current method and filter only rows with non-missing values
    imc_column <- paste0("imc", cell_type)
    method_column <- paste0(method, cell_type)
    imc_data <- allimmune_cellular[, imc_column]
    method_data <- allimmune_cellular[, method_column]
    
    # Calculate Pearson correlation using weightedCorr
    Pcorr <- weightedCorr(imc_data, method_data, method = "Pearson")
    Scorr <- weightedCorr(imc_data, method_data, method = "Spearman")
    wPcorr <- weightedCorr(imc_data, method_data, method = "Pearson", weights = allimmune_cellular$total_objects)
    wScorr <- weightedCorr(imc_data, method_data, method = "Spearman", weights = allimmune_cellular$total_objects)
    rmse <- rmse(imc_data, method_data)
    data <- data.frame(imc_data, method_data)
    model <- lm(imc_data ~ method_data, data = data)
    pR2 <- R2(model, type = "pearson")
    adjpR2 <- pR2["R.squared.adj"]
    sR2 <- R2(model, type = "spearman")
    adjsR2 <- sR2["R.squared.adj"]
    theilu <- TheilU(imc_data, method_data, type = 2, na.rm = FALSE)
    meanae <- mae(imc_data, method_data)
    medianae <- mdae(imc_data, method_data)
    bias <- bias(imc_data, method_data)
    kld <- KLD(imc_data, method_data)
    
    # Store results in the dataframe
    results_cellularity <- rbind(results_cellularity, c(method, cell_type, "Moderate", Pcorr, Scorr, wPcorr, wScorr, rmse, adjpR2, adjsR2, theilu, meanae, medianae, bias, kld[["mean.sum.KLD"]]))
  }}

# Nested loops to iterate over methods and immune cell types in High cellularity
allimmune_cellular <- allimmune %>%
  filter(Cellularity == "High")
for (method in methods) {
  print(method)
  imc_cells <- grep("^imc", names(allimmune_cellular), value = TRUE)
  method_cells <- grep(method, names(allimmune_cellular), value = TRUE)
  shared_cells <- intersect(sub("imc", "", imc_cells), sub(method, "", method_cells)) # Extracting immune cell names shared between 'imc' and 'abis'
  print(shared_cells)
  
  for (cell_type in shared_cells) {
    # Extract columns for 'imc' and current method and filter only rows with non-missing values
    imc_column <- paste0("imc", cell_type)
    method_column <- paste0(method, cell_type)
    imc_data <- allimmune_cellular[, imc_column]
    method_data <- allimmune_cellular[, method_column]
    
    # Calculate Pearson correlation using weightedCorr
    Pcorr <- weightedCorr(imc_data, method_data, method = "Pearson")
    Scorr <- weightedCorr(imc_data, method_data, method = "Spearman")
    wPcorr <- weightedCorr(imc_data, method_data, method = "Pearson", weights = allimmune_cellular$total_objects)
    wScorr <- weightedCorr(imc_data, method_data, method = "Spearman", weights = allimmune_cellular$total_objects)
    rmse <- rmse(imc_data, method_data)
    data <- data.frame(imc_data, method_data)
    model <- lm(imc_data ~ method_data, data = data)
    pR2 <- R2(model, type = "pearson")
    adjpR2 <- pR2["R.squared.adj"]
    sR2 <- R2(model, type = "spearman")
    adjsR2 <- sR2["R.squared.adj"]
    theilu <- TheilU(imc_data, method_data, type = 2, na.rm = FALSE)
    meanae <- mae(imc_data, method_data)
    medianae <- mdae(imc_data, method_data)
    bias <- bias(imc_data, method_data)
    kld <- KLD(imc_data, method_data)
    
    # Store results in the dataframe
    results_cellularity <- rbind(results_cellularity, c(method, cell_type, "High", Pcorr, Scorr, wPcorr, wScorr, rmse, adjpR2, adjsR2, theilu, meanae, medianae, bias, kld[["mean.sum.KLD"]]))
  }}

###Kass results
imc_kassandra <- allimmune %>%
  merge(kassandra, by = "metabric_id") #kass only has 440 samples, put back in
  
print(method)
imc_cells <- grep("^imc", names(imc_kassandra), value = TRUE)
method_cells <- grep("^kassandra", names(imc_kassandra), value = TRUE)
shared_cells <- intersect(sub("imc", "", imc_cells), sub("kassandra", "", method_cells)) # Extracting immune cell names shared between 'imc' and 'abis'
print(shared_cells)

#loops to iterate over kassandra and immune cell types in High cellularity
imc_kassandra_cellular <- imc_kassandra %>%
  filter(Cellularity == "High")
for (cell_type in shared_cells) {
  # Extract columns for 'imc' and current method and filter only rows with non-missing values
  imc_column <- paste0("imc", cell_type)
  method_column <- paste0("kassandra", cell_type)
  imc_data <- imc_kassandra_cellular[, imc_column]
  method_data <- imc_kassandra_cellular[, method_column]
  
  # Calculate Pearson correlation using weightedCorr
  Pcorr <- weightedCorr(imc_data, method_data, method = "Pearson")
  Scorr <- weightedCorr(imc_data, method_data, method = "Spearman")
  wPcorr <- weightedCorr(imc_data, method_data, method = "Pearson", weights = imc_kassandra_cellular$total_objects)
  wScorr <- weightedCorr(imc_data, method_data, method = "Spearman", weights = imc_kassandra_cellular$total_objects)
  rmse <- rmse(imc_data, method_data)
  data <- data.frame(imc_data, method_data)
  model <- lm(imc_data ~ method_data, data = data)
  pR2 <- R2(model, type = "pearson")
  adjpR2 <- pR2["R.squared.adj"]
  sR2 <- R2(model, type = "spearman")
  adjsR2 <- sR2["R.squared.adj"]
  theilu <- TheilU(imc_data, method_data, type = 2, na.rm = FALSE)
  meanae <- mae(imc_data, method_data)
  medianae <- mdae(imc_data, method_data)
  bias <- bias(imc_data, method_data)
  kld <- KLD(imc_data, method_data)
  
  # Store results in the dataframe
  method <- "kassandra"
  results_cellularity <- rbind(results_cellularity, c(method, cell_type, "High", Pcorr, Scorr, wPcorr, wScorr, rmse, adjpR2, adjsR2, theilu, meanae, medianae, bias, kld[["mean.sum.KLD"]]))
}

#loops to iterate over kassandra and immune cell types in Moderate cellularity
imc_kassandra_cellular <- imc_kassandra %>%
  filter(Cellularity == "Moderate")
for (cell_type in shared_cells) {
  # Extract columns for 'imc' and current method and filter only rows with non-missing values
  imc_column <- paste0("imc", cell_type)
  method_column <- paste0("kassandra", cell_type)
  imc_data <- imc_kassandra_cellular[, imc_column]
  method_data <- imc_kassandra_cellular[, method_column]
  
  # Calculate Pearson correlation using weightedCorr
  Pcorr <- weightedCorr(imc_data, method_data, method = "Pearson")
  Scorr <- weightedCorr(imc_data, method_data, method = "Spearman")
  wPcorr <- weightedCorr(imc_data, method_data, method = "Pearson", weights = imc_kassandra_cellular$total_objects)
  wScorr <- weightedCorr(imc_data, method_data, method = "Spearman", weights = imc_kassandra_cellular$total_objects)
  rmse <- rmse(imc_data, method_data)
  data <- data.frame(imc_data, method_data)
  model <- lm(imc_data ~ method_data, data = data)
  pR2 <- R2(model, type = "pearson")
  adjpR2 <- pR2["R.squared.adj"]
  sR2 <- R2(model, type = "spearman")
  adjsR2 <- sR2["R.squared.adj"]
  theilu <- TheilU(imc_data, method_data, type = 2, na.rm = FALSE)
  meanae <- mae(imc_data, method_data)
  medianae <- mdae(imc_data, method_data)
  bias <- bias(imc_data, method_data)
  kld <- KLD(imc_data, method_data)
  
  # Store results in the dataframe
  method <- "kassandra"
  results_cellularity <- rbind(results_cellularity, c(method, cell_type, "Moderate", Pcorr, Scorr, wPcorr, wScorr, rmse, adjpR2, adjsR2, theilu, meanae, medianae, bias, kld[["mean.sum.KLD"]]))
}

#loops to iterate over kassandra and immune cell types in Low cellularity
imc_kassandra_cellular <- imc_kassandra %>%
  filter(Cellularity == "Low") %>%
  select(-imc_NK_Cells) # "Remove' NK_Cells as it is not present in the imc data among the kassandra samples
imc_cells <- grep("^imc", names(imc_kassandra_cellular), value = TRUE)
method_cells <- grep("^kassandra", names(imc_kassandra_cellular), value = TRUE)
shared_cells <- intersect(sub("imc", "", imc_cells), sub("kassandra", "", method_cells)) # Extracting immune cell names shared between 'imc' and method
for (cell_type in shared_cells) {
  # Extract columns for 'imc' and current method and filter only rows with non-missing values
  imc_column <- paste0("imc", cell_type)
  method_column <- paste0("kassandra", cell_type)
  imc_data <- imc_kassandra_cellular[, imc_column]
  method_data <- imc_kassandra_cellular[, method_column]
  
  # Calculate Pearson correlation using weightedCorr
  Pcorr <- weightedCorr(imc_data, method_data, method = "Pearson")
  Scorr <- weightedCorr(imc_data, method_data, method = "Spearman")
  wPcorr <- weightedCorr(imc_data, method_data, method = "Pearson", weights = imc_kassandra_cellular$total_objects)
  wScorr <- weightedCorr(imc_data, method_data, method = "Spearman", weights = imc_kassandra_cellular$total_objects)
  rmse <- rmse(imc_data, method_data)
  data <- data.frame(imc_data, method_data)
  model <- lm(imc_data ~ method_data, data = data)
  pR2 <- R2(model, type = "pearson")
  adjpR2 <- pR2["R.squared.adj"]
  sR2 <- R2(model, type = "spearman")
  adjsR2 <- sR2["R.squared.adj"]
  theilu <- TheilU(imc_data, method_data, type = 2, na.rm = FALSE)
  meanae <- mae(imc_data, method_data)
  medianae <- mdae(imc_data, method_data)
  bias <- bias(imc_data, method_data)
  kld <- KLD(imc_data, method_data)
  
  # Store results in the dataframe
  method <- "kassandra"
  results_cellularity <- rbind(results_cellularity, c(method, cell_type, "Low", Pcorr, Scorr, wPcorr, wScorr, rmse, adjpR2, adjsR2, theilu, meanae, medianae, bias, kld[["mean.sum.KLD"]]))
}

#make columns 3-13 of results numeric
results_cellularity[,4:14] <- sapply(results_cellularity[,4:14], as.numeric)

###ploting results in pheatmap, , stratified by cellularity####
#remove epithelial and abis for now
#results <- results[results$Immune_Cell != "Epithelial",]
results_cellularity <- results_cellularity[results$Method != "abis",]

#list of R correlations methods

methods_list <- c("Unweighted Pearson", "Unweighted Spearman", "Weighted Pearson", 
                  "Weighted Spearman", "adjusted pearson R2", "adjusted spearman R2")
cellularity <- c("High", "Moderate", "Low")

# Loop over each R correlations method
for (cellular in cellularity) {
for (method in methods_list) {
  #make heatmap matrix
  mat <- results_cellularity %>%
    filter(Cellularity == cellular) %>%
    
    select(Method, Immune_Cell, !!method) %>%
    pivot_wider(names_from = Immune_Cell, values_from = !!method)
  
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$Method
  mat <- mat[,-1]
  mat <- as.matrix(mat)
  
  mat <- t(mat)
  
  # Sort rows and columns by average
  row_order <- order(rowMedians(mat, na.rm=TRUE), decreasing = TRUE)
  col_order <- order(colMedians(mat, na.rm=TRUE), decreasing = TRUE)
  
  # Create a new heatmap with transposed and sorted data
  pheatmap(
    mat[row_order, col_order],
    color = colorRampPalette(c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"))(100),
    na_col = "grey",
    name = "Value",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    display_numbers = T,
    fontsize_number = 12,
    main = paste(method, cellular, sep = ", Cellularity "),
    fontsize = 14,
    angle_col = 45
  )
}}

# list of "lower is better" methods

methods_list <- c("RMSE", "Theil's U Index", "Mean Absolute Error", 
                  "Median Absolute Error", "Bias")

# Loop over each method
for (cellular in cellularity) {
for (method in methods_list) {
  #make heatmap matrix
  mat <- results_cellularity %>%
    filter(Cellularity == cellular) %>%
    select(Method, Immune_Cell, !!method) %>%
    pivot_wider(names_from = Immune_Cell, values_from = !!method)
  
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$Method
  mat <- mat[,-1]
  mat <- as.matrix(mat)
  
  mat <- t(mat)
  
  # Sort rows and columns by average
  row_order <- order(rowMedians(mat, na.rm=TRUE), decreasing = TRUE)
  col_order <- order(colMedians(mat, na.rm=TRUE), decreasing = TRUE)

  col_fun = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"))
  
  # Create a new heatmap with transposed and sorted data
  mat_new <- mat[row_order, col_order]
  
  plot <- Heatmap(
    mat_new,
    col = col_fun,
    na_col = "grey",
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.3f", mat_new[i, j]), x, y, gp = gpar(fontsize = 10))},
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = paste(cellular, "Cellularity"),
    rect_gp = gpar(col = "white", lwd = 2),
    name = method
  )
  print(plot)
}}

#plotting MAE heatmap (stratified by cellularity)####
#rename methods and cells
results_cellularity$Method <- case_when(
  results_cellularity$Method == "abis" ~ "Abis",
  results_cellularity$Method == "instaprism" ~ "InstaPrism",
  results_cellularity$Method == "cibersort" ~ "Cibersort",
  results_cellularity$Method == "epic" ~ "Epic",
  results_cellularity$Method == "mcp" ~ "MCPCounter",
  results_cellularity$Method == "quantiseq" ~ "Quantiseq",
  results_cellularity$Method == "scaden" ~ "Scaden",
  results_cellularity$Method == "timer" ~ "TIMER",
  results_cellularity$Method == "xcell" ~ "xCell",
  results_cellularity$Method == "consensus" ~ "ConsensusTME",
  results_cellularity$Method == "kassandra" ~ "Kassandra",
  TRUE ~ results_cellularity$Method
)
results_cellularity$Immune_Cell <- case_when(
  results_cellularity$Immune_Cell == "_B_cells" ~ "B Cells",
  results_cellularity$Immune_Cell == "_Granulocytes" ~ "Granulocytes",
  results_cellularity$Immune_Cell == "_NK_Cells" ~ "NK Cells",
  results_cellularity$Immune_Cell == "_CD8_T_Cells" ~ "CD8 T Cells",
  results_cellularity$Immune_Cell == "_Endothelial" ~ "Endothelial",
  results_cellularity$Immune_Cell == "_Fibroblasts" ~ "Fibroblasts",
  results_cellularity$Immune_Cell == "_Macrophages" ~ "Macrophages",
  results_cellularity$Immune_Cell == "_CD4_T_Cells" ~ "CD4 T Cells",
  results_cellularity$Immune_Cell == "_Tregs" ~ "Tregs",
  TRUE ~ results_cellularity$Immune_Cell
)

method_order <- c("Abis", "TIMER", "Scaden", "Kassandra", "ConsensusTME", "MCPCounter", "Epic", "xCell", "Cibersort", "Quantiseq", "InstaPrism")
results_cellularity$Method <- factor(results_cellularity$Method, levels = method_order)

# Function to create heatmap
create_heatmap <- function(data, title) {
  mat <- data %>%
    select(Method, Immune_Cell, "Median Absolute Error") %>%
    pivot_wider(names_from = Immune_Cell, values_from = "Median Absolute Error")
  
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$Method
  mat <- mat[,-1]
  mat <- as.matrix(mat)
  
  mat <- t(mat)
  
  # Sort rows and columns by average
  row_order <- order(rowMedians(mat, na.rm=TRUE), decreasing = TRUE)
  #col_order <- order(colMedians(mat, na.rm=TRUE), decreasing = TRUE)
  
  # Create a new heatmap with transposed and sorted data
  col_fun = colorRamp2(c(0, 0.65), c("#ffffff", "#DB162F"))
  
  mat_new <- mat[row_order, method_order]
  
  Heatmap(
    mat_new,
    col = col_fun,
    na_col = "grey",
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.3f", mat_new[i, j]), x, y, gp = gpar(fontsize = 10))},
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title = title,
    rect_gp = gpar(col = "white", lwd = 2),
    name = "Median Absolute Error"
  )
}

# Create heatmaps for each cellularity level
heatmap_low <- create_heatmap(results_cellularity %>% filter(Cellularity == "Low"), "Low Cellularity")
heatmap_mod <- create_heatmap(results_cellularity %>% filter(Cellularity == "Moderate"), "Moderate Cellularity")
heatmap_high <- create_heatmap(results_cellularity %>% filter(Cellularity == "High"), "High Cellularity")

heatmap_low %v% heatmap_mod %v% heatmap_high
#plotting RMSE (stratified by cellularity)####
#rename methods and cells
results_cellularity$Method <- case_when(
  results_cellularity$Method == "abis" ~ "Abis",
  results_cellularity$Method == "instaprism" ~ "InstaPrism",
  results_cellularity$Method == "cibersort" ~ "Cibersort",
  results_cellularity$Method == "epic" ~ "Epic",
  results_cellularity$Method == "mcp" ~ "MCPCounter",
  results_cellularity$Method == "quantiseq" ~ "Quantiseq",
  results_cellularity$Method == "scaden" ~ "Scaden",
  results_cellularity$Method == "timer" ~ "TIMER",
  results_cellularity$Method == "xcell" ~ "xCell",
  results_cellularity$Method == "consensus" ~ "ConsensusTME",
  results_cellularity$Method == "kassandra" ~ "Kassandra",
  TRUE ~ results_cellularity$Method
)
results_cellularity$Immune_Cell <- case_when(
  results_cellularity$Immune_Cell == "_B_cells" ~ "B Cells",
  results_cellularity$Immune_Cell == "_Granulocytes" ~ "Granulocytes",
  results_cellularity$Immune_Cell == "_NK_Cells" ~ "NK Cells",
  results_cellularity$Immune_Cell == "_CD8_T_Cells" ~ "CD8 T Cells",
  results_cellularity$Immune_Cell == "_Endothelial" ~ "Endothelial",
  results_cellularity$Immune_Cell == "_Fibroblasts" ~ "Fibroblasts",
  results_cellularity$Immune_Cell == "_Macrophages" ~ "Macrophages",
  results_cellularity$Immune_Cell == "_CD4_T_Cells" ~ "CD4 T Cells",
  results_cellularity$Immune_Cell == "_Tregs" ~ "Tregs",
  TRUE ~ results_cellularity$Immune_Cell
)

method_order <- c("Abis", "TIMER", "Epic", "Kassandra", "ConsensusTME", "MCPCounter", "xCell", "Scaden", "Cibersort", "InstaPrism", "Quantiseq")
results_cellularity$Method <- factor(results_cellularity$Method, levels = method_order)

# Function to create heatmap
create_heatmap <- function(data, title) {
  mat <- data %>%
    select(Method, Immune_Cell, "RMSE") %>%
    pivot_wider(names_from = Immune_Cell, values_from = "RMSE")
  
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$Method
  mat <- mat[,-1]
  mat <- as.matrix(mat)
  
  mat <- t(mat)
  
  # Sort rows and columns by average
  row_order <- order(rowMedians(mat, na.rm=TRUE), decreasing = TRUE)
  #col_order <- order(colMedians(mat, na.rm=TRUE), decreasing = TRUE)
  
  # Create a new heatmap with transposed and sorted data
  col_fun = colorRamp2(c(0, 0.65), c("#ffffff", "#DB162F"))
  
  mat_new <- mat[row_order, method_order]
  
  Heatmap(
    mat_new,
    col = col_fun,
    na_col = "grey",
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.3f", mat_new[i, j]), x, y, gp = gpar(fontsize = 10))},
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title = title,
    rect_gp = gpar(col = "white", lwd = 2),
    name = "RMSE"
  )
}

# Create heatmaps for each cellularity level
heatmap_low <- create_heatmap(results_cellularity %>% filter(Cellularity == "Low"), "Low Cellularity")
heatmap_mod <- create_heatmap(results_cellularity %>% filter(Cellularity == "Moderate"), "Moderate Cellularity")
heatmap_high <- create_heatmap(results_cellularity %>% filter(Cellularity == "High"), "High Cellularity")

heatmap_low %v% heatmap_mod %v% heatmap_high

#plotting spearman (stratified by cellularity)####
#rename methods and cells
results_cellularity$Method <- case_when(
  results_cellularity$Method == "abis" ~ "Abis",
  results_cellularity$Method == "instaprism" ~ "InstaPrism",
  results_cellularity$Method == "cibersort" ~ "Cibersort",
  results_cellularity$Method == "epic" ~ "Epic",
  results_cellularity$Method == "mcp" ~ "MCPCounter",
  results_cellularity$Method == "quantiseq" ~ "Quantiseq",
  results_cellularity$Method == "scaden" ~ "Scaden",
  results_cellularity$Method == "timer" ~ "TIMER",
  results_cellularity$Method == "xcell" ~ "xCell",
  results_cellularity$Method == "consensus" ~ "ConsensusTME",
  results_cellularity$Method == "kassandra" ~ "Kassandra",
  TRUE ~ results_cellularity$Method
)
results_cellularity$Immune_Cell <- case_when(
  results_cellularity$Immune_Cell == "_B_cells" ~ "B Cells",
  results_cellularity$Immune_Cell == "_Granulocytes" ~ "Granulocytes",
  results_cellularity$Immune_Cell == "_NK_Cells" ~ "NK Cells",
  results_cellularity$Immune_Cell == "_CD8_T_Cells" ~ "CD8 T Cells",
  results_cellularity$Immune_Cell == "_Endothelial" ~ "Endothelial",
  results_cellularity$Immune_Cell == "_Fibroblasts" ~ "Fibroblasts",
  results_cellularity$Immune_Cell == "_Macrophages" ~ "Macrophages",
  results_cellularity$Immune_Cell == "_CD4_T_Cells" ~ "CD4 T Cells",
  results_cellularity$Immune_Cell == "_Tregs" ~ "Tregs",
  TRUE ~ results_cellularity$Immune_Cell
)

method_order <- c("Cibersort", "InstaPrism", "Scaden", "Abis", "TIMER", "Kassandra", "ConsensusTME", "Quantiseq", "MCPCounter", "xCell", "Epic") 
results_cellularity$Method <- factor(results_cellularity$Method, levels = method_order)

# Function to create heatmap
create_heatmap <- function(data, title) {
  mat <- data %>%
    select(Method, Immune_Cell, "Unweighted Spearman") %>%
    pivot_wider(names_from = Immune_Cell, values_from = "Unweighted Spearman")
  
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$Method
  mat <- mat[,-1]
  mat <- as.matrix(mat)
  
  mat <- t(mat)
  
  # Sort rows and columns by average
  row_order <- order(rowMedians(mat, na.rm=TRUE))
  col_order <- order(colMedians(mat, na.rm=TRUE), decreasing = TRUE)
  print(col_order)
  # Create a new heatmap with transposed and sorted data
  col_fun = colorRamp2(c(-0.8, 0, 0.6), c("#DB162F", "#ffffff", "#0E6BA8"))
  
  mat_new <- mat[row_order, method_order]
  
  Heatmap(
    mat_new,
    col = col_fun,
    na_col = "grey",
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.3f", mat_new[i, j]), x, y, gp = gpar(fontsize = 10))},
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title = title,
    rect_gp = gpar(col = "white", lwd = 2),
    name = "Unweighted Spearman"
  )
}

# Create heatmaps for each cellularity level
heatmap_low <- create_heatmap(results_cellularity %>% filter(Cellularity == "Low"), "Low Cellularity")
heatmap_mod <- create_heatmap(results_cellularity %>% filter(Cellularity == "Moderate"), "Moderate Cellularity")
heatmap_high <- create_heatmap(results_cellularity %>% filter(Cellularity == "High"), "High Cellularity")

heatmap_low %v% heatmap_mod %v% heatmap_high

###benchmarking results, stratified by pam50####
methods <- c("abis", "instaprism", "cibersort", "epic", "mcp", "quantiseq", "scaden", "timer", "xcell", "consensus")
results_pam <- data.frame(Method = character(), CellType = character(), cellularity = character(), Pearson_R = numeric(), Spearman_R = numeric(), rmse = numeric(), adjpR2 = numeric(), adjsR2 = numeric(), theilU = numeric(), meanae = numeric(), medianae = numeric(), bias = numeric(), kld = numeric())
pam50_classes <- c("Normal", "LumA", "LumB", "Her2", "Basal")

for (pam_class in pam50_classes) {
  # Filter the dataset for the current PAM50 classification
  allimmune_pam <- allimmune %>%
    filter(Pam50Subtype == pam_class)
  
  for (method in methods) {
    print(paste("Method:", method, "PAM50 Class:", pam_class))
    imc_cells <- grep("^imc", names(allimmune_pam), value = TRUE)
    method_cells <- grep(method, names(allimmune_pam), value = TRUE)
    shared_cells <- intersect(sub("imc", "", imc_cells), sub(method, "", method_cells)) # Extracting immune cell names shared between 'imc' and 'abis'
    print(shared_cells)
    
    for (cell_type in shared_cells) {
      # Extract columns for 'imc' and current method and filter only rows with non-missing values
      imc_column <- paste0("imc", cell_type)
      method_column <- paste0(method, cell_type)
      imc_data <- allimmune_pam[, imc_column]
      method_data <- allimmune_pam[, method_column]
      
      # Calculate Pearson correlation using weightedCorr
      Pcorr <- weightedCorr(imc_data, method_data, method = "Pearson")
      Scorr <- weightedCorr(imc_data, method_data, method = "Spearman")
      wPcorr <- weightedCorr(imc_data, method_data, method = "Pearson", weights = allimmune_pam$total_objects)
      wScorr <- weightedCorr(imc_data, method_data, method = "Spearman", weights = allimmune_pam$total_objects)
      rmse <- rmse(imc_data, method_data)
      data <- data.frame(imc_data, method_data)
      model <- lm(imc_data ~ method_data, data = data)
      pR2 <- R2(model, type = "pearson")
      adjpR2 <- pR2["R.squared.adj"]
      sR2 <- R2(model, type = "spearman")
      adjsR2 <- sR2["R.squared.adj"]
      theilu <- TheilU(imc_data, method_data, type = 2, na.rm = FALSE)
      meanae <- mae(imc_data, method_data)
      medianae <- mdae(imc_data, method_data)
      bias <- bias(imc_data, method_data)
      kld <- KLD(imc_data, method_data)
      
      # Store results in the dataframe
      results_pam <- rbind(results_pam, c(method, cell_type, pam_class, Pcorr, Scorr, wPcorr, wScorr, rmse, adjpR2, adjsR2, theilu, meanae, medianae, bias, kld[["mean.sum.KLD"]]))
    }
  }
}

colnames(results_pam) <- c("Method", "Immune_Cell", "Pam50", "Unweighted Pearson", "Unweighted Spearman", "Weighted Pearson", "Weighted Spearman", "RMSE", "adjusted pearson R2", "adjusted spearman R2", "Theil's U Index", "Mean Absolute Error", "Median Absolute Error", "Bias", "KLD")

###Kass results
imc_kassandra <- allimmune %>%
  merge(kassandra, by = "metabric_id") #kass only has 440 samples, put back in

print(method)
imc_cells <- grep("^imc", names(imc_kassandra), value = TRUE)
method_cells <- grep("^kassandra", names(imc_kassandra), value = TRUE)
shared_cells <- intersect(sub("imc", "", imc_cells), sub("kassandra", "", method_cells)) # Extracting immune cell names shared between 'imc' and 'abis'
print(shared_cells)
pam50_classes <- c("Normal", "LumA", "LumB", "Her2", "Basal")

#loops to iterate over kassandra and immune cell types in High cellularity
for (pam_class in pam50_classes) {
  # Filter the dataset for the current PAM50 classification
  imc_kassandra_pam <- imc_kassandra %>%
    filter(Pam50Subtype == pam_class)
  
  for (cell_type in shared_cells) {
    # Extract columns for 'imc' and current method and filter only rows with non-missing values
    imc_column <- paste0("imc", cell_type)
    method_column <- paste0("kassandra", cell_type)
    imc_data <- imc_kassandra_pam[, imc_column]
    method_data <- imc_kassandra_pam[, method_column]
    
    # Calculate Pearson correlation using weightedCorr
    Pcorr <- weightedCorr(imc_data, method_data, method = "Pearson")
    Scorr <- weightedCorr(imc_data, method_data, method = "Spearman")
    wPcorr <- weightedCorr(imc_data, method_data, method = "Pearson", weights = imc_kassandra_pam$total_objects)
    wScorr <- weightedCorr(imc_data, method_data, method = "Spearman", weights = imc_kassandra_pam$total_objects)
    rmse <- rmse(imc_data, method_data)
    data <- data.frame(imc_data, method_data)
    model <- lm(imc_data ~ method_data, data = data)
    pR2 <- R2(model, type = "pearson")
    adjpR2 <- pR2["R.squared.adj"]
    sR2 <- R2(model, type = "spearman")
    adjsR2 <- sR2["R.squared.adj"]
    theilu <- TheilU(imc_data, method_data, type = 2, na.rm = FALSE)
    meanae <- mae(imc_data, method_data)
    medianae <- mdae(imc_data, method_data)
    bias <- bias(imc_data, method_data)
    kld <- KLD(imc_data, method_data)
    
    # Store results in the dataframe
    method <- "kassandra"
    results_pam <- rbind(results_pam, c(method, cell_type, pam_class, Pcorr, Scorr, wPcorr, wScorr, rmse, adjpR2, adjsR2, theilu, meanae, medianae, bias, kld[["mean.sum.KLD"]]))
  }
}

#plotting MAE (stratified by pam50)####
#rename methods and cells
results_pam$Method <- case_when(
  results_pam$Method == "abis" ~ "Abis",
  results_pam$Method == "instaprism" ~ "InstaPrism",
  results_pam$Method == "cibersort" ~ "Cibersort",
  results_pam$Method == "epic" ~ "Epic",
  results_pam$Method == "mcp" ~ "MCPCounter",
  results_pam$Method == "quantiseq" ~ "Quantiseq",
  results_pam$Method == "scaden" ~ "Scaden",
  results_pam$Method == "timer" ~ "TIMER",
  results_pam$Method == "xcell" ~ "xCell",
  results_pam$Method == "consensus" ~ "ConsensusTME",
  results_pam$Method == "kassandra" ~ "Kassandra",
  TRUE ~ results_pam$Method
)
results_pam$Immune_Cell <- case_when(
  results_pam$Immune_Cell == "_B_cells" ~ "B Cells",
  results_pam$Immune_Cell == "_Granulocytes" ~ "Granulocytes",
  results_pam$Immune_Cell == "_NK_Cells" ~ "NK Cells",
  results_pam$Immune_Cell == "_CD8_T_Cells" ~ "CD8 T Cells",
  results_pam$Immune_Cell == "_Endothelial" ~ "Endothelial",
  results_pam$Immune_Cell == "_Fibroblasts" ~ "Fibroblasts",
  results_pam$Immune_Cell == "_Macrophages" ~ "Macrophages",
  results_pam$Immune_Cell == "_CD4_T_Cells" ~ "CD4 T Cells",
  results_pam$Immune_Cell == "_Tregs" ~ "Tregs",
  TRUE ~ results_pam$Immune_Cell
)

method_order <- c("Abis", "TIMER", "Scaden", "Kassandra", "ConsensusTME", "MCPCounter", "Epic", "xCell", "Cibersort", "Quantiseq", "InstaPrism")
results_pam$Method <- factor(results_pam$Method, levels = method_order)

# Function to create heatmap
create_heatmap <- function(data, title) {
  mat <- data %>%
    select(Method, Immune_Cell, "Median Absolute Error") %>%
    pivot_wider(names_from = Immune_Cell, values_from = "Median Absolute Error")
  
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$Method
  mat <- mat[,-1]
  mat <- as.matrix(mat)
  
  mat <- t(mat)
  
  # Sort rows and columns by average
  mat <- apply(mat, 2, as.numeric)
  row_order <- order(rowMedians(mat, na.rm=TRUE), decreasing = TRUE)
  #col_order <- order(colMedians(mat, na.rm=TRUE), decreasing = TRUE)
  
  # Create a new heatmap with transposed and sorted data
  col_fun = colorRamp2(c(0, 0.65), c("#ffffff", "#DB162F"))
  
  mat_new <- mat[row_order, method_order]
  
  Heatmap(
    mat_new,
    col = col_fun,
    na_col = "grey",
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.3f", mat_new[i, j]), x, y, gp = gpar(fontsize = 10))},
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title = title,
    rect_gp = gpar(col = "white", lwd = 2),
    name = "Median Absolute Error"
  )
}

# Create heatmaps for each cellularity level
heatmap_luma <- create_heatmap(results_pam %>% filter(Pam50 == "LumA"), "LumA")
heatmap_lumb <- create_heatmap(results_pam %>% filter(Pam50 == "LumB"), "LumB")
heatmap_basal <- create_heatmap(results_pam %>% filter(Pam50 == "Basal"), "Basal")
heatmap_her <- create_heatmap(results_pam %>% filter(Pam50 == "Her2"), "Her2")
heatmap_normal <- create_heatmap(results_pam %>% filter(Pam50 == "Normal"), "Normal")

heatmap_luma %v% heatmap_lumb %v% heatmap_basal %v% heatmap_her %v% heatmap_normal

#plotting RMSE (stratified by pam50)####
#rename methods and cells
results_pam$Method <- case_when(
  results_pam$Method == "abis" ~ "Abis",
  results_pam$Method == "instaprism" ~ "InstaPrism",
  results_pam$Method == "cibersort" ~ "Cibersort",
  results_pam$Method == "epic" ~ "Epic",
  results_pam$Method == "mcp" ~ "MCPCounter",
  results_pam$Method == "quantiseq" ~ "Quantiseq",
  results_pam$Method == "scaden" ~ "Scaden",
  results_pam$Method == "timer" ~ "TIMER",
  results_pam$Method == "xcell" ~ "xCell",
  results_pam$Method == "consensus" ~ "ConsensusTME",
  results_pam$Method == "kassandra" ~ "Kassandra",
  TRUE ~ results_pam$Method
)
results_pam$Immune_Cell <- case_when(
  results_pam$Immune_Cell == "_B_cells" ~ "B Cells",
  results_pam$Immune_Cell == "_Granulocytes" ~ "Granulocytes",
  results_pam$Immune_Cell == "_NK_Cells" ~ "NK Cells",
  results_pam$Immune_Cell == "_CD8_T_Cells" ~ "CD8 T Cells",
  results_pam$Immune_Cell == "_Endothelial" ~ "Endothelial",
  results_pam$Immune_Cell == "_Fibroblasts" ~ "Fibroblasts",
  results_pam$Immune_Cell == "_Macrophages" ~ "Macrophages",
  results_pam$Immune_Cell == "_CD4_T_Cells" ~ "CD4 T Cells",
  results_pam$Immune_Cell == "_Tregs" ~ "Tregs",
  TRUE ~ results_pam$Immune_Cell
)

method_order <- c("Abis", "TIMER", "Kassandra", "Scaden", "Epic", "ConsensusTME", "Cibersort", "MCPCounter", "InstaPrism", "xCell", "Quantiseq")
results_pam$Method <- factor(results_pam$Method, levels = method_order)

# Function to create heatmap
create_heatmap <- function(data, title) {
  mat <- data %>%
    select(Method, Immune_Cell, "RMSE") %>%
    pivot_wider(names_from = Immune_Cell, values_from = "RMSE")
  
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$Method
  mat <- mat[,-1]
  mat <- as.matrix(mat)
  
  mat <- t(mat)
  
  # Sort rows and columns by average
  mat <- apply(mat, 2, as.numeric)
  row_order <- order(rowMedians(mat, na.rm=TRUE), decreasing = TRUE)
  #col_order <- order(colMedians(mat, na.rm=TRUE), decreasing = TRUE)
  
  # Create a new heatmap with transposed and sorted data
  col_fun = colorRamp2(c(0, 0.65), c("#ffffff", "#DB162F"))
  
  mat_new <- mat[row_order, method_order]
  
  Heatmap(
    mat_new,
    col = col_fun,
    na_col = "grey",
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.3f", mat_new[i, j]), x, y, gp = gpar(fontsize = 10))},
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title = title,
    rect_gp = gpar(col = "white", lwd = 2),
    name = "RMSE"
  )
}

# Create heatmaps for each cellularity level
heatmap_luma <- create_heatmap(results_pam %>% filter(Pam50 == "LumA"), "LumA")
heatmap_lumb <- create_heatmap(results_pam %>% filter(Pam50 == "LumB"), "LumB")
heatmap_basal <- create_heatmap(results_pam %>% filter(Pam50 == "Basal"), "Basal")
heatmap_her <- create_heatmap(results_pam %>% filter(Pam50 == "Her2"), "Her2")
heatmap_normal <- create_heatmap(results_pam %>% filter(Pam50 == "Normal"), "Normal")

heatmap_luma %v% heatmap_lumb %v% heatmap_basal %v% heatmap_her %v% heatmap_normal



#plotting MAE boxplots (stratified by cellularity)####
#rename methods and cells
results_cellularity$Method <- case_when(
  results_cellularity$Method == "abis" ~ "Abis",
  results_cellularity$Method == "instaprism" ~ "InstaPrism",
  results_cellularity$Method == "cibersort" ~ "Cibersort",
  results_cellularity$Method == "epic" ~ "Epic",
  results_cellularity$Method == "mcp" ~ "MCPCounter",
  results_cellularity$Method == "quantiseq" ~ "Quantiseq",
  results_cellularity$Method == "scaden" ~ "Scaden",
  results_cellularity$Method == "timer" ~ "TIMER",
  results_cellularity$Method == "xcell" ~ "xCell",
  results_cellularity$Method == "consensus" ~ "ConsensusTME",
  results_cellularity$Method == "kassandra" ~ "Kassandra",
  TRUE ~ results_cellularity$Method
)
results_cellularity$Immune_Cell <- case_when(
  results_cellularity$Immune_Cell == "_B_cells" ~ "B Cells",
  results_cellularity$Immune_Cell == "_Granulocytes" ~ "Granulocytes",
  results_cellularity$Immune_Cell == "_NK_Cells" ~ "NK Cells",
  results_cellularity$Immune_Cell == "_CD8_T_Cells" ~ "CD8 T Cells",
  results_cellularity$Immune_Cell == "_Endothelial" ~ "Endothelial",
  results_cellularity$Immune_Cell == "_Fibroblasts" ~ "Fibroblasts",
  results_cellularity$Immune_Cell == "_Macrophages" ~ "Macrophages",
  results_cellularity$Immune_Cell == "_CD4_T_Cells" ~ "CD4 T Cells",
  results_cellularity$Immune_Cell == "_Tregs" ~ "Tregs",
  TRUE ~ results_cellularity$Immune_Cell
)

results_cellularity$Method <- factor(results_cellularity$Method, levels = method_order)
results_cellularity <- results_cellularity %>%
  group_by(Cellularity) %>%
  mutate(Method = reorder(Method, `Median Absolute Error`, median))

results_cellularity$Cellularity <- factor(results_cellularity$Cellularity, 
                                          levels = c("Low", "Moderate", "High"))

ggplot(results_cellularity, aes(x = Method, y = `Median Absolute Error`, fill = Cellularity)) +
  geom_boxplot() + 
  geom_point(color = "black", alpha = 0.5, position = position_dodge(width = 0.75)) + 
  theme_bw() +
  facet_wrap(~Cellularity, scales = "free_y") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Increase size of x-axis text
    axis.text.y = element_text(size = 12), # Increase size of y-axis text
    axis.title = element_text(size = 14), # Increase size of axis titles
    plot.title = element_text(size = 16, face = "bold"), # Increase size of the plot title
    strip.text = element_text(size = 14) # Increase size of facet labels
  ) +
  guides(fill = "none") + # Remove the legend for the fill aesthetic
  labs(title = "Median Absolute Error by Method and Cellularity", 
       x = "Method", 
       y = "Median Absolute Error") +
  scale_fill_manual(values = c("Low" = "#A0C4FF", "Moderate" = "#CAFFBF", "High" = "#FFADAD")) + # Set custom colors for Cellularity
  guides(fill = "none") + # Remove the legend for the fill aesthetic
  labs(title = NULL, 
       x = "Method", 
       y = "Median Absolute Error")

#A0C4FF
#CAFFBF
#FFADAD

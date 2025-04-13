#Libraries####
library(data.table)
library(dplyr)
library(tidyr)
library(writexl)
library(readxl)
library(ggplot2)
library(GGally)
library(plyr)
library(hrbrthemes)

###data prep####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
imc <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 1)
instaprism <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 2)
abis <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 3)
consensus <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 4)
epic <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 5)
estimate <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 6)
kassandra <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 7)
mcp <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 8)
quantiseq <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 9)
scaden <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 10)
timer <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 11)
xcell <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 12)
digpath <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 13)
methylayer <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 14)
ascat <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 15)
bcps <- read_xlsx("Tissue Level Characterization.xlsx", sheet = 16)


setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
clinical <- fread('metabric_clinical.txt')
clinical <- clinical %>%
  rename_with(~"metabric_id", .cols = "METABRIC.ID") %>%
  select(metabric_id, Cellularity, iC10)

allimmune <- merge(imc, instaprism, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, abis, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, consensus, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, epic, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, estimate, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, kassandra, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, mcp, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, quantiseq, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, scaden, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, timer, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, xcell, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, digpath, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, methylayer, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, ascat, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, bcps, by = "metabric_id", all = TRUE)
allimmune <- merge(allimmune, clinical, by = "metabric_id", all = TRUE)

na_indices <- which(is.na(allimmune$Cellularity))
allimmune$Cellularity[na_indices] <- ifelse(allimmune$IMC_epithelial[na_indices] < 0.33, "low",
                                       ifelse(allimmune$IMC_epithelial[na_indices] >= 0.33 & allimmune$IMC_epithelial[na_indices] <= 0.67, "moderate",
                                              ifelse(allimmune$IMC_epithelial[na_indices] > 0.67, "high", NA)))

#Parallel coordinate plots####
immune_cols <- grep("_immune", names(allimmune), value = TRUE)
immune_data <- allimmune[, c(immune_cols, "Cellularity")]
ggparcoord(immune_data, columns = 1:13, groupColumn = 14, scale = "center", scaleSummary = "mean",  missing = "median", boxplot = TRUE, title = "Immune Parallel Coordinates Plot")

epi_cols <- grep("_epithelial", names(allimmune), value = TRUE)
epi_data <- allimmune[, c(epi_cols, "Cellularity")]
ggparcoord(epi_data, columns = 1:6, groupColumn =7, scale = "center", scaleSummary = "mean",  missing = "median", boxplot = TRUE, title = "Epithelial Parallel Coordinates Plot")

stroma_cols <- grep("_stroma", names(allimmune), value = TRUE)
stroma_data <- allimmune[, c(stroma_cols, "Cellularity")]
ggparcoord(stroma_data, columns = 1:10, groupColumn =11, scale = "center", scaleSummary = "mean",  missing = "median", boxplot = TRUE, title = "Stroma Parallel Coordinates Plot")

#tissue-level landscape plot####
tissue_data <- melt(allimmune)
tissue_data <- select(tissue_data, metabric_id, variable, value, Cellularity)
tissue_data <- separate(tissue_data, variable, into = c("method", "tissue"), sep = "_")
gg <- ggplot(tissue_data, aes(x = Cellularity, y = value, fill = tissue)) +
  geom_violin() +
  facet_wrap(~method, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_ft() +
  scale_fill_ft() +
  theme_ipsum(grid="XY", axis="xy")+
  labs(title = "METABRIC Deconvolution", x = "Cellularity", y = "Fraction/Score")


flush_ticks(gg)
                      
#GGpairs plots####
#immune plots
immune_cols <- grep("_immune", names(allimmune), value = TRUE)
immune_data <- allimmune[, c(immune_cols, "Cellularity")]
immune_data$instaprism_immune[immune_data$instaprism_immune == 0] <- NA
immune_data$scaden_immune[immune_data$scaden_immune == 0] <- NA
immune_data <- immune_data[!is.na(immune_data$IMC_immune), ]
setDT(immune_data)

#cols_to_transform <- names(immune_data)[1:13] # Column names to be transformed
#for (col in cols_to_transform) {
#  immune_data[, (col) := ifelse(is.na(get(col)), median(get(col), na.rm = TRUE), get(col)), by = Cellularity]
#} # Loop through each column and swap NAs with the median of the columnm grouped by Cellularity

ggpairs(immune_data, columns = c(1:6,8:14),
        ggplot2::aes(colour = Cellularity),
        title = paste("Immune"), 
        axisLabels = "show", 
        columnLabels = c("IMC", "InstaPrism", "Abis", "Consensus", "EPIC", "ESTIMATE", "MCPCounter", "Quantiseq", "Scaden", "TIMER", "XCell", "Digital Pathology", "Methylayer"), 
        theme = ggplot2::theme_minimal() +
          ggplot2::theme(text = ggplot2::element_text(size = 14)))


#to compare variation in data modality
ggpairs(immune_data, columns = c(1,8,13,14),
        ggplot2::aes(colour = Cellularity),
        title = paste("Immune"), 
        axisLabels = "show", 
        columnLabels = c("IMC", "MCPCounter", "Digital Pathology", "Methylayer"))

#epithelial plots
epi_cols <- grep("_epithelial", names(allimmune), value = TRUE)
epi_data <- allimmune[, c(epi_cols, "Cellularity")]
epi_data$instaprism_epithelial[epi_data$instaprism_epithelial == 0] <- NA
epi_data <- epi_data[!is.na(epi_data$IMC_epithelial), ]
setDT(epi_data)

#cols_to_transform <- names(epi_data)[1:8] # Column names to be transformed
#for (col in cols_to_transform) {
#  epi_data[, (col) := ifelse(is.na(get(col)), median(get(col), na.rm = TRUE), get(col)), by = Cellularity]
#} # Loop through each column and swap NAs with the median of the columnm grouped by Cellularity

ggpairs(epi_data, columns = c(1:4,6:8),
        ggplot2::aes(colour = Cellularity),
        title = paste("Epithelial"), 
        axisLabels = "show", 
        columnLabels = c("IMC", "InstaPrism", "EPIC", "ESTIMATE", "Digital Pathology", "ASCAT", "BCPS"))

#to compare variation in data modality
ggpairs(epi_data, columns = c(1,2,6,7),
        ggplot2::aes(colour = Cellularity),
        title = paste("Epithelial"), 
        axisLabels = "show", 
        columnLabels = c("IMC", "InstaPrism",  "Digital Pathology", "ASCAT"))


#stroma plots
stroma_cols <- grep("_stroma", names(allimmune), value = TRUE)
stroma_data <- allimmune[, c(stroma_cols, "Cellularity")]
stroma_data$instaprism_stroma[stroma_data$instaprism_stroma == 0] <- NA
stroma_data <- stroma_data[!is.na(stroma_data$IMC_stroma), ]
setDT(stroma_data)

#cols_to_transform <- names(stroma_data)[1:10] # Column names to be transformed
#for (col in cols_to_transform) {
#  stroma_data[, (col) := ifelse(is.na(get(col)), median(get(col), na.rm = TRUE), get(col)), by = Cellularity]
#} # Loop through each column and swap NAs with the median of the columnm grouped by Cellularity

ggpairs(stroma_data, columns = c(1:5,7:10),
        ggplot2::aes(colour = Cellularity),
        title = paste("Stroma"), 
        axisLabels = "show", 
        columnLabels = c("IMC", "InstaPrism", "Consensus", "EPIC", "ESTIMATE", "MCPCounter", "XCell", "Digital Pathology", "Methylayer"))

#to compare variation in data modality
ggpairs(stroma_data, columns = c(1,2,9,10),
        ggplot2::aes(colour = Cellularity),
        title = paste("Stroma"), 
        axisLabels = "show", 
        columnLabels = c("IMC", "InstaPrism", "Digital Pathology", "Methylayer"))

#data moadlity plots####
immune_cols <- grep("_immune", names(allimmune), value = TRUE)
immune_data <- allimmune[, c(immune_cols, "Cellularity")]
immune_data$instaprism_immune[immune_data$instaprism_immune == 0] <- NA
immune_data$scaden_immune[immune_data$scaden_immune == 0] <- NA
immune_data <- immune_data[!is.na(immune_data$IMC_immune), ]
setDT(immune_data)


#Epithelial tisue by cellularity####
epi_cols <- grep("_epithelial", names(allimmune), value = TRUE)
epi_data <- allimmune[, c(epi_cols, "Cellularity")]
epi_data$instaprism_epithelial[epi_data$instaprism_epithelial == 0] <- NA
epi_data <- epi_data[!is.na(epi_data$IMC_epithelial), ]
setDT(epi_data)

epi_data_long <- melt(epi_data, id.vars = "Cellularity", variable.name = "Method", value.name = "Value")
epi_data_long$Cellularity <- factor(epi_data_long$Cellularity, levels = c("low", "moderate", "high"))
epi_data_long$Method <- gsub("_epithelial", "", epi_data_long$Method)

epi_data_long$Method <- gsub("instaprism", "Instaprism", epi_data_long$Method)
epi_data_long$Method <- gsub("epic", "EPIC", epi_data_long$Method)
epi_data_long$Method <- gsub("estimate", "Estimate", epi_data_long$Method)
epi_data_long$Method <- gsub("kassandra", "Kassandra", epi_data_long$Method)
epi_data_long$Method <- gsub("ascat", "ASCAT", epi_data_long$Method)
epi_data_long$Method <- gsub("bcps", "BCPS", epi_data_long$Method)
epi_data_long$Method <- gsub("digpath", "Digital Pathology", epi_data_long$Method)

method_order <- c("IMC", "Instaprism", "EPIC", "Estimate", "Kassandra", "ASCAT", "BCPS", "Digital Pathology")
epi_data_long$Method <- factor(epi_data_long$Method, levels = method_order)

ggplot(epi_data_long, aes(x = Method, y = Value, fill = Cellularity)) + 
  geom_boxplot() + 
  theme_minimal() + 
  facet_wrap(~Method, scales = "free", ncol = 9) +
  theme(axis.text.x = element_text(size = 12),
        strip.text = element_blank()) + 
  labs(title = "Epithelial Tissue by Cellularity", x = "", y = "Tumor Purity Score/Fraction")+
  scale_fill_manual(values = c("low" = "#A0C4FF", "moderate" = "#CAFFBF", "high" = "#FFADAD"))

#Epithelial tisue by cellularity, force range between 0-1####
epi_cols <- grep("_epithelial", names(allimmune), value = TRUE)
epi_data <- allimmune[, c(epi_cols, "Cellularity")]
epi_data$instaprism_epithelial[epi_data$instaprism_epithelial == 0] <- NA
epi_data <- epi_data[!is.na(epi_data$IMC_epithelial), ]
setDT(epi_data)

cols_to_scale <- grep("_epithelial$", names(epi_data), value = TRUE)

# Apply scaling to force a range of 0-1, 
epi_data[, (cols_to_scale) := lapply(.SD, function(x) {
  if (is.numeric(x)) {
    scaled_x <- scale(x, center = min(x, na.rm = TRUE), scale = max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    return(as.numeric(scaled_x))  # Ensure the result is numeric
  }
  return(x)  # Leave non-numeric columns unchanged
}), .SDcols = cols_to_scale]


epi_data_long <- melt(epi_data, id.vars = "Cellularity", variable.name = "Method", value.name = "Value")
epi_data_long$Cellularity <- factor(epi_data_long$Cellularity, levels = c("low", "moderate", "high"))
epi_data_long$Method <- gsub("_epithelial", "", epi_data_long$Method)

epi_data_long$Method <- gsub("instaprism", "Instaprism", epi_data_long$Method)
epi_data_long$Method <- gsub("epic", "EPIC", epi_data_long$Method)
epi_data_long$Method <- gsub("estimate", "Estimate", epi_data_long$Method)
epi_data_long$Method <- gsub("kassandra", "Kassandra", epi_data_long$Method)
epi_data_long$Method <- gsub("ascat", "ASCAT", epi_data_long$Method)
epi_data_long$Method <- gsub("bcps", "BCPS", epi_data_long$Method)
epi_data_long$Method <- gsub("digpath", "Digital Pathology", epi_data_long$Method)

method_order <- c("IMC", "Instaprism", "EPIC", "Estimate", "Kassandra", "ASCAT", "BCPS", "Digital Pathology")
epi_data_long$Method <- factor(epi_data_long$Method, levels = method_order)

ggplot(epi_data_long, aes(x = Method, y = Value, fill = Cellularity)) + 
  geom_boxplot() + 
  theme_minimal() + 
  facet_wrap(~Method, scales = "free_x", ncol = 9) +
  theme(axis.text.x = element_text(size = 12),
        strip.text = element_blank()) + 
  labs(title = "Epithelial Tissue by Cellularity", x = "", y = "Tumor Purity Score/Fraction")+
  scale_fill_manual(values = c("low" = "#A0C4FF", "moderate" = "#CAFFBF", "high" = "#FFADAD"))

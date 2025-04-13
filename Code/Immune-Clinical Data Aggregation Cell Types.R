#libraries####
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(mice)

#HatzisD cleaning####
#import data
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
HatzisD <- fread("Hatzis-Discovery_clinical.tsv")
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
HatzisD_immune <- read_xlsx("hatzis-discovery InstaPrism results tpm.xlsx", sheet = 4)

#rename and select variables
HatzisD <- HatzisD %>%
  mutate(drfs_even_time_days = round(drfs_even_time_years * 365))
HatzisD <- HatzisD %>%
  dplyr::select(Sample_ID, age_years, clinical_nodal_status, clinical_t_stage, 
         er_status_ihc, her2_status, grade, pam50_class, 
         pathologic_response_pcr_rd, pathologic_response_rcb_class, 
         source, drfs_even_time_days, drfs_1_event_0_censored, intclust) %>%
  dplyr::rename(patient_id = Sample_ID, age = age_years, node_status = clinical_nodal_status, 
         size = clinical_t_stage, ER_status = er_status_ihc, HER2_status = her2_status, 
         grade = grade, pam50 = pam50_class, pcr = pathologic_response_pcr_rd, 
         rcb = pathologic_response_rcb_class, source = source, 
         MFS_time = drfs_even_time_days, MFS_event = drfs_1_event_0_censored)%>%
  mutate(treatment = "taxane-anthracycline")%>%
  mutate(treatment_catagory = "chemotherapy")%>%
  mutate(study = "HatzisDiscovery")
HatzisD <- HatzisD %>%
  mutate(
    sample_id = NA,
    OS_event = NA,
    OS_time = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    primary = "true",
    metastases_location = NA)
HatzisD <- HatzisD %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )
HatzisD_merged <- merge(HatzisD, HatzisD_immune, by = "patient_id", all.x = TRUE)

#HatzisV cleaning####
#import data
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
HatzisV <- fread("Hatzis-Validation_clinical.tsv")
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
HatzisV_immune <- read_xlsx("hatzis-validation InstaPrism results tpm.xlsx", sheet = 4)

#rename and select variables
HatzisV <- HatzisV %>%
  mutate(drfs_even_time_days = round(drfs_even_time_years * 365))
HatzisV <- HatzisV %>%
  dplyr::select(Sample_ID, age_years, clinical_nodal_status, clinical_t_stage, 
         er_status_ihc, her2_status, grade, pam50_class, 
         pathologic_response_pcr_rd, pathologic_response_rcb_class, source, 
         drfs_even_time_days, drfs_1_event_0_censored, intclust) %>%
  dplyr::rename(patient_id = Sample_ID, age = age_years, node_status = clinical_nodal_status, 
         size = clinical_t_stage, ER_status = er_status_ihc, HER2_status = her2_status, 
         grade = grade, pam50 = pam50_class, pcr = pathologic_response_pcr_rd, 
         rcb = pathologic_response_rcb_class, source = source, 
         MFS_time = drfs_even_time_days, MFS_event = drfs_1_event_0_censored)%>%
  mutate(treatment = "taxane-anthracycline")%>%
  mutate(treatment_catagory = "chemotherapy")%>%
  mutate(study = "HatzisValidation")
HatzisV <- HatzisV %>%
  mutate(
    sample_id = NA,
    OS_event = NA,
    OS_time = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    primary = "true",
    metastases_location = NA)
HatzisV <- HatzisV %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )
HatzisV_merged <- merge(HatzisV, HatzisV_immune, by = "patient_id", all.x = TRUE)

#metabric cleaning####
#import data
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
metabric <- fread("metabric_clinical.txt")
metabric2 <- fread("metabric_clinical_cbioportal.tsv")
metabric2 <- dplyr::rename(metabric2, METABRIC.ID = "Patient ID")
metabric <- merge(metabric, metabric2, by = "METABRIC.ID", all.x = TRUE)
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
metabric_immune <- read_xlsx("metabric InstaPrism results tpm.xlsx", sheet = 4)

#rename and select variables
metabric <- metabric %>%
  mutate(
    node_status = case_when(
      Lymph.Nodes.Positive %in% 0 ~ "N0",
      Lymph.Nodes.Positive %in% 1:3 ~ "N1",
      Lymph.Nodes.Positive %in% 4:9 ~ "N2",
      Lymph.Nodes.Positive >= 10 ~ "N3",
      TRUE ~ NA_character_  # Include this line to handle other cases, if any
    )
  )
metabric <- metabric %>%
  mutate(RFS_event = 1 * (DeathBreast == 1 | LR == 1 | DR == 1)) %>%
  mutate(RFS_time = apply(., 1, function(row) min(row[c("TLR", "TDR", "T")], na.rm = TRUE))) %>%
  dplyr::select(METABRIC.ID, Age.At.Diagnosis, ER.Status, `HER2 Status`, iC10, Grade, Size, node_status, Pam50Subtype, RFS_time, RFS_event, T, DeathBreast, `Tumor Stage`) %>% 
  mutate(study = "METABRIC")%>%
  dplyr::rename(patient_id = METABRIC.ID, age = Age.At.Diagnosis, ER_status = ER.Status, HER2_status = `HER2 Status`, 
         grade = Grade, pam50 = Pam50Subtype, size = `Tumor Stage`, BCSS_time = T, BCSS_event = DeathBreast, intclust = iC10)%>%
  mutate(
    sample_id = NA,
    OS_event = NA,
    OS_time = NA,
    MFS_event = NA,
    MFS_time = NA,
    primary = "true",
    treatment_catagory = NA,
    treatment = NA,
    pcr = NA,
    rcb = NA,
    metastases_location = NA)
metabric <- metabric %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )
metabric_merged <- merge(metabric, metabric_immune, by = "patient_id", all.x = TRUE)

#SCANB cleaning####
#import data
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
scanb <- read_xlsx("SCANB Clinical.xlsx")
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
scanb_immune <- read_xlsx("SCANB InstaPrism results.xlsx", sheet = 4)
names(scanb_immune)[1] <- "sample_id"
names(scanb)[44] <- "age"

#rename and select variables
scanb <- scanb %>%
  dplyr::select(Patient, GEX.assay, age, ER_add_mclust, HER2_add_mclust, LN.spec, Size.mm, NHG, pT, DRFi_days, DRFi_event, OS_days, OS_event, RFi_days, RFi_event, BCFi_days, BCFi_event, NCN.PAM50, intclust) %>%
  dplyr::rename(patient_id = Patient, sample_id = GEX.assay, ER_status = ER_add_mclust, HER2_status = HER2_add_mclust, grade = NHG, node_status = LN.spec, size = pT, MFS_time = DRFi_days, pam50 = NCN.PAM50,
         MFS_event = DRFi_event, OS_time = OS_days, OS_event = OS_event, RFS_time = RFi_days, RFS_event = RFi_event, BCSS_time = BCFi_days, BCSS_event = BCFi_event)%>%
  mutate(study = "SCANB")
scanb <- scanb %>%
  mutate(
    primary = "true",
    treatment_catagory = NA,
    treatment = NA,
    pcr = NA,
    rcb = NA,
    metastases_location = NA)
scanb <- scanb %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )
scanb_merged <- merge(scanb, scanb_immune, by = "sample_id", all.x = TRUE)

#SMC cleaning####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
smc_patient <- fread("smc_clinical_patient.txt")
smc_sample <- fread("smc_clinical_sample.txt")
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
smc_immune <- read_xlsx("SMC InstaPrism results.xlsx", sheet = 4)

smc <- merge(smc_patient, smc_sample, by = "PATIENT_ID", all.x = TRUE)
smc <- smc %>%
  mutate(ER_status = ifelse(SUBTYPE_CONSENSUS == "TN", "Negative", "Positive")) %>%
  mutate(HER2_status = ifelse(SUBTYPE_CONSENSUS == "HER2+" | SUBTYPE_CONSENSUS == "ER+HER2+" | SUBTYPE_CONSENSUS == "TN", "Positive", "Negative"))%>%
  dplyr::select(PATIENT_ID, AGE, INITIAL_TNM_STAGE, PAM50_SUBTYPE, ER_status, HER2_status, intclust) %>%
  dplyr::rename(patient_id = PATIENT_ID, age = AGE, pam50 = PAM50_SUBTYPE)%>%
  mutate(study = "SMC")
smc <- smc %>%
  mutate(
    sample_id = NA,
    grade = NA,
    size = NA,
    node_status = NA,
    OS_event = NA,
    OS_time = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    MFS_event = NA,
    MFS_time = NA,
    primary = "true",
    treatment_catagory = NA,
    treatment = NA,
    pcr = NA,
    rcb = NA,
    metastases_location = NA)
smc <- smc %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )
smc_merged <- merge(smc, smc_immune, by = "patient_id", all.x = TRUE)



#TCGA cleaning####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
tcga_patient <- fread("tcga_clinical_patient.txt")
tcga_sample <- fread("tcga_clinical_sample.txt")
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
tcga_immune <- read_xlsx("TCGA InstaPrism results.xlsx", sheet = 4)

tcga <- merge(tcga_patient, tcga_sample, by = "PATIENT_ID", all.x = TRUE)
tcga <- tcga %>%
  mutate(OS_DAYS = round(OS_MONTHS * 30))%>%
  mutate(DSS_DAYS = round(DSS_MONTHS * 30))%>%
  mutate(DFS_DAYS = round(DFS_MONTHS * 30))

tcga <- tcga %>%
  dplyr::select(PATIENT_ID, SAMPLE_ID, AGE, PATH_N_STAGE, PATH_T_STAGE, RACE, OS_STATUS, OS_DAYS, DSS_STATUS, DSS_DAYS, DFS_STATUS, 
         DFS_DAYS, SUBTYPE, intclust, ER_status, HER2_status)%>%
  dplyr::rename(patient_id = PATIENT_ID, sample_id = SAMPLE_ID, age = AGE, pam50 = SUBTYPE, node_status = PATH_N_STAGE, size = PATH_T_STAGE, race = RACE,
         OS_event = OS_STATUS, OS_time = OS_DAYS, BCSS_event = DSS_STATUS, BCSS_time = DSS_DAYS, RFS_event = DFS_STATUS, 
         RFS_time = DFS_DAYS)%>%
  mutate(study = "TCGA")
tcga <- tcga %>%
  mutate(
    grade = NA,
    MFS_event = NA,
    MFS_time = NA,
    intclust = NA,
    primary = "true",
    treatment_catagory = NA,
    treatment = NA,
    pcr = NA,
    rcb = NA,
    metastases_location = NA)
tcga <- tcga %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )
tcga_immune$patient_id <- gsub("-01", "", tcga_immune$patient_id)
tcga_merged <- merge(tcga, tcga_immune, by = "patient_id", all.x = TRUE)

#MBC cleaning####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
mbc_patient <- fread("mbc_clinical_patient.txt")
mbc_sample <- fread("mbc_clinical_sample.txt")
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
mbc_immune <- read_xlsx("MBC InstaPrism results.xlsx", sheet = 4)

mbc <- merge(mbc_patient, mbc_sample, by = "PATIENT_ID", all.x = TRUE)
mbc <- mbc %>%
  dplyr::select(PATIENT_ID, SAMPLE_ID, CALC_MET_SETTING, BX_LOCATION, BX_ER, BX_HER2OVERALL, intclust, pam50)%>%
  dplyr::rename(patient_id = PATIENT_ID, ER_status = BX_ER, HER2_status = BX_HER2OVERALL, sample_id = SAMPLE_ID, metastases_location = BX_LOCATION)%>%
  mutate(study = "MBC") %>%
  mutate(primary = ifelse(metastases_location == "BREAST" & CALC_MET_SETTING == "NO_METASTATIC_DISEASE_PRESENT", "true", "false"))%>%
  filter(!grepl("BLOOD", sample_id))

mbc <- mbc %>%
  mutate(
    grade = NA,
    size = NA,
    node_status = NA,
    age = NA,
    OS_event = NA,
    OS_time = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    MFS_event = NA,
    MFS_time = NA,
    treatment_catagory = NA,
    treatment = NA,
    pcr = NA,
    rcb = NA)
mbc <- mbc %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )

mbc_immune$sample_id <- gsub("\\.", "-", mbc_immune$sample_id)
mbc_merged <- merge(mbc, mbc_immune, by = "sample_id", all.x = TRUE)
column_names <- names(mbc_merged)
new_column_order <- c(column_names[2], column_names[-2])
mbc_merged <- mbc_merged[, ..new_column_order, with = FALSE] 
mbc_merged <- mbc_merged %>%
  filter(!(grepl("BREAST", metastases_location) & primary == "false"))

#ISPY cleaning####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
ispy <- read_xlsx("ispy2_clinical.xlsx")
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
ispy_immune <- read_xlsx("ispy InstaPrism results tpm.xlsx", sheet = 4)

ispy <- ispy %>%
  dplyr::select("Patient Identifier", Arm, pCR, HR, HER2, PAM50.Subtype, intclust) %>%
  dplyr::rename(patient_id = "Patient Identifier", ER_status = HR, HER2_status = HER2, pam50 = PAM50.Subtype, 
         pcr = pCR, treatment = Arm)%>%
  mutate(treatment_catagory = "chemotherapy")%>%
  mutate(study = "ispy2")
ispy <- ispy %>%
  mutate(
    sample_id = NA,
    age = NA,
    grade = NA,
    size = NA,
    node_status = NA,
    OS_event = NA,
    OS_time = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    MFS_event = NA,
    MFS_time = NA,
    rcb = NA,
    primary = "true",
    metastases_location = NA)
ispy <- ispy %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )

ispy_immune$patient_id <- gsub("X", "", ispy_immune$patient_id)
ispy_merged <- merge(ispy, ispy_immune, by = "patient_id")

#Derouane cleaning####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
derouane_immune <- read_xlsx("Derouane InstaPrism results.xlsx", sheet = 4)
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets/Subtypes")
derouane <- read_xlsx("Derouane_subtype.xlsx")

derouane <- derouane %>%
  mutate(
    sample_id = NA,
    study = "Derouane",
    ER_status = "neg",
    HER2_status = "neg",
    age = NA,
    grade = NA,
    size = NA,
    node_status = NA,
    OS_event = NA,
    OS_time = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    MFS_event = NA,
    MFS_time = NA,
    treatment_catagory = NA,
    treatment = NA,
    pcr = NA,
    rcb = NA,
    primary = "true",
    metastases_location = NA)
derouane <- derouane %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )

derouane_merged <- merge(derouane, derouane_immune, by = "patient_id")

#calgb cleaning####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
calgb_immune <- read_xlsx("CALGB InstaPrism results.xlsx", sheet = 4)
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets/Subtypes")
calgb <- read_xlsx("CALGB_subtypes.xlsx")

calgb <- calgb %>%
  mutate(
    sample_id = NA,
    study = "CALGB",
    ER_status = "neg",
    HER2_status = "neg",
    age = NA,
    grade = NA,
    size = NA,
    node_status = NA,
    OS_event = NA,
    OS_time = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    MFS_event = NA,
    MFS_time = NA,
    treatment_catagory = NA,
    treatment = NA,
    pcr = NA,
    rcb = NA,
    primary = "true",
    metastases_location = NA)
calgb <- calgb %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )

calgb_immune$patient_id <- gsub("//.", "-", calgb$patient_id)
calgb_merged <- merge(calgb, calgb_immune, by = "patient_id")

#matador cleaning####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
matador_immune <- read_xlsx("Matador InstaPrism results.xlsx", sheet = 4)
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
matador <- read_xlsx("matador_subtypes.xlsx")

matador <- matador %>%
  mutate(
    sample_id = NA,
    study = "Matador",
    age = NA,
    grade = NA,
    size = NA,
    node_status = NA,
    OS_event = NA,
    OS_time = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    MFS_event = NA,
    MFS_time = NA,
    treatment_catagory = NA,
    treatment = NA,
    pcr = NA,
    rcb = NA,
    primary = "true",
    metastases_location = NA)
matador <- matador %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )

matador_merged <- merge(matador, matador_immune, by = "patient_id")

#AURORA cleaning####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
aurora_patient <- read_excel("aurora_clinical.xlsx", sheet = 3)
aurora_sample <- read_excel("aurora_clinical.xlsx", sheet = 2)
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
aurora_immune <- read_xlsx("Aurora InstaPrism results.xlsx", sheet = 4)

aurora <- merge(aurora_patient, aurora_sample, by = "patient_id", all.x = TRUE)
aurora <- aurora %>%
  dplyr::select(patient_id, "PAM50 Call", intclust, "Short barcode", "Anatomic Site", "Inferred ER from RNAseq", "Inferred HER2 from RNAseq", "Anatomic Site Simplified")%>%
  dplyr::rename(sample_id = "Short barcode", pam50 = "PAM50 Call", ER_status = "Inferred ER from RNAseq", HER2_status = "Inferred HER2 from RNAseq", metastases_location = "Anatomic Site Simplified")%>%
  mutate(study = "aurora") %>%
  mutate(primary = ifelse(metastases_location == "Breast", "true", "false"))
aurora <- aurora %>%
  mutate(
    grade = NA,
    size = NA,
    node_status = NA,
    age = NA,
    OS_event = NA,
    OS_time = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    MFS_event = NA,
    MFS_time = NA,
    treatment_catagory = NA,
    treatment = NA,
    pcr = NA,
    rcb = NA)
aurora <- aurora %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )

aurora_merged <- merge(aurora, aurora_immune, by = "sample_id", all.x = TRUE)
column_names <- names(aurora_merged)
new_column_order <- c(column_names[2], column_names[-2])
aurora_merged <- aurora_merged[, new_column_order]

#TransNEO cleaning####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
transneo <- read_xlsx("TransNeo_clinical.xlsx", sheet = 1)
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
transneo_immune <- read_xlsx("TransNeo InstaPrism results.xlsx", sheet = 4)

transneo <- transneo %>%
  dplyr::select("Donor.ID", Age, "T.stage", "ER.status", "HER2.status", "Grade.pre.NAT", "NAT.regimen", RCB.category, pCR.RD, PAM50, iC10) %>%
  dplyr::rename(patient_id = "Donor.ID", age = Age, grade = "Grade.pre.NAT", size = "T.stage", ER_status = "ER.status", HER2_status = "HER2.status", pam50 = PAM50, intclust = iC10, treatment = NAT.regimen, pcr = pCR.RD, rcb = RCB.category)%>%
  mutate(treatment_catagory = "chemotherapy")%>%
  mutate(study = "transneo")
transneo <- transneo %>%
  mutate(
    sample_id = NA,
    node_status = NA,
    OS_event = NA,
    OS_time = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    MFS_event = NA,
    MFS_time = NA,
    primary = "true",
    metastases_location = NA)
transneo <- transneo %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )

transneo_merged <- merge(transneo, transneo_immune, by = "patient_id")

#Newton cleaning####
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
newton <- read_xlsx("Newton_clinical.xlsx", sheet = 1)
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
newton_immune <- read_xlsx("Newton InstaPrism results.xlsx", sheet = 4)

newton <- newton %>%
  dplyr::select("patient_id", AgeInt, ER_status1, HER2_status, Stage1, Grade1, VitalStatus, ic10, pam50, OS_days) %>%
  dplyr::rename(age = AgeInt, grade = Grade1, size = Stage1, ER_status = ER_status1, intclust = ic10, OS_event = VitalStatus, OS_time = OS_days)%>%
  mutate(study = "newton")
newton <- newton %>%
  mutate(
    sample_id = NA,
    node_status = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    MFS_event = NA,
    MFS_time = NA,
    primary = "true",
    metastases_location = NA,
    treatment_catagory = NA,
    treatment = NA,
    pcr = NA,
    rcb = NA)
newton <- newton %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )

newton_merged <- merge(newton, newton_immune, by = "patient_id")


#combine into dataframe and standardize variables####
combined_standard <- rbind(HatzisD_merged, HatzisV_merged, metabric_merged, scanb_merged, smc_merged, tcga_merged, mbc_merged, ispy_merged, derouane_merged, calgb_merged, matador_merged, aurora_merged, transneo_merged, newton_merged)
combined_standard <- dplyr::rename(combined_standard, stage = size)

#remove duplicate patients
combined_standard <- combined_standard %>%
  filter(primary == "true") %>%
  distinct(patient_id, study, .keep_all = TRUE) %>%
  bind_rows(combined_standard %>% filter(primary == "false"))

#change ER status (pos vs neg)
combined_standard$ER_status <- gsub("\\bP\\b|\\bPositive\\b|\\bPOSITIVE\\b|\\b1\\b|\\bPOS\\b", "pos", combined_standard$ER_status)
combined_standard$ER_status <- gsub("\\bN\\b|\\bNegative\\b|\\b0\\b|\\bNEGATIVE\\b|\\bnegA\\b|\\bNEG\\b", "neg", combined_standard$ER_status)
combined_standard$ER_status <- gsub("\\bI\\b|\\bIndeterminate\\b|\\bNA\\b|\\/\\/\\[Not Evaluated\\/\\/\\]|\\#N\\/A\\#|\\bTESTING_PERFORMED_ON_DIFFERENT_SAMPLE\\b|\\bNOT_FOUND_IN_RECORD\\b", NA, combined_standard$ER_status)
combined_standard$ER_status <- gsub("\\[Not Evaluated\\]|#neg/A|888", NA, combined_standard$ER_status)
unique(combined_standard$ER_status)
#change HER2 status (pos vs neg)
combined_standard$HER2_status <- gsub("\\bP\\b|\\bPositive\\b|\\bPOSITIVE\\b|\\b1\\b|\\bPOS\\b", "pos", combined_standard$HER2_status)
combined_standard$HER2_status <- gsub("\\bN\\b|\\bNegative\\b|\\b0\\b|\\bNEGATIVE\\b|\\bnegA\\b|\\bNEG\\b", "neg", combined_standard$HER2_status)
combined_standard$HER2_status <- gsub("\\bI\\b|\\bIndeterminate\\b|\\bNA\\b|\\/\\/\\[Not Evaluated\\/\\/\\]|\\#N\\/A\\#|\\bTESTING_PERFORMED_ON_DIFFERENT_SAMPLE\\b|\\bNOT_FOUND_IN_RECORD\\b|\\#neg\\/A", NA, combined_standard$HER2_status)
combined_standard$HER2_status <- gsub("CopyNum Not Available|EQUIVOCAL", NA, combined_standard$HER2_status)
unique(combined_standard$HER2_status)

#change grade status (1-3)
combined_standard$grade <- gsub("4\\=Indeterminate|888", NA, combined_standard$grade)
unique(combined_standard$grade)

#change size/stage status (0-4)
combined_standard$stage <- gsub("T1|T1c|T1b|T1a|T1mi|T1C|T1B|T1A", 1, combined_standard$stage)
combined_standard$stage <- gsub("T2|T2B|T2A", 2, combined_standard$stage)
combined_standard$stage <- gsub("T3|T3A", 3, combined_standard$stage)
combined_standard$stage <- gsub("T4|T4B|T4D", 4, combined_standard$stage)
combined_standard$stage <- gsub("0|TX|T0|Tis|888", NA, combined_standard$stage)
unique(combined_standard$stage)

#change pathological node status (0-3)
combined_standard$node_status <- gsub("N0|N0 \\(I\\+\\)|N0 \\(I\\-\\)|N0 \\(MOL\\+\\)", "0", combined_standard$node_status)
combined_standard$node_status <- gsub("N1|1to3|N1A|N1MI|N1B|N1C|SubMicroMet", 1, combined_standard$node_status)
combined_standard$node_status <- gsub("N2|N2A", 2, combined_standard$node_status)
combined_standard$node_status <- gsub("N3|N3C|N3A|N3B", 3, combined_standard$node_status)
combined_standard$node_status <- gsub("4toX|NX|N3A|N3B", NA, combined_standard$node_status)
unique(combined_standard$node_status)

#change age into categories (5-year ranges)
combined_standard$age <- cut(
  combined_standard$age,
  breaks = c(20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100.01, Inf),
  labels = c("21-25", "26-30", "31-35", "36-40", "41-45", "46-50", "51-55", "56-60", "61-65", "66-70", "71-75", "76-80", "81-85", "86-90", "91-95", "96-100", "100+"),
  include.lowest = TRUE,
  right = FALSE
)

#standardize pam50
combined_standard$pam50 <- gsub("LuminalA|BRCA_LumA", "LumA", combined_standard$pam50)
combined_standard$pam50 <- gsub("LuminalB|BRCA_LumB", "LumB", combined_standard$pam50)
combined_standard$pam50 <- gsub("BRCA_Normal", "Normal", combined_standard$pam50)
combined_standard$pam50 <- gsub("BRCA_Basal", "Basal", combined_standard$pam50)
combined_standard$pam50 <- gsub("BRCA_Her2", "Her2", combined_standard$pam50)
combined_standard$pam50 <- gsub("unclassified|NA|Unk", NA, combined_standard$pam50)
unique(combined_standard$pam50)

#standardize intclust
combined_standard <- combined_standard %>%
  mutate(intclust = case_when(
    intclust == 4 & ER_status == "pos" ~ "4ER+",
    intclust == 4 & ER_status == "neg" ~ "4ER-",
    TRUE ~ as.character(intclust)  # Keeping other values of intclust as is
  ))
unique(combined_standard$intclust)

#standardize OS_event (1 = event, 0 = no event)
combined_standard$OS_event <- gsub("1:DECEASED|Deceased", "1", combined_standard$OS_event)
combined_standard$OS_event <- gsub("0:LIVING|Alive", "0", combined_standard$OS_event)
unique(combined_standard$OS_event)

#standardize RFS_event (1 = event, 0 = no event)
combined_standard$RFS_event <- gsub("1:Recurred/Progressed", "1", combined_standard$RFS_event)
combined_standard$RFS_event <- gsub("0:DiseaseFree", "0", combined_standard$RFS_event)
unique(combined_standard$RFS_event)

#standardize BCSS_event (1 = event, 0 = no event)
combined_standard$BCSS_event <- gsub("1:DEAD WITH TUMOR", "1", combined_standard$BCSS_event)
combined_standard$BCSS_event <- gsub("0:ALIVE OR DEAD TUMOR FREE", "0", combined_standard$BCSS_event)
unique(combined_standard$BCSS_event)

#standardize pcr (1 = true, 0 = false)
combined_standard$pcr <- gsub("pCR", "1", combined_standard$pcr)
combined_standard$pcr <- gsub("RD", "0", combined_standard$pcr)
combined_standard$pcr <- gsub("NA", NA, combined_standard$pcr)
unique(combined_standard$pcr)

#standardize rcb (1 = true, 0 = false)
combined_standard$rcb <- gsub("NA", NA, combined_standard$rcb)
unique(combined_standard$rcb)

#standardize metastases location
combined_standard$metastases_location <- gsub("Lung |LUNG|Left lung", "Lung", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("Liver |LIVER|Left liver|Right liver", "Liver", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("Brain |BRAIN", "Brain", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("CHEST WALL|ChestWall", "Chest", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("BREAST", "Breast", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("BONE|Bone ", "Bone", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("AXILLARY LYMPH NODE|LYMPH NODE DISTANT|LYMPH NODE REGIONAL|Lymph node", "LymphNode", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("ASCENDING COLON", "Rectum", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("AXILLA", "Axilla", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("FALLOPIAN TUBE", "FallopianTube", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("SOFT TISSUE|Soft tissue|Soft Tissue", "SoftTissue", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("BREAST SKIN|Breast SKIN", "Skin", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("OMENTUM", "Peritoneum", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("PLEURAL FLUID", "Pleura", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("Bone MARROW", "BoneMarrow", combined_standard$metastases_location)
combined_standard$metastases_location <- gsub("GASTRIC BODY", "Stomach", combined_standard$metastases_location)
unique(combined_standard$metastases_location)

#standardize rna expression platform
combined_standard <- combined_standard %>%
  mutate(platform = ifelse(study %in% c("HatzisDiscovery", "HatzisValidation", "METABRIC", "Zhang_U133", "Zhang_U133p2", "ispy2"), "microarray", "rnaseq"))

#Revise treatment_catagory column
combined_standard <- combined_standard %>%
  dplyr::select(-treatment_catagory)
combined_standard$treatment_catagory <- combined_standard$treatment
combined_standard$treatment_catagory <- gsub("T-DM1 \\+ Pertuzumab", "immuno", combined_standard$treatment_catagory)
combined_standard$treatment_catagory <- gsub("Paclitaxel \\+ Ganitumab|Paclitaxel \\+ Pembrolizumab|Paclitaxel \\+ Pertuzumab \\+ Trastuzumab|Paclitaxel \\+ Trastuzumab|Paclitaxel \\+ AMG 386 \\+ Trastuzumab|Paclitaxel \\+ MK-2206 \\+ Trastuzumab|FEC-T \\+ Trastuzumab|EC-T \\+ Trastuzumab|TC \\+ Trastuzumab|FEC-T \\+ Trastuzumab \\+ Pertuzumab|EC \\+ Trastuzumab|T \\+ Pertuzumab \\+ Trastuzumab|T-FEC \\+ Trastuzumab \\+ Pertuzumab|T-FEC \\+ Trastuzumab|TC \\+ Pertuzumab \\+ Trastuzumab|P \\+ Trastuzumab", 
                                             "chemoimmuno", combined_standard$treatment_catagory)
combined_standard$treatment_catagory <- gsub("taxane-anthracycline|Paclitaxel \\+ Ganetespib|Paclitaxel \\+ MK-2206|Paclitaxel|Paclitaxel \\+ Neratinib|Paclitaxel \\+ AMG 386|Paclitaxel \\+ ABT 888 \\+ Carboplatin|T-FEC|FEC-T|TC|T|T-EC|T-Carboplatin|EC-T|P-FEC|P-EC|P-Carboplatin", "chemo", combined_standard$treatment_catagory)
unique(combined_standard$treatment_catagory)
unique(combined_standard$treatment)

#remove samples without microenviornment results
combined_standard <- combined_standard %>%
  filter(!is.na(CAFs))

#export to xlsx
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization")
write.table(combined_standard, "Clinical Associations with Cell Types.txt", row.names = FALSE)

#archive code (Zhang)####
#Zhang_U133 cleaning
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
Zhang_U133 <- fread("Zhang-U133_clinical.tsv")
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
Zhang_U133_immune <- read_xlsx("Zhang-U133 InstaPrism results.xlsx", sheet = 3)

Zhang_U133 <- Zhang_U133 %>%
  dplyr::select(Sample_ID, distant_metastasis, pam50, intclust) %>%
  dplyr::rename(patient_id = Sample_ID, metastases_location = distant_metastasis) %>%
  mutate(metastases_location = gsub("metastasis", "", metastases_location)) %>%
  mutate(study = "Zhang_U133")
Zhang_U133 <- Zhang_U133 %>%
  mutate(
    sample_id = NA,
    ER_status = NA,
    HER2_status = NA,
    grade = NA,
    size = NA,
    node_status = NA,
    age = NA,
    OS_event = NA,
    OS_time = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    MFS_event = NA,
    MFS_time = NA,
    primary = "false",
    treatment_catagory = NA,
    treatment = NA,
    pcr = NA,
    rcb = NA)
Zhang_U133 <- Zhang_U133 %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )

Zhang_U133_merged <- merge(Zhang_U133, Zhang_U133_immune, by = "patient_id", all.x = TRUE)



#Zhang_U133p cleaning ()
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Clinical Datasets")
Zhang_U133p <- fread("Zhang-U133plus2_clinical.tsv")
setwd("~/Library/CloudStorage/OneDrive-CRUKCambridgeInstitute/Immune Characterization/Instaprism Results")
Zhang_U133p_immune <- read_xlsx("Zhang-U133plus2 InstaPrism results.xlsx", sheet = 3)

Zhang_U133p <- Zhang_U133p %>%
  dplyr::select(Sample_ID, distant_metastasis, pam50, intclust) %>%
  dplyr::rename(patient_id = Sample_ID, metastases_location = distant_metastasis) %>%
  mutate(metastases_location = gsub("metastasis", "", metastases_location)) %>%
  mutate(primary = "no")%>%
  mutate(study = "Zhang_U133p2")
Zhang_U133p <- Zhang_U133p %>%
  mutate(
    sample_id = NA,
    ER_status = NA,
    HER2_status = NA,
    grade = NA,
    size = NA,
    node_status = NA,
    age = NA,
    OS_event = NA,
    OS_time = NA,
    RFS_event = NA,
    RFS_time = NA,
    BCSS_event = NA,
    BCSS_time = NA,
    MFS_event = NA,
    MFS_time = NA,
    primary = "false",
    treatment_catagory = NA,
    treatment = NA,
    pcr = NA,
    rcb = NA)
Zhang_U133p <- Zhang_U133p %>%
  dplyr::select(
    patient_id, sample_id, study, ER_status, HER2_status, grade, size, node_status,
    age, intclust, pam50, OS_event, OS_time, RFS_event, RFS_time, BCSS_event, BCSS_time,
    MFS_event, MFS_time, treatment_catagory, treatment, pcr, rcb, primary, metastases_location
  )

Zhang_U133p_merged <- merge(Zhang_U133p, Zhang_U133p_immune, by = "patient_id", all.x = TRUE)




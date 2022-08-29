.libPaths( c( .libPaths(), "/hpc/users/jayarp02/R/lib/4.0.3/") )
.libPaths()

#BiocManager::install("transplantr", lib="/hpc/users/jayarp02/R/lib/4.0.3/")

library("dplyr")
library("stringr")
library("lubridate")
library("transplantr")

setwd(getwd())

scr_data = read.csv("Creatinine_data.csv", header=TRUE, stringsAsFactors = FALSE)
length(unique(scr_data$Subject_ID))

##convert to numeric
scr_data$cr = ifelse(scr_data$LabComponentEpicId == "2991",as.numeric(as.character(scr_data$Value)),0)
scr_data = scr_data[!is.na(scr_data$cr),]

##somalogic has everything we need
load("somalogic_after_aki_encounter_clinical_data_with_protein.RData")

##this wll give you final_df dataframe
merged_1 = inner_join(scr_data, final_df, by = "Subject_ID")

####merge patient class from clinical rds - renamed for public codeuse
clinicalcov = readRDS("Biobank_clinical_data_table.RDS")

#clincov = data.frame(clinicalcov$Subject_ID,clinicalcov$Event_Date,clinicalcov$Encounter_Patient_Class)
#head(clincov)
clincov = clinicalcov[c("Subject_ID","Encounter_Patient_Class")]
table(clincov$Encounter_Patient_Class)

clincov = distinct(clincov)
table(clincov$Encounter_Patient_Class)
clincov = clincov[clincov$Encounter_Patient_Class=="OUTPATIENT",]
head(clincov)
dim(clincov)
#head(merged_1)

print("merging two dfs")
merged = inner_join(clincov, merged_1, by = "Subject_ID")
head(merged)

print("done merging two dfs")

merged$sex_char <- ifelse(merged$sex=="Male", "M","F")
merged$sex_val = ifelse(merged$sex=="Male", 1.000, 0.742)

table(merged$age)
table(merged$cr)
str(merged)

merged$eth = "non-black"

merged$units = "US"

print("starting ckd-epi")

### use non-black for everyone - workaround for race-free calculations

pheno_df_2 <- merged %>%
  mutate(
	eGFR = ckd_epi(creat = cr, age = age, sex = sex_char, eth = eth, units = units),
	eGFR2 = 175 * ((cr)**(-1.154)) * ((age)**(-0.203)) * (sex_val)
	)

### re-compute eGFR - because some values come up with Inf
pheno_df_2$cr <- ifelse(pheno_df_2$cr==0, pheno_df_2$ValueNumeric, pheno_df_2$cr)
pheno_df_2$eGFR <- ckd_epi(creat = pheno_df_2$cr, age = pheno_df_2$age, sex = pheno_df_2$sex_char, eth = pheno_df_2$eth, units = pheno_df_2$units)

print("done with ckd-epi")
length(unique(pheno_df_2$Subject_ID))


saveRDS(pheno_df_2, file="pheno_df_egfr_merged_alllabs.RDS")


setwd("/Users/iparanjpe/Documents/mount_sinai_research/covid19/biobank_somalogic/")
.libPaths( c( .libPaths(), "/Users/jayaramanp/workspace/R/lib/4.04/") )
.libPaths()

library(limma)
library(gtools)
# install.packages("gtools")

save(proteomic, file="data/preprocessed_data/somalogic_after_aki_oct17_2022.RData")


####Load data
source("code/limma_after_aki/load_data.R")

##### 
## Load data frame after being saved
#####


# ### restrict to only proteins significant main effect
# load( file="results/logisitic/all_aki/somalogic_after_aki_limma_all_proteins_only_aki.RData")
# 
# ttable$entrez <- protein_metadata$EntrezGeneID[match(rownames(ttable), protein_metadata$TargetFullName)]
# ttable <- ttable[abs(ttable$FC)>1.2,]
# 
# ttable <- ttable[ttable$adj.P.Val<0.05,]


### change outliers to NA
# final_df[["HBE1"]]
load("somalogic_after_aki_encounter_clinical_data_with_protein (1).RData")
load("somalogic_after_aki_oct17_2022.RData")
load("Somalogics_patients_clinCharacterists.rData")


final_df$max_aki_outcome<- ifelse(final_df$max_aki_stage>1,1,0 )
# final_df <- final_df[!is.na(final_df$max_wbc),]
vars <- colnames(proteomic)
vars <- vars[vars!="sample"]

temp <-  final_df[,colnames(final_df) %in% vars]
z_scores <- scale(temp)

# temp[abs(z_scores)>3] <- NA
table(is.na(temp))
# temp <- scale(temp)
temp <- log2(temp )
# temp <- temp[,rownames(ttable)]
temp[1:4,1:10]


# final_df <- test
library(limma)
# aki_var <- factor(aki_var)
# sofa_binary <-
# prop.table(table(sofa, aki_var), margin = 2)
# severe_aki <- ifelse(max_aki_stage>1,1,0)
# sofa_binary <- sofa>=1

# table(sofa_binary)
# table(final_df$ventilation)
##### PJ testing!
#table(somalogic_patients_RRT_notDup$Uncorrected_Blood_Sample)
#table(final_df$Uncorrected_Blood_Sample)
### --done PJ testing
table(somalogic_patients_RRT_notDup$COMORBID_DIABETES)
table(somalogic_patients_RRT_notDup$COMORBID_HEART_FAILURE)
table(somalogic_patients_RRT_notDup$MED_METHYLPREDNISOLONE)
table(somalogic_patients_RRT_notDup$MED_PREDNISONE)
table(somalogic_patients_RRT_notDup$MED_DEXAMETHASONE)
table(somalogic_patients_RRT_notDup$COMORBID_DIABETES, somalogic_patients_RRT_notDup$COMORBID_HEART_FAILURE)
new_final_df = merge(final_df, somalogic_patients_RRT_notDup, by = 'Uncorrected_Blood_Sample')

final_df = new_final_df
library(tidyverse)
##PJ test 18th october 2022
#med_final_df <- final_df %>% dplyr:: select(starts_with("MED_"))
#colnames(med_final_df)

### Impute WHO ORdinal score
final_df$WHO_Ordinal[is.na(final_df$WHO_Ordinal)] <- median(final_df$WHO_Ordinal, na.rm = T)
# design <- model.matrix(~0+age+ckd+sex+min_platelet+max_bun+modified_sofa*aki, data = final_df[!is.na(final_df$min_platelet),])
# design <- model.matrix(~0+age+ckd+sex+min_platelet+max_bun+modified_sofa*aki, data = final_df[!is.na(final_df$min_platelet) &!is.na(final_df$max_bun),])

# design <- model.matrix(~0+age+ckd+sex+modified_sofa*aki, data = final_df)
# design <- model.matrix(~0+age+ckd+sex+ventilation+severe_aki, data = final_df)
# design <- model.matrix(~0+age+ckd+sex+ventilation+max_wbc+max_aki_outcome, data = final_df[!is.na(final_df$max_wbc),])

# design <- model.matrix(~0+age+ckd+sex+WHO_Ordinal+severe_aki, data = final_df)
#head(design)
### delta_cr
#design <- model.matrix(~age+ckd+sex+ventilation+delta_cr, data = final_df)
design <- model.matrix(~age+ckd+sex+ventilation+
                         #COMORBID_HEART_FAILURE+
                         #COMORBID_DIABETES+
                       #MED_METHYLPREDNISOLONE+MED_DEXAMETHASONE+MED_PREDNISONE+ 
                         max_aki_outcome, data = final_df)

head(design)
fit <-lmFit(t(temp), design)

# fit <-lmFit(t(temp[!is.na(final_df$min_platelet),]), design)
# fit <-lmFit(t(temp[!is.na(final_df$min_platelet) &!is.na(final_df$max_bun),]), design)

cont.matrix <- makeContrasts(contrasts="max_aki_outcome", levels=design)
fit2<- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
# fit2 <- eBayes(fit)

FDR_cutoff <- 0.05
ttable <- topTable(fit2, n= Inf, adjust.method = "bonferroni", coef = "max_aki_outcome", p.value =1)
# ttable <- topTreat(fit2, n= Inf, adjust.method = "BH", coef = "aki")

table(ttable$adj.P.Val<0.05)
ttable$FC <- logratio2foldchange(ttable$logFC, base = 2)
# ttable$FC  <- 2^ttable$FC
#rownames(ttable)[!rownames(ttable) %in% protein_metadata$TargetFullName]
#ttable$gene_name <- protein_metadata$EntrezGeneSymbol[match(rownames(ttable), protein_metadata$TargetFullName)]
ttable["Platelet factor 4",]
ttable[ttable$adj.P.Val<FDR_cutoff,]
old_ttable <- ttable

head(ttable)


save(ttable, file="somalogic_adjust_ventilation_aki_18thOct22_mergedcols.RData")

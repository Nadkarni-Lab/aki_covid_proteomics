####
## Find proteins associated with in-hospital AKI
####

setwd("/Users/iparanjpe/Documents/mount_sinai_research/covid19/biobank_somalogic/")
library(limma)

####Load data
source("code/limma_after_aki/load_data.R")

### Create design matrix for limma
### Model 1
### severe_aki = AKI stage 2 or 3
design <- model.matrix(~age+ckd+sex+ventilation+severe_aki, data = final_df)

### Model 2
### Max change in creatinine
design <- model.matrix(~age+ckd+sex+ventilation+delta_cr, data = final_df)



#### Run Limma for each model 
head(design)
fit <-lmFit(t(temp), design)
fit2 <- eBayes(fit)

FDR_cutoff <- 0.05
ttable <- topTable(fit2, n= Inf, adjust.method = "bonferroni", coef = "delta_cr", p.value =1)
ttable$FC <- logratio2foldchange(ttable$logFC, base = 2)
rownames(ttable)[!rownames(ttable) %in% protein_metadata$TargetFullName]
ttable$gene_name <- protein_metadata$EntrezGeneSymbol[match(rownames(ttable), protein_metadata$TargetFullName)]



save(ttable, file="results/somalogic.RData")

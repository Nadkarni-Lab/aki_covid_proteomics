setwd("/Users/iparanjpe/Documents/mount_sinai_research/covid19/biobank_somalogic/")
#BiocManager::install('PCAtools')
library(PCAtools)

####Load data
source("code/limma_after_aki/load_data.R")
vars <- colnames(proteomic)
vars <- vars[vars!="sample"]
temp <-  final_df[,colnames(final_df) %in% vars]

# temp <- scale(temp)
temp <- log2(temp )



# metadata <- data.frame(aki= as.numeric(aki_var), sex= sex, )
metadata <- final_df[,c("sex","age","ckd","aki","ventilation", "max_aki_stage")]
metadata$max_aki_outcome <- ifelse(metadata$max_aki_stage >1, 1, 0)
metadata$max_aki_outcome <- as.factor(metadata$max_aki_outcome)

metadata$aki <- as.factor(metadata$aki)
metadata$sex <- ifelse(metadata$sex=="Male",1,0)
rownames(metadata) <- colnames(final_df_t)


### Regress out covariates
resdesign <- model.matrix(~age+ckd+sex+ventilation
                          #+COMORBID_HEART_FAILURE
                          #+COMORBID_DIABETES
                          #+MED_METHYLPREDNISOLONE+MED_DEXAMETHASONE+MED_PREDNISONE
                          , final_df)



vobjDream = t(temp)
fitmm_res <-lmFit(vobjDream, resdesign)


# Fit Empirical Bayes for moderated t-statistics
fitDupCor_res <- eBayes( fitmm_res )

print("names fitDupCor")
names(fitDupCor_res)


print("fitting residual")
res = residuals(fitDupCor_res, vobjDream)
dim(res)
res[1:4,1:4]

### Run pca
pca <- PCAtools::pca(res, metadata = metadata,scale = T )
colkey = c("0" = "gold", "1"="blue")
pairsplot(pca,colby = "max_aki_outcome" , legendPosition = "none",colkey=colkey,
          components = c("PC1","PC2","PC3","PC4"))



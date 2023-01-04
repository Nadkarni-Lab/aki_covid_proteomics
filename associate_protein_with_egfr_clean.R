# library("tidyverse")
library("stringr")
library("lubridate")

library("ggplot2")
library(plotly)
setwd("/Users/iparanjpe/Documents/mount_sinai_research/covid19/biobank_somalogic/")

### Load post discharge eGFR data
load("data/clarity/processed_data_frames/pheno_df_egfr_merged_final (2).RData")

### Load clinical and protein expression data
load("data/preprocessed_data/somalogic_after_aki_encounter_clinical_data_with_protein.RData")

### Subset to egfr after discharge
pheno <- pheno[pheno$date> (pheno$chosen_somalogic_date + (pheno$hospital_stay- pheno$somalogic_hospital_day)),]

pheno<- pheno %>% group_by( Subject_ID, date) %>% 
  slice(1) %>%
  as.data.frame()

pheno <- pheno[pheno$eGFR<150,]

pheno$days_since_somalogic <- as.numeric(pheno$date-  pheno$chosen_somalogic_date)

pheno$days_since_discharge <- as.numeric(pheno$date- (pheno$chosen_somalogic_date + (pheno$hospital_stay- pheno$somalogic_hospital_day)))




#########  Plot  #######################
protein_name <- "Trefoil factor 3"
cutoff <- quantile(log2(final_df[[protein_name]]), na.rm = T, probs = c(.33,.67,1))
cutoff <- c(-Inf, cutoff)
pheno$protein_cat <- cut(log2(pheno[[protein_name]]), cutoff, labels=c("low","intermediate","high"))
summary(pheno$protein_cat)

ggplot(data = pheno, aes(x = days_since_somalogic, y =eGFR, group=protein_cat, color= factor(protein_cat)))+
  stat_smooth()+
  theme_minimal() +
  ylab("eGFR") +
  xlab("Days since protein measurement")+
  labs(color=paste0(protein_name, " level"))+
  scale_color_brewer(palette="Set1")+
  theme(axis.title = element_text(size = 0),
        plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"),
        axis.text = element_text(size = 35, color="black", face="bold"),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 15),
        legend.position = "none",
        title =element_text(size=20, face='bold'),
        # plot.title = element_text(size = 30,hjust = 0.5), 
  )


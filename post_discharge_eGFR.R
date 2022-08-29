########
#### Find proteins associated with post-discharge eGFR measurements
########

.libPaths()
library("tidyverse")
library("stringr")
library("lubridate")
library("transplantr")
library("ggplot2")
setwd(getwd())


#########  ANALYZE  #######################

### Load eGFR measurements
pheno <- readRDS("pheno_df_egfr_merged_alllabs.RDS")

## Load clinical covariates
clinicalcov = readRDS("Biobank_clinical_data_table_by_RNA_sample.RDS")

clincov = data.frame(clinicalcov$Subject_ID,clinicalcov$Event_Date,clinicalcov$Encounter_Patient_Class)

### Change date format
pheno$date <-lubridate::as_date(pheno$Event_Date)

### Select somalogic date
pheno = pheno%>% group_by(Subject_ID) %>% filter(date>=chosen_somalogic_date)

#### Plot graph
pdf("plot_egfr.pdf", height = 100, width = 500)
p <- ggplot(data = pheno, aes(x = date, y = eGFR, group = Subject_ID, col = Subject_ID))
p + geom_line()


ggplot(pheno, aes(x=Subject_ID, y=eGFR, col = Subject_ID)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
                outlier.size=4) + 
  geom_jitter(shape=16, position=position_jitter(0.2))


pheno_summary = pheno %>% group_by(Subject_ID) %>% 
		summarize(
			Subject_ID = unique(Subject_ID),
			num_egfr = n(),
			age = unique(age),
			sex = unique(sex),
			death = unique(death),
			median_egfr = median(eGFR, na.rm = TRUE),
			max_egfr = max(eGFR, na.rm=TRUE),
			min_egfr = min(eGFR, na.rm = TRUE),
			mean_egfr = mean(eGFR, na.rm = TRUE),
			var_egfr = var(eGFR, na.rm = TRUE),
			history_ckd = unique(ckd),
			history_aki = unique(aki)
			

)

ggplot(pheno_summary, aes(x=Subject_ID, y=num_egfr, fill=Subject_ID)) + 
	geom_bar(stat="identity")+ facet_grid(sex~.)

dev.off()
saveRDS(pheno_summary, file = "pheno_summary.RDS")

save(pheno, file = "pheno_df_egfr_merged_final.RData")

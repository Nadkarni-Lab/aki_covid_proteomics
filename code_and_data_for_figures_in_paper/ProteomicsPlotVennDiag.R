.libPaths( c( .libPaths(), "/Users/jayaramanp/workspace/R/lib/4.04/") )
.libPaths()
library(devtools)
#install.packages("sf", lib="/Users/jayaramanp/workspace/R/lib/4.04/", dependencies = TRUE)
#install.packages("ggVennDiagram", lib="/Users/jayaramanp/workspace/R/lib/4.04/")
#install.packages("nVennR", lib="/Users/jayaramanp/workspace/R/lib/4.04/")
library(tidyverse)
library(ggVennDiagram)
library(EnhancedVolcano)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
#d <- loadRData("~/blah/ricardo.RData")

original <- loadRData("somalogic_adjust_ventilation_aki_18thOct22_mergedcols.RData")
model_HF <- loadRData("somalogic_adjust_ventilation_aki_18thOct22_mergedcols_HF.RData")
model_DM <- loadRData("somalogic_adjust_ventilation_aki_18thOct22_mergedcols_DM.RData")
model_DM_HF <- loadRData("somalogic_adjust_ventilation_aki_18thOct22_mergedcols_DM_HF.RData")
model_MEDS <- loadRData("somalogic_adjust_ventilation_aki_18thOct22_mergedcols_MEDS.RData")
model_DM_MEDS <- loadRData("somalogic_adjust_ventilation_aki_18thOct22_mergedcols_DM_MEDS.RData")
model_HF_MEDS <- loadRData("somalogic_adjust_ventilation_aki_18thOct22_mergedcols_HF_MEDS.RData")
model_HF_DM_MEDS <- loadRData("somalogic_adjust_ventilation_aki_18thOct22_mergedcols_HF_DM_MEDS.RData")

rownames(original[original$adj.P.Val<0.05,])

pdf("volcanoPlot.pdf")
EnhancedVolcano(toptable = original,
                x = "logFC",
                y = "adj.P.Val",
                lab = rownames(original),
                xlim = c(-3, +3),
                ylim = c(0,50),
                pCutoff = 5e-02,  
                FCcutoff = 1.3, 
                labSize = 2.0,
                pointSize = 0.8,
                title = "AKI in COVID vs noAKI Proteomics \n (fold change cutoff = 1.5, p-value cutoff = 5e-02)"
)
dev.off()

set.seed(20210419)
x <- list(
  base_model=rownames(original[original$adj.P.Val<0.05,]),
  DM=rownames(model_DM[model_DM$adj.P.Val<0.05,]),
  HF=rownames(model_HF[model_HF$adj.P.Val<0.05,]),
  MEDS=rownames(model_MEDS[model_MEDS$adj.P.Val<0.05,])
  #DM_HF=rownames(model_DM_HF[model_DM_HF$adj.P.Val<0.05,]),
  #DM_MEDS=rownames(model_DM_MEDS[model_DM_MEDS$adj.P.Val<0.05,]),
  #HF_MEDS=rownames(model_HF_MEDS[model_HF_MEDS$adj.P.Val<0.05,]),
  #DM_HF_MEDS=rownames(model_HF_DM_MEDS[model_HF_DM_MEDS$adj.P.Val<0.05,])
  )



library(ggplot2)
library(nVennR)
myNV <- plotVenn(x, showPlot = T)
showSVG(myNV, outFile="nVennR.svg", opacity=0.2)

pdf("venndiag_DM_HF_MEDS.pdf")
ggVennDiagram(x) + scale_fill_gradient(low="blue",high = "red")
region_data = process_region_data(Venn(x))
dev.off()
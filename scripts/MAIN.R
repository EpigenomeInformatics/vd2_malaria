# nolint start
####################################################################################
# Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
# Script: MAIN.R
# Author: MIJM
# Last Update: 27.08.26
# Description: Full analysis for ATACseq data
####################################################################################

source("/icbb/projects/mmaran/SCRIPTS/CCB/malariaScripts/utils.R")

plotDir <- "/icbb/projects/mmaran/malaria/2025/Manuscript/LONG/BatchCorDM/" #paste0(getwd(), "/Manuscript/LONG/BatchCorDM/")


################
# Load object
################

OBJECT <- "long" #c("long", "regional", )

data_ <- loadOBJECT_(OBJECT)
data_ <- removeSamples_(data_)


#TODO:
1. exploratory

2. chromvAR



GRP <- "EirGrp"
GRP1 <- "high"
GRP2 <- "low"
CUTOFF <- 0.05

group1 <- getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP1] #%>% droplevels()
group2 <- getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP2] #%>% droplevels()

cv <- "Altius" #c("JASPAR", "Altius")
GR <- ".peaks.cons"

if(cv =="JASPAR"){
  ### chromVAR with JASPAR motifs
  cvRes <- getChromVarDev(data_, GR)  
}else{
  ### ChromVAR Motif Family
  cvRes <- chromVAR_Altius(data_, GR)
}
devZ <- chromVAR::deviationScores(cvRes)


CORRECTION <- "donor"
ann <- getSampleAnnot(data_)
ann <- ann[colnames(devZ),]
identical(rownames(ann), colnames(devZ))

design <- model.matrix(~ EirGrp, data = ann)
devZ <- limma::removeBatchEffect(devZ,
                           batch      = ann$donor,
                           design     = design)



mycolors <- c("darkred", "darkblue", "grey")
names(mycolors) <- c(GRP1, GRP2, "NO")

# nolint end

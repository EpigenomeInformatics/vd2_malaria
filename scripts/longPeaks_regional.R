# nolint start
####################################################################################
# Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
# Script: longPeaks_regional
# Author: MIJM
# Last Update: 27.08.26
# Description: Aggregating Regional data with Longitudinal peakset 
####################################################################################

source("/icbb/projects/mmaran/malariaScripts/utils.R")

# ### load datasets
long_ <- loadOBJECT_("long")

PATH = "/icbb/projects/mmaran/malaria/pilot/run_ChrAccR_v202004/ChrAccR_run/data/dsATAC_filtered"
regional_ <- loadDsAcc(PATH)
regional_@sampleAnnot <- getSampleAnnot(regional_)[, 1:14]

#remove US and CD4 samples
temp <- getSampleAnnot(regional_)[,"cellType"] %in% c("CD4") | getSampleAnnot(regional_)[,"region"] %in% c("US")
regional_ <- removeSamples(regional_, temp)
rm(temp)


# call peaks on long data: load long info prior
coord_ <- cleanCoord_(long_, ".peaks.cons")


regional_aggregated <- regionAggregation(
  regional_,
  coord_,
  type= "longPeaks",
  signal = "insertions",
  # aggrFun = "sum",
  dropEmpty = TRUE,
  bySample = TRUE
  # chunkSize = 5000L
)

saveDsAcc(Pvd2, "/icbb/projects/mmaran/malaria/pilot/longPeaksYear0.5/", forceDiskDump = FALSE, updateDiskRef = TRUE)

# nolint end

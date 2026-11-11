# nolint start
####################################################################################
#' Project:     Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
#' Script:      longPeaks_regional.R
#' Author:      MIJM
#' Last Update: 2026-11-02
#' Description: Re-quantify regional cohort VD2 samples (Uganda only) over the
#'              longitudinal consensus peak set and save the aggregated dsATAC object
#'              used as the "regional" object in loadOBJECT_().
#' Depends:     scripts/ATACseq/utils.R
#' Inputs:      loadOBJECT_("long"); pilot dsATAC_filtered (cluster path)
#' Outputs:     dsATAC saved to PATH_regional_longPeakYear0.5 (cluster path)
#' Usage:       Rscript scripts/ATACseq/longPeaks_regional.R  # run once, before REGIONAL_chromVAR.R
####################################################################################

source("./scripts/ATACseq/utils.R")

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

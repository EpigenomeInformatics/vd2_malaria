# nolint start
####################################################################################
#' Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
#' Script: utils.R
#' Author: MIJM
#' Last Update: 27.08.26
#' Description: ATACseq's All libraries, utils, and helper functions. Load conda atac
#' Functions:
#'   - loadOBJECT_: load dsATAC oject
#'      @return dsATAC object
#'   - removeSamples_:
#'      @param 
#'      @return cleaned dsATAC obejct
#'  - cleanCoord_
#'      @param
#'      @return 
#'  - takeResLimma

####################################################################################


# 1. Load libraries
set.seed(123)

suppressPackageStartupMessages({
  library(ChrAccR)
  library(factoextra)
  library(FactoMineR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(DESeq2)
  library(ggrepel)
  library(LOLA)
  library(ggrastr) 
  library(circlize)
})

options(ggrepel.max.overlaps = Inf)

plotDir <- "results/ATACseq/"

# 2. Global Variables & Paths
PATH_long_0.5 <- "/icbb/projects/mmaran/malaria/dsATAC/longSaturationAnalysis/vd2.year0.5/data/dsATAC_filtered"
PATH_regional_longPeakYear0.5 <- "/icbb/projects/mmaran/malaria/pilot/longPeaksYear0.5/"
PATH_ANNOT <- "/icbb/projects/mmaran/malaria/2025/comprehensiveSampleAnnot.csv"

RESPONDERS <- c(
  "VD2_2016_d3340_88",  "VD2_2016_d3381_164", "VD2_2016_d3401_216",
  "VD2_2016_d3278_60",  "VD2_2016_d3384_192", "VD2_2016_d3280_168",
  "VD2_2016_d3369_112", "VD2_2013_d3278_52",  "VD2_2013_d3340_80", 
  "VD2_2014_d3381_160", "VD2_2013_d3401_208", "VD2_2013_d3280_44", 
  "VD2_2013_d3369_104", "VD2_2014_d3369_108", "VD2_2014_d3398_200"
)

# 3. FUNC: loadOBJECT_

loadOBJECT_ <- function(obj_type) {
  if (obj_type == "long") {
    long <- loadDsAcc(PATH_long_0.5) 
    
    # Sample Swap
    swap_samples <- function(count_matrix) {
      as.data.frame(count_matrix) %>% 
        rename(VD2_2015_d3092_248 = "VD2_2016_d3345_244", 
               VD2_2016_d3345_244 = "VD2_2015_d3092_248") %>%
        as.matrix()
    }
    
    for (assay in c(".peaks.cons", "tiling200bp")) {
      long@counts[[assay]] <- swap_samples(long@counts[[assay]])
      long@sampleAnnot <- long@sampleAnnot[colnames(long@counts[[assay]]), ]
    }
    
    # Comprehensive Annotation
    cols_to_factor <- c("EirGrp", "LAMP", "malariacat", "anyparasitemia",  
                        "eirB", "priormal12B", "priormal12Grp", "YearGrp")
    
    comp_annot <- read.csv(PATH_ANNOT, row.names = 1) %>%
      mutate(priormal12Grp = ifelse(clin_priormalincidence_last12months < 2, "low", "high")) %>%
      mutate(across(all_of(cols_to_factor), as.factor))
    
    stopifnot("Sample IDs do not match exact row order" = identical(rownames(comp_annot), rownames(long@sampleAnnot)))
    long@sampleAnnot[cols_to_factor] <- comp_annot[cols_to_factor]
    
    annot_df <- getSampleAnnot(long) %>% 
      mutate(
        responders = as.factor(ifelse(sampleId %in% RESPONDERS, "yes", "no")),
        alt_aeir = as.factor(paste0(responders, "_", EirGrp))
      )
    
    long <- addSampleAnnotCol(long, "alt_eir", annot_df$alt_aeir)
    long <- addSampleAnnotCol(long, "responders", annot_df$responders)
    
    return(long)
    
  } else if (obj_type == "regional") {
    regional <- loadDsAcc(PATH_regional_longPeakYear0.5)
    regional@sampleAnnot <- getSampleAnnot(regional)[, 1:14]
    annot_pilot <- getSampleAnnot(regional) %>% mutate(malariaPrevYear =  ifelse(malariaPrevYear == 0, "Negative", "Positive"))
    regional <- addSampleAnnotCol(regional, "malariaPrevYear_", annot_pilot$malariaPrevYear)
    return(regional)
  } else {
    stop("obj_type must be either 'long' or 'regional'")
  }
}

# 4. FUNC: removeSamples_

removeSamples_ <- function(data_obj, type = "EirGrp") {
  annot <- getSampleAnnot(data_obj)
  lev <- levels(as.factor(annot[[type]]))

    # Total: 46 samples
    # >1: 42 samples
    # aeir TP: 33 samples
    # Responders: 13 samples

  if (type == "responders") {
    # Remove non-responders
    to_remove <- annot$responders == "no"
    data_obj <- removeSamples(data_obj, to_remove)
    
  } else {
    # Keep only donors in BOTH groups
    valid_donors <- annot %>%
      group_by(donor) %>%
      filter(any(.data[[type]] == lev[1]) & any(.data[[type]] == lev[2])) %>%
      pull(donor) %>% 
      unique()
    to_remove <- !annot$donor %in% valid_donors
    data_obj <- removeSamples(data_obj, to_remove)
  }
    # Drop donor levels
    data_obj <- addSampleAnnotCol(data_obj, "donor", droplevels(getSampleAnnot(data_obj)$donor))
  
  return(data_obj)
}

# 5. FUNC: cleanCoord_

cleanCoord_ <- function(OBJECT, RT= ".peaks.cons"){
    getCoord(long_, ".peaks.cons")
    tssGr <- ChrAccRAnnotationHg38::getGeneAnnotation(anno="gencode_coding", type="tssGr")
    peakTssDist <- GenomicRanges::distanceToNearest(coord_, tssGr, ignore.strand=TRUE)
    mcols(coord_)$nearestGeneDist <- S4Vectors::mcols(peakTssDist)[,"distance"]
    mcols(coord_)$nearestGene <- S4Vectors::elementMetadata(tssGr)[S4Vectors::subjectHits(peakTssDist),"gene_name"]
    colnames(data.frame(coord_))
    names(coord_) <- sub("\\..*?_peak", "_peak", names(coord_))
    return(coord_)
}


# 6. Colors for donors
colors <- c("black", "red", "blue", "green3", "gold", "purple", "deeppink", "gray50", "cyan", "magenta", "darkgreen", "navy", "darkred", "deepskyblue",  
"darkviolet", "yellowgreen", "chocolate", "turquoise4", "orangered", "sienna")  

# 7. FUNC: takeResLimma

takeResLimma_ <- function(cm, annot, corCol, compCol="EirGrp") {

    annot <- as.data.frame(annot)

    stopifnot(all(colnames(cm) == rownames(annot)))

    f <- as.formula(paste("~", compCol))
    design <- model.matrix(f, data = annot)
    res_ <- limma::removeBatchEffect(as.matrix(cm),
                              batch      = annot[[corCol]],
                              design     = design)

    
    return(res_)
}

# nolint end

# nolint start
####################################################################################
#' Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
#' Script: chromVAR.R
#' Author: MIJM
#' Last Update: 27.08.26
#' Description: chormVAR analysis with JASPAR2018 or Altius motif families
#' Functions:
#'   - 
#'      @return dsATAC object
#'  - chromVAR altius
####################################################################################
source("/icbb/projects/mmaran/SCRIPTS/CCB/malariaScripts/utils.R")


# 1. chromVAR Altius
# code snippet from F.Muller

chromVAR_Altius <- function(data_, GR){
  print("Loading chromVAR Motifs")

  cvMot <- readRDS("/scratch/icbb/mmaran/malaria/tfMotifClusters_hg38_uniqueClust.rds")

  OBJECT <- data_
  # RT <- "peaks"

  gr <- getCoord(OBJECT, GR)
  countSe <- getCountsSE(OBJECT, GR, naIsZero=TRUE)
  motifClusterOccM <- as(do.call(cbind, lapply(cvMot$clusterOcc, FUN=function(x){
    overlapsAny(gr, x, ignore.strand=TRUE)
  })), "sparseMatrix")
  mmObj_mc <- SummarizedExperiment::SummarizedExperiment(
    assays=list(motifMatches=motifClusterOccM), 
      rowRanges=gr, 
      colData=S4Vectors::DataFrame(cvMot$clusterAnnot))

  genomeObj <- getGenomeObject(getGenome(OBJECT))
  genome(countSe) <- S4Vectors::metadata(genomeObj)$genome #override inconsistent naming of genome versions (e.g. hg38 and GRCh38)
  countSe <- chromVAR::addGCBias(countSe, genome=genomeObj)
  ridx <- safeMatrixStats(assay(countSe), "rowSums") > 0 & safeMatrixStats(assay(mmObj_mc), "rowSums") > 0 # only consider regions where there is an actual motif match and counts

  cvRes <- chromVAR::computeDeviations(object=countSe[ridx,], annotations=mmObj_mc[ridx,])
  return(cvRes)
}

#' TODO: Wilcoxon Test
#'  @return diff_acc

zDiff <- function(data_, devZ, GRP1 = GRP1, GRP2 = GRP2){

    group1 <- getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP1] #%>% droplevels()
    group2 <- getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP2] #%>% droplevels()

    group1 <- as.matrix(devZ[,group1])
    group2 <- as.matrix(devZ[,group2])

    results <- lapply(1:nrow(group1), function(i) {
    wilcox.test(group1[i, ], group2[i, ], alternative = "two.sided", paired = FALSE)
    })

    diff_acc <- data.frame(pval = sapply(results, function(res) res$p.value))
    diff_acc <- diff_acc %>% mutate(p_value_adjusted = p.adjust(pval, method = "fdr"),
                Significant = p_value_adjusted < 0.05) 
    rownames(diff_acc) <- rownames(group1)
    return(diff_acc)

}



# filter
# High 35, low 34 Param
# 25, 29          nonParam

annot_zDiff <- function(diff_acc, data_, GRP1 = GRP1, GRP2 = GRP2, CUTOFF = 1, padj = 0.05){
    
    group1 <- getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP1] #%>% droplevels()
    group2 <- getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP2] #%>% droplevels()

    group1 <- as.data.frame(devZ[,group1])
    group2 <- as.data.frame(devZ[,group2])

    group1_ <- rowMeans(group1)
    group2_ <- rowMeans(group2)
    grp <- data.frame(group1_, group2_)

    grp$zDiff <- grp$group1_ - grp$group2_
    grp$padj <- diff_acc$p_value_adjusted
    grp$name <- rownames(grp)
    names(grp) <- c(GRP1, GRP2, "zDiff", "padj", "name")

    grp$DiffAcc <- "NO"
    grp$DiffAcc[grp$zDiff > CUTOFF & grp$padj < padj] <- GRP1
    grp$DiffAcc[grp$zDiff < -(CUTOFF) & grp$padj < padj] <- GRP2

    return(grp)
}


# nolint end

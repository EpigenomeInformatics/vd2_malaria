# nolint start
####################################################################################
# Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
# Script: RNAseq.R
# Author: MIJM
# Last Update: 12 Dec 2025
#' Description: RNAseq Analysis
# Combined invivo and invitro
# Common Normalization before biotype. Same background set for all analyses
####################################################################################

source("./scripts/RNAseq/RNA_utils.R")

#' Main analysis
#' Per analysis: > 0 samples per group per donor 
#' dgeMaster: Imported from RNA_utils.R

dge <- dgeMaster 
groupCol <- "aeir_group" 

if(groupCol == "invitro"){
  print("INVITRO")
  targetSamples <- inVitroSamples

} else {
  tempAnnot <- sampleAnnot
  tempAnnot <- droplevels(tempAnnot[!tempAnnot[[groupCol]] == "", ])

  grp1 <- levels(factor(tempAnnot[[groupCol]]))[1]
  grp2 <- levels(factor(tempAnnot[[groupCol]]))[2]

  filteredDf <- tempAnnot %>%
    group_by(donors) %>% 
    filter(any(.data[[groupCol]] == grp1) & any(.data[[groupCol]] == grp2)) %>%      
    data.frame()
  
  targetSamples <- allSamples[allSamples %in% filteredDf$sampleId]
}

tempAnnot <- tempAnnot %>% filter(sampleId %in% targetSamples)

dge$samples$group <- factor(sampleAnnot[[groupCol]], levels = groupLevels[[groupCol]])
dge <- dge[, rownames(dge$samples) %in% targetSamples]
dge$samples$donors <- droplevels(dge$samples$donors)
table(dge$samples$donors)

# DEG analysis

if(groupCol == "invitro"){
  designMatrix <- model.matrix(~group, data = dge$samples)
} else {
  designMatrix <- model.matrix(~group + donors, data = dge$samples)
}

y <- estimateDisp(dge, designMatrix, robust = TRUE)
fit <- glmFit(y, designMatrix, robust = TRUE)
lrt <- glmLRT(fit, coef = 2)

topGenes <- as.data.frame(topTags(lrt, n = nrow(lrt)))
log10Fdr <- -log10(topGenes$FDR)
log10Fdr[log10Fdr == 0] <- 1e-10
topGenes$log10FDR <- log10Fdr

isDe1 <- decideTests(lrt, p.value = 0.05, lfc = 0.6)

print(paste(groupCol, "summary:"))
print(summary(isDe1))

lev <- levels(dge$samples$group)
degDf <- topGenes
degDf$diffExp <- "NO"
degDf$diffExp[degDf$FDR <= 0.05 & degDf$logFC >= 0.6] <- lev[2]                                                           
degDf$diffExp[degDf$FDR <= 0.05 & degDf$logFC <= -0.6] <- lev[1]
# write.csv(degDf, paste0(plotDir, "RNA_vd2.", groupCol ,".csv"))

geneLabels <- degDf[degDf$name %in% genesOfInterest, ]

#nolint end
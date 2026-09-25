# nolint start
####################################################################################
# Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
# Script: RNAseq.R
# Author: MIJM
# Last Update: 25Sept26
#' Description: Saving data for P2G analysis
####################################################################################

source("./scripts/RNAseq/RNAseq.R")

# For P2G
# Using log-converted counts

dataDir <- "./data/RNA"


logCounts <- cpm(dge, log = TRUE, prior.count = 2)
saveRDS(logCounts, file.path(dataDir, "vd2LogCounts.rds"))

geneInfo <- geneInfo[geneInfo$genes %in% dgeMaster$genes$genes, ]

geneInfo <- geneInfo %>% 
  select(seqnames, start, end, width, strand, type, gene_type, name, genesVersion, genes)

saveRDS(logCounts, file.path(dataDir, "aeirLogCounts.rds"))
saveRDS(geneInfo, file.path(dataDir, "rnaGR.rds"))

saveRDS(sampleAnnot, file.path(dataDir, "rnaAnnot.rds"))

# nolint end

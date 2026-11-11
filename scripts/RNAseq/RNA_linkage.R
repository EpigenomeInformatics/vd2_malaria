# nolint start
####################################################################################
#' Project:     Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
#' Script:      RNA_linkage.R
#' Author:      MIJM
#' Last Update: 2026-11-02
#' Description: Save log-CPM counts, gene coordinates and sample annotation as inputs
#'              for peak-to-gene (P2G) linkage.
#' Depends:     scripts/RNAseq/RNAseq.R
#' Inputs:      dge, dgeMaster, geneInfo, sampleAnnot (from RNAseq.R)
#' Outputs:     data/RNA/vd2LogCounts.rds, data/RNA/aeirLogCounts.rds,
#'              data/RNA/rnaGR.rds, data/RNA/rnaAnnot.rds
#' Usage:       Rscript scripts/RNAseq/RNA_linkage.R  # from repository root
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

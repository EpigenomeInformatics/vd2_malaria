# nolint start
####################################################################################
# Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
# Script: RNAseq.R
# Author: MIJM
# Last Update: 12 Dec 2025
#' Description: RNAseq Analysis
# Combined invivo and invitro
# Common Normalization before biotype
# Same background set for all analyses
# 7 Jan 2026
# ADDED GENCODE to remove sex chr
####################################################################################

source("/icbb/projects/mmaran/SCRIPTS/GitHubRepos/vd2_malaria/scripts/RNAseq/RNAseq.R")


# for P2G

#use log converted counts

logCounts <- cpm(DGE, log = TRUE, prior.count = 2)
dim(logCounts)
saveRDS(logCounts, paste0(plotDir, "vd2LogCounts.rds"))

genes_ <- genes_[genes_$genes %in% dge_master$genes$genes,]

summary(genes_$genes %in% dge_master$genes$genes)
summary(dge_master$genes$genes %in% genes_$genes)

genes_ <- genes_ %>% select(seqnames, start, end, width, strand, type, gene_type, name, genes_version, genes)
head(genes_)

rnagr <- rnagr %>% 
            rename(chromosome_name = "seqnames",
                  start_position = "start",
                  end_position = "end") %>%
            mutate(strand = ifelse(strand %in% "-1", "-", "+"),
                    seqnames = paste0("chr", seqnames)) %>% 
            select(-genes)
head(rnagr)

saveRDS(logCounts, paste0(plotDir, "aeirLogCounts.rds"))
saveRDS(genes_, paste0(plotDir, "rnaGR.rds"))
saveRDS(annot, paste0(plotDir, "rnaAnnot.rds"))

# nolint end

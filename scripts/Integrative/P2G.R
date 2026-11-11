# nolint start
####################################################################################
#' Project:     Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
#' Script:      P2G.R
#' Author:      MIJM
#' Last Update: 2026-11-11
#' Description: Peak-to-gene (P2G) linkage. Pearson correlation between donor-corrected
#'              ATAC accessibility of differential peaks and RNA expression of DEGs, for
#'              every peak within ±250 kb of a gene's TSS. ATAC: log2 CPM of consensus
#'              peaks, donors with >= 2 samples; RNA: log CPM. Both limma-corrected
#'              for donor.
#' Depends:     scripts/ATACseq/utils.R (loadOBJECT_, removeSamples_, takeResLimma_)
#' Inputs:      RNA_vd2.aeir_group.csv, aeirLogCounts.rds, rnaGR.rds, rnaAnnot.rds (rnaPath)
#'              ./data/LONG_vd2.EirGrp-AltiusSubset.diffacc.csv
#' Outputs:     results/integrative/p2g-Correlation.csv
#'              Session objects for P2G_plotting.R: data_ (responder samples, with
#'              donors column), atacDiffAcc (full diffAcc table), correlation_
#' Usage:       source("./scripts/Integrative/P2G.R")  # from repository root;
#'              then source("./scripts/Integrative/P2G_plotting.R") in the same session
####################################################################################

source("./scripts/ATACseq/utils.R")

library(data.table)


rnaPath <- "/icbb/projects/mmaran/malaria/2025/Manuscript/RNA/COMBINED/noSexChr"
out_dir <- file.path(getwd(), "results/integrative")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

FLANK <- 250000  # bp either side of each gene's TSS


################
# 1. RNA: DEGs, log CPM, donor-corrected
################

rnaDEG <- read.csv(file.path(rnaPath, "RNA_vd2.aeir_group.csv"), row.names = 1) %>%
  filter(diffExp != "NO")

rnacm <- readRDS(file.path(rnaPath, "aeirLogCounts.rds")) %>% as.matrix()
rnacm <- rnacm[rownames(rnacm) %in% rownames(rnaDEG), ]

rnagr <- readRDS(file.path(rnaPath, "rnaGR.rds"))
rnagr <- rnagr[rownames(rnagr) %in% rownames(rnaDEG), ]

rnaAnnot <- readRDS(file.path(rnaPath, "rnaAnnot.rds")) %>%
  filter(sampleId %in% colnames(rnacm)) %>%
  arrange(match(sampleId, colnames(rnacm)))
rownames(rnaAnnot) <- rnaAnnot$sampleId

# ±250 kb window around each gene's TSS
rnagr <- as(rnagr, "GRanges")
rnagr <- flank(rnagr, width = FLANK, start = TRUE, both = TRUE)
rnagr <- trim(rnagr, use.names = TRUE)

rnacm_ <- takeResLimma_(rnacm, rnaAnnot, "donors", "aeir_group")


################
# 2. ATAC: differential peaks, log2 CPM, donor-corrected
################

data_ <- loadOBJECT_("long")

# Drop donors with < 2 samples
annot <- getSampleAnnot(data_)
singleDonors <- names(which(table(annot$donor) < 2))
atac_ <- removeSamples(data_, annot$donor %in% singleDonors)
atac_ <- addSampleAnnotCol(atac_, "donors", droplevels(getSampleAnnot(atac_)$donor))
summary(getSampleAnnot(atac_)$donors)

data_cpm <- transformCounts(atac_, method = "CPM", region = ".peaks.cons")
data_cpm <- transformCounts(data_cpm, method = "log2", region = ".peaks.cons")
atacgr <- getCoord(data_cpm, ".peaks.cons")
ataccm <- getCounts(data_cpm, ".peaks.cons")
rownames(ataccm) <- atacgr$name

ataccm_ <- takeResLimma_(ataccm, getSampleAnnot(atac_), "donors", "EirGrp")

atacDiffAcc <- read.csv("./data/LONG_vd2.EirGrp-AltiusSubset.diffacc.csv", row.names = 1)
diffPeaks <- rownames(atacDiffAcc)[!atacDiffAcc$diffaccess %in% "NO"]
table(atacDiffAcc$diffaccess)

ataccm_ <- ataccm_[rownames(ataccm_) %in% diffPeaks, ]
atacgr <- atacgr[atacgr$name %in% diffPeaks]


################
# 3. Overlap: differential peaks within a DEG window, shared samples
################

atacgr1 <- subsetByOverlaps(atacgr, rnagr, type = "within")
rnagr1 <- subsetByOverlaps(rnagr, atacgr1, type = "any")

sharedSamples <- intersect(colnames(rnacm_), colnames(ataccm_))
ataccm1 <- ataccm_[rownames(ataccm_) %in% atacgr1$name, sharedSamples]
rnacm1 <- rnacm_[rownames(rnacm_) %in% names(rnagr1), sharedSamples]

stopifnot(identical(colnames(ataccm1), colnames(rnacm1)))
dim(ataccm1)  # 7711
dim(rnacm1)   # 7798


################
# 4. Pearson correlation per peak-gene pair
################

calc_row_cor <- function(x, y) {
  x_c <- x - rowMeans(x)
  y_c <- y - rowMeans(y)
  rowSums(x_c * y_c) / (sqrt(rowSums(x_c^2) * rowSums(y_c^2)) + 1e-16)
}

hits <- findOverlaps(rnagr1, atacgr1)
rna_idx  <- queryHits(hits)
atac_idx <- subjectHits(hits)

rna_ids_vec  <- rnagr1$genes_version[rna_idx]
atac_ids_vec <- atacgr1$name[atac_idx]

rna_mat  <- as.matrix(rnacm1[rna_ids_vec, , drop = FALSE])
atac_mat <- as.matrix(ataccm1[atac_ids_vec, , drop = FALSE])

correlation_ <- data.table(
  pearson           = calc_row_cor(rna_mat, atac_mat),
  gene              = rna_ids_vec,
  seqnames          = as.character(seqnames(rnagr1)[rna_idx]),
  genename          = rnagr1$name[rna_idx],
  rna.start         = start(rnagr1)[rna_idx],
  rna.end           = end(rnagr1)[rna_idx],
  peakname          = atac_ids_vec,
  atac.start        = start(atacgr1)[atac_idx],
  atac.end          = end(atacgr1)[atac_idx],
  nearest_gene      = atacgr1$.nearest_gene_name[atac_idx],
  nearest_gene_dist = atacgr1$.nearest_gene_dist[atac_idx]
)

fwrite(correlation_, file.path(out_dir, "p2g-Correlation.csv"))
print("CORRELATION COMPLETE")


# nolint end

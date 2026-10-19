# nolint start
####################################################################################
#' Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
#' Script: activity.R
#' Author: MIJM
#' Last Update: 27.08.26
#' Description: computing activity scores for ATACseq data and plotting with RNAseq data
#' ATACseq samples are donor-paired and adjusted (with limma::removeBatchEffects) responder samples identified from chromVAR differentials by aEIR groups. Mean difference between the scores are calculated and plotted with RNAseq data as scatterplots.
# responders only; limma correction applied

####################################################################################



library(rtracklayer)


GOI_Master <- c("IL7R", "LYN", "TNF", "CD8A", "IFNG", "GZMK", "PLCG2", "KIR2DL4")


out_dir <- file.path(getwd(), "./results/integrative")
dir.create(out_dir, showWarnings = TRUE) 


aggrMethod <- "mean" # "sum", "mean", "median"
GR <- c("tiling200bp", ".peaks.cons")


atacDiffAcc <- read.csv("./data/LONG_data_.EirGrp-AltiusSubset.diffacc.csv", row.names=1)


linkage <- if(exists("correlation_")) {
    print("Correlation present. Using loaded data")
    correlation_
} else {
    print("Loading P2G linkages")
    read.csv(paste0(plotDir,"p2g-Correlation.csv"), row.names=1)  
}



data_cpm <- transformCounts(data_, method="CPM", regionTypes= GR[1])
data_cpm <- transformCounts(data_cpm, method="log2", region = GR[1]) 
annot <- getSampleAnnot(data_) 
ataccm <- getCounts(data_cpm, GR[1])
identical(rownames(annot), colnames(ataccm))

data_batch <- data_cpm

data_batch@counts$tiling200bp <-  takeResLimma_(ataccm, annot, "donors", "EirGrp")
head(data_batch@counts$tiling200bp)

print(paste("aggregation method:", aggrMethod))
merged <- mergeSamples(data_batch, "EirGrp", countAggrFun= aggrMethod) # "sum", "mean"

getCounts(merged, GR[1]) %>% head()

exportCountTracks(merged, GR[1], out_dir, formats = "bed")
exportCountTracks(data_batch, GR[1], out_dir, formats = "bed")

coord_ <- getCoord(data_, GR[2])
tmp <- mcols(coord_)$".nearest_gene_name"
names(coord_) <- mcols(coord_)$name
mcols(coord_) <- NULL
mcols(coord_)$nearest_gene <- tmp

# export(coord_, file.path(out_dir, "consensus.bed"))



tmp <- atacDiffAcc %>% filter(!diffaccess %in% "NO") %>% rownames(.)
peaks_df <- as.data.frame(coord_)
peaks_df <- peaks_df %>%
  mutate(
    rgb_color = ifelse(rownames(.) %in% tmp, "0,255,0", "128,128,128")
  )

# 0-indexing
bed9 <- data.frame(
  chrom = peaks_df$seqnames,
  chromStart = peaks_df$start - 1, 
  chromEnd = peaks_df$end,
  name = rownames(peaks_df),
  score = 0,                   
  strand = ".",               
  thickStart = peaks_df$start - 1, 
  thickEnd = peaks_df$end,        
  itemRgb = peaks_df$rgb_color 
)

write.table(bed9, file.path(out_dir, "consensus_colored.bed"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


source("/icbb/projects/mmaran/SCRIPTS/CCB/p2g/2025/gencode.r")
gencode <- getAnnotGrl.gencode("gencode.v46")
names(gencode)
flat_gencode <- unlist(gencode, use.names = FALSE)


toKeepBiotype <- c("protein_coding",
      "antisense", 
      "bidirectional_promoter_lncRNA", 
      "lincRNA", 
      "lncRNA",
      "macro_lncRNA", 
      "sense_intronic", 
      "sense_overlapping",
      "TR_C_gene",
      "TR_D_gene" ,
      "TR_J_gene",                    
      "TR_V_gene"
    )


gr <- flat_gencode[mcols(flat_gencode)$gene_type %in% toKeepBiotype]
table(mcols(gr)$gene_type)
names(gr) <- mcols(gr)$gene_name
rtracklayer::export(gr, file.path(out_dir, "ref.gtf"), format = "GTF")


gr <- gencode[[1]]
gr <- gr[mcols(gr)$gene_type %in% toKeepBiotype]
table(mcols(gr)$gene_type)
names(gr) <- mcols(gr)$gene_name


# Looping GOI

for (goi in GOI_Master) {
  message(sprintf("Processing gene: %s", goi))
  
  # 1. Base linkages for the current gene
  linkage_ <- linkage %>% 
    as.data.frame() %>% 
    filter(genename == goi)
  
  diffAcc_ <- atacDiffAcc %>% 
    filter(diffaccess != "NO", rownames(.) %in% linkage_$peakname)
  
  # 2. Process vlines
  vline <- coord_[names(coord_) %in% rownames(diffAcc_)]
  
  if(length(vline) > 0){
    vline <- c(
      GRanges(seqnames(vline), IRanges(start(vline), start(vline) + 1)),
      GRanges(seqnames(vline), IRanges(end(vline), end(vline) + 1))
    )
    rtracklayer::export(vline, file.path(out_dir, paste0(goi, "_vlines.bed")), format = "BED")
  } else {
    message(sprintf("No vlines to export for %s.", goi))
  }
  
  # 3. Gene Body and TSS coordinates
  goi_gr <- gr[names(gr) == goi]
  if(length(goi_gr) == 0) {
      message(sprintf("Gene %s not found in reference. Skipping.", goi))
      next
  }
  
  goi_tss <- promoters(goi_gr, upstream = 0, downstream = 1) 
  
  # 4.  Linkage Processing
  shared_linkage <- linkage %>% 
    filter(pearson > 0.3, genename == goi) %>%
    arrange(pearson) %>%
    mutate(
      chrom1 = seqnames,
      chrom2 = seqnames,
      atac.mid.start = round((atac.start + atac.end) / 2), 
      atac.mid.end = atac.mid.start + 1
    )
  
  # 5a. Create P2G Links (Using TSS coordinates)
  p2g_links <- shared_linkage %>%
    mutate(
      gene.start = start(goi_tss)[1],
      gene.end = end(goi_tss)[1]
    ) %>%
    dplyr::select(chrom1, atac.mid.start, atac.mid.end, chrom2, gene.start, gene.end)
    
  write.table(p2g_links, file.path(out_dir, paste0(goi, "_p2g.links")), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # 5b. Create View Links (Using full Gene Body coordinates)
  view_links <- shared_linkage %>%
    mutate(
      gene.start = start(goi_gr)[1],
      gene.end = end(goi_gr)[1]
    ) %>%
    dplyr::select(chrom1, atac.mid.start, atac.mid.end, chrom2, gene.start, gene.end)
    
  write.table(view_links, file.path(out_dir, paste0(goi, "_view.links")), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}


goi_gr <- gr[names(gr) %in% GOI_Master]
  if(length(goi_gr) == 0) {
      message(sprintf("Gene %s not found in reference. Skipping.", goi))
      next
  }
goi_tss <- promoters(goi_gr, upstream = 0, downstream = 1) 
mcols(goi_tss) <- NULL
# mcols(goi_tss)$gene_name <- names(goi_tss)

write.table(as.data.frame(goi_tss), file.path(out_dir, "tss.csv"), row.names = TRUE, col.names = FALSE, quote = FALSE, sep = ",")

# nolint end



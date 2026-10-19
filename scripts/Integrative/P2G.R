# Created 2025
# Author: Midhuna I JM
# Description: Pygenome
# Last Modified: 12Mar2026
# Load 2 tabs; one with conda pygenometracks and the other with atac. Inital R commands in R; rest in pygenometracks
# responders only; limma correction applied



```{r= remove donor specific changes}

takeResLimma <- function(cm, annot, corCol, compCol="EirGrp") {

    annot <- as.data.frame(annot)

    stopifnot(all(colnames(cm) == rownames(annot)))

    f <- as.formula(paste("~", compCol))
    design <- model.matrix(f, data = annot)
    res_ <- limma::removeBatchEffect(as.matrix(cm),
                              batch      = annot[[corCol]],
                              design     = design)

    
    return(res_)
}

```

```{r=}
library(rtracklayer)
library(dplyr)

# GOI_Master <- c("CD44", "GZMK", "IFNG", "IL7R", "LYN", "PLCG2", "NCR1", "TNF", "CD8A", "HAVCR2", "KLRC4", "KIR2DL1", "KIR2DL4", "KIR3DL1", "KIR3DL2", "CD244")
GOI_Master <- c("IL7R", "LYN", "TNF", "CD8A", "IFNG", "GZMK", "PLCG2", "KIR2DL4")


out_dir <- file.path(getwd(), "2026/p2g/pyGenome/2Sep26")
dir.create(out_dir, showWarnings = TRUE) 


aggrMethod <- "mean" # "sum", "mean", "median"
```

```{r=}

atacDiffAcc <- read.csv("/icbb/projects/mmaran/malaria/2025/Manuscript/LONG/LONG_vd2.EirGrp-AltiusSubset.diffacc.csv", row.names=1)
table(atacDiffAcc$diffaccess)
head(atacDiffAcc)

linkage <- if(exists("correlation_")) {
    print("Correlation present. Using loaded data")
    correlation_
} else {
    print("Loading P2G linkages")
    read.csv(paste0(plotDir,"p2g-Correlation.csv"), row.names=1)  
}

```

```{r=}

GR <- c("tiling200bp", ".peaks.cons")

data_cpm <- transformCounts(vd2, method="CPM", regionTypes= GR[1])
data_cpm <- transformCounts(data_cpm, method="log2", region = GR[1]) 
# data_batch <- transformCounts(data_cpm, method="batchCorrect", batch= getSampleAnnot(vd2)$donors, regionTypes= GR[1])


annot <- getSampleAnnot(vd2) 
ataccm <- getCounts(data_cpm, GR[1])
identical(rownames(annot), colnames(ataccm))

data_batch <- data_cpm

if (!exists("takeResLimma", mode = "function")) {
  stop("NO LIMMA FUNCTION, CHECK!")
  data_batch@counts$tiling200bp <-  takeResLimma(ataccm, annot, "donors", "EirGrp")
}


head(data_batch@counts$tiling200bp)

# tmp <- data_batch@counts$.peaks.cons
# tmp2 <- getCounts(data_batch, GR[2])
# identical(as.integer(tmp), as.integer(tmp2))


print(paste("aggregation method:", aggrMethod))
merged <- mergeSamples(data_batch, "EirGrp", countAggrFun= aggrMethod) # "sum", "mean"

getCounts(merged, GR[1]) %>% head()

exportCountTracks(merged, GR[1], out_dir, formats = "bed")
exportCountTracks(data_batch, GR[1], out_dir, formats = "bed")

coord_ <- getCoord(vd2, GR[2])
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
  score = 0,                   # Placeholder
  strand = ".",                # Placeholder 
  thickStart = peaks_df$start - 1, # Placeholder 
  thickEnd = peaks_df$end,         # Placeholder
  itemRgb = peaks_df$rgb_color 
)

write.table(bed9, file.path(out_dir, "consensus_colored.bed"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

```


```{r=}

# gr <- readRDS(file.path(getwd(), "2026/p2g/pyGenome/processFiles/ref.rds"))
# table(mcols(gr)$gene_type)
# export(gr, file.path(getwd(), "2026/p2g/pyGenome/ref.bed"))


# features_to_export <- c("transcript", "exon", "CDS", "UTR")
# gencode_export <- sort(
#   do.call(c, unname(as.list(gencode[features_to_export])))
# )

# head(mcols(gencode_export)[, c("gene_id", "transcript_id", "gene_name", "type")])

# # Export
# rtracklayer::export(
#   gencode_export,
#   file.path(out_dir, "reference.gtf"),
#   format = "GTF"
# )
# message("GENCODE GTF exported: ", length(gencode_export), " features.")


# source("/icbb/projects/mmaran/DEEP/SolveRD/2025/FUNCTIONS/GENCODE.R")
# gencode <- getAnnotGrl.gencode_("gencode.v27")

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
# goi_gencode <- flat_gencode[flat_gencode$gene_name %in% GOI_Master]
rtracklayer::export(gr, file.path(out_dir, "ref.gtf"), format = "GTF")


gr <- gencode[[1]]
gr <- gr[mcols(gr)$gene_type %in% toKeepBiotype]
table(mcols(gr)$gene_type)
names(gr) <- mcols(gr)$gene_name

```

```{r= Looping GOI}

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
  
  # 3. Get both Gene Body and TSS coordinates
  goi_gr <- gr[names(gr) == goi]
  if(length(goi_gr) == 0) {
      message(sprintf("Gene %s not found in reference. Skipping.", goi))
      next
  }
  
  goi_tss <- promoters(goi_gr, upstream = 0, downstream = 1) 
  
  # 4. Shared Linkage Processing
  # Filter and do the math once per gene to save time
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

# write.table(as.data.frame(goi_tss), file.path(out_dir, paste0("tss.csv")), row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(goi_tss), file.path(out_dir, "tss.csv"), row.names = TRUE, col.names = FALSE, quote = FALSE, sep = ",")

```



```{r= BASH; pygenome}

#!/bin/bash

# Convert into bedGraph
for f in *.bed; do
  if [ -f "$f" ]; then
    awk '{printf "%s\t%d\t%d\t%2.3f\n", $1,$2,$3,$5}' "$f" > "${f%.bed}.bedgraph"
  fi
done


for p2g_file in *_p2g.links; do
  
  filename=$(basename "$p2g_file")
  gene="${filename%_p2g.links}"
  view_file="${gene}_view.links"

  echo "========================================"
  echo "Generating PyGenomeTracks for: $gene"

  if [ ! -s "$p2g_file" ]; then
    echo "No linkages found for $gene. Skipping."
    continue
  fi

  chr=$(head -n 1 "$view_file" | awk '{print $1}')


  # start=$(awk -F',' -v g="$gene" '$1 == g || $2 == g { s=$3; print (s - 250000 < 0 ? 0 : s - 250000) }' tss.csv)
  # end=$(awk -F',' -v g="$gene" '$1 == g || $2 == g { print $3 + 250000 }' tss.csv)

  # start=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$p2g_file" | sort -n | head -n 1 | awk '{print $1-10000}')
  # end=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$p2g_file" | sort -n | tail -n 1 | awk '{print $1+10000}')

  start=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$view_file" | sort -n | head -n 1 | awk '{print $1-10000}')
  end=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$view_file" | sort -n | tail -n 1 | awk '{print $1+10000}')
  
  # : "${start:=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$view_file" | sort -n | head -n 1 | awk '{print $1-10000}')}"
  # : "${end:=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$view_file" | sort -n | tail -n 1 | awk '{print $1+10000}')}"

  echo "Plotting region: ${chr}:${start}-${end}"

  max_y=$(awk -v c="$chr" -v s="$start" -v e="$end" '
    $1 == c && $3 >= s && $2 <= e {
      if ($4 > max) max = $4
    }
    END {
      if (max == "") max = 1
      print int(max * 1.1) + 1
    }
  ' *.bedgraph)

  echo "Calculated shared max y-axis: $max_y"

  sed -e "s/MAX_Y/$max_y/g" \
      -e "s|GENE_P2G_LINKS|$p2g_file|g" \
      -e "s|GENE_VIEW_LINKS|$view_file|g" \
      -e "s|GENE_VLINES|${gene}_vlines.bed|g" ./MAIN.ini > "MAIN_${gene}.ini"

  pyGenomeTracks --tracks "MAIN_${gene}.ini" --region "${chr}:${start}-${end}" -o "./plots/mean_full_${gene}.pdf"

unset start end

done


```


```{r= zoomed in}

ls *_p2g.links


p2g_file="PLCG2_p2g.links"
start=81770000
end=81800000

filename=$(basename "$p2g_file")
gene="${filename%_p2g.links}"
view_file="${gene}_view.links"

echo "========================================"
echo "Generating PyGenomeTracks for: $gene"

if [ ! -s "$p2g_file" ]; then
  echo "No linkages found for $gene. Skipping."
  continue
fi

chr=$(head -n 1 "$view_file" | awk '{print $1}')


: "${start:=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$view_file" | sort -n | head -n 1 | awk '{print $1-10000}')}"
: "${end:=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$view_file" | sort -n | tail -n 1 | awk '{print $1+10000}')}"

echo "Plotting region: ${chr}:${start}-${end}"

max_y=$(awk -v c="$chr" -v s="$start" -v e="$end" '
  $1 == c && $3 >= s && $2 <= e {
    if ($4 > max) max = $4
  }
  END {
    if (max == "") max = 1
    print int(max * 1.1) + 1
  }
' *.bedgraph)

echo "Calculated shared max y-axis: $max_y"

sed -e "s/MAX_Y/$max_y/g" \
    -e "s|GENE_P2G_LINKS|$p2g_file|g" \
    -e "s|GENE_VIEW_LINKS|$view_file|g" \
    -e "s|GENE_VLINES|${gene}_vlines.bed|g" ./MAIN.ini > "MAIN_${gene}.ini"

pyGenomeTracks --tracks "MAIN_${gene}.ini" --region "${chr}:${start}-${end}" -o "./plots/mean_view_reduced_${gene}.pdf"

unset start end



```
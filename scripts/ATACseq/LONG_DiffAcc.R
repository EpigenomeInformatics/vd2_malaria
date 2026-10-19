# nolint start
####################################################################################
# Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
# Script: MAIN.R
# Author: MIJM
# Last Update: 27.08.26
# Description: Full analysis for Longitudinal ATACseq data
####################################################################################

source("./scripts/ATACseq/utils.R")

################
# Load object
################

OBJECT <- "long" #c("long", "regional", )
data_ <- loadOBJECT_(OBJECT)



coord_ <- getCoord(data_, ".peaks.cons")
head(coord_)
tssGr <- ChrAccRAnnotationHg38::getGeneAnnotation(anno="gencode_coding", type="tssGr")
peakTssDist <- GenomicRanges::distanceToNearest(coord_, tssGr, ignore.strand=TRUE)
mcols(coord_)$nearestGeneDist <- S4Vectors::mcols(peakTssDist)[,"distance"]
mcols(coord_)$nearestGene <- S4Vectors::elementMetadata(tssGr)[S4Vectors::subjectHits(peakTssDist),"gene_name"]
colnames(data.frame(coord_))



compList <- c("y2013v14", "y2013v16", "YearGrp", "EirGrp", "priormal12Grp", "anyparasitemia", "LAMP", "priormal12B", "eirB" )

DAG <- list()

for (GROUP in compList) {
  print(paste("computing:", GROUP))
  
  NAME <- GROUP
  dataRun <- data_  

  if(GROUP %in% "y2013v14"){
      dataRun <- addSampleAnnotCol(dataRun, GROUP, getSampleAnnot(dataRun)$year)
      toDel <- !getSampleAnnot(dataRun)[[GROUP]] %in% c(2013, 2014)
      print(summary(toDel))
      dataRun <- removeSamples(dataRun, toDel) 
  }else if (GROUP %in% "y2013v16") {
    dataRun <- addSampleAnnotCol(dataRun, GROUP, getSampleAnnot(dataRun)$year)
    toDel <- !getSampleAnnot(dataRun)[[GROUP]] %in% c(2013, 2016)
    print(summary(toDel))
    dataRun <- removeSamples(dataRun, toDel) 
  }

  annot_ <- getSampleAnnot(dataRun)
  
  annot_ <- droplevels(annot_[!annot_[[GROUP]] == "", ])
  print(table(annot_[[GROUP]]))
  
  group_levels <- levels(factor(annot_[[GROUP]]))
  if (length(group_levels) < 2) {
    warning(paste("Skipping", GROUP, "- not enough levels for comparison."))
    next
  }
  
  grp1 <- group_levels[1]
  grp2 <- group_levels[2]
  

  paired_sample_ids <- annot_ %>%
    group_by(donor) %>%
    filter(any(.data[[GROUP]] == grp1) & any(.data[[GROUP]] == grp2)) %>%
    pull(sampleId)
  
  toDel <- !getSampleAnnot(dataRun)$sampleId %in% paired_sample_ids
  print(summary(toDel))
  
  dataRun <- removeSamples(dataRun, toDel)
  
  dataRun <- addSampleAnnotCol(dataRun, "donors", droplevels(getSampleAnnot(dataRun)$donor))
  
  cat("\n--- samples per group:", GROUP, "---\n")
  print(table(getSampleAnnot(dataRun)[[GROUP]]))
  cat("\n--- Donors:", GROUP, "---\n")
  print(table(getSampleAnnot(dataRun)$donors))
  
  dataRun <- addSampleAnnotCol(dataRun, GROUP, droplevels(as.factor(getSampleAnnot(dataRun)[[GROUP]])))

  
  diffAcc <- getDiffAcc(dataRun, ".peaks.cons", GROUP, grp1Name = grp1, grp2Name = grp2, adjustCols = c("donors"))
  
  rownames(diffAcc) <- sub("\\..*?_peak", "_peak", rownames(diffAcc))

  diffAcc <- diffAcc %>% 
    mutate(
        Significant = padj <= 0.05 & abs(log2FoldChange) > 0.6,
        diffaccess = case_when(
            log2FoldChange > 0.6 & padj <= 0.05 ~ grp1,
            log2FoldChange < -0.6 & padj <= 0.05 ~ grp2,
            TRUE ~ "NO"),
        nearest_genes = mcols(coord_)$nearestGene,
        nearestGeneDist = as.numeric(as.character(mcols(coord_)$nearestGeneDist)),
        relation = case_when(
            abs(nearestGeneDist) <= 2000 ~ "proximal",
            abs(nearestGeneDist) > 2000 & abs(nearestGeneDist) <= 250000 ~ "distal",
            abs(nearestGeneDist) > 250000 ~ "NO",
            TRUE ~ NA_character_
        ) )

    cat("\n--- Summary for:", GROUP, "---\n")
    print(table(diffAcc$diffaccess))

    cat("\n--- relation to nearest genes:", GROUP, "---\n")
    print(table(diffAcc$relation))

    DAG[[NAME]] <- diffAcc
    rm(diffAcc)
}

names(DAG)



GOI <- read.csv("./data/malaria_GOI.txt", header=F) %>% pull(V1)
label_ <- diffAcc %>% filter(!diffaccess %in% c("NO", "Not differential")) %>% 
              filter(nearest_genes %in% GOI | nearest_genes %in% c("LYN", "HAVCR1", "JUN", "FOS", "BATF", "IL7R", "EOMES")) %>%
              filter(!relation %in% "NO")

table(label_$relation)
label_$relationColor <- ifelse(label_$relation %in% "proximal", "darkgreen", "black")

label_ <- label_ %>% group_by(diffaccess, nearest_genes, relation) %>% slice_max(order_by = abs(log2FoldChange), n = 1) %>% ungroup()


diffAcc$category <- diffAcc$diffaccess
diffAcc$category[diffAcc$category %in% c("NO", "Not differential")] <- "Not differential"

# Counts
high_n <- sum(diffAcc$diffaccess == "high")
low_n  <- sum(diffAcc$diffaccess == "low")
no_n   <- sum(diffAcc$diffaccess == "NO")


# jpeg("image.jpg", width=16, height=14, units="in", res=600)
pdf(paste0(plotDir, "vd2.diffAcc.Responder_EirGrp.relation.pdf"), width=16,height=12 )
ggplot(diffAcc, aes(x = log2BaseMean, y = log2FoldChange)) +
geom_bin2d(data = subset(diffAcc, category == "Not differential"), bins = 500) +
scale_fill_gradient(low = "gray90", high = "gray10", guide = "none") +
geom_point(data = subset(diffAcc, category != "Not differential"), aes(color = category), size = 2) +
scale_color_manual(
    values = c("high" = "darkred", "low" = "darkblue", "Not differential" = "grey"),
    labels = c(paste0("high = ", high_n),
               paste0("low = ", low_n),
               paste0("NO = ", no_n))) +

geom_hline(yintercept = 0, color = "black", lwd = 1) +
geom_label_repel(data = label_,
                   aes(label = nearest_genes),
                   fill = "white",
                   color = label_$relationColor, #"black",
                   size = 7,
                   fontface = "bold",
                   box.padding = 0.5,
                   point.padding = 0.3,
                   segment.color = "black",
                   segment.size = 0.3,
                   max.overlaps = Inf) +

labs(x = "Mean Accessibility",
       y = "Log2FC",
       title = "Long: EirGrp- differentials, paired samples",
       color = "Differential access") +

guides(color = guide_legend(override.aes = list(size = 3))) +
theme_classic(base_size = 18)
dev.off()




# nolint end

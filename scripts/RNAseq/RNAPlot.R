# nolint start

####################################################################################
# Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
# Script: RNAPlot.R
# Author: MIJM
# Last Update: 24Sept26
#' Description: RNAseq Analysis
# Combined invivo and invitro
# Common Normalization before biotype
# Same background set for all analyses
####################################################################################

source("./scripts/RNAseq/RNAseq.R")
source("./scripts/RNAseq/RNA_utils.R")

# plotting

############
# 1. MA plot
############


library(ggrepel)

highN <- sum(degDf$diffExp == lev[2])
lowN  <- sum(degDf$diffExp == lev[1])
noN   <- sum(degDf$diffExp == "NO")
degDf$diffExp[degDf$diffExp == "NO"] <- "Not differential"

pdf(paste0(plotDir, "RNA-vd2.MA.", groupCol, "label.pdf"), width = 12, height = 10)
ggplot(degDf, aes(x = logCPM, y = logFC)) +
  geom_bin2d(data = subset(degDf, diffExp == "Not differential"), bins = 500) +
  scale_fill_gradient(low = "gray90", high = "gray10", guide = "none") +
  geom_point(data = subset(degDf, diffExp != "Not differential"), aes(color = diffExp), size = 2) +
  scale_color_manual(
    values = setNames(c("darkred", "darkblue", "gray60"), c(lev[2], lev[1], "Not differential")),
    breaks = c(lev[2], lev[1], "Not differential"),
    labels = c(
      paste0(lev[2], " = ", highN),
      paste0(lev[1], " = ", lowN),
      paste0("Not Differential = ", noN)
    ),
    drop = FALSE
  ) +
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_label_repel(
    data = geneLabels,
    aes(label = name),
    fill = "white",
    color = "black",
    fontface = "bold",
    size = 6,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = 20,
    arrow = arrow(length = unit(0.015, "npc")), 
    min.segment.length = 0
  ) +
  labs(
    x = "Mean Expression (log2-CPM)",
    y = "LogFC",
    title = paste("RNA:allsamples:", groupCol),
    color = "DEG"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic(base_size = 18)
dev.off()



############
# 2. Heatmap
############

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

subDeg <- degDf[!degDf$diffExp %in% c("NO", "Not differential"), ]
logCounts <- cpm(dge, log = TRUE, prior.count = 2)

heatmapDesign <- model.matrix(~ group, data = dge$samples)
logCounts <- removeBatchEffect(
  logCounts,
  batch = dge$samples$donors, 
  design = heatmapDesign
)

heatMat <- logCounts[rownames(subDeg), ]

stopifnot(identical(rownames(subDeg), rownames(heatMat)))
rownames(heatMat) <- subDeg$name

annotColDf <- data.frame(
  "year" = tempAnnot$year, 
  "log10aEIR" = log10(tempAnnot$aeir + 1), 
  "aEIRGrp" = tempAnnot$aeir_group 
)
rownames(annotColDf) <- tempAnnot$sampleId

heatAnnoColors <- list(
  log10aEIR = colorRamp2(breaks = c(min(annotColDf$log10aEIR, na.rm = TRUE), max(annotColDf$log10aEIR, na.rm = TRUE)), colors = c("green", "red")),
  year = colorRamp2(breaks = c(2013, 2016), colors = c("brown", "gold")),
  aEIRGrp = c("low" = "darkblue", "high" = "darkred")
)

# Order annotations by aEIR and re-align matrix
annotColDf <- annotColDf[order(-annotColDf$log10aEIR), ]
heatMat <- heatMat[, rownames(annotColDf)]

stopifnot(identical(rownames(annotColDf), colnames(heatMat)))

# Within-group Kmeans clustering 
group1Samples <- annotColDf %>% subset(aEIRGrp == "high") %>% rownames()
group2Samples <- annotColDf %>% subset(aEIRGrp == "low") %>% rownames()

highExp <- as.data.frame(heatMat[, group1Samples])
lowExp <- as.data.frame(heatMat[, group2Samples])

highClust <- as.data.frame(kmeans(t(highExp), centers = 2)$cluster) %>% setNames("cluster") 
lowClust <- as.data.frame(kmeans(t(lowExp), centers = 2)$cluster) %>% setNames("cluster") 

highClust <- highClust %>% arrange(-cluster)
lowClust <- lowClust %>% arrange(cluster)

combClust <- rbind(highClust, lowClust)

annotColDf <- annotColDf[match(rownames(combClust), rownames(annotColDf)), ]
heatMat <- heatMat[, rownames(annotColDf)]

# Plotting heatmap
numRows <- nrow(heatMat)
# calcHeight <- 20 + (numRows * 0.15) 

topGenesLabels <- subDeg %>% 
  group_by(diffExp) %>% 
  mutate(orderVal = logFC * -log10(PValue)) %>% 
  group_by(diffExp) %>% 
  arrange(-orderVal) %>% 
  slice_head(n = 20) %>% 
  pull(name)

# Inherited genesOfInterest from RNA_utils.R instead of re-typing the array
heatmapLabels <- unique(c(genesOfInterest, topGenesLabels)) 

labelIndices <- which(rownames(heatMat) %in% heatmapLabels)
labelsToShow <- rownames(heatMat)[labelIndices]

rowAnno <- rowAnnotation(
  genes = anno_mark(at = labelIndices, labels = labelsToShow, labels_gp = gpar(fontsize = 15))
)

colFun <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

pdf(paste0(plotDir, "RNA-heatmap.vd2.", groupCol, ".pdf"), width = 10, height = 12) 
Heatmap(
  t(scale(t(heatMat))), 
  top_annotation = HeatmapAnnotation(df = annotColDf, col = heatAnnoColors), 
  cluster_rows = TRUE, 
  cluster_columns = FALSE, 
  right_annotation = rowAnno, 
  show_row_names = FALSE,     
  show_column_names = FALSE,
  col = colFun
)
dev.off()

#nolint end
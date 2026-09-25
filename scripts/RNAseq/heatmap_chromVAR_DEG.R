# nolint start
####################################################################################
#' Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
#' Script: heatmap_chromVAR_DEG.R
#' Author: MIJM
#' Last Update: 27.08.26
#' Description: Plotting ATACseq chromVAR heatmap over DEGs

####################################################################################

suppressPackageStartupMessages(c(
library(circlize),
library(ComplexHeatmap)
))

# order from chromVAR heatmap
order_ <- c("VD2_2013_d3278_52", "VD2_2013_d3340_80", "VD2_2014_d3381_160",
"VD2_2013_d3401_208", "VD2_2013_d3280_44", "VD2_2013_d3369_104",
"VD2_2014_d3369_108", "VD2_2014_d3398_200", "VD2_2014_d3340_84", 
"VD2_2013_d3381_156", "VD2_2014_d3401_212", "VD2_2013_d3375_144",
"VD2_2013_d3203_224", "VD2_2014_d3203_228", "VD2_2014_d3278_56", 
"VD2_2014_d3384_188", "VD2_2013_d3398_196", "VD2_2013_d3257_40", 
"VD2_2013_d3326_68", "VD2_2014_d3375_148", "VD2_2016_d3375_152",
"VD2_2016_d3203_232", "VD2_2013_d3384_184", "VD2_2015_d3398_204",
"VD2_2015_d3257_48", "VD2_2014_d3326_72", "VD2_2016_d3340_88", 
"VD2_2016_d3381_164", "VD2_2016_d3401_216", "VD2_2016_d3278_60", 
"VD2_2016_d3384_192", "VD2_2016_d3280_168", "VD2_2016_d3369_112")

responderSamples <- c(
"VD2_2016_d3340_88" ,
"VD2_2016_d3381_164",
"VD2_2016_d3401_216",
"VD2_2016_d3278_60" ,
"VD2_2016_d3384_192",
"VD2_2016_d3280_168",
"VD2_2016_d3369_112",

"VD2_2013_d3278_52" ,
"VD2_2013_d3340_80" ,
"VD2_2014_d3381_160",
"VD2_2013_d3401_208",
"VD2_2013_d3280_44" ,
"VD2_2013_d3369_104",
"VD2_2014_d3369_108",
"VD2_2014_d3398_200"
)


devZ <- readRDS("./data/RNA/devZmotifFamilies.rds")
motifAnnotDf <- read.csv("./data/RNA/annotcol4motifs.csv", row.names = 1)

devZ <- devZ[, order_]
if (!all(order_ == colnames(devZ))) stop("Columns must be identical")

logCounts <- cpm(dge, log = TRUE, prior.count = 2)
heatmapDesign <- model.matrix(~ group, data = dge$samples)
logCounts <- removeBatchEffect(
  logCounts,
  batch  = dge$samples$donors,
  design = heatmapDesign
)

subDeg <- degDf[!degDf$diffExp %in% c("NO", "Not differential"), ]
logCounts <- logCounts[rownames(subDeg), ]

subDeg <- subDeg %>% arrange(match(rownames(subDeg), rownames(logCounts)))
stopifnot(identical(rownames(subDeg), rownames(logCounts)))
rownames(logCounts) <- subDeg$name

degAligned <- matrix(
  NA, 
  nrow = nrow(logCounts), 
  ncol = ncol(devZ), 
  dimnames = list(rownames(logCounts), colnames(devZ))
)

commonCols <- intersect(colnames(devZ), colnames(logCounts))
degAligned[, commonCols] <- logCounts[, commonCols]

if (!all(colnames(degAligned) == colnames(devZ))) stop("Columns must be identical")
stopifnot(identical(colnames(degAligned), colnames(devZ))) 

motifAnnotDf <- motifAnnotDf %>% select(year, log10_aeir, EirGrp)

motifAnnotDf$responder <- ifelse(rownames(motifAnnotDf) %in% responderSamples, "Res", "NonRes")
motifAnnotDf$year <- as.numeric(motifAnnotDf$year)

responderColors <- c("goldenrod3", "olivedrab3") 
names(responderColors) <- unique(motifAnnotDf$responder)

heatAnnoColors <- list(
  log10_aeir = colorRamp2(breaks = c(min(motifAnnotDf$log10_aeir, na.rm = TRUE), max(motifAnnotDf$log10_aeir, na.rm = TRUE)), colors = c("green", "red")),
  year       = colorRamp2(breaks = c(2013, 2016), colors = c("brown", "gold")),
  EirGrp     = c("low" = "green", "high" = "red"),
  responder  = responderColors
)

colFunChromVar <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
colFunDeg      <- colorRamp2(c(-3, 0, 3), c("darkblue", "white", "darkred"))

# Plotting

motifAnnotDf <- motifAnnotDf[colnames(devZ), ]
stopifnot(
  identical(rownames(motifAnnotDf), colnames(devZ)),
  identical(rownames(motifAnnotDf), colnames(degAligned))
)

colAnno <- HeatmapAnnotation(
  df = motifAnnotDf,
  col = heatAnnoColors,
  annotation_legend_param = list(
    labels_gp = gpar(fontsize = 10) 
  )
)

ht1 <- Heatmap(
  devZ, 
  name = "chromVAR", 
  col = colFunChromVar,
  top_annotation = colAnno,
  cluster_columns = FALSE,
  show_row_names = FALSE,     
  show_column_names = TRUE
)

htTemp <- draw(ht1)
rowOrder1 <- row_order(htTemp)

ht2 <- Heatmap(
  t(scale(t(degAligned))), 
  name = "DEG", 
  col = colFunDeg,
  na_col = "white",       
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = FALSE,     
  show_column_names = FALSE
)

htList <- ht1 %v% ht2

pdf(paste0(plotDir, "RNA-vd2.chromVAR-DEG.", groupCol, ".pdf"), width = 12, height = 16)
draw(
  htList, 
  annotation_legend_side = "right", 
  heatmap_legend_side = "right"
)
dev.off()
# nolint end
# nolint start
####################################################################################
#' Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
#' Script: heatmap_chromVAR_DEG.R
#' Author: MIJM
#' Last Update: 27.08.26
#' Description: Plotting ATACseq chromVAR heatmap over DEGs

####################################################################################


library(circlize)
library(ComplexHeatmap)

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

colors_ <- c(
  "black","red","blue","green3","gold","purple","deeppink","gray50","cyan",
  "darkred","darkviolet","chocolate","orange",
  "navy","forestgreen","magenta","turquoise4","orangered","sienna",
  "dodgerblue3","firebrick","goldenrod3","darkorange3","brown","olivedrab3"
)


devZ <- readRDS("/icbb/projects/mmaran/malaria/devZmotifFamilies.rds")
annot_col <- read.csv("/icbb/projects/mmaran/malaria/annotcol4motifs.csv", row.names=1)
devZ <- devZ[,order_]
if(!all(order_ == colnames(devZ))) stop("Columns must be identical")

# aEIR DGE
logCounts <- cpm(DGE, log = TRUE, prior.count = 2)
design <- model.matrix(~ group, data = DGE$samples)
logCounts <- removeBatchEffect(logCounts,
                           batch      = DGE$samples$donor,
                           design     = design)

sub <- DEG[!DEG$diffExp %in% "NO",]
logCounts <- logCounts[rownames(sub),]


summary(rownames(logCounts) %in% rownames(sub))
summary(rownames(sub) %in% rownames(logCounts))
sub <- sub %>% arrange(match(rownames(sub), rownames(logCounts)))
identical(rownames(sub), rownames(logCounts))
rownames(logCounts) <- sub$name

degAligned <- matrix(NA, 
                       nrow = nrow(logCounts), 
                       ncol = ncol(devZ), 
                       dimnames = list(rownames(logCounts), colnames(devZ)))

common_cols <- intersect(colnames(devZ), colnames(logCounts))
degAligned[, common_cols] <- logCounts[, common_cols]

if(!all(colnames(degAligned) == colnames(devZ))) stop("Columns must be identical")
identical(colnames(degAligned), colnames(devZ)) 




annot_col <- annot_col %>% select(year, log10_aeir, EirGrp)
annot_col$responder <- ifelse(rownames(annot_col) %in% responderSamples, "Res", "NonRes")
annot_col$year <- as.numeric(annot_col$year)

responder_ <- c("goldenrod3", "olivedrab3") 
names(responder_) <- unique(annot_col$responder)

annoCol<-list(
    log10_aeir = colorRamp2(breaks = c(min(annot_col$log10_aeir), max(annot_col$log10_aeir)), colors = c("green", "red")),
    year = colorRamp2(breaks = c(2013, 2016), colors = c("brown", "gold")),
    EirGrp = c("low" = "green", "high" = "red"),
    responder = responder_
)



h1_ = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
h2_ = colorRamp2(c(-3, 0, 3), c("darkblue", "white", "darkred"))


# plotting

annot_col <- annot_col[colnames(devZ),]
identical(rownames(annot_col), colnames(devZ))
identical(rownames(annot_col), colnames(degAligned))

colAnn <- HeatmapAnnotation(
  df = annot_col,
  col = annoCol,
  annotation_legend_param = list(
    labels_gp = gpar(fontsize = 10) 
  ))


h1 <- Heatmap(devZ, 
            name = "chromVAR", 
            col = h1_,
            top_annotation = colAnn,
            cluster_columns = FALSE,
            show_row_names = FALSE,     
            show_column_names = TRUE)

ht = draw(h1)
ro_1 <- row_order(ht)

h2 <- Heatmap(t(scale(t(degAligned))), 
              name = "DEG", 
        col = h2_,
            na_col = "white",       
            cluster_columns = FALSE,
            cluster_rows = TRUE,
            show_row_names = FALSE,     
            show_column_names = FALSE)

ht_list = h1 %v% h2

jpeg("image.jpg", width=12, height=16, units="in", res=600)
# pdf(paste0(plotDir, "vd2.chromVAR-DEG.", GROUP, ".pdf"), width=12, height=16)
draw(ht_list, 
     annotation_legend_side = "right", 
     heatmap_legend_side = "right")
dev.off()

# nolint end
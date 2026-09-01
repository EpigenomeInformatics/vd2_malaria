# nolint start
####################################################################################
#' Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
#' Script: activity.R
#' Author: MIJM
#' Last Update: 27.08.26
#' Description: computing activity scores for ATACseq data and plotting with RNAseq data
#' ATACseq samples are donor-paired and adjusted (with limma::removeBatchEffects) responder samples identified from chromVAR differentials by aEIR groups. Mean difference between the scores are calculated and plotted with RNAseq data as scatterplots.

####################################################################################

source("/icbb/projects/mmaran/SCRIPTS/CCB/malariaScripts/utils.R")

# load data:
rnaPath <- "/icbb/projects/mmaran/malaria/2025/Manuscript/RNA/COMBINED/noSexChr"
rnaDEG <- read.csv(file.path(rnaPath, "RNA_vd2.aeir_group.csv"), row.names = 1)
rnaAnnot <- readRDS(file.path(rnaPath, "rnaAnnot.rds")) 
diffAcc <- read.csv("/icbb/projects/mmaran/malaria/2025/Manuscript/LONG/LONG_vd2.EirGrp-AltiusSubset.diffacc.csv", row.names =1) 


OBJ <- "vd2"
RT <- ".peaks.cons"
normCounts <- "weighted"  # "weighted", "total", "none"

GRP <- "EirGrp"
GRP1 <- "high"
GRP2 <- "low"
CUTOFF <- 0.05


# Taking only TPs of responders (13 samples)

OBJECT <- "long" #c("long", "regional", )

data_ <- loadOBJECT_(OBJECT)
data_ <- removeSamples_(data_, "responders") #TODO:fix this



# calculate gene activities using RBF

RBF <- getRBFGeneActivities(get(OBJ), 
                RT, 
                tssGr= NULL,
                maxDist=100000L,
                sigma=10000,
                minWeight=0.25,
                binarize=FALSE,
                useGeneBody=TRUE,
                normCounts = normCounts)  

activity <- assays(RBF)$gA

# CPM and limma correction
activity <- activity * 1e6
annot <- getSampleAnnot(vd2)
stopifnot(identical(rownames(annot), colnames(activity)))
activity <- takeResLimma(activity, annot, "donors", "EirGrp")

actGR <- rowRanges(RBF)
mcols(actGR) <- mcols(actGR) %>% as.data.frame() %>% select(gene_id, gene_type, gene_name)
rownames(activity) <- mcols(actGR)$gene_id 


# prep data:
rownames(activity) <- gsub("\\..*$", "", rownames(activity))
rownames(rnaDEG) <- gsub("\\..*$", "", rownames(rnaDEG))

rnaDEG <- rnaDEG[rownames(rnaDEG) %in% rownames(activity),]
atac_ <- activity[rownames(activity) %in% rownames(rnaDEG),]
rnaDEG <- rnaDEG[rownames(atac_),]

mcols(actGR)$gene <- gsub("\\..*$", "", mcols(actGR)$gene_id)
actGR <- actGR[mcols(actGR)$gene %in% rownames(atac_)]
#checks
stopifnot(identical(rownames(atac_), rownames(rnaDEG)))
stopifnot(identical(rownames(atac_), mcols(actGR)$gene))
stopifnot(identical(nrow(atac_), nrow(rnaDEG)))


# Calculate delta(mean) for activity

group1_ <- getSampleAnnot(vd2)$sampleId[getSampleAnnot(vd2)[[GRP]] %in% GRP1] 
group2_ <- getSampleAnnot(vd2)$sampleId[getSampleAnnot(vd2)[[GRP]] %in% GRP2] 

group1 <- as.matrix(atac_[,group1_])
group2 <- as.matrix(atac_[,group2_])

results <- lapply(1:nrow(group1), function(i) {
  wilcox.test(group1[i, ], group2[i, ], alternative = "two.sided", paired = FALSE)
})

diff_acc <- data.frame(pval = sapply(results, function(res) res$p.value))
head(diff_acc)
diff_acc <- diff_acc %>% mutate(p_value_adjusted = p.adjust(pval, method = "fdr"),
              Significant = p_value_adjusted < 0.05) 
rownames(diff_acc) <- rownames(group1)

group1 <- as.data.frame(atac_[,group1_])
group2 <- as.data.frame(atac_[,group2_])

group1_ <- rowMeans(group1)
group2_ <- rowMeans(group2)
atacMean <- data.frame(group1_, group2_)
atacMean$delta <- log2(atacMean$group1_+1) - log2(atacMean$group2_+1)

names(atacMean) <- c(GRP1, GRP2, "atacDelta")

atac_Diff <- cbind(atacMean, diff_acc)

atac_Diff <- atac_Diff %>% 
    mutate(
      diffaccess = case_when(
        atacDelta > 0.6 & p_value_adjusted <= 0.05 ~ "high",
        atacDelta < -0.6 & p_value_adjusted <= 0.05 ~ "low",
        TRUE ~ "NO" 
      )
    )


#scatterplot of activity vs RNAseq logFC
tab <- data.frame(activity = atac_Diff$atacDelta, rna = rnaDEG$logFC, DEG = rnaDEG$diffExp, DiffAcc = atac_Diff$diffaccess)

tab <- tab %>% 
    mutate(
      combDiff = case_when(
        DEG %in% "high" & DiffAcc %in% "high" ~ "high",
        DEG %in% "low" & DiffAcc %in% "low" ~ "low",
        DEG %in% "high" & DiffAcc %in% "low" ~ "high_DEG,low_Activity",
        DEG %in% "low" & DiffAcc %in% "high" ~ "low_DEG,high_Activity",
        DEG %in% "high" & DiffAcc %in% "NO" ~ "high_DEG",
        DEG %in% "low" & DiffAcc %in% "NO" ~ "low_DEG",   
        DEG %in% "NO" & DiffAcc %in% "high" ~ "high_Activity",
        DEG %in% "NO" & DiffAcc %in% "low" ~ "low_Activity",        
        TRUE ~ "NO" 
      )
    )

table(tab$combDiff)
tab$names <- as.character(names(actGR))


mycolors <- c("darkred", "red1", "pink", "black", "blue", "purple", "turquoise", "orange", "grey")
names(mycolors) <- c("high", "high_DEG", "high_Activity", "low", "low_DEG", "low_Activity", "low_DEG,high_Activity", "high_DEG,low_Activity", "NO")    


pdf(paste0(plotDir, "/activity.pdf"), width=14, height=12)
# jpeg(paste0(getwd(), "/image.jpg"), width=12, height=12, units="in", res=600)
  ggplot(tab, aes(x = activity, y = rna)) +
  rasterise(geom_bin2d(data = subset(tab, combDiff == "NO"), bins = 400), dpi = 600) +
    scale_fill_gradient(low = "gray90", high = "gray10", guide = "none") +
  geom_point(data = subset(tab, combDiff != "NO"), aes(col = combDiff), size =3) +
    theme_classic() +
    geom_vline(xintercept = c(0), col = "black") +
    geom_hline(yintercept = c(0), col = "black") +
    ggtitle(paste("Activity vs Exp. aEIR", "Cor =", round(cor(tab$activity, tab$rna), 3))) +
    scale_colour_manual(values = mycolors) +
    geom_label_repel(data = subset(tab, combDiff %in%  c("high", "low")), aes(label = subset(tab, combDiff %in%  c("high", "low"))$names), size =6, max.overlaps = 30) +
    theme(text = element_text(size = 20))
dev.off()

# nolint end
# nolint start
####################################################################################
# Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
# Script: MAIN.R
# Author: MIJM
# Last Update: 27.08.26
# Description: Full analysis for Longitudinal ATACseq data
####################################################################################

source("./scripts/ATACseq/utils.R")
source("./scripts/ATACseq/chromVAR.R")


################
# Load object
################

OBJECT <- "long" #c("long", "regional", )
GRP <- "EirGrp"
GRP1 <- "high"
GRP2 <- "low"
CUTOFF <- 1

data_ <- loadOBJECT_(OBJECT)
data_ <- removeSamples_(data_, GRP) #c("responders", "EirGrp", )
# data_ <- removeSamples_(data_, "responders") #c("responders", "EirGrp", )
data_ <- addSampleAnnotCol(data_, "donors", droplevels(getSampleAnnot(data_)[,"donor"]))


################
# chromVAR
################

group1 <- getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP1] #%>% droplevels()
group2 <- getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP2] #%>% droplevels()

cv <- "Altius" #c("JASPAR", "Altius")
GR <- ".peaks.cons"

if(cv =="JASPAR"){
  ### chromVAR with JASPAR motifs
  cvRes <- getChromVarDev(data_, GR, motifs="jaspar")
}else{
  ### ChromVAR Motif Family
  cvRes <- chromVAR_Altius_(data_, GR)
}
devZ <- chromVAR::deviationScores(cvRes)


CORRECTION <- "donor"
ann <- getSampleAnnot(data_)
ann <- ann[colnames(devZ),]
identical(rownames(ann), colnames(devZ))

devZ <- takeResLimma_(devZ, ann, CORRECTION, GRP) 
diff_acc <- zDiff_(data_, devZ, GRP1 = GRP1, GRP2 = GRP2)
head(diff_acc)
summary(diff_acc$Significant)

chromDiffAcc <- annot_zDiff_(diff_acc, data_, GRP1 = GRP1, GRP2 = GRP2, CUTOFF = 1, padj = 0.05)
head(chromDiffAcc)
table(chromDiffAcc$DiffAcc)


# plot volcano plot 
mycolors <- c("darkred", "darkblue", "grey")
names(mycolors) <- c(GRP1, GRP2, "NO")

label_ <- chromDiffAcc %>% subset(DiffAcc != "NO")

# jpeg("image.jpg", width=12,height=12, units="in", res=600)
pdf(paste0(plotDir, "vd2.zDiff.EirGrp.pdf"), width=12,height=12 )
ggplot(data=chromDiffAcc, aes(x=zDiff, y=-log10(padj), col=DiffAcc)) + 
geom_point() + theme_classic() + 
geom_vline(xintercept=c(-CUTOFF, CUTOFF), col="red") +
geom_hline(yintercept=-log10(0.05), col="red") +
scale_colour_manual(values = mycolors) +
geom_text_repel(data=label_, aes(label=name), col="black", size = 6) +
theme(text = element_text(size = 20))
dev.off()

devZ <- devZ %>% subset(rownames(.) %in% label_$name)


group1 <- getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP1] #%>% droplevels()
group2 <- getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP2] #%>% droplevels()

group1 <- as.data.frame(devZ[,group1])
group2 <- as.data.frame(devZ[,group2])

group1_clust <- as.data.frame(kmeans(t(group1), centers = 2)$cluster) %>%
                setNames("cluster") 
                
group2_clust <- as.data.frame(kmeans(t(group2), centers = 2)$cluster) %>%
                setNames("cluster") 


group1_clust <- group1_clust %>% arrange(cluster)
group2_clust <- group2_clust %>% arrange(-cluster)

comb <- rbind(group1_clust, group2_clust)


annot_col <- data.frame("year" = as.numeric(getSampleAnnot(data_)$year), 
                        "log10_aeir" = log10(getSampleAnnot(data_)$aeir+1), 
                        "EirGrp" = getSampleAnnot(data_)$EirGrp
                        )
rownames(annot_col) <- getSampleAnnot(data_)$sampleId

colors_ <- c(
  "black","red","blue","green3","gold","purple","deeppink","gray50","cyan",
  "darkred","darkviolet","chocolate","orange",
  "navy","forestgreen","magenta","turquoise4","orangered","sienna",
  "dodgerblue3","firebrick","goldenrod3","darkorange3","brown","olivedrab3"
)

colors_ <- colors_[1:length(unique(annot_col$donor))]
names(colors_) <- unique(annot_col$donor)


annoCol<-list(
    EirGrp = c("low" = "green", "high" = "red"),
    log10_aeir = colorRampPalette( colors = c("green", "red"))(100),
    year = colorRampPalette( colors = c("brown", "gold"))(100)
)


annot_col <- annot_col[match(rownames(comb), rownames(annot_col)),]
head(annot_col)
identical(rownames(annot_col), rownames(comb))

annot_col <- annot_col %>%
                mutate(cluster = as.factor(comb$cluster))

devZ <- devZ[,rownames(annot_col)]
identical(rownames(annot_col), colnames(devZ))


# jpeg("./image.jpg", width=12,height=12, units="in", res=600)
pdf(paste0(plotDir, "vd2.chromVAR.EirGrp.Wilcox.pdf"), width=12,height=14 )
pheatmap(devZ, scale ='none', 
    cluster_rows = T,
    cluster_cols = F, 
    annotation_col = annot_col, 
    annotation_colors = annoCol, 
    fontsize = 13, 
    fontsize_row = 13, 
    fontsize_col = 7, 
    angle_col = 270, 
    border_color= NA, 
    show_rownames = T, 
    show_colnames = T, 
    fontfamily = "sans",
    ,breaks = seq(-5, 5, length.out = 100)
    )
dev.off()


# nolint end

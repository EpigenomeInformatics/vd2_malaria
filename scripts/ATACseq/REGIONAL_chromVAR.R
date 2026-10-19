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
library(limma)

################
# Load object
################

OBJECT <- "regional" #c("long", "regional", )
GRP <- "region"
GRP1 <- "Tororo"
GRP2 <- "Jinja"
CUTOFF <- 1

data_ <- loadOBJECT_(OBJECT)

################
# chromVAR
################

group1 <- getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP1] #%>% droplevels()
group2 <- getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP2] #%>% droplevels()

cv <- "JASPAR" #c("JASPAR", "Altius")
GR <- ".peaks.cons"

if(cv =="JASPAR"){
  ### chromVAR with JASPAR motifs
  cvRes <- getChromVarDev(data_, GR, motifs="jaspar")
}else{
  ### ChromVAR Motif Family
  cvRes <- chromVAR_Altius_(data_, GR)
}
devZ <- chromVAR::deviationScores(cvRes)


# Limma

annot <- getSampleAnnot(data_)

design <- model.matrix(~ 0 + region, data = getSampleAnnot(data_))
colnames(design) <- c("Jinja", "Tororo") 
design

corfit <- duplicateCorrelation(devZ, design)

fit <- lmFit(devZ, design, correlation = corfit$consensus)
fit

contrast.matrix <- makeContrasts(TvJ = Tororo - Jinja, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

top <- topTable(fit2, coef = "TvJ", number = Inf)

TOP <- top %>% filter(adj.P.Val < 0.05) %>% nrow()
grp <- top %>%
  mutate(
    rank_pi = row_number(adj.P.Val),             
    direction = case_when(
      rank_pi <= TOP & logFC > 0 ~ "Tororo",
      rank_pi <= TOP & logFC < 0 ~ "Jinja",
      TRUE ~ "NO"
    )
  ) %>%
  select(-rank_pi)
table(grp$direction)



label_ <- grp %>% subset(direction != "NO") %>% mutate(name = sub(".*_", "", rownames(.)))

mycolors <- c("darkred", "darkblue", "grey")
names(mycolors) <- c(GRP1, GRP2, "NO")

# jpeg("image.jpg", width=12,height=12, units="in", res=600)
pdf(paste0(plotDir, "vd2.zDiff.regional.pdf"), width=12,height=12 )
ggplot(data=grp, aes(x=logFC, y=-log10(adj.P.Val), col=direction)) + 
geom_point() + theme_classic() + 
geom_vline(xintercept=c(-CUTOFF, CUTOFF), col="red") +
geom_hline(yintercept=-log10(0.05), col="red") +
scale_colour_manual(values = mycolors) +
geom_text_repel(data=label_, aes(label=name), col="black", size = 4) +
theme(text = element_text(size = 20))
dev.off()


rownames(devZ) <- sub(".*_", "", rownames(devZ))
devZ <- devZ %>% subset(rownames(.) %in% label_$name)
dim(devZ)
rownames(devZ)

devZ <- devZ[,c(getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP2], "vd2_Tororo_3345", "vd2_Tororo_3348", "vd2_Tororo_3289","vd2_Tororo_3339",  "vd2_Tororo_3349")]


colnames(getSampleAnnot(data_))
dim(getSampleAnnot(data_))

Pvd2_ <- getSampleAnnot(data_)
annot_col <- Pvd2_ %>% select(region,  malariaPrevYear, parasiteStatus)

annot_col <- annot_col[c(getSampleAnnot(data_)$sampleId[getSampleAnnot(data_)[[GRP]] %in% GRP2], "vd2_Tororo_3345", "vd2_Tororo_3348", "vd2_Tororo_3289","vd2_Tororo_3339",  "vd2_Tororo_3349"),]


annoCol<-list(
    region=c("Jinja" =  "darkblue", "Tororo" = "darkred"), 
    parasiteStatus= c("Negative" = "green", "Positive" = "red")
)

dim(devZ)


pdf(paste0(plotDir, "vd2.chromVAR.jaspar.", GRP ,".pdf"), width=12,height=14 )
# jpeg("image.jpg", width=18, height= 20,units="in", res= 600)
pheatmap(devZ, 
  scale ='none', 
  cluster_rows = T, 
  cluster_cols = F, 
  annotation_col = annot_col, 
  annotation_colors = annoCol, 
  fontsize = 18, 
  fontsize_row = 15, 
#   fontsize_col = 13, 
  angle_col = 270, 
  border_color= NA, 
  show_rownames = T, 
  show_colnames = F,  
  fontfamily = "sans"
  ,breaks = seq(-5, 5, length.out = 100)
)
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

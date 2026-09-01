# nolint start
####################################################################################
# Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
# Script: plotting.R
# Author: MIJM
# Last Update: 27.08.26
# Description: ATACseq plotting
####################################################################################

source("/icbb/projects/mmaran/malariaScripts/utils.R")


# 1. TSS-Nfrag

data_ <- loadDsAcc("/icbb/projects/mmaran/malaria/dsATAC/data_filtered_aeir")
toKeep <- getSampleAnnot(data_) %>% filter(cellType %in% c("CD4", "VD2")) %>% pull(sampleId)
data_ <- removeSamples(data_, !getSampleAnnot(data_)$sampleId %in% toKeep)

annot <- getSampleAnnot(data_) %>%
            mutate(nFragments = getFragmentNum(data_),
                    tssEnrichment = getTssEnrichmentBatch(data_)$tssEnrichment.smoothed)

# labelling the poor QC samples
label_ <- annot %>% filter(sampleId %in% c("VD2_2013_d3089_21", "VD2_2014_d3089_25", "VD2_2013_d3371_128", "VD2_2016_d3326_76"))

jpeg(paste0(getwd(), "/image.jpg"), width=12, height=12, units="in", res=600)
# pdf(paste0(plotDir, "qc.vd2-cd4.pdf"), width=11, height=11) 
ggplot(annot, aes(x=nFragments, y=tssEnrichment, color = cellType)) +
geom_point(size = 5) +   theme_classic() + 
    labs(x = "Number of Fragments", y = "TSS Enrichment") +
    theme(text = element_text(size = 20)) +
    geom_label_repel(data = label_, aes(label = sampleId), size = 3)+
    geom_hline(yintercept=3, col="red") +
    geom_vline(xintercept=1000000, col="red") +
    scale_colour_manual(values = colors_)
dev.off()



# 2. PCA

# TODO: include long and pilot full as well
OBJECT <- "long" #c(regional_reduced, long_reduced, long, regional) 
data_ <- loadOBJECT_(OBJECT)


RT <- ".peaks.cons"
data_qnorm <- transformCounts(data_, method="quantile", region = RT) 
data_log <- transformCounts(data_qnorm, method="log2", region = RT) 
peaks_gen <- data.frame(getCoord(data_log, RT))
cm <- getCounts(data_log, RT)
rownames(cm) <- peaks_gen$name
cm_prcomp <- muRtools::getDimRedCoords.pca(t(cm), components = c(1, 2))
annot <- getSampleAnnot(data_)
colnames(annot)


if(OBJECT == "long"){
    colors_ <- colors[1:length(unique(getSampleAnnot(data_)$donor))]
    names(colors_) <- unique(getSampleAnnot(data_)$donor)

    # # color by responder
    # colors <- c("yellowgreen", "sienna")  
    # colors_ <- colors[1:length(unique(getSampleAnnot(data_)$responder))]
    # names(colors_) <- unique(getSampleAnnot(data_)$responder)
} else{
    print("not long data; using default colors")
}

# Define:
colorCol_ <- "donor"
shapeCol_ <- "EirGrp"

jpeg("./image.jpg",width=10,height=11,units="in",res=600)
# pdf(paste0(plotDir, "PCA.vd2.donors.pdf"), width=11, height=11) 
    muRtools::getDimRedPlot(cm_prcomp, 
        annot = getSampleAnnot(data_), 
        colScheme = if (OBJECT == "long") colors_ else NULL, 
        colorCol = colorCol_, 
        shapeCol= shapeCol_, 
        ptSize = 5, addLabels = F) +
    theme_classic() +
    theme(text = element_text(size = 15))
dev.off()

# 3. 

# nolint end

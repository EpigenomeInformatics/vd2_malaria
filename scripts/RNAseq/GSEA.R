# nolint start
####################################################################################
# Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
# Script: GSEA.R
# Author: MIJM
# Last Update: 24Sept26
#' Description: RNAseq Analysis

####################################################################################

suppressPackageStartupMessages(c(
library(biomaRt),
library(ggplot2),
library(dplyr),
library(clusterProfiler),
library(org.Hs.eg.db),
library(enrichplot),
library(msigdbr)))

plotDir <- "./results/RNAseq/"
groupCol <- "aeir_group"

if(!exists("degDf")){
  print(paste("loading DEG:", groupCol))
  degDf <- read.csv(paste0("./data/RNA_vd2.", groupCol ,".csv"), row.names = 1)
} else{
  print("using already exiting degDf")
}

degDf <- degDf %>% dplyr::filter(!is.na(logFC))

geneList <- degDf %>% 
  mutate(rank = -log10(PValue) * sign(logFC)) %>% 
  arrange(desc(rank)) %>% 
  dplyr::select(genes, rank) %>%
  tibble::deframe()


# RELEVANT PATHWAYS

immuneKeywords <- c(
  "immune", "cytokine", "lymphocyte", "leukocyte", 
  "interferon", "interleukin", "chemokine", "antigen", 
  "T cell", "B cell", "myeloid", "defense response",
  "inflammatory", "parasite", "macrophage"
)

gseGo <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db, 
  ont          = "BP",        
  keyType      = "ENSEMBL",  
  minGSSize    = 10,          
  maxGSSize    = 500,        
  pvalueCutoff = 0.05,
  eps          = 0,
  verbose      = TRUE,
  seed         = TRUE      
)

immuneResults <- gseGo@result %>%
  filter(grepl(paste(immuneKeywords, collapse = "|"), Description, ignore.case = TRUE)) %>%
  filter(p.adjust <= 0.05) 

gseImmune <- gseGo
gseImmune@result <- immuneResults

pdf(paste0(plotDir, "RNA-vd2.dotPlot.GSEA_GO.", groupCol, ".pdf"), width = 14, height = 16)
dotplot(gseImmune, showCategory = 20, split = ".sign", orderBy = "x") + 
  facet_grid(. ~ .sign, labeller = as_labeller(c("activated" = "High aEIR", "suppressed" = "Low aEIR")))
dev.off()

gseImmune <- pairwise_termsim(gseImmune)

# emapplot
pdf(paste0(plotDir, "RNA-vd2.emaPlot.GSEA_GO.", groupCol, ".pdf"), width = 14, height = 16)
# jpeg("image.jpg", width=14, height=16, units="in", res=600)
enrichplot::emapplot(
  gseImmune, 
  showCategory = 30,
  layout = "kk",          
  cex_category = 1.2,    
  cex_label_category = 1 
)
dev.off()
# nolint end
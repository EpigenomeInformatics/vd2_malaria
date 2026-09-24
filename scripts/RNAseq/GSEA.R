# nolint start
####################################################################################
# Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
# Script: GSEA.R
# Author: MIJM
# Last Update: 24Sept26
#' Description: RNAseq Analysis

####################################################################################


library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(biomaRt)
library(enrichplot)
library(msigdbr)

plotDir <- "/icbb/projects/mmaran/malaria/2025/Manuscript/RNA/COMBINED/noSexChr/"


GROUP <- "aeir_group"
DEG  <- read.csv(paste0(plotDir, "RNA_vd2.", GROUP ,".csv"), row.names = 1)
head(DEG)
dim(DEG)
summary(is.na(DEG$logFC))

DEG %>% count(genes) %>%  
  dplyr::filter(n > 1)

DEG <- DEG %>% dplyr::filter(!is.na(logFC))
geneList <- DEG %>% 
                mutate(rank = -log10(PValue) * sign(logFC)) %>% 
                arrange(desc(rank)) %>% 
                dplyr::select(genes, rank) %>%
                tibble::deframe()

head(geneList)
tail(geneList)


# RELEVANT PATHWAYS

keywords <- c("immune", "cytokine", "lymphocyte", "leukocyte", 
              "interferon", "interleukin", "chemokine", "antigen", 
              "T cell", "B cell", "myeloid", "defense response",
              "inflammatory", "parasite", "macrophage")

gse_go <- gseGO(
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

immune_results <- gse_go@result %>%
  filter(grepl(paste(keywords, collapse = "|"), Description, ignore.case = TRUE)) %>%
  filter(p.adjust <= 0.05) 

gse_immune <- gse_go
gse_immune@result <- immune_results

# pdf(paste0(plotDir, "vd2.dotPlot.GSEA_GO.pval-signLogFC", GROUP, ".pdf"), width=14, height=16)
jpeg("image.jpg", width=14, height=16, units="in", res=600)
  dotplot(gse_immune, showCategory=20, split=".sign", orderBy = "x") + 
  facet_grid(.~.sign,  labeller = as_labeller(c("activated" = "High aEIR", "suppressed" = "Low aEIR")))
dev.off()

gse_immune <- pairwise_termsim(gse_immune)

# emapplot
pdf(paste0(plotDir, "vd2.emaPlot.GSEA_GO.pval-signLogFC", GROUP, ".pdf"), width=14, height=16)
# jpeg("image.jpg", width=14, height=16, units="in", res=600)
    enrichplot::emapplot(gse_immune, 
               showCategory = 30,
               layout = "kk",          
               cex_category = 1.2,    
               cex_label_category = 1 
)
dev.off()


# nolint end
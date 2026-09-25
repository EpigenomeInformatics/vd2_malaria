# nolint start
####################################################################################
#' Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
#' Script: utils.R
#' Author: MIJM
#' Last Update: 27.08.26
#'  Description: RNAseq libraries, utils, helper functions, and initial data processing.

####################################################################################

# 0. Libraries and Setup
set.seed(123)

suppressPackageStartupMessages(c(
library(edgeR),
library(tximport),
library(biomaRt),
library(ggplot2),
library(dplyr),
library(stringi),
library(rhdf5),
library(ggrepel)
))
options(ggrepel.max.overlaps = Inf)

plotDir <- "./results/RNAseq/"
setwd("/icbb/projects/mmaran/SCRIPTS/GitHubRepos/vd2_malaria/")

########### GENCODEv27 #################

toKeepBiotype <- c(
  "protein_coding", "antisense", "bidirectional_promoter_lncRNA", 
  "lincRNA", "lncRNA", "macro_lncRNA", "sense_intronic", 
  "sense_overlapping", "TR_C_gene", "TR_D_gene", 
  "TR_J_gene", "TR_V_gene"
)

groupLevels <- list(year_group = c("post", "pre"), 
                aeir_group = c("low", "high" ), 
                priormal12Grp = c("low", "high"), 
                y13v14 = c(2014, 2013), 
                y13v16 = c(2016, 2013), 
                eirB = c(0,1), 
                priormal12 = c(0,1), 
                LAMP = c("Negative", "Positive"), 
                anyparasitemia = c(0, 1),
                invitro = c("Unstim", "Stim"))

genesOfInterest <- c(
  "HAVCR2", "KIR3DL2", "KIR2DL4", "PLCG2", "CD8A", "KIR2DL1", "LYN", 
  "KIR3DL1", "IL10", "CD244", "AKAP5", "CEP78", "SLCO4C1", "CD44", 
  "TNF", "IFNG", "IL7R", "GZMK", "IL23R"
)

# 1. tximport
filePath <- "/icbb/projects/mmaran/malaria/rna/rna_seq_v202205/snakeRNA/output/quant"

inVivoSamples <- list.files(file.path(filePath,"kallisto/vd2"))
inVitroSamples <- list.files(file.path(filePath, "invitro/vd2"))
allSamples <- c(inVivoSamples, inVitroSamples)

inVivoFiles <- file.path(filePath,"kallisto/vd2", inVivoSamples, "abundance.h5")
inVitroFiles <- file.path(filePath, "invitro/vd2", inVitroSamples, "abundance.h5")
allFiles <- c(inVivoFiles, inVitroFiles)

stopifnot(identical(allSamples, gsub(".*/([^/]+)/abundance\\.h5", "\\1", allFiles)))
names(allFiles) <- allSamples
stopifnot(all(file.exists(allFiles)))

abundanceTemp <- read.table("/icbb/projects/mmaran/malaria/rna/abundanceforTximport.tsv", header = TRUE, sep = '\t')
splitValues <- strsplit(as.character(abundanceTemp$target_id), "|", fixed = TRUE)
splitDf <- data.frame(matrix(unlist(splitValues), ncol = length(splitValues[[1]]), byrow = TRUE))
splitDf <- splitDf %>% select(X1, X2, X6, X7, X8)
names(splitDf) <- c("transcript", "genes", "name", "id", "type")
splitDf$id <- NULL
tx2Gene <- splitDf[, 1:2]

vd2 <- tximport(allFiles, type = "kallisto", tx2gene = tx2Gene, txOut = FALSE, ignoreAfterBar = TRUE)


sampleAnnot <- read.csv("./data/RNA/RNAcomprehensiveSampleAnnot.csv", row.names = 1)

# 3. Obtain gene information from abundance file

geneInfo <- splitDf[!duplicated(splitDf$genes), ]
rownames(geneInfo) <- geneInfo$genes
geneInfo <- geneInfo[rownames(vd2$counts), ]

stopifnot(identical(rownames(vd2$counts), rownames(geneInfo)))
sampleAnnot <- sampleAnnot[order(allSamples), ]

stopifnot(
  identical(as.character(allSamples), sampleAnnot$sampleId),
  identical(sampleAnnot$sampleId, colnames(vd2$counts))
)

geneInfo$genes <- gsub("\\..*", "", geneInfo$genes)

# 4. Filter and Normalize

normMat <- as.matrix(vd2$length)
normMat <- normMat / exp(rowMeans(log(normMat)))
normCts <- vd2$counts / normMat
effLib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, effLib, "*")
normMat <- log(normMat)

dgeMaster <- DGEList(counts = round(vd2$counts), genes = geneInfo)
dgeMaster <- scaleOffset(dgeMaster, normMat)

groupCol <- "aeir_group"
tempAnnot <- sampleAnnot
tempAnnot <- droplevels(tempAnnot[!tempAnnot[[groupCol]] == "", ])

grp1 <- levels(factor(tempAnnot[[groupCol]]))[1]
grp2 <- levels(factor(tempAnnot[[groupCol]]))[2]

filteredDf <- tempAnnot %>%
  group_by(donors) %>%
  filter(any(.data[[groupCol]] == grp1) & any(.data[[groupCol]] == grp2)) %>%      
  data.frame()

targetSamples <- allSamples[allSamples %in% filteredDf$sampleId]
tempAnnot <- tempAnnot %>% mutate(aeirTp = ifelse(sampleId %in% targetSamples, aeir_group, "NO"))
keepGenes <- filterByExpr(dgeMaster, group = tempAnnot$aeirTp)
dgeMaster <- dgeMaster[keepGenes, ]

# 5. GENCODE annotation

gencodeGr <- readRDS("./data/gencodeV27_genes.rds")
gencodeGr <- gencodeGr[[1]]
gencodeGr <- keepStandardChromosomes(gencodeGr, pruning.mode = "coarse")
gencodeGr <- gencodeGr[!seqnames(gencodeGr) %in% c("chrM", "chrX", "chrY")]
gencodeGr <- dropSeqlevels(gencodeGr, c("chrM", "chrX", "chrY"))

mcols(gencodeGr) <- mcols(gencodeGr)[, names(mcols(gencodeGr)) %in% c("type", "gene_id", "gene_type", "gene_name")]
mcols(gencodeGr)$genes <- gsub("\\..*", "", mcols(gencodeGr)$gene_id)

geneInfo <- geneInfo[geneInfo$genes %in% mcols(gencodeGr)$genes, ]
geneInfo$genesVersion <- rownames(geneInfo)

geneInfo <- geneInfo[geneInfo$type %in% toKeepBiotype, ]
geneInfo$type <- NULL

geneInfo <- merge(as.data.frame(gencodeGr), geneInfo, by = "genes") 
rownames(geneInfo) <- geneInfo$genesVersion

dgeMaster <- dgeMaster[dgeMaster$genes$genes %in% geneInfo$genes, ]
dgeMaster$samples$donors <- as.factor(sampleAnnot[["donors"]])

# nolint end
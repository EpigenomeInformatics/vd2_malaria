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




db <- readRDS("/icbb/projects/mmaran/malaria/lolaDatabase.rds")

coord_ <- getCoord(data_, ".peaks.cons")
names(coord_) <- coord_$name
identical(names(coord_), rownames(diffAcc))

high <- diffAcc %>% subset(diffaccess %in% "high")
dim(high)

low <- diffAcc %>% subset(diffaccess %in% "low")
dim(low)

lolaRes_high <- LOLA::runLOLA(coord_[rownames(high)], coord_, db, cores=8)
lolaRes_low <- LOLA::runLOLA(coord_[rownames(low)], coord_, db, cores=8)


table(lolaRes_high$collection)

COLLECTION <- "TF_motif_clusters"
high_ <- lolaRes_high %>% mutate(direction = "High aEIR") %>% 
                        filter(collection %in% COLLECTION) %>%
                        # mutate(qValue= qvalue::qvalue(p= 10^(-pValueLog))$qvalues) %>%
                        mutate(log2OR = log2(oddsRatio+0.001),
                        significant = ifelse(qValue < 0.05 & log2OR >= 0.3, TRUE, FALSE)) %>%
                        arrange(-log2OR) %>%
                        head(50)

low_ <- lolaRes_low %>% mutate(direction = "Low aEIR") %>% 
                        filter(collection %in% COLLECTION) %>%
                        # mutate(qValue= qvalue::qvalue(p= 10^(-pValueLog))$qvalues) %>%
                        mutate(log2OR = log2(oddsRatio+0.001),
                        significant = ifelse(qValue < 0.05 & log2OR >= 0.3, TRUE, FALSE)) 

table(high_$significant)
table(low_$significant)

high_ %>% filter(significant %in% TRUE) %>% select(description, log2OR, qValue) #%>% pull(abs(log2OR)) %>% summary()
low_ %>% filter(significant %in% TRUE) %>% select(description, log2OR, qValue) #%>% pull(abs(log2OR)) %>% summary()

high_ %>% filter(significant %in% TRUE) %>% select(description)
summary(low_$description %in% high_$description)


# Common TFs
tmp1 <- high_ %>%
  filter(significant %in% TRUE) %>% 
  pull(description) 

tmp2 <- low_ %>%
  filter(significant %in% TRUE) %>% 
  pull(description) 

summary(tmp1 %in% tmp2)
summary(tmp2 %in% tmp1)
common <- tmp1[tmp1 %in% tmp2]


# Subsetting
high_ <- high_ %>% 
          filter(significant %in% TRUE,
                  !description %in% common) 
low_ <- low_ %>% filter(significant %in% TRUE,
                  !description  %in% common) 


# plotting

df2p <- rbind(high_, low_) %>% 
  filter(significant == TRUE) %>%
  arrange(direction, -log2OR) %>% 
  mutate(target = factor(description, levels = unique(description))) 
  


aesObj <- aes(x = target, y = log2OR, fill = qValue)
BASE_SIZE <- 20

pp <- ggplot(df2p) + aesObj + 
      geom_col(position = "dodge") + 
      scale_x_discrete(name = "Transcription Factor") + 
      ylab("Log2 Odds Ratio") +
      theme_bw(base_size = BASE_SIZE) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE_SIZE * 1.1),
        axis.title = element_text(face = "bold", size = BASE_SIZE),
        strip.text = element_text(face = "bold", size = BASE_SIZE * 1.1),
        legend.text = element_text(size = BASE_SIZE * 0.9),
        legend.title = element_text(size = BASE_SIZE)
      )


if ("collection" %in% colnames(df2p) && length(unique(df2p$collection)) > 1){
  pp <- pp + facet_grid(. ~ direction + collection, scales = "free", space = "free")
} else {
  pp <- pp + facet_grid( . ~ direction, scales = "free_x", space = "free")
}



pdf(paste0(plotDir, "LOLA_Long.pdf"), width=24, height= 14)
pp
dev.off()


# nolint end

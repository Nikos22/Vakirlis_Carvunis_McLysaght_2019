library(dplyr)
library(ggplot2)
library(mltools)
library(gridExtra)
library(reshape2)

setwd("~/Documents/research/ortho_buster/")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

df <- read.csv("scripts_submission/final_scripts/Supplementary_Table_1.csv")

df$species <- df$speciesShort
df <- dplyr::select(df, -speciesShort)

#
df <- filter(df, species!="Caenorhabditis_elegans")
#

####################

df <- arrange(df, div)
df$species <- factor(df$species, ordered = TRUE, levels = df$species)
df$tag <- factor(df$tag, levels=c("yeast", "fruitfly", "human"), ordered = TRUE)

hist_fn <- ggplot() +
  geom_bar(data = df, 
              aes(x=species, 
                  y=total,
                  fill=tag),
              stat="identity") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme(text = element_text(size=16),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("") +
  ylab("Total # of genes in conserved micro-synteny regions") +
  facet_wrap(~tag, scales="free", dir = "v", drop=TRUE)

ggsave("scripts_submission/final_scripts/SupplementaryFigure1.png", height=10, width=5)

corPred_fn <- ggplot() +
  geom_jitter(data = df, 
           aes(x=div, 
               y=(foundPred)/(foundPred+found_elsewhere),
               colour=tag),
           size=2.4,
           alpha=0.6) +
  scale_colour_manual(values = cbPalette) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines")) +
  xlab("Time since divergence (my)") +
  ylab("Proportion of matches found opposite")

ggsave("scripts_submission/final_scripts/Figure3B.png", width=6, height=4.5)

No_of_comparisons_with_high_pred_found <- length(with(df, table(which(foundPred/(foundPred+found_elsewhere) > 0.9))))

#####

cerCor <- cor.test(with(filter(df, tag=="yeast"), div), with(filter(df, tag=="yeast"), no_match/total), method="pearson")
drosCor <- cor.test(with(filter(df, tag=="fruitfly"), div), with(filter(df, tag=="fruitfly"), no_match/total), method="pearson")
verCor <- cor.test(with(filter(df, tag=="human"), div), with(filter(df, tag=="human"), no_match/total), method="pearson")


div_fn <- ggplot() +
  geom_jitter(data = df, 
              aes(x=div, 
                  y=no_match/total,
                  colour=tag),
              size=3,
              alpha=0.7) +
  geom_smooth(data = df, 
              aes(x=div, 
                  y=no_match/total,
                  colour=tag),
              method="lm") +
  xlim(0,600) + 
  theme(text = element_text(size=16),
        axis.text = element_text(size=16, colour="black"),
        axis.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines")) +
  xlab("Time since divergence (my)") +
  ylab("Missed homology proportion") +
  scale_colour_manual(values = cbPalette)

ggsave(plot = div_fn, "scripts_submission/final_scripts/Figure4A.pdf", width=8, height=6)
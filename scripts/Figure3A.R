library(dplyr)
library(ggplot2)
library(mltools)
library(gridExtra)
library(tidyr)

setwd("~/Documents/research/Vakirlis_Carvunis_McLysaght_2019/")

################

#
# load data
#

df <- read.csv("Figure_3-source_data_1.csv")
df$tag <- factor(df$tag, levels=c("yeast", "fruitfly", "human"), ordered = TRUE)

#
df <- filter(df, species!="Caenorhabditis_elegans")
#

df_tbl_min_ev <- group_by(df, by=speciesShort) %>% slice(which.min(evalue))

#
# plot undetectable homology percentages 
#


fn_pl = ggplot() +
  geom_line(data = df,
            aes(x=-log(evalue,10),
                y=((not_found)/(total)),
                colour=div,
                group=speciesShort),
            size=1) +
  geom_text(data = df_tbl_min_ev,
            aes(x=28,
                y=((not_found)/(total)),
                label=speciesShort),
            hjust=-0.2,
            vjust=1) +
  xlim(0,36) +
  ylab("Putative missed homology proportion") +
  xlab("-log(E-value)") +
  scale_color_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 400) +
  facet_grid(~tag) + 
  theme(axis.text = element_text(size=16, colour="black"),
        axis.title = element_text(size=16),
        strip.text.x = element_text(size=16),
        legend.position = "none",
        panel.spacing = unit(1, "lines"))

#####

df <- read.csv("Figure_3-source_data_2.csv")
df$tag <- factor(df$tag, levels=c("yeast", "fruitfly", "human"), ordered = TRUE)
#
df <- filter(df, species!="Caenorhabditis_elegans")
#
df_tbl_min_ev <- group_by(df, by=speciesShort) %>% slice(which.min(evalue))

#
# plot false homology percentages 
#

fp_pl = ggplot() +
  geom_line(data = df,
            aes(x=-log(evalue,10),
                y=found/total,
                colour=div,
                group=speciesShort),
            size=1) +
  ylab("False homology proportion") +
  xlab("-log(E-value)") +
  scale_x_continuous(breaks = seq(0,16, by=2), limits = c(0,15)) +
  scale_color_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 400) +
  facet_grid(~tag) +
  labs(colour = "Time since\ndivergence (my)") +
  facet_grid(~tag) + 
  theme(axis.text = element_text(size=16, colour="black"),
        axis.title = element_text(size=16),
        strip.text.x = element_text(size=16),
        legend.position = "bottom",
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        panel.spacing = unit(1, "lines"),
        legend.key.width = unit(2, "lines"))

#####

pl <- grid.arrange(fn_pl, fp_pl,  ncol=1, heights = c(2,2.4))

ggsave(plot = pl, "figures/Figure3A.pdf", width=10, height = 10)

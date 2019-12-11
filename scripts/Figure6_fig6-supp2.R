library(dplyr)
library(ggplot2)
library(mltools)
library(gridExtra)
library(reshape2)

setwd("~/Documents/research/Vakirlis_Carvunis_McLysaght_2019/")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### Load and prepare data tables ####

sp1 <- read.csv("Figure3-supplement_1.csv",
                row.names=NULL,
                stringsAsFactors = FALSE)

vert_div <- read.table("divergence_times/vert_div.txt", row.names=1)
dros_div <- read.table("divergence_times/dros_div.txt", row.names=1)
cer_div <- read.table("divergence_times/cer_div.txt", row.names=1)

#
dros_div <- dros_div[rownames(dros_div)!="Caenorhabditis_elegans",]
#

cer_div$species <- rownames(cer_div)
dros_div$species <- rownames(dros_div)
vert_div$species <- rownames(vert_div)

cer_div_d <- cer_div
cer_div_d$species <- rownames(cer_div_d)
cer_div_d <- select(cer_div_d, "V2", "species")
cer_div_d$total_genes_checked <- 5997
colnames(cer_div_d) <- c("Divergence_time", "species", "total_genes_checked")
cer_div_d <- cer_div_d[order(cer_div_d$Divergence_time),]

vert_div_d <- vert_div
vert_div_d$species <- rownames(vert_div_d)
vert_div_d <- select(vert_div_d, "V2", "species")
vert_div_d$total_genes_checked <- 19892
colnames(vert_div_d) <- c("Divergence_time", "species", "total_genes_checked")
vert_div_d <- vert_div_d[order(vert_div_d$Divergence_time),]

dros_div_d <- dros_div
dros_div_d$species <- rownames(dros_div_d)
dros_div_d <- select(dros_div_d, "V5", "species")
dros_div_d$total_genes_checked <- 13929
colnames(dros_div_d) <- c("Divergence_time", "species", "total_genes_checked")
dros_div_d <- dros_div_d[order(dros_div_d$Divergence_time),]

cer_sim <- read.csv("synt_simil_tables/yeast_similarity_df.csv",
                    stringsAsFactors = FALSE)
vert_sim <- read.csv("synt_simil_tables/human_similarity_df.csv",
                    stringsAsFactors = FALSE)
dros_sim <- read.csv("synt_simil_tables/fruitfly_similarity_df.csv",
                    stringsAsFactors = FALSE)
cer_synt <- read.csv("synt_simil_tables/yeast_synteny_df.csv",
                    stringsAsFactors = FALSE)
vert_synt <- read.csv("synt_simil_tables/human_synteny_df.csv",
                     stringsAsFactors = FALSE)
dros_synt <- read.csv("synt_simil_tables/fruitfly_synteny_df.csv",
                     stringsAsFactors = FALSE)

cer_sim <- select(cer_sim, -Gene_ID)
cer_orphan <- melt(cer_sim,
                   id.vars = c("Gene_name"),
                   measure.vars=colnames(cer_sim)[2:ncol(cer_sim)])
cer_orphan <- filter(cer_orphan, value==FALSE) %>%
  select(-value)
colnames(cer_orphan) <- c("Gene_focal", "species")

cer_synt <- select(cer_synt, -Gene_ID)
cer_gp_oo <- melt(cer_synt,
                   id.vars = c("Gene_name"),
                   measure.vars=colnames(cer_synt)[2:ncol(cer_synt)])
cer_gp_oo <- filter(cer_gp_oo, value==TRUE) %>%
  select(-value)
colnames(cer_gp_oo) <- c("Gene_focal", "species")
###
vert_sim <- select(vert_sim, -Gene_ID)
vert_orphan <- melt(vert_sim,
                   id.vars = c("Gene_name"),
                   measure.vars=colnames(vert_sim)[2:ncol(vert_sim)])
vert_orphan <- filter(vert_orphan, value==FALSE) %>%
  select(-value)
colnames(vert_orphan) <- c("Gene_focal", "species")

vert_synt <- select(vert_synt, -Gene_ID)
vert_gp_oo <- melt(vert_synt,
                  id.vars = c("Gene_name"),
                  measure.vars=colnames(vert_synt)[2:ncol(vert_synt)])
vert_gp_oo <- filter(vert_gp_oo, value==TRUE) %>%
  select(-value)
colnames(vert_gp_oo) <- c("Gene_focal", "species")
###
dros_sim <- select(dros_sim, -Gene_ID)
dros_orphan <- melt(dros_sim,
                    id.vars = c("Gene_name"),
                    measure.vars=colnames(dros_sim)[2:ncol(dros_sim)])
dros_orphan <- filter(dros_orphan, value==FALSE) %>%
  select(-value)
colnames(dros_orphan) <- c("Gene_focal", "species")

dros_synt <- select(dros_synt, -Gene_ID)
dros_gp_oo <- melt(dros_synt,
                   id.vars = c("Gene_name"),
                   measure.vars=colnames(dros_synt)[2:ncol(dros_synt)])
dros_gp_oo <- filter(dros_gp_oo, value==TRUE) %>%
  select(-value)
colnames(dros_gp_oo) <- c("Gene_focal", "species")


#### Calculate phylogeny-based proportions ####


cer_div$species <- rownames(cer_div)
cer_div <- cer_div[order(cer_div$V2),]
pcts_cer_all <- c()
for (i in unique(cer_div$V2))
{
  species_to_get <- cer_div[which(cer_div$V2 >= i), "species"]
  cer_orphan_p <- cer_orphan %>%
    group_by(Gene_focal) %>%
    filter(all(species_to_get %in% species)) %>%
    summarise(no_sp = length(unique(species)))
  pcts_cer_all <- c(pcts_cer_all,length(unique((cer_orphan_p$Gene_focal)))/cer_div_d$total_genes_checked[1])
}


pcts_cer_cons <- c()
pcts_cer_cons_num <- c()

for (i in unique(cer_div$V2))
{
  species_to_get <- cer_div[which(cer_div$V2 >= i), "species"]
  genes_with_further <- unique(cer_gp_oo[which(cer_gp_oo$species %in% species_to_get), "Gene_focal"])
  
  cer_orphan_p <- cer_orphan %>%
    group_by(Gene_focal) %>%
    filter(all(species_to_get %in% species)) %>%
    summarise(no_sp = length(unique(species)))
  
  pcts_cer_cons <- c(pcts_cer_cons,length(unique(which(genes_with_further %in% cer_orphan_p$Gene_focal)))/length(unique(genes_with_further)))
  pcts_cer_cons_num <- c(pcts_cer_cons_num, length(unique(genes_with_further)))
}

vert_div$species <- rownames(vert_div)
vert_div <- vert_div[order(vert_div$V2),]

pcts_vert_all <- c()
for (i in unique(vert_div$V2))
{
  species_to_get <- vert_div[which(vert_div$V2 >= i), "species"]
  vert_orphan_p <- vert_orphan %>%
    group_by(Gene_focal) %>%
    filter(all(species_to_get %in% species)) %>%
    summarise(no_sp = length(unique(species)))
  pcts_vert_all <- c(pcts_vert_all,length(unique((vert_orphan_p$Gene_focal)))/vert_div_d$total_genes_checked[1])
}


pcts_vert_cons <- c()
pcts_vert_cons_num <- c()
for (i in unique(vert_div$V2))
{
  
  species_to_get <- vert_div[which(vert_div$V2 >= i), "species"]
  genes_with_further <- unique(vert_gp_oo[which(vert_gp_oo$species %in% species_to_get), "Gene_focal"])

  vert_orphan_p <- vert_orphan %>%
    group_by(Gene_focal) %>%
    filter(all(species_to_get %in% species)) %>%
    summarise(no_sp = length(unique(species)))
  
  pcts_vert_cons <- c(pcts_vert_cons,length(unique(which(genes_with_further %in% vert_orphan_p$Gene_focal)))/length(unique(genes_with_further)))
  pcts_vert_cons_num <- c(pcts_vert_cons_num, length(unique(genes_with_further)))
}


dros_div <- dros_div[order(dros_div$V5),]
dros_div$species <- rownames(dros_div)
pcts_dros_all <- c()
for (i in unique(dros_div$V5))
{
  species_to_get <- dros_div[which(dros_div$V5 >= i), "species"]
  dros_orphan_p <- dros_orphan %>%
    group_by(Gene_focal) %>%
    filter(all(species_to_get %in% species)) %>%
    summarise(no_sp = length(unique(species)))
  pcts_dros_all <- c(pcts_dros_all,length(unique((dros_orphan_p$Gene_focal)))/dros_div_d$total_genes_checked[1])
}


pcts_dros_cons <- c()
pcts_dros_cons_num <- c()
for (i in unique(dros_div$V5))
{
  
  species_to_get <- dros_div[which(dros_div$V5 >= i), "species"]
  genes_with_further <- unique(dros_gp_oo[which(dros_gp_oo$species %in% species_to_get), "Gene_focal"])
  
  dros_orphan_p <- dros_orphan %>%
    group_by(Gene_focal) %>%
    filter(all(species_to_get %in% species)) %>%
    summarise(no_sp = length(unique(species)))
  
  pcts_dros_cons <- c(pcts_dros_cons,length(unique(which(genes_with_further %in% dros_orphan_p$Gene_focal)))/length(unique(genes_with_further)))
  pcts_dros_cons_num <- c(pcts_dros_cons_num, length(unique(genes_with_further)))
}


##### create data frames and plot pairwise proportions first ####

cer_df_hist_synt <- sp1 %>%
  filter(tag=="yeast") %>%
  mutate(pct = (no_match/total),
         sep = sqrt(pct*(1-pct)/total)) %>%
  arrange(div)


cer_df_hist_all <- as.data.frame(table(cer_orphan$species)) %>%
  mutate(pct = ((Freq)/cer_div_d$total_genes_checked),
         sep = sqrt(pct*(1-pct)/cer_div_d$total_genes_checked)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = cer_div_d[match(species, cer_div_d$species), "Divergence_time"]) %>%
  arrange(div)

cer_div_jitter <- jitter(cer_df_hist_synt$div, factor=5)


ver_df_hist_synt <- sp1 %>%
  filter(tag=="human") %>%
  mutate(pct = (no_match/total),
         sep = sqrt((pct*(1-pct))/total)) %>%
  arrange(div)


ver_df_hist_all <- as.data.frame(table(vert_orphan$species)) %>%
  mutate(pct = ((Freq)/vert_div_d$total_genes_checked),
         sep = sqrt((pct*(1-pct))/vert_div_d$total_genes_checked)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = vert_div_d[match(species, vert_div_d$species), "Divergence_time"]) %>%
  arrange(div)



dros_df_hist_synt <- sp1 %>%
  filter(tag=="fruitfly") %>%
  mutate(pct = (no_match/total),
         sep = sqrt((pct*(1-pct))/total)) %>%
  arrange(div)

dros_orphan <- filter(dros_orphan , species != "Solenopsis_invicta", species != "Caenorhabditis_elegans")
dros_df_hist_all <- as.data.frame(table(dros_orphan$species)) %>%
  mutate(pct = ((Freq)/dros_div_d$total_genes_checked),
         sep = sqrt((pct*(1-pct))/dros_div_d$total_genes_checked)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = dros_div_d[match(species, dros_div_d$species), "Divergence_time"]) %>%
  arrange(div)

dros_div_jitter <- jitter(dros_df_hist_synt$div, factor=5)



#######

ver_df_hist_synt$tag <- "With undetected homologues"
ver_df_hist_synt_h <- select(ver_df_hist_synt, pct, sep, div, tag, species)
ver_df_hist_all$tag <- "All genes"
ver_df_hist_all_h <- select(ver_df_hist_all, pct, sep, div, tag, species)
ver_df_c <- bind_rows(ver_df_hist_synt_h, ver_df_hist_all_h)

ver_df_c$species <- as.character(ver_df_c$species)
ver_df_c$species <- paste(substr(ver_df_c$species, 1, 1), substr(unlist(strsplit(ver_df_c$species, split = "_"))[2*(1:length(ver_df_c$species))],1,3), sep='')

ver_df_c$species <- factor(ver_df_c$species, levels = unique(ver_df_c$species), ordered = TRUE)

ver_hist_sep_pl <- ggplot(ver_df_c,
                          aes(x=species, 
                              y=pct,
                              ymin=pct-sep,
                              ymax=pct+sep,
                              alpha=tag)) + 
  geom_bar(stat="identity",
           position = position_dodge(), 
           fill="#56B4E9") +
  geom_errorbar(position=position_dodge(width=0.9),
                width=0.4) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("Phylostrata, time since divergence") +
  xlab("")+
  ylab("")


dros_df_hist_synt$tag <- "With undetected homologues"
dros_df_hist_synt_h <- select(dros_df_hist_synt, pct, sep, div, tag, species)
dros_df_hist_all$tag <- "All genes"
dros_df_hist_all_h <- select(dros_df_hist_all, pct, sep, div, tag, species)
dros_df_c <- bind_rows(dros_df_hist_synt_h, dros_df_hist_all_h)

dros_df_c$species <- as.character(dros_df_c$species)
dros_df_c$species <- paste(substr(dros_df_c$species, 1, 1), substr(unlist(strsplit(dros_df_c$species, split = "_"))[2*(1:length(dros_df_c$species))],1,3), sep='')

dros_df_c$species <- factor(dros_df_c$species, levels = unique(dros_df_c$species), ordered = TRUE)

dros_hist_sep_pl <- ggplot(dros_df_c,
                          aes(x=species, 
                              y=pct,
                              ymin=pct-sep,
                              ymax=pct+sep,
                              alpha=tag)) + 
  geom_bar(stat="identity",
           position = position_dodge(), 
           fill="#E69F00") +
  geom_errorbar(position=position_dodge(width=0.9),
                width=0.4) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("Phylostrata, time since divergence") +
  xlab("") +
  ylab("")


cer_df_hist_synt$tag <- "With undetected homologues"
cer_df_hist_synt_h <- select(cer_df_hist_synt, pct, sep, div, tag, species)
cer_df_hist_all$tag <- "All genes"
cer_df_hist_all_h <- select(cer_df_hist_all, pct, sep, div, tag, species)
cer_df_c <- bind_rows(cer_df_hist_synt_h, cer_df_hist_all_h)

cer_df_c$species <- as.character(cer_df_c$species)
cer_df_c$species <- paste(substr(cer_df_c$species, 1, 1), substr(unlist(strsplit(cer_df_c$species, split = "_"))[2*(1:length(cer_df_c$species))],1,3), sep='')

cer_df_c$species <- factor(cer_df_c$species, levels = cer_df_hist_synt$speciesShort, ordered = TRUE)

cer_hist_sep_pl <- ggplot(cer_df_c,
                           aes(x=species, 
                               y=pct,
                               ymin=pct-sep,
                               ymax=pct+sep,
                               alpha=tag)) + 
  geom_bar(stat="identity",
           position = position_dodge(), 
           fill="#999999") +
  geom_errorbar(position=position_dodge(width=0.9),
                width=0.4) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("Phylostrata, time since divergence") +
  xlab("") +
  ylab("")


####

ver_df_hist_synt$tag <- "With undetected homologues"
ver_df_hist_synt_h <- select(ver_df_hist_synt, pct, sep, div, tag, species)
ver_df_hist_all$tag <- "All genes"
ver_df_hist_all_h <- select(ver_df_hist_all, pct, sep, div, tag, species)
ver_df_hist_synt_h$pct <- ver_df_hist_synt_h$pct/ver_df_hist_all_h$pct
ver_df_hist_all_h$pct <- 1

ver_df_c <- bind_rows(ver_df_hist_all_h, ver_df_hist_synt_h)

ver_df_c$species <- as.character(ver_df_c$species)
ver_df_c$species <- paste(substr(ver_df_c$species, 1, 1), substr(unlist(strsplit(ver_df_c$species, split = "_"))[2*(1:length(ver_df_c$species))],1,3), sep='')

ver_df_c$species <- factor(ver_df_c$species, levels = ver_df_hist_synt$speciesShort, ordered = TRUE)

hist_ver <- ggplot() +
  geom_bar(data = ver_df_c, 
           aes(x=species, 
               y=pct,
               fill=tag),
           stat="identity",
           position = position_identity(), 
           colour="black") +
  scale_fill_manual(values = c("white", "#56B4E9")) +
  geom_hline(aes(yintercept=mean(ver_df_c[ver_df_c$tag=="With undetected homologues", "pct"])),
             colour="red",
             size=2) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("") +
  ylab("")


dros_df_hist_synt$tag <- "With undetected homologues"
dros_df_hist_synt_h <- select(dros_df_hist_synt, pct, sep, div, tag, species)
dros_df_hist_all$tag <- "All genes"
dros_df_hist_all_h <- select(dros_df_hist_all, pct, sep, div, tag, species)
dros_df_hist_synt_h$pct <- dros_df_hist_synt_h$pct/dros_df_hist_all_h$pct
dros_df_hist_all_h$pct <- 1

dros_df_c <- bind_rows(dros_df_hist_all_h, dros_df_hist_synt_h)

dros_df_c$species <- as.character(dros_df_c$species)
dros_df_c$species <- paste(substr(dros_df_c$species, 1, 1), substr(unlist(strsplit(dros_df_c$species, split = "_"))[2*(1:length(dros_df_c$species))],1,3), sep='')

dros_df_c$species <- factor(dros_df_c$species, levels = dros_df_hist_synt$speciesShort, ordered = TRUE)

hist_dros <- ggplot() +
  geom_bar(data = dros_df_c, 
           aes(x=species, 
               y=pct,
               fill=tag),
           stat="identity",
           position = position_identity(), 
           colour="black") +
  scale_fill_manual(values = c("white", "#E69F00")) +
  geom_hline(aes(yintercept=mean(dros_df_c[dros_df_c$tag=="With undetected homologues", "pct"])),
             colour="red",
             size=2) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("") +
  ylab("")



cer_df_hist_synt$tag <- "With undetected homologues"
cer_df_hist_synt_h <- select(cer_df_hist_synt, pct, sep, div, tag, species)
cer_df_hist_all$tag <- "All genes"
cer_df_hist_all_h <- select(cer_df_hist_all, pct, sep, div, tag, species)
cer_df_hist_synt_h$pct <- cer_df_hist_synt_h$pct/cer_df_hist_all_h$pct
cer_df_hist_all_h$pct <- 1

cer_df_c <- bind_rows(cer_df_hist_all_h, cer_df_hist_synt_h)

cer_df_c$species <- as.character(cer_df_c$species)
cer_df_c$species <- paste(substr(cer_df_c$species, 1, 1), substr(unlist(strsplit(cer_df_c$species, split = "_"))[2*(1:length(cer_df_c$species))],1,3), sep='')

cer_df_c$species <- factor(cer_df_c$species, levels = unique(cer_df_c$species), ordered = TRUE)

hist_cer <- ggplot() +
  geom_bar(data = cer_df_c, 
           aes(x=species, 
               y=pct,
               fill=tag),
           stat="identity",
           position = position_identity(), 
           colour="black") +
  scale_fill_manual(values = c("white", "#999999")) +
  geom_hline(aes(yintercept=mean(cer_df_c[cer_df_c$tag=="With undetected homologues", "pct"])),
             colour="red",
             size=2) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("") +
  ylab("")


complete_fig <- grid.arrange(cer_hist_sep_pl,
                             dros_hist_sep_pl,
                             ver_hist_sep_pl,
                             hist_cer,
                             hist_dros,
                             hist_ver,
                             heights=c(0.5,0.5),
                             ncol=3,
                             nrow=2)

ggsave(plot = complete_fig, "figures/Figure6.pdf", width=13, height = 8)

# stats

all_3_pw <- bind_rows(ver_df_hist_synt_h, dros_df_hist_synt_h, cer_df_hist_synt_h)
ov_min_pw <- all_3_pw[which.min(all_3_pw$pct),]
ov_max_pw <- all_3_pw[which.max(all_3_pw$pct),]
avg_pw <- mean(all_3_pw$pct)


##### create data frames and plot phylogeny based proportions ####

ver_phylo_df_p1 <- data.frame(pct_found = pcts_vert_cons/pcts_vert_all,
                              tag = "With predicted undetectable homologues",
                              div = unique(vert_div$V2),
                              synt_missed = pcts_vert_cons,
                              all_missed = pcts_vert_all)

psynt <- data.frame(pct = pcts_vert_cons,
                 tag = "With undetected homologues",
                 div = as.character(unique(vert_div$V2)),
                 sep = sqrt(pcts_vert_cons*(1-pcts_vert_cons)/pcts_vert_cons_num))

pall <- data.frame(pct = pcts_vert_all,
                    tag = "all genes",
                    div = as.character(unique(vert_div$V2)),
                    sep = sqrt(pcts_vert_all*(1-pcts_vert_all)/vert_div_d$total_genes_checked[1]))

ver_df_c <- bind_rows(psynt, pall)

ver_df_c$div <- factor(ver_df_c$div, levels = unique(ver_df_c$div), ordered = TRUE)

ver_hist_sep_pl <- ggplot(ver_df_c,
                          aes(x=div, 
                              y=pct,
                              ymin=pct-sep,
                              ymax=pct+sep,
                              alpha=tag)) + 
  geom_bar(stat="identity",
           position = position_dodge(), 
           fill="#56B4E9") +
  geom_errorbar(position=position_dodge(width=0.9),
                width=0.4) +
  scale_x_discrete(labels=c(0,as.character(ver_df_c$div)[1:(length(ver_df_c$div)-1)])) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("Phylostrata, time since divergence") +
  xlab("") +
  ylab("")


dros_phylo_df_p1 <- data.frame(pct_found = pcts_dros_cons/pcts_dros_all,
                              tag = "With predicted undetectable homologues",
                              div = unique(dros_div$V5),
                              synt_missed = pcts_dros_cons,
                              all_missed = pcts_dros_all)

psynt <- data.frame(pct = pcts_dros_cons,
                    tag = "With undetected homologues",
                    div = as.character(unique(dros_div$V5)),
                    sep = sqrt(pcts_dros_cons*(1-pcts_dros_cons)/pcts_dros_cons_num))

pall <- data.frame(pct = pcts_dros_all,
                   tag = "all genes",
                   div = as.character(unique(dros_div$V5)),
                   sep = sqrt(pcts_dros_all*(1-pcts_dros_all)/dros_div_d$total_genes_checked[1]))

dros_df_c <- bind_rows(psynt, pall)

dros_df_c$div <- factor(dros_df_c$div, levels = unique(dros_df_c$div), ordered = TRUE)


dros_hist_sep_pl <- ggplot(dros_df_c,
                           aes(x=div, 
                               y=pct,
                               ymin=pct-sep,
                               ymax=pct+sep,
                               alpha=tag)) + 
  geom_bar(stat="identity",
           position = position_dodge(), 
           fill="#E69F00") +
  geom_errorbar(position=position_dodge(width=0.9),
                width=0.4) +
  scale_x_discrete(labels=c(0,as.character(dros_df_c$div)[1:(length(dros_df_c$div)-1)])) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("Phylostrata, time since divergence") +
  xlab("") +
  ylab("")


cer_phylo_df_p1 <- data.frame(pct_found = pcts_cer_cons/pcts_cer_all,
                               tag = "With predicted undetectable homologues",
                               div = unique(cer_div$V2),
                               synt_missed = pcts_cer_cons,
                               all_missed = pcts_cer_all)

psynt <- data.frame(pct = pcts_cer_cons,
                    tag = "With undetected homologues",
                    div = as.character(unique(cer_div$V2)),
                    sep = sqrt(pcts_cer_cons*(1-pcts_cer_cons)/pcts_cer_cons_num))

pall <- data.frame(pct = pcts_cer_all,
                   tag = "all genes",
                   div = as.character(unique(cer_div$V2)),
                   sep = sqrt(pcts_cer_all*(1-pcts_cer_all)/cer_div_d$total_genes_checked[1]))

cer_df_c <- bind_rows(psynt, pall)

cer_df_c$div <- factor(cer_df_c$div, levels = unique(cer_df_c$div), ordered = TRUE)
cer_hist_sep_pl <- ggplot(cer_df_c,
                          aes(x=div, 
                              y=pct,
                              ymin=pct-sep,
                              ymax=pct+sep,
                              alpha=tag)) + 
  geom_bar(stat="identity",
           position = position_dodge(), 
           fill="#999999") +
  geom_errorbar(position=position_dodge(width=0.9),
                width=0.4) +
  scale_x_discrete(labels=c(0,as.character(cer_df_c$div)[1:(length(cer_df_c$div)-1)])) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("Phylostrata, time since divergence") +
  xlab("") +
  ylab("")


####

ver_df_c <- bind_rows(mutate(ver_phylo_df_p1, pct_found=1, tag="all genes"),ver_phylo_df_p1)

ver_df_c$div <- factor(ver_df_c$div, levels = unique(ver_df_c$div), ordered = TRUE)

hist_ver <- ggplot() +
  geom_bar(data = ver_df_c, 
           aes(x=div, 
               y=pct_found,
               fill=tag),
           stat="identity",
           position = position_identity(), 
           colour="black") +
  scale_fill_manual(values = c("white", "#56B4E9")) +
  xlab("") +
  scale_x_discrete(labels=c(0,as.character(ver_df_c$div)[1:(length(ver_df_c$div)-1)])) +
  geom_hline(aes(yintercept=mean(ver_df_c[ver_df_c$tag=="With predicted undetectable homologues", "pct_found"])),
             colour="red",
             size=2) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("") +
  ylab("")


dros_df_c <- bind_rows(mutate(dros_phylo_df_p1, pct_found=1, tag="all genes"),dros_phylo_df_p1)

dros_df_c$div <- factor(dros_df_c$div, levels = unique(dros_df_c$div), ordered = TRUE)

hist_dros <- ggplot() +
  geom_bar(data = dros_df_c, 
           aes(x=div, 
               y=pct_found,
               fill=tag),
           stat="identity",
           position = position_identity(), 
           colour="black") +
  scale_fill_manual(values = c("white", "#E69F00")) +
  xlab("") +
  scale_x_discrete(labels=c(0,as.character(dros_df_c$div)[1:(length(dros_df_c$div)-1)])) +
  geom_hline(aes(yintercept=mean(dros_df_c[dros_df_c$tag=="With predicted undetectable homologues", "pct_found"])),
             colour="red",
             size=2) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("") +
  ylab("")



cer_df_c <- bind_rows(mutate(cer_phylo_df_p1, pct_found=1, tag="all genes"),cer_phylo_df_p1)

cer_df_c$div <- factor(cer_df_c$div, levels = unique(cer_df_c$div), ordered = TRUE)

hist_cer <- ggplot() +
  geom_bar(data = cer_df_c, 
           aes(x=div, 
               y=pct_found,
               fill=tag),
           stat="identity",
           position = position_identity(), 
           colour="black") +
  scale_fill_manual(values = c("white", "#999999")) +
  scale_x_discrete(labels=c(0,as.character(cer_df_c$div)[1:(length(cer_df_c$div)-1)])) +
  geom_hline(aes(yintercept=mean(cer_df_c[cer_df_c$tag=="With predicted undetectable homologues", "pct_found"])),
             colour="red",
             size=2) +
  xlab("") +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("") +
  ylab("")



#####

complete_fig <- grid.arrange(cer_hist_sep_pl,
                             dros_hist_sep_pl,
                             ver_hist_sep_pl,
                             hist_cer,
                             hist_dros,
                             hist_ver,
                             heights=c(0.5, 0.5),
                             ncol=3,
                             nrow=2)

ggsave(plot = complete_fig, "figures/Figure6-supplement_2.pdf", width=13, height = 7)

# stats

all_3_ph <- bind_rows(ver_phylo_df_p1, dros_phylo_df_p1, cer_phylo_df_p1)
ov_min_ph <- all_3_ph[which.min(all_3_ph$pct),]
ov_max_ph <- all_3_ph[which.max(all_3_ph$pct),]
avg_ph <- mean(all_3_ph$pct)

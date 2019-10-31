library(dplyr)
library(ggplot2)
library(mltools)
library(gridExtra)
library(reshape2)

setwd("~/Documents/research/ortho_buster/")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cer_res <- read.table("scripts_submission/no_residues/yeast_residues.txt", row.names=1)
dros_res <- read.table("scripts_submission/no_residues/droso_residues.txt", row.names=1)
ver_res <- read.table("scripts_submission/no_residues/vert_residues.txt", row.names=1)

vert_div <- read.table("scripts_submission/divergence_times/vert_div.txt", row.names=1)
dros_div <- read.table("scripts_submission/divergence_times/dros_div.txt", row.names=1)
cer_div <- read.table("scripts_submission/divergence_times/cer_div.txt", row.names=1)

#
dros_div <- dros_div[rownames(dros_div)!="Caenorhabditis_elegans",]
#

cer_div$species <- rownames(cer_div)
dros_div$species <- rownames(dros_div)
vert_div$species <- rownames(vert_div)

####################
#### dataframes ####
####################


dros_orphan <- read.table("scripts_submission/all_genes_blast/droso_all_genes_det_final.txt",
                          header=TRUE,
                          stringsAsFactors = FALSE)

dros_gp_tbl <- read.table("droso/normal/full_results_gene_pairs_001_25_9.txt",
                        as.is=TRUE,
                        header=TRUE,
                        stringsAsFactors = FALSE,
                        fill=TRUE)


dros_gp_oo <- read.table("scripts_submission/oo_dfs/full_results_slim_onlyOrth_dros.txt",
                         header=TRUE,
                         stringsAsFactors = FALSE,
                         row.names=NULL,
                         fill=TRUE)

dros_gp_tbl <- select(dros_gp_tbl, species, Gene_focal)
dros_orphan <- select(dros_orphan, species, Gene_focal)
dros_gp_oo <- select(dros_gp_oo, species, Gene_focal)

dros_gp_tbl <- unique(dros_gp_tbl)
dros_orphan <- unique(dros_orphan)
dros_gp_oo <- unique(dros_gp_oo)

dros_gp_tbl <- filter(dros_gp_tbl,
                     paste(species, Gene_focal) %in% paste(dros_orphan$species, dros_orphan$Gene_focal))

dros_gp_oo$Found_simil <- TRUE
dros_gp_oo[which(paste(dros_gp_oo$species, dros_gp_oo$Gene_focal) %in% paste(dros_gp_tbl$species, dros_gp_tbl$Gene_focal)), "Found_simil"] <- FALSE

dros_orphan <- dros_orphan[which(!paste(dros_orphan$species, dros_orphan$Gene_focal) %in% paste(dros_gp_oo[which(dros_gp_oo$Found_simil), "species"], dros_gp_oo[which(dros_gp_oo$Found_simil), "Gene_focal"])), ]


dros_div_d <- dros_div
dros_div_d$species <- rownames(dros_div_d)
dros_div_d <- select(dros_div_d, "V5", "species")
dros_div_d$total_genes_checked <- 13929
colnames(dros_div_d) <- c("Divergence_time", "species", "total_genes_checked")


#### yeast


cer_orphan <- read.table("scripts_submission/all_genes_blast/cer_all_genes_det_final.txt",
                          header=TRUE,
                          stringsAsFactors = FALSE)

cer_gp_tbl <- read.table("yeast/normal/cerevisiae_full_results_gene_pairs_001_25_9_slim.txt",
                        as.is=TRUE,
                        stringsAsFactors = FALSE,
                        header=TRUE)


cer_gp_oo <- read.table("scripts_submission/oo_dfs/cerevisiae_full_results_onlyOrth_slim.txt",
                        header=TRUE,
                        stringsAsFactors = FALSE,
                        row.names=NULL,
                        fill=TRUE)

cer_gp_tbl <- select(cer_gp_tbl, species, Gene_focal)
cer_orphan <- select(cer_orphan, species, Gene_focal)
cer_gp_oo <- select(cer_gp_oo, species, Gene_focal)

cer_gp_tbl <- unique(cer_gp_tbl)
cer_orphan <- unique(cer_orphan)
cer_gp_oo <- unique(cer_gp_oo)

cer_gp_tbl <- filter(cer_gp_tbl,
                     paste(species, Gene_focal) %in% paste(cer_orphan$species, cer_orphan$Gene_focal))

cer_gp_oo$Found_simil <- TRUE
cer_gp_oo[which(paste(cer_gp_oo$species, cer_gp_oo$Gene_focal) %in% paste(cer_gp_tbl$species, cer_gp_tbl$Gene_focal)), "Found_simil"] <- FALSE


cer_div_d <- cer_div
cer_div_d$species <- rownames(cer_div_d)
cer_div_d <- select(cer_div_d, "V2", "species")
cer_div_d$total_genes_checked <- 5997
colnames(cer_div_d) <- c("Divergence_time", "species", "total_genes_checked")

cer_df <- cer_gp_oo %>%
  group_by(species) %>%
  summarise(missed = length(which(Found_simil==FALSE))/length(Found_simil)) %>%
  mutate(div = cer_div_d[species, "Divergence_time"])


#### vertebrates


vert_orphan <- read.table("scripts_submission/all_genes_blast/vertebrates_forOrphan_all_-_final.txt",
                         header=TRUE,
                         stringsAsFactors = FALSE,
                         fill=TRUE)

vert_gp_tbl <- read.table("vertebrates/normal/all_gene_pair_results_slim_27_9_19.txt",
                       stringsAsFactors = FALSE,
                       header=TRUE,
                       sep='\t',
                       fill=TRUE)

vert_gp_oo <- read.table("scripts_submission/oo_dfs/vertebrates_full_results_onlyOrth_11_8.txt",
                         header=TRUE,
                         stringsAsFactors = FALSE,
                         row.names=NULL,
                         fill=TRUE)


vert_gp_tbl <- select(vert_gp_tbl, species, Gene_focal)
vert_orphan <- select(vert_orphan, species, Gene_focal)
vert_gp_oo <- select(vert_gp_oo, species, Gene_focal)

vert_gp_tbl <- unique(vert_gp_tbl)
vert_orphan <- unique(vert_orphan)
vert_gp_oo <- unique(vert_gp_oo)

# no match in oprhans
vert_gp_tbl <- filter(vert_gp_tbl,
                     paste(species, Gene_focal) %in% paste(vert_orphan$species, vert_orphan$Gene_focal))

vert_gp_oo$Found_simil <- TRUE
vert_gp_oo[which(paste(vert_gp_oo$species, vert_gp_oo$Gene_focal) %in% paste(vert_gp_tbl$species, vert_gp_tbl$Gene_focal)), "Found_simil"] <- FALSE


vert_div_d <- vert_div
vert_div_d$species <- rownames(vert_div_d)
vert_div_d <- select(vert_div_d, "V2", "species")
vert_div_d$total_genes_checked <- 19892
colnames(vert_div_d) <- c("Divergence_time", "species", "total_genes_checked")


vert_df <- vert_gp_oo %>%
  group_by(species) %>%
  summarise(missed = length(which(Found_simil==FALSE))/length(Found_simil)) %>%
  mutate(div = vert_div_d[species, "Divergence_time"])


##### prepare final dataframes to be loaded by next final script #####
##### if we add a list of all genes then we can have everything in one master table
##### and then manipulate it at will

vert_gp_oo$tag <- "human"
cer_gp_oo$tag <- "yeast"
dros_gp_oo$tag <- "fruitfly"

vert_all_genes_list <- read.table("scripts_submission/full_lists_genes/all_human_genes_checked_compl.txt",
                                  stringsAsFactors = FALSE,
                                  header=TRUE)

cer_all_genes_list <- read.table("scripts_submission/full_lists_genes/all_yeast_genes_checked_compl.txt",
                                  stringsAsFactors = FALSE,
                                  header=TRUE)

dros_all_genes_list <- read.table("scripts_submission/full_lists_genes/all_fruitfly_genes_checked_compl.txt",
                                  stringsAsFactors = FALSE,
                                  sep='\t',
                                  header=TRUE)
##
cer_simil_df <- cer_all_genes_list
for (sp in rownames(cer_div))
{
  cer_simil_df[,sp] <- TRUE
  cer_simil_df[which(paste(cer_simil_df$Gene_name, sp) %in% paste(cer_orphan$Gene_focal, cer_orphan$species)), sp] <- FALSE
}

cer_synt_df <- cer_all_genes_list
for (sp in rownames(cer_div))
{
  cer_synt_df[,sp] <- FALSE
  cer_synt_df[which(paste(cer_synt_df$Gene_name, sp) %in% paste(cer_gp_oo$Gene_focal, cer_gp_oo$species)), sp] <- TRUE
}
##
dros_simil_df <- dros_all_genes_list
for (sp in rownames(dros_div))
{
  dros_simil_df[,sp] <- TRUE
  dros_simil_df[which(paste(dros_simil_df$Gene_name, sp) %in% paste(dros_orphan$Gene_focal, dros_orphan$species)), sp] <- FALSE
}

dros_synt_df <- dros_all_genes_list
for (sp in rownames(dros_div))
{
  dros_synt_df[,sp] <- FALSE
  dros_synt_df[which(paste(dros_synt_df$Gene_name, sp) %in% paste(dros_gp_oo$Gene_focal, dros_gp_oo$species)), sp] <- TRUE
}
##
vert_simil_df <- vert_all_genes_list
for (sp in rownames(vert_div))
{
  vert_simil_df[,sp] <- TRUE
  vert_simil_df[which(paste(vert_simil_df$Gene_name, sp) %in% paste(vert_orphan$Gene_focal, vert_orphan$species)), sp] <- FALSE
}

vert_synt_df <- vert_all_genes_list
for (sp in rownames(vert_div))
{
  vert_synt_df[,sp] <- FALSE
  vert_synt_df[which(paste(vert_synt_df$Gene_name, sp) %in% paste(vert_gp_oo$Gene_focal, vert_gp_oo$species)), sp] <- TRUE
}

#
#dros_simil_df <- select(dros_simil_df, -Caenorhabditis_elegans)
#

write.csv(vert_simil_df, "scripts_submission/final_scripts/synt_simil_tables/human_similarity_df.csv", row.names = FALSE)
write.csv(cer_simil_df, "scripts_submission/final_scripts/synt_simil_tables/yeast_similarity_df.csv", row.names = FALSE)
write.csv(dros_simil_df, "scripts_submission/final_scripts/synt_simil_tables/fruitfly_similarity_df.csv", row.names = FALSE)

#
#dros_synt_df <- select(dros_synt_df, -Caenorhabditis_elegans)
#

write.csv(vert_synt_df, "scripts_submission/final_scripts/synt_simil_tables/human_synteny_df.csv", row.names = FALSE)
write.csv(cer_synt_df, "scripts_submission/final_scripts/synt_simil_tables/yeast_synteny_df.csv", row.names = FALSE)
write.csv(dros_synt_df, "scripts_submission/final_scripts/synt_simil_tables/fruitfly_synteny_df.csv", row.names = FALSE)


##### get the numbers to answer reviewer #2 comment #####


species_to_get <- cer_div[which(cer_div$V2 >= 15), "species"]
cer_orphan_p <- cer_orphan %>%
  group_by(Gene_focal) %>%
  filter(all(species_to_get %in% species)) %>%
  summarise(no_sp = length(unique(species)))
rel_orph_cer <- unique((cer_orphan_p$Gene_focal))

table(rel_orph_cer %in% cer_gp_oo$Gene_focal)

species_to_get <- vert_div[which(vert_div$V2 >= 3), "species"]
vert_orphan_p <- vert_orphan %>%
  group_by(Gene_focal) %>%
  filter(all(species_to_get %in% species)) %>%
  summarise(no_sp = length(unique(species)))
rel_orph_vert <- unique((vert_orphan_p$Gene_focal))

table(rel_orph_vert %in% vert_gp_oo$Gene_focal)



##### Figure 4 out of final data

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

#### phylostrat for cons synt

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

####

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

#####

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

sp1 <- read.csv("scripts_submission/final_scripts/Supplementary_Table_1.csv",
                row.names=NULL,
                stringsAsFactors = FALSE)

#####

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

cer_all_pl <- ggplot() + 
  geom_point(data = cer_df_hist_synt, 
             aes(x=cer_div_jitter, 
                 y=pct),
             size=3,
             colour="#999999") +
  geom_smooth(data = cer_df_hist_synt, 
              aes(x=cer_div_jitter, 
                  y=pct),
              method="lm",
              colour="#999999") +
  geom_point(data = cer_df_hist_all, 
             aes(x=cer_div_jitter, 
                 y=pct),
             size=3,
             colour="#999999",
             shape=2) +
  # geom_text(data = cer_df_tbl_opt,
  #           aes(x=cer_div_jitter,
  #               y=(no_match-found_tblastn)/total,
  #               label=species),
  #           hjust=-0.2,
  #           vjust=1) +
  geom_linerange(data = cer_df_hist_synt,
                 aes(x = cer_div_jitter,
                     ymin = pct,
                     ymax = with(cer_df_hist_all, pct)),
                 linetype = 3, 
                 alpha=0.5) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 19, hjust=0.5),
        axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("Time since divergence (my)") + ylab("Proportion without similarity\n(out of total genes)") +
  xlab("") + ylab("") +
  xlim(0,130) +
  ggtitle("yeast")


ver_df_hist_synt <- sp1 %>%
  filter(tag=="human") %>%
  mutate(pct = (no_match/total),
         sep = sqrt(pct*(1-pct)/total)) %>%
  arrange(div)


ver_df_hist_all <- as.data.frame(table(vert_orphan$species)) %>%
  mutate(pct = ((Freq)/vert_div_d$total_genes_checked),
         sep = sqrt(pct*(1-pct)/vert_div_d$total_genes_checked)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = vert_div_d[match(species, vert_div_d$species), "Divergence_time"]) %>%
  arrange(div)

ver_div_jitter <- jitter(ver_df_hist_synt$div, factor=5)



ver_all_pl <- ggplot() + 
  geom_point(data = ver_df_hist_synt, 
             aes(x=ver_div_jitter, 
                 y=pct),
             size=3,
             colour="#56B4E9") +
  geom_smooth(data = ver_df_hist_synt, 
              aes(x=ver_div_jitter, 
                  y=pct),
              method="lm",
              colour="#56B4E9") +
  geom_point(data = ver_df_hist_all, 
             aes(x=ver_div_jitter, 
                 y=pct),
             size=3,
             colour="#56B4E9",
             shape=2) +
  # geom_text(data = cer_df_tbl_opt,
  #           aes(x=cer_div_jitter,
  #               y=(no_match-found_tblastn)/total,
  #               label=species),
  #           hjust=-0.2,
  #           vjust=1) +
  geom_linerange(data = ver_df_hist_synt,
                 aes(x = ver_div_jitter,
                     ymin = pct,
                     ymax = with(ver_df_hist_all, pct)),
                 linetype = 3, 
                 alpha=0.5) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 19, hjust=0.5),
        axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("Time since divergence (my)") +
  xlab("") +
  ylab("")+
  xlim(0,500) +
  ggtitle("human")



dros_df_hist_synt <- sp1 %>%
  filter(tag=="fruitfly") %>%
  mutate(pct = (no_match/total),
         sep = sqrt(pct*(1-pct)/total)) %>%
  arrange(div)

dros_orphan <- filter(dros_orphan , species != "Solenopsis_invicta", species != "Caenorhabditis_elegans")
dros_df_hist_all <- as.data.frame(table(dros_orphan$species)) %>%
  mutate(pct = ((Freq)/dros_div_d$total_genes_checked),
         sep = sqrt(pct*(1-pct)/dros_div_d$total_genes_checked)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = dros_div_d[match(species, dros_div_d$species), "Divergence_time"]) %>%
  arrange(div)

dros_div_jitter <- jitter(dros_df_hist_synt$div, factor=5)


dros_all_pl <- ggplot() + 
  geom_point(data = dros_df_hist_synt, 
             aes(x=dros_div_jitter, 
                 y=pct),
             size=3,
             colour="#E69F00") +
  geom_smooth(data = dros_df_hist_synt, 
              aes(x=dros_div_jitter, 
                  y=pct),
              method="lm",
              colour="#E69F00") +
  geom_point(data = dros_df_hist_all, 
             aes(x=dros_div_jitter, 
                 y=pct),
             size=3,
             colour="#E69F00",
             shape=2) +
  # geom_text(data = cer_df_tbl_opt,
  #           aes(x=cer_div_jitter,
  #               y=(no_match-found_tblastn)/total,
  #               label=species),
  #           hjust=-0.2,
  #           vjust=1) +
  geom_linerange(data = dros_df_hist_synt,
                 aes(x = dros_div_jitter,
                     ymin = pct,
                     ymax = with(dros_df_hist_all, pct)),
                 linetype = 3, 
                 alpha=0.5) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 19, hjust=0.5),
        axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("Time since divergence (my)") +
  xlab("") +
  ylab("") +
  xlim(0,400) +
  ggtitle("fruitfly")



pl_all <- grid.arrange(
  cer_all_pl,
  dros_all_pl,
  ver_all_pl,
  nrow=1)

ggsave(plot = pl_all, "scripts_submission/Figure_5C.pdf", width=13, height = 5)



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

cer_df_c$species <- factor(cer_df_c$species, levels = unique(cer_df_c$species), ordered = TRUE)

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

pl_all <- grid.arrange(
  cer_hist_sep_pl,
  dros_hist_sep_pl,
  ver_hist_sep_pl,
  nrow=1)

ggsave(plot = pl_all, "scripts_submission/Figure_5B1.pdf", width=13, height = 3)


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

ver_df_c$species <- factor(ver_df_c$species, levels = unique(ver_df_c$species), ordered = TRUE)

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

dros_df_c$species <- factor(dros_df_c$species, levels = unique(dros_df_c$species), ordered = TRUE)

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


pl_all <- grid.arrange(
  hist_cer,
  hist_dros,
  hist_ver,
  nrow=1)

ggsave(plot = pl_all, "scripts_submission/Figure_5B2.pdf", width=13, height = 3)

complete_fig <- grid.arrange(cer_all_pl,
                             dros_all_pl,
                             ver_all_pl,
                             cer_hist_sep_pl,
                             dros_hist_sep_pl,
                             ver_hist_sep_pl,
                             hist_cer,
                             hist_dros,
                             hist_ver,
                             heights=c(2.8/6,1.6/6,1.6/6),
                             ncol=3,
                             nrow=3)

ggsave(plot = complete_fig, "scripts_submission/final_scripts/Fig4_compl.pdf", width=13, height = 10)

##### stats for pairwise ######

all_3 <- bind_rows(ver_df_hist_synt_h, dros_df_hist_synt_h, cer_df_hist_synt_h)
ov_min <- all_3[which.min(all_3$pct),]
ov_max <- all_3[which.max(all_3$pct),]
avg <- mean(all_3$pct)

###############

##################
###################
#######################

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

pl_sep_all <- grid.arrange(
  cer_hist_sep_pl,
  dros_hist_sep_pl,
  ver_hist_sep_pl,
  nrow=1)

ggsave(plot = pl_sep_all, "scripts_submission/Figure_supp_5B1.pdf", width=13, height = 3)


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


pl_all <- grid.arrange(
  hist_cer,
  hist_dros,
  hist_ver,
  nrow=1)

ggsave(plot = pl_all, "scripts_submission/supp_Figure_5B2.pdf", width=13, height = 3)


#####

cer_all_pl <- ggplot() + 
  geom_point(data = cer_phylo_df_p1, 
             aes(x=div, 
                 y=synt_missed),
             size=3,
             colour="#999999") +
  geom_smooth(data = cer_phylo_df_p1, 
              aes(x=div, 
                  y=synt_missed),
              method="lm",
              colour="#999999") +
  geom_point(data = cer_phylo_df_p1, 
             aes(x=div, 
                 y=all_missed),
             size=3,
             shape=2,
             colour="#999999") +
  # geom_text(data = cer_df_tbl_opt,
  #           aes(x=cer_div_jitter,
  #               y=(no_match-found_tblastn)/total,
  #               label=species),
  #           hjust=-0.2,
  #           vjust=1) +
  geom_linerange(data = cer_phylo_df_p1,
                 aes(x = div,
                     ymin = synt_missed,
                     ymax = all_missed),
                 linetype = 3, 
                 alpha=0.5) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 19, hjust=0.5),
        axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("Time since divergence (my)") + ylab("Proportion without similarity\n(out of total genes)") +
  xlab("") + ylab("") +
  xlim(0,130) +
  ggtitle("yeast")


dros_all_pl <- ggplot() + 
  geom_point(data = dros_phylo_df_p1, 
             aes(x=div, 
                 y=synt_missed),
             size=3,
             colour="#E69F00") +
  geom_smooth(data = dros_phylo_df_p1, 
              aes(x=div, 
                  y=synt_missed),
              method="lm",
              colour="#E69F00") +
  geom_point(data = dros_phylo_df_p1, 
             aes(x=div, 
                 y=all_missed),
             size=3,
             shape=2,
             colour="#E69F00") +
  # geom_text(data = dros_df_tbl_opt,
  #           aes(x=dros_div_jitter,
  #               y=(no_match-found_tblastn)/total,
  #               label=species),
  #           hjust=-0.2,
  #           vjust=1) +
  geom_linerange(data = dros_phylo_df_p1,
                 aes(x = div,
                     ymin = synt_missed,
                     ymax = all_missed),
                 linetype = 3, 
                 alpha=0.5) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 19, hjust=0.5),
        axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("Time since divergence (my)") +
  ylab("") +
  xlab("") +
  xlim(0,400) +
  ggtitle("fruitfly")



vert_all_pl <- ggplot() + 
  geom_point(data = ver_phylo_df_p1, 
             aes(x=div, 
                 y=synt_missed),
             size=3,
             colour="#56B4E9") +
  geom_smooth(data = ver_phylo_df_p1, 
              aes(x=div, 
                  y=synt_missed),
              method="lm",
              colour="#56B4E9") +
  geom_point(data = ver_phylo_df_p1, 
             aes(x=div, 
                 y=all_missed),
             size=3,
             shape=2,
             colour="#56B4E9") +
  # geom_text(data = vert_df_tbl_opt,
  #           aes(x=vert_div_jitter,
  #               y=(no_match-found_tblastn)/total,
  #               label=species),
  #           hjust=-0.2,
  #           vjust=1) +
  geom_linerange(data = ver_phylo_df_p1,
                 aes(x = div,
                     ymin = synt_missed,
                     ymax = all_missed),
                 linetype = 3, 
                 alpha=0.5) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 19, hjust=0.5),
        axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "none") +
  xlab("Time since divergence (my)") +
  xlab("") +
  ylab("")+
  xlim(0,500) +
  ggtitle("human")

pl_all <- grid.arrange(
  cer_all_pl,
  dros_all_pl,
  vert_all_pl,
  nrow=1)

ggsave(plot = pl_all, "scripts_submission/Figure_5C_alt.pdf", width=13, height = 5)

complete_fig <- grid.arrange(cer_all_pl,
                             dros_all_pl,
                             vert_all_pl,
                             cer_hist_sep_pl,
                             dros_hist_sep_pl,
                             ver_hist_sep_pl,
                             hist_cer,
                             hist_dros,
                             hist_ver,
                             heights=c(2.8/6,1.6/6,1.6/6),
                             ncol=3,
                             nrow=3)

ggsave(plot = complete_fig, "scripts_submission/final_scripts/Supp_fig6_compl.pdf", width=13, height = 10)

##### stats for phylogen ######

all_3 <- bind_rows(ver_phylo_df_p1, dros_phylo_df_p1, cer_phylo_df_p1)
ov_min <- all_3[which.min(all_3$pct),]
ov_max <- all_3[which.max(all_3$pct),]
avg <- mean(all_3$pct)

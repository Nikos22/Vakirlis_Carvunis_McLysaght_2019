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

cer_gp_oo_old <- cer_gp_oo
cer_gp_oo_old$species <- as.character(cer_gp_oo_old$species)

cer_gp_oo <- read.table("synt_relaxed_data/yeast_synt_relaxed_genes.txt",
                        stringsAsFactors = FALSE,
                        header=TRUE)
cer_gp_oo_str <- read.table("synt_relaxed_data/yeast_synt_stringent_genes.txt",
                        stringsAsFactors = FALSE,
                        header=TRUE)

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

vert_gp_oo_old <- vert_gp_oo
vert_gp_oo_old$species <- as.character(vert_gp_oo_old$species)

vert_gp_oo <- read.table("synt_relaxed_data/human_synt_relaxed_genes.txt",
                        stringsAsFactors = FALSE,
                        header=TRUE)
vert_gp_oo_str <- read.table("synt_relaxed_data/human_synt_stringent_genes.txt",
                            stringsAsFactors = FALSE,
                            header=TRUE)


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

dros_gp_oo_old <- dros_gp_oo
dros_gp_oo_old$species <- as.character(dros_gp_oo_old$species)

dros_gp_oo <- read.table("synt_relaxed_data/fly_synt_relaxed_genes.txt",
                        stringsAsFactors = FALSE,
                        header=TRUE)
dros_gp_oo_str <- read.table("synt_relaxed_data/fly_synt_stringent_genes.txt",
                            stringsAsFactors = FALSE,
                            header=TRUE)

########################################################
#### Calculate phylogeny-based proportions, relaxed ####
########################################################

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

##########################################################
#### Calculate phylogeny-based proportions, stringent ####
##########################################################

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


pcts_cer_cons_str <- c()
pcts_cer_cons_num_str <- c()

for (i in unique(cer_div$V2))
{
  species_to_get <- cer_div[which(cer_div$V2 >= i), "species"]
  genes_with_further <- unique(cer_gp_oo_str[which(cer_gp_oo_str$species %in% species_to_get), "Gene_focal"])
  
  cer_orphan_p <- cer_orphan %>%
    group_by(Gene_focal) %>%
    filter(all(species_to_get %in% species)) %>%
    summarise(no_sp = length(unique(species)))
  
  pcts_cer_cons_str <- c(pcts_cer_cons_str,length(unique(which(genes_with_further %in% cer_orphan_p$Gene_focal)))/length(unique(genes_with_further)))
  pcts_cer_cons_num_str <- c(pcts_cer_cons_num_str, length(unique(genes_with_further)))
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


pcts_vert_cons_str <- c()
pcts_vert_cons_num_str <- c()
for (i in unique(vert_div$V2))
{
  
  species_to_get <- vert_div[which(vert_div$V2 >= i), "species"]
  genes_with_further <- unique(vert_gp_oo_str[which(vert_gp_oo_str$species %in% species_to_get), "Gene_focal"])
  
  vert_orphan_p <- vert_orphan %>%
    group_by(Gene_focal) %>%
    filter(all(species_to_get %in% species)) %>%
    summarise(no_sp = length(unique(species)))
  
  pcts_vert_cons_str <- c(pcts_vert_cons_str,length(unique(which(genes_with_further %in% vert_orphan_p$Gene_focal)))/length(unique(genes_with_further)))
  pcts_vert_cons_num_str <- c(pcts_vert_cons_num_str, length(unique(genes_with_further)))
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


pcts_dros_cons_str <- c()
pcts_dros_cons_num_str <- c()
for (i in unique(dros_div$V5))
{
  
  species_to_get <- dros_div[which(dros_div$V5 >= i), "species"]
  genes_with_further <- unique(dros_gp_oo_str[which(dros_gp_oo_str$species %in% species_to_get), "Gene_focal"])
  
  dros_orphan_p <- dros_orphan %>%
    group_by(Gene_focal) %>%
    filter(all(species_to_get %in% species)) %>%
    summarise(no_sp = length(unique(species)))
  
  pcts_dros_cons_str <- c(pcts_dros_cons_str,length(unique(which(genes_with_further %in% dros_orphan_p$Gene_focal)))/length(unique(genes_with_further)))
  pcts_dros_cons_num_str <- c(pcts_dros_cons_num_str, length(unique(genes_with_further)))
}


##########################################################
#### Calculate phylogeny-based proportions, current ####
##########################################################

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


pcts_cer_cons_old <- c()
pcts_cer_cons_num_old <- c()

for (i in unique(cer_div$V2))
{
  species_to_get <- cer_div[which(cer_div$V2 >= i), "species"]
  genes_with_further <- unique(cer_gp_oo_old[which(cer_gp_oo_old$species %in% species_to_get), "Gene_focal"])
  
  cer_orphan_p <- cer_orphan %>%
    group_by(Gene_focal) %>%
    filter(all(species_to_get %in% species)) %>%
    summarise(no_sp = length(unique(species)))
  
  pcts_cer_cons_old <- c(pcts_cer_cons_old,length(unique(which(genes_with_further %in% cer_orphan_p$Gene_focal)))/length(unique(genes_with_further)))
  pcts_cer_cons_num_old <- c(pcts_cer_cons_num_old, length(unique(genes_with_further)))
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


pcts_vert_cons_old <- c()
pcts_vert_cons_num_old <- c()
for (i in unique(vert_div$V2))
{
  
  species_to_get <- vert_div[which(vert_div$V2 >= i), "species"]
  genes_with_further <- unique(vert_gp_oo_old[which(vert_gp_oo_old$species %in% species_to_get), "Gene_focal"])
  
  vert_orphan_p <- vert_orphan %>%
    group_by(Gene_focal) %>%
    filter(all(species_to_get %in% species)) %>%
    summarise(no_sp = length(unique(species)))
  
  pcts_vert_cons_old <- c(pcts_vert_cons_old,length(unique(which(genes_with_further %in% vert_orphan_p$Gene_focal)))/length(unique(genes_with_further)))
  pcts_vert_cons_num_old <- c(pcts_vert_cons_num_old, length(unique(genes_with_further)))
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


pcts_dros_cons_old <- c()
pcts_dros_cons_num_old <- c()
for (i in unique(dros_div$V5))
{
  
  species_to_get <- dros_div[which(dros_div$V5 >= i), "species"]
  genes_with_further <- unique(dros_gp_oo_old[which(dros_gp_oo_old$species %in% species_to_get), "Gene_focal"])
  
  dros_orphan_p <- dros_orphan %>%
    group_by(Gene_focal) %>%
    filter(all(species_to_get %in% species)) %>%
    summarise(no_sp = length(unique(species)))
  
  pcts_dros_cons_old <- c(pcts_dros_cons_old,length(unique(which(genes_with_further %in% dros_orphan_p$Gene_focal)))/length(unique(genes_with_further)))
  pcts_dros_cons_num_old <- c(pcts_dros_cons_num_old, length(unique(genes_with_further)))
}

##########################################################################
##### create data frames and plot pairwise proportions first, relaxed ####
##########################################################################


cer_orphan_synt <- cer_orphan %>% filter(paste(Gene_focal, species) %in% paste(cer_gp_oo$Gene_focal, cer_gp_oo$species))
cer_orphan_synt$species <- as.character(cer_orphan_synt$species)

cer_df_hist_synt <- as.data.frame(table(cer_orphan_synt$species)) %>%
  mutate(pct = ((Freq)/as.data.frame(table(cer_gp_oo$species))$Freq),
         sep = sqrt(pct*(1-pct)/as.data.frame(table(cer_gp_oo$species))$Freq)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = cer_div_d[match(species, cer_div_d$species), "Divergence_time"]) %>%
  arrange(div)



cer_df_hist_all <- as.data.frame(table(cer_orphan$species)) %>%
  mutate(pct = ((Freq)/cer_div_d$total_genes_checked),
         sep = sqrt(pct*(1-pct)/cer_div_d$total_genes_checked)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = cer_div_d[match(species, cer_div_d$species), "Divergence_time"]) %>%
  arrange(div)


cer_df_hist_all <- cer_df_hist_all[match(cer_df_hist_synt$species, cer_df_hist_all$species), ]

cer_div_jitter <- jitter(cer_df_hist_synt$div, factor=5)




vert_orphan_synt <- vert_orphan %>% filter(paste(Gene_focal, species) %in% paste(vert_gp_oo$Gene_focal, vert_gp_oo$species))
vert_orphan_synt$species <- as.character(vert_orphan_synt$species)

ver_df_hist_synt <- as.data.frame(table(vert_orphan_synt$species)) %>%
  mutate(pct = ((Freq)/as.data.frame(table(vert_gp_oo$species))$Freq),
         sep = sqrt(pct*(1-pct)/as.data.frame(table(vert_gp_oo$species))$Freq)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = vert_div_d[match(species, vert_div_d$species), "Divergence_time"]) %>%
  arrange(div)


ver_df_hist_all <- as.data.frame(table(vert_orphan$species)) %>%
  mutate(pct = ((Freq)/vert_div_d$total_genes_checked),
         sep = sqrt((pct*(1-pct))/vert_div_d$total_genes_checked)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = vert_div_d[match(species, vert_div_d$species), "Divergence_time"]) %>%
  arrange(div)

ver_df_hist_all <- ver_df_hist_all[match(ver_df_hist_synt$species, ver_df_hist_all$species), ]


dros_orphan <- filter(dros_orphan , species != "Solenopsis_invicta", species != "Caenorhabditis_elegans")
dros_gp_oo <- filter(dros_gp_oo, species != "Caenorhabditis_elegans")
dros_orphan_synt <- dros_orphan %>% filter(paste(Gene_focal, species) %in% paste(dros_gp_oo$Gene_focal, dros_gp_oo$species))
dros_orphan_synt$species <- as.character(dros_orphan_synt$species)


dros_df_hist_synt <- as.data.frame(table(dros_orphan_synt$species)) %>%
  mutate(pct = ((Freq)/as.data.frame(table(dros_gp_oo$species))$Freq),
         sep = sqrt(pct*(1-pct)/as.data.frame(table(dros_gp_oo$species))$Freq)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = dros_div_d[match(species, dros_div_d$species), "Divergence_time"]) %>%
  arrange(div)

dros_df_hist_all <- as.data.frame(table(dros_orphan$species)) %>%
  mutate(pct = ((Freq)/dros_div_d$total_genes_checked),
         sep = sqrt((pct*(1-pct))/dros_div_d$total_genes_checked)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = dros_div_d[match(species, dros_div_d$species), "Divergence_time"]) %>%
  arrange(div)

dros_df_hist_all <- dros_df_hist_all[match(dros_df_hist_synt$species, dros_df_hist_all$species), ]
dros_div_jitter <- jitter(dros_df_hist_synt$div, factor=5)


############################################################################
##### create data frames and plot pairwise proportions first, stringent ####
############################################################################


cer_orphan_synt <- cer_orphan %>% filter(paste(Gene_focal, species) %in% paste(cer_gp_oo_str$Gene_focal, cer_gp_oo_str$species))
cer_orphan_synt$species <- as.character(cer_orphan_synt$species)
temp <- data.frame(Var1=unique(cer_gp_oo_str$species)[which(!unique(cer_gp_oo_str$species) %in% unique(cer_orphan_synt$species))], Freq=0)

cer_df_hist_synt_str <- rbind(as.data.frame(table(cer_orphan_synt$species)), temp) %>%
  arrange(as.character(Var1)) %>%
  mutate(pct = ((Freq)/as.data.frame(table(cer_gp_oo_str$species))$Freq),
         sep = sqrt(pct*(1-pct)/as.data.frame(table(cer_gp_oo_str$species))$Freq)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = cer_div_d[match(species, cer_div_d$species), "Divergence_time"]) %>%
  arrange(div)


vert_orphan_synt <- vert_orphan %>% filter(paste(Gene_focal, species) %in% paste(vert_gp_oo_str$Gene_focal, vert_gp_oo_str$species))
vert_orphan_synt$species <- as.character(vert_orphan_synt$species)
temp <- data.frame(Var1=unique(vert_gp_oo_str$species)[which(!unique(vert_gp_oo_str$species) %in% unique(vert_orphan_synt$species))], Freq=0)


ver_df_hist_synt_str <- rbind(as.data.frame(table(vert_orphan_synt$species)), temp) %>%
  arrange(as.character(Var1)) %>%
  mutate(pct = ((Freq)/as.data.frame(table(vert_gp_oo_str$species))$Freq),
         sep = sqrt(pct*(1-pct)/as.data.frame(table(vert_gp_oo_str$species))$Freq)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = vert_div_d[match(species, vert_div_d$species), "Divergence_time"]) %>%
  arrange(div)



dros_orphan <- filter(dros_orphan , species != "Solenopsis_invicta", species != "Caenorhabditis_elegans")
dros_gp_oo_str <- filter(dros_gp_oo_str, species != "Caenorhabditis_elegans")
dros_orphan_synt <- dros_orphan %>% filter(paste(Gene_focal, species) %in% paste(dros_gp_oo_str$Gene_focal, dros_gp_oo_str$species))
dros_orphan_synt$species <- as.character(dros_orphan_synt$species)
temp <- data.frame(Var1=unique(dros_gp_oo_str$species)[which(!unique(dros_gp_oo_str$species) %in% unique(dros_orphan_synt$species))], Freq=0)


dros_df_hist_synt_str <- rbind(as.data.frame(table(dros_orphan_synt$species)), temp) %>%
  arrange(as.character(Var1)) %>%
  mutate(pct = ((Freq)/as.data.frame(table(dros_gp_oo_str$species))$Freq),
         sep = sqrt(pct*(1-pct)/as.data.frame(table(dros_gp_oo_str$species))$Freq)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = dros_div_d[match(species, dros_div_d$species), "Divergence_time"]) %>%
  arrange(div)


############################################################################
##### create data frames and plot pairwise proportions first, current ####
############################################################################


cer_orphan_synt <- cer_orphan %>% filter(paste(Gene_focal, species) %in% paste(cer_gp_oo_old$Gene_focal, cer_gp_oo_old$species))
cer_orphan_synt$species <- as.character(cer_orphan_synt$species)

cer_df_hist_synt_old <- as.data.frame(table(cer_orphan_synt$species)) %>%
  mutate(pct = ((Freq)/as.data.frame(table(cer_gp_oo_old$species))$Freq),
         sep = sqrt(pct*(1-pct)/as.data.frame(table(cer_gp_oo_old$species))$Freq)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = cer_div_d[match(species, cer_div_d$species), "Divergence_time"]) %>%
  arrange(div)


vert_orphan_synt <- vert_orphan %>% filter(paste(Gene_focal, species) %in% paste(vert_gp_oo_old$Gene_focal, vert_gp_oo_old$species))
vert_orphan_synt$species <- as.character(vert_orphan_synt$species)


ver_df_hist_synt_old <- as.data.frame(table(vert_orphan_synt$species)) %>%
  mutate(pct = ((Freq)/as.data.frame(table(vert_gp_oo_old$species))$Freq),
         sep = sqrt(pct*(1-pct)/as.data.frame(table(vert_gp_oo_old$species))$Freq)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = vert_div_d[match(species, vert_div_d$species), "Divergence_time"]) %>%
  arrange(div)



dros_orphan <- filter(dros_orphan , species != "Solenopsis_invicta", species != "Caenorhabditis_elegans")
dros_gp_oo_old <- filter(dros_gp_oo_old, species != "Caenorhabditis_elegans")
dros_orphan_synt <- dros_orphan %>% filter(paste(Gene_focal, species) %in% paste(dros_gp_oo_old$Gene_focal, dros_gp_oo_old$species))
dros_orphan_synt$species <- as.character(dros_orphan_synt$species)
temp <- data.frame(Var1=unique(dros_gp_oo_old$species)[which(!unique(dros_gp_oo_old$species) %in% unique(dros_orphan_synt$species))], Freq=0)


dros_df_hist_synt_old <- rbind(as.data.frame(table(dros_orphan_synt$species)), temp) %>%
  arrange(as.character(Var1)) %>%
  mutate(pct = ((Freq)/as.data.frame(table(dros_gp_oo_old$species))$Freq),
         sep = sqrt(pct*(1-pct)/as.data.frame(table(dros_gp_oo_old$species))$Freq)) %>%
  dplyr::rename(species=Var1) %>%
  mutate(div = dros_div_d[match(species, dros_div_d$species), "Divergence_time"]) %>%
  arrange(div)

#######

### relaxed ###
ver_df_hist_synt$tag <- "With undetected homologues"
ver_df_hist_synt_h <- select(ver_df_hist_synt, pct, sep, div, tag, species)
ver_df_hist_synt_h$pct <- ver_df_hist_synt_h$pct/ver_df_hist_all$pct
ver_df_hist_synt_h$species <- as.character(ver_df_hist_synt_h$species)
ver_df_hist_synt_h$species <- paste(substr(ver_df_hist_synt_h$species, 1, 1), substr(unlist(strsplit(ver_df_hist_synt_h$species, split = "_"))[2*(1:length(ver_df_hist_synt_h$species))],1,3), sep='')
ver_df_hist_synt_h$species <- factor(ver_df_hist_synt_h$species, levels = sp1[which(sp1$tag=="human"), "speciesShort"], ordered = TRUE)
ver_df_hist_synt_h$type="relaxed"

### stringent ###
ver_df_hist_synt_str$tag <- "With undetected homologues"
ver_df_hist_synt_str_h <- select(ver_df_hist_synt_str, pct, sep, div, tag, species)
ver_df_hist_synt_str_h$pct <- ver_df_hist_synt_str_h$pct/ver_df_hist_all$pct
ver_df_hist_synt_str_h$species <- as.character(ver_df_hist_synt_str_h$species)
ver_df_hist_synt_str_h$species <- paste(substr(ver_df_hist_synt_str_h$species, 1, 1), substr(unlist(strsplit(ver_df_hist_synt_str_h$species, split = "_"))[2*(1:length(ver_df_hist_synt_str_h$species))],1,3), sep='')
ver_df_hist_synt_str_h$species <- factor(ver_df_hist_synt_str_h$species, levels = sp1[which(sp1$tag=="human"), "speciesShort"], ordered = TRUE)
ver_df_hist_synt_str_h$type="stringent"

### current ###
ver_df_hist_synt_old$tag <- "With undetected homologues"
ver_df_hist_synt_old_h <- select(ver_df_hist_synt_old, pct, sep, div, tag, species)
ver_df_hist_synt_old_h$pct <- ver_df_hist_synt_old_h$pct/ver_df_hist_all$pct
ver_df_hist_synt_old_h$species <- as.character(ver_df_hist_synt_old_h$species)
ver_df_hist_synt_old_h$species <- paste(substr(ver_df_hist_synt_old_h$species, 1, 1), substr(unlist(strsplit(ver_df_hist_synt_old_h$species, split = "_"))[2*(1:length(ver_df_hist_synt_old_h$species))],1,3), sep='')
ver_df_hist_synt_old_h$species <- factor(ver_df_hist_synt_old_h$species, levels = sp1[which(sp1$tag=="human"), "speciesShort"], ordered = TRUE)

ver_df_hist_synt_old_h$type="current"

ver_df_c <- bind_rows(ver_df_hist_synt_old_h,
                      ver_df_hist_synt_str_h,
                      ver_df_hist_synt_h)
ver_df_c$type <- factor(ver_df_c$type, levels=c("relaxed", "current", "stringent"), ordered = TRUE)


hist_ver <- ggplot() +
  geom_bar(data = ver_df_c, 
           aes(x=species, 
               y=pct,
               fill=type,
               alpha=type),
           stat="identity",
           position = position_dodge(), 
           colour="black") +
  scale_fill_manual(values = c("skyblue1", "dodgerblue2", "blue3")) +
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


### relaxed ###
dros_df_hist_synt$tag <- "With undetected homologues"
dros_df_hist_synt_h <- select(dros_df_hist_synt, pct, sep, div, tag, species)
dros_df_hist_synt_h$pct <- dros_df_hist_synt_h$pct/dros_df_hist_all$pct
dros_df_hist_synt_h$species <- as.character(dros_df_hist_synt_h$species)
dros_df_hist_synt_h$species <- paste(substr(dros_df_hist_synt_h$species, 1, 1), substr(unlist(strsplit(dros_df_hist_synt_h$species, split = "_"))[2*(1:length(dros_df_hist_synt_h$species))],1,3), sep='')
dros_df_hist_synt_h$species <- factor(dros_df_hist_synt_h$species, levels = sp1[which(sp1$tag=="fruitfly"), "speciesShort"], ordered = TRUE)
dros_df_hist_synt_h$type="relaxed"

### stringent ###
dros_df_hist_synt_str$tag <- "With undetected homologues"
dros_df_hist_synt_str_h <- select(dros_df_hist_synt_str, pct, sep, div, tag, species)
dros_df_hist_synt_str_h$pct <- dros_df_hist_synt_str_h$pct/dros_df_hist_all$pct
dros_df_hist_synt_str_h$species <- as.character(dros_df_hist_synt_str_h$species)
dros_df_hist_synt_str_h$species <- paste(substr(dros_df_hist_synt_str_h$species, 1, 1), substr(unlist(strsplit(dros_df_hist_synt_str_h$species, split = "_"))[2*(1:length(dros_df_hist_synt_str_h$species))],1,3), sep='')
dros_df_hist_synt_str_h$species <- factor(dros_df_hist_synt_str_h$species, levels = sp1[which(sp1$tag=="fruitfly"), "speciesShort"], ordered = TRUE)
dros_df_hist_synt_str_h$type="stringent"

### current ###
dros_df_hist_synt_old$tag <- "With undetected homologues"
dros_df_hist_synt_old_h <- select(dros_df_hist_synt_old, pct, sep, div, tag, species)
dros_df_hist_synt_old_h$pct <- dros_df_hist_synt_old_h$pct/dros_df_hist_all$pct
dros_df_hist_synt_old_h$species <- as.character(dros_df_hist_synt_old_h$species)
dros_df_hist_synt_old_h$species <- paste(substr(dros_df_hist_synt_old_h$species, 1, 1), substr(unlist(strsplit(dros_df_hist_synt_old_h$species, split = "_"))[2*(1:length(dros_df_hist_synt_old_h$species))],1,3), sep='')
dros_df_hist_synt_old_h$species <- factor(dros_df_hist_synt_old_h$species, levels = sp1[which(sp1$tag=="fruitfly"), "speciesShort"], ordered = TRUE)

dros_df_hist_synt_old_h$type="current"

dros_df_c <- bind_rows(dros_df_hist_synt_old_h,
                      dros_df_hist_synt_str_h,
                      dros_df_hist_synt_h)
dros_df_c$type <- factor(dros_df_c$type, levels=c("relaxed", "current", "stringent"), ordered = TRUE)



hist_dros <- ggplot() +
  geom_bar(data = dros_df_c, 
           aes(x=species, 
               y=pct,
               fill=type,
               alpha=type),
           stat="identity",
           position = position_dodge(), 
           colour="black") +
  scale_fill_manual(values = c("tan1", "darkorange", "darkorange4")) +
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


###
### relaxed ###
cer_df_hist_synt$tag <- "With undetected homologues"
cer_df_hist_synt_h <- select(cer_df_hist_synt, pct, sep, div, tag, species)
cer_df_hist_synt_h$pct <- cer_df_hist_synt_h$pct/cer_df_hist_all$pct
cer_df_hist_synt_h$species <- as.character(cer_df_hist_synt_h$species)
cer_df_hist_synt_h$species <- paste(substr(cer_df_hist_synt_h$species, 1, 1), substr(unlist(strsplit(cer_df_hist_synt_h$species, split = "_"))[2*(1:length(cer_df_hist_synt_h$species))],1,3), sep='')
cer_df_hist_synt_h$species <- factor(cer_df_hist_synt_h$species, levels = sp1[which(sp1$tag=="yeast"), "speciesShort"], ordered = TRUE)
cer_df_hist_synt_h$type="relaxed"

### stringent ###
cer_df_hist_synt_str$tag <- "With undetected homologues"
cer_df_hist_synt_str_h <- select(cer_df_hist_synt_str, pct, sep, div, tag, species)
cer_df_hist_synt_str_h$pct <- cer_df_hist_synt_str_h$pct/cer_df_hist_all$pct
cer_df_hist_synt_str_h$species <- as.character(cer_df_hist_synt_str_h$species)
cer_df_hist_synt_str_h$species <- paste(substr(cer_df_hist_synt_str_h$species, 1, 1), substr(unlist(strsplit(cer_df_hist_synt_str_h$species, split = "_"))[2*(1:length(cer_df_hist_synt_str_h$species))],1,3), sep='')
cer_df_hist_synt_str_h$species <- factor(cer_df_hist_synt_str_h$species, levels = sp1[which(sp1$tag=="yeast"), "speciesShort"], ordered = TRUE)
cer_df_hist_synt_str_h$type="stringent"

### current ###
cer_df_hist_synt_old$tag <- "With undetected homologues"
cer_df_hist_synt_old_h <- select(cer_df_hist_synt_old, pct, sep, div, tag, species)
cer_df_hist_synt_old_h$pct <- cer_df_hist_synt_old_h$pct/cer_df_hist_all$pct
cer_df_hist_synt_old_h$species <- as.character(cer_df_hist_synt_old_h$species)
cer_df_hist_synt_old_h$species <- paste(substr(cer_df_hist_synt_old_h$species, 1, 1), substr(unlist(strsplit(cer_df_hist_synt_old_h$species, split = "_"))[2*(1:length(cer_df_hist_synt_old_h$species))],1,3), sep='')
cer_df_hist_synt_old_h$species <- factor(cer_df_hist_synt_old_h$species, levels = sp1[which(sp1$tag=="yeast"), "speciesShort"], ordered = TRUE)

cer_df_hist_synt_old_h$type="current"

cer_df_c <- bind_rows(cer_df_hist_synt_old_h,
                       cer_df_hist_synt_str_h,
                       cer_df_hist_synt_h)
cer_df_c$type <- factor(cer_df_c$type, levels=c("relaxed", "current", "stringent"), ordered = TRUE)


hist_cer <- ggplot() +
  geom_bar(data = cer_df_c, 
           aes(x=species, 
               y=pct,
               fill=type),
           stat="identity",
           position = position_dodge(), 
           colour="black") +
  scale_fill_manual(values = c("grey78", "grey41", "grey16")) +
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


# comparing relaxed to current

summary(c(cer_df_hist_synt_old_h$pct, dros_df_hist_synt_old_h$pct, ver_df_hist_synt_old_h$pct))
summary(c(cer_df_hist_synt_h$pct, dros_df_hist_synt_h$pct, ver_df_hist_synt_h$pct))

avg_change <- mean(abs(c(cer_df_hist_synt_old_h$pct, dros_df_hist_synt_old_h$pct, ver_df_hist_synt_old_h$pct)-
                        c(cer_df_hist_synt_h$pct, dros_df_hist_synt_h$pct, ver_df_hist_synt_h$pct)))
                       

sp_with_more_300_in_string <- rbind(as.data.frame(table(cer_gp_oo_str$species)),
                                    as.data.frame(table(dros_gp_oo_str$species)),
                                    as.data.frame(table(vert_gp_oo_str$species))) %>%
                              filter(Freq>200) %>%
                              select(Var1) %>%
                              rename(species=Var1)

sp_with_more_300_in_string$species <- as.character(sp_with_more_300_in_string$species)
sp_with_more_300_in_string$species <- paste(substr(sp_with_more_300_in_string$species, 1, 1), substr(unlist(strsplit(sp_with_more_300_in_string$species, split = "_"))[2*(1:length(sp_with_more_300_in_string$species))],1,3), sep='')


all_cur <- bind_rows(cer_df_hist_synt_old_h,
                     dros_df_hist_synt_old_h,
                     ver_df_hist_synt_old_h) %>%
          filter(species %in% sp_with_more_300_in_string$species)

all_str <- bind_rows(cer_df_hist_synt_str_h,
                     dros_df_hist_synt_str_h,
                     ver_df_hist_synt_str_h) %>%
  filter(species %in% sp_with_more_300_in_string$species)

summary(all_str$pct)
summary(all_cur$pct)


##### create data frames and plot phylogeny based proportions ####

ver_df_c2 <- data.frame(pct_found = pcts_vert_cons/pcts_vert_all,
                              tag = "With predicted undetectable homologues",
                              div = unique(vert_div$V2))
ver_df_c2$div <- factor(ver_df_c2$div, levels = unique(ver_df_c2$div), ordered = TRUE)
ver_df_c2$type="relaxed"

ver_df_c2_str <- data.frame(pct_found = pcts_vert_cons_str/pcts_vert_all,
                        tag = "With predicted undetectable homologues",
                        div = unique(vert_div$V2))
ver_df_c2_str$div <- factor(ver_df_c2_str$div, levels = unique(ver_df_c2_str$div), ordered = TRUE)
ver_df_c2_str$type="stringent"

ver_df_c2_old <- data.frame(pct_found = pcts_vert_cons_old/pcts_vert_all,
                            tag = "With predicted undetectable homologues",
                            div = unique(vert_div$V2))
ver_df_c2_old$div <- factor(ver_df_c2_old$div, levels = unique(ver_df_c2_old$div), ordered = TRUE)
ver_df_c2_old$type="current"

ver_df_c2_tot <- bind_rows(ver_df_c2, ver_df_c2_old, ver_df_c2_str)
ver_df_c2_tot$type <- factor(ver_df_c2_tot$type, levels=c("relaxed", "current", "stringent"), ordered = TRUE)




cer_df_c2 <- data.frame(pct_found = pcts_cer_cons/pcts_cer_all,
                        tag = "With predicted undetectable homologues",
                        div = unique(cer_div$V2))
cer_df_c2$div <- factor(cer_df_c2$div, levels = unique(cer_df_c2$div), ordered = TRUE)
cer_df_c2$type="relaxed"

cer_df_c2_str <- data.frame(pct_found = pcts_cer_cons_str/pcts_cer_all,
                            tag = "With predicted undetectable homologues",
                            div = unique(cer_div$V2))
cer_df_c2_str$div <- factor(cer_df_c2_str$div, levels = unique(cer_df_c2_str$div), ordered = TRUE)
cer_df_c2_str$type="stringent"

cer_df_c2_old <- data.frame(pct_found = pcts_cer_cons_old/pcts_cer_all,
                            tag = "With predicted undetectable homologues",
                            div = unique(cer_div$V2))
cer_df_c2_old$div <- factor(cer_df_c2_old$div, levels = unique(cer_df_c2_old$div), ordered = TRUE)
cer_df_c2_old$type="current"

cer_df_c2_tot <- bind_rows(cer_df_c2, cer_df_c2_old, cer_df_c2_str)
cer_df_c2_tot$type <- factor(cer_df_c2_tot$type, levels=c("relaxed", "current", "stringent"), ordered = TRUE)



dros_df_c2 <- data.frame(pct_found = pcts_dros_cons/pcts_dros_all,
                        tag = "With predicted undetectable homologues",
                        div = unique(dros_div$V5))
dros_df_c2$div <- factor(dros_df_c2$div, levels = unique(dros_df_c2$div), ordered = TRUE)
dros_df_c2$type="relaxed"

dros_df_c2_str <- data.frame(pct_found = pcts_dros_cons_str/pcts_dros_all,
                            tag = "With predicted undetectable homologues",
                            div = unique(dros_div$V5))
dros_df_c2_str$div <- factor(dros_df_c2_str$div, levels = unique(dros_df_c2_str$div), ordered = TRUE)
dros_df_c2_str$type="stringent"

dros_df_c2_old <- data.frame(pct_found = pcts_dros_cons_old/pcts_dros_all,
                            tag = "With predicted undetectable homologues",
                            div = unique(dros_div$V5))
dros_df_c2_old$div <- factor(dros_df_c2_old$div, levels = unique(dros_df_c2_old$div), ordered = TRUE)
dros_df_c2_old$type="current"

dros_df_c2_tot <- bind_rows(dros_df_c2, dros_df_c2_old, dros_df_c2_str)
dros_df_c2_tot$type <- factor(dros_df_c2_tot$type, levels=c("relaxed", "current", "stringent"), ordered = TRUE)


####


hist_ver_phyl <- ggplot() +
  geom_bar(data = ver_df_c2_tot, 
           aes(x=div, 
               y=pct_found,
               fill=type,
               alpha=type),
           stat="identity",
           position = position_dodge(), 
           colour="black") +
  scale_fill_manual(values = c("skyblue1", "dodgerblue2", "blue3")) +
  xlab("") +
  scale_x_discrete(labels=c(0,as.character(ver_df_c2$div)[1:(length(ver_df_c2$div)-1)])) +
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


hist_dros_phyl <- ggplot() +
  geom_bar(data = dros_df_c2_tot, 
           aes(x=div, 
               y=pct_found,
               fill=type,
               alpha=type),
           stat="identity",
           position = position_dodge(), 
           colour="black") +
  scale_fill_manual(values = c("tan1", "darkorange", "darkorange4")) +
  xlab("") +
  scale_x_discrete(labels=c(0,as.character(dros_df_c2$div)[1:(length(dros_df_c2$div)-1)])) +
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



hist_cer_phyl <- ggplot() +
  geom_bar(data = cer_df_c2_tot, 
           aes(x=div, 
               y=pct_found,
               fill=type,
               alpha=type),
           stat="identity",
           position = position_dodge(), 
           colour="black") +
  scale_fill_manual(values = c("grey78", "grey41", "grey16")) +
  xlab("") +
  scale_x_discrete(labels=c(0,as.character(cer_df_c2$div)[1:(length(cer_df_c2$div)-1)])) +
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

complete_fig <- grid.arrange(hist_cer,
                             hist_cer_phyl,
                             hist_dros,
                             hist_dros_phyl,
                             hist_ver,
                             hist_ver_phyl,
                             heights=c(0.3, 0.3, 0.3),
                             widths=c(0.65, 0.35),
                             ncol=2,
                             nrow=3)

ggsave(plot = complete_fig, "figures/Figure6-supplement_3A.pdf", width=9, height = 10)

# stats

summary(c(cer_df_c2_old$pct, dros_df_c2_old$pct, ver_df_c2_old$pct))
summary(c(cer_df_c2$pct, dros_df_c2$pct, ver_df_c2$pct))

summary(c(cer_df_c2_str$pct, dros_df_c2_str$pct[1:5], ver_df_c2_str$pct))
summary(c(cer_df_c2_old$pct, dros_df_c2_old$pct[1:5], ver_df_c2_old$pct))



#####

cer_cur <- as.data.frame(table(cer_gp_oo_old$species))
cer_cur$type="current"
cer_rel <- as.data.frame(table(cer_gp_oo$species))
cer_rel$type="relaxed"
cer_str <- as.data.frame(table(cer_gp_oo_str$species))
cer_str$type="stringent"

cer_df <- bind_rows(cer_cur, cer_rel, cer_str)
cer_df$type <- factor(cer_df$type, levels=c("relaxed", "current", "stringent"), ordered=TRUE)
cer_df$Var1 <- as.character(cer_df$Var1)
cer_df$Var1 <- paste(substr(cer_df$Var1, 1, 1), substr(unlist(strsplit(cer_df$Var1, split = "_"))[2*(1:length(cer_df$Var1))],1,3), sep='')
cer_df$Var1 <- factor(cer_df$Var1, levels = sp1[which(sp1$tag=="yeast"), "speciesShort"], ordered = TRUE)
cer_df$tag <- "yeast"
cer_df <- arrange(cer_df, Var1)

dros_cur <- as.data.frame(table(dros_gp_oo_old$species))
dros_cur$type="current"
dros_rel <- as.data.frame(table(dros_gp_oo$species))
dros_rel$type="relaxed"
dros_str <- as.data.frame(table(dros_gp_oo_str$species))
dros_str$type="stringent"

dros_df <- bind_rows(dros_cur, dros_rel, dros_str)
dros_df$type <- factor(dros_df$type, levels=c("relaxed", "current", "stringent"), ordered=TRUE)
dros_df$Var1 <- as.character(dros_df$Var1)
dros_df$Var1 <- paste(substr(dros_df$Var1, 1, 1), substr(unlist(strsplit(dros_df$Var1, split = "_"))[2*(1:length(dros_df$Var1))],1,3), sep='')
dros_df$Var1 <- factor(dros_df$Var1, levels = sp1[which(sp1$tag=="fruitfly"), "speciesShort"], ordered = TRUE)
dros_df$tag <- "fly"
dros_df <- arrange(dros_df, Var1)

vert_cur <- as.data.frame(table(vert_gp_oo_old$species))
vert_cur$type="current"
vert_rel <- as.data.frame(table(vert_gp_oo$species))
vert_rel$type="relaxed"
vert_str <- as.data.frame(table(vert_gp_oo_str$species))
vert_str$type="stringent"

vert_df <- bind_rows(vert_cur, vert_rel, vert_str)
vert_df$type <- factor(vert_df$type, levels=c("relaxed", "current", "stringent"), ordered=TRUE)
vert_df$Var1 <- as.character(vert_df$Var1)
vert_df$Var1 <- paste(substr(vert_df$Var1, 1, 1), substr(unlist(strsplit(vert_df$Var1, split = "_"))[2*(1:length(vert_df$Var1))],1,3), sep='')
vert_df$Var1 <- factor(vert_df$Var1, levels = sp1[which(sp1$tag=="human"), "speciesShort"], ordered = TRUE)
vert_df$tag <- "human"
vert_df <- arrange(vert_df, Var1)

all_df <- bind_rows(cer_df, dros_df, vert_df)
all_df$tag <- factor(all_df$tag, levels=c("yeast", "fly", "human"), ordered=TRUE)

all_df$Var1 <- factor(all_df$Var1, levels=unique(all_df$Var1), ordered=TRUE)


hist_plot <- ggplot() +
  geom_bar(data = all_df, 
           aes(x=Var1, 
               y=Freq,
               fill=tag,
               alpha=type),
           position="dodge",
           stat="identity") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  scale_alpha_manual(values=c(0.3,0.7,1)) +
  theme_bw() +
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
  facet_wrap(~tag, scales="free", dir = "h", drop=TRUE)

ggsave(plot = hist_plot, "figures/Figure6-supplement_3B.pdf", width=12, height = 4)

library(dplyr)
library(ggplot2)
library(mltools)
library(reshape2)
library(Hmisc)
library(tidyr)
library(grid)
library(gridExtra)
library(effsize)

setwd("~/Documents/research/Vakirlis_Carvunis_McLysaght_2019/")


##############
vert_div <- read.table("divergence_times/vert_div.txt", row.names=1)
vert_div$species <- rownames(vert_div)
vert_div <- vert_div[order(vert_div$V2),]

dros_div <- read.table("divergence_times/dros_div.txt", row.names=1)
dros_div <- dros_div[order(dros_div$V5),]

cer_div <- read.table("divergence_times/cer_div.txt", row.names=1)
cer_div$species <- rownames(cer_div)
cer_div <- cer_div[order(cer_div$V2),]

cer_gp_tbl <- read.csv("all_gene_pairs/cer_genePairs_df.csv",
         as.is=TRUE,
         stringsAsFactors = FALSE)

dros_gp_tbl <- read.csv("all_gene_pairs/dros_genePairs_df.csv",
                       as.is=TRUE,
                       stringsAsFactors = FALSE)

#
dros_gp_tbl <- filter(dros_gp_tbl, species!="Caenorhabditis_elegans")
#

vert_gp_tbl <- read.csv("all_gene_pairs/vert_genePairs_df.csv",
                        as.is=TRUE,
                        stringsAsFactors = FALSE)

vert_gp_oo <- read.table("oo_dfs/human_all_gp.txt",
                         header=TRUE,
                         stringsAsFactors = FALSE,
                         row.names=NULL,
                         fill=TRUE)

cer_gp_oo <- read.table("oo_dfs/yeast_all_gp.txt",
                         header=TRUE,
                         stringsAsFactors = FALSE,
                         row.names=NULL,
                         fill=TRUE)

dros_gp_oo <- read.table("oo_dfs/fruitfly_all_gp.txt",
                        header=TRUE,
                        stringsAsFactors = FALSE,
                        row.names=NULL,
                        fill=TRUE)


cer_nr <- read.table("nr_results/yeast_001_prot_vs_nr_noSC_clean.out",
                     stringsAsFactors = FALSE,
                     fill=TRUE,
                     row.names=NULL,
                     header=FALSE,
                     sep='\t')

cer_nr <- filter(cer_nr, V2 != "", V3 == "") %>%
  mutate(V2=gsub(" ", "_", V2)) %>%
  dplyr::select(-V3)

dros_nr <- read.table("nr_results/dros_001_prot_vs_nr_noDM_clean_noU.out",
                      stringsAsFactors = FALSE,
                      fill=TRUE,
                      row.names=NULL,
                      header=FALSE,
                      sep='\t')

dros_nr <- filter(dros_nr, V2 != "") %>%
  mutate(V2=gsub(" ", "_", V2))

vert_nr <- read.table("nr_results/human_0001_prot_vs_nr_noHS_clean_noU.out",
                      stringsAsFactors = FALSE,
                      fill=TRUE,
                      row.names=NULL,
                      header=FALSE,
                      sep='\t')

vert_nr <- filter(vert_nr, V2 != "", V3 == "") %>%
  mutate(V2=gsub(" ", "_", V2)) %>%
  dplyr::select(-V3)

#####


dros_gp_tbl <- dros_gp_tbl  %>%
  filter(!(lc_ortho>50 & lc_focal>50)) %>%
  ungroup()

dros_gp_tbl <- as.data.frame(dros_gp_tbl)

cer_gp_tbl <- cer_gp_tbl  %>%
  filter(!(lc_ortho>50 & lc_focal>50)) %>%
  filter(Gene_focal != "BSC4")

cer_gp_tbl <- as.data.frame(cer_gp_tbl)

vert_gp_tbl <- vert_gp_tbl  %>%
  filter(!(lc_ortho>50 & lc_focal>50)) %>%
  ungroup()

vert_gp_tbl <- as.data.frame(vert_gp_tbl)

####### GC and CDS correlation in human

CDS_vert_cor <- cor.test(vert_gp_tbl$len_focal, vert_gp_tbl$len_ortho, method="spearman")
GC_vert_cor <- cor.test(vert_gp_tbl$GC_focal, vert_gp_tbl$GC_ortho, method="spearman")

####### Protein properties

Exp_cer_cor <- cor.test(cer_gp_tbl$Epct.1_focal, cer_gp_tbl$Epct.1_ortho, method="spearman")
Coil_cer_cor <- cor.test(cer_gp_tbl$Cpct_focal, cer_gp_tbl$Cpct_ortho, method="spearman")
Coil_vert_cor <- cor.test(vert_gp_tbl$Cpct_focal, vert_gp_tbl$Cpct_ortho, method="spearman")
Hel_vert_cor <- cor.test(vert_gp_tbl$Hpct_focal, vert_gp_tbl$Hpct_ortho, method="spearman")




##########
## pfam results
##########

cer_focal_pfam <- read.table("Pfam_search_raw_data/yeast_focal_vs_pfam.txt",
                             header=TRUE,
                             stringsAsFactors = FALSE,
                             row.names=NULL,
                             fill=TRUE)
colnames(cer_focal_pfam) <- c(c("Gene_focal"), colnames(cer_focal_pfam)[2:15])
cer_ortho_pfam <- read.table("Pfam_search_raw_data/yeast_target_vs_pfam.txt",
                             header=TRUE,
                             stringsAsFactors = FALSE,
                             row.names=NULL,
                             fill=TRUE)
colnames(cer_ortho_pfam) <- c(c("Gene_ortho"), colnames(cer_ortho_pfam)[2:15])

cer_gp_pfam <- left_join(cer_gp_tbl, cer_focal_pfam, by="Gene_focal")

#number with pfam match in focal
focal <- filter(cer_gp_pfam, E.value <0.001)
cer_focal_prop <- length(unique(focal$Gene_focal))/length(unique(cer_gp_tbl$Gene_focal))
cer_focal_pfam_raw <- length(unique(focal$Gene_focal))
cer_focal_pfam_all <- length(unique(cer_gp_tbl$Gene_focal))

#number with pfam match in focal
ortho <-left_join(cer_gp_pfam, cer_ortho_pfam, by="Gene_ortho") %>%
  filter(E.value.y <0.001)
cer_ortho_prop <- length(unique(ortho$Gene_ortho))/length(unique(cer_gp_tbl$Gene_ortho))
cer_ortho_pfam_raw <- length(unique(ortho$Gene_ortho))
cer_ortho_pfam_all <- length(unique(cer_gp_tbl$Gene_ortho))

cer_gp_pfam <- left_join(cer_gp_pfam, cer_ortho_pfam, by="Gene_ortho") %>%
  filter(E.value.x <0.001 & E.value.y < 0.001)
cer_gp_pfam_slim <- unique(cer_gp_pfam[,c("Gene_focal", "Gene_ortho", "species", "hmm_name.x", "hmm_name.y")])

cer_gp_pfam_slim$same = FALSE
cer_gp_pfam_slim[which(cer_gp_pfam_slim$hmm_name.x==cer_gp_pfam_slim$hmm_name.y), "same"] = TRUE
cer_gp_pfam_slim <- dplyr::select(cer_gp_pfam_slim, -hmm_name.x, -hmm_name.y)
cer_gp_pfam_slim <- unique(cer_gp_pfam_slim)

cer_same_cluster <- unique(cer_gp_tbl[which(cer_gp_tbl$Same_Cluster=="Same_ortho_clusterTrue"), c("Gene_focal", "Gene_ortho", "species")])
cer_same_cluster$same <- TRUE

cer_gp_pfam_slim <- unique(bind_rows(cer_gp_pfam_slim, cer_same_cluster))

cer_gp_pfam_slim$focal <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "len_focal"])
cer_gp_pfam_slim$ortho <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "len_ortho"])
cer_gp_pfam_slim$var <- "CDS length"

cer_main_df <- cer_gp_pfam_slim

cer_gp_pfam_slim$focal <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "disoPCT_focal"])
cer_gp_pfam_slim$ortho <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "disoPCT_ortho"])
cer_gp_pfam_slim$var <- "ISD pct"

cer_main_df <- rbind(cer_main_df, cer_gp_pfam_slim)


cer_gp_pfam_slim$focal <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "lc_focal"])
cer_gp_pfam_slim$ortho <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "lc_ortho"])
cer_gp_pfam_slim$var <- "LowComp pct"

cer_main_df <- rbind(cer_main_df, cer_gp_pfam_slim)

cer_gp_pfam_slim$focal <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "GC_focal"])
cer_gp_pfam_slim$ortho <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "GC_ortho"])
cer_gp_pfam_slim$var <- "GC pct"

cer_main_df <- rbind(cer_main_df, cer_gp_pfam_slim)

cer_gp_pfam_slim$focal <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "tm_focal"])
cer_gp_pfam_slim$ortho <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "tm_ortho"])
cer_gp_pfam_slim$var <- "TM pct"

cer_main_df <- rbind(cer_main_df, cer_gp_pfam_slim)
cer_gp_pfam_slim$focal <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "Hpct_focal"])
cer_gp_pfam_slim$ortho <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "Hpct_ortho"])
cer_gp_pfam_slim$var <- "Helix pct"

cer_main_df <- rbind(cer_main_df, cer_gp_pfam_slim)
cer_gp_pfam_slim$focal <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "Epct_focal"])
cer_gp_pfam_slim$ortho <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "Epct_ortho"])
cer_gp_pfam_slim$var <- "Strand pct"

cer_main_df <- rbind(cer_main_df, cer_gp_pfam_slim)
cer_gp_pfam_slim$focal <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "Bpct_focal"])
cer_gp_pfam_slim$ortho <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "Bpct_ortho"])
cer_gp_pfam_slim$var <- "Buried pct"


cer_main_df <- rbind(cer_main_df, cer_gp_pfam_slim)
cer_gp_pfam_slim$focal <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "Epct.1_focal"])
cer_gp_pfam_slim$ortho <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "Epct.1_ortho"])
cer_gp_pfam_slim$var <- "Exposed pct"


cer_main_df <- rbind(cer_main_df, cer_gp_pfam_slim)
cer_gp_pfam_slim$focal <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "Cpct_focal"])
cer_gp_pfam_slim$ortho <- with(cer_gp_pfam_slim, cer_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(cer_gp_tbl$Gene_focal, cer_gp_tbl$Gene_ortho)), "Cpct_ortho"])
cer_gp_pfam_slim$var <- "Coil pct"

cer_main_df <- rbind(cer_main_df, cer_gp_pfam_slim)


#####

dros_focal_pfam <- read.table("Pfam_search_raw_data/fruitfly_focal_vs_pfam.txt",
                              header=TRUE,
                              stringsAsFactors = FALSE,
                              row.names=NULL,
                              fill=TRUE)
colnames(dros_focal_pfam) <- c(c("Gene_focal"), colnames(dros_focal_pfam)[2:15])
dros_ortho_pfam <- read.table("Pfam_search_raw_data/fruitfly_target_vs_pfam.txt",
                              header=TRUE,
                              stringsAsFactors = FALSE,
                              row.names=NULL,
                              fill=TRUE)
colnames(dros_ortho_pfam) <- c(c("Gene_ortho"), colnames(dros_ortho_pfam)[2:15])

dros_gp_pfam <- left_join(dros_gp_tbl, dros_focal_pfam, by="Gene_focal")

#number with pfam match in focal
focal <- filter(dros_gp_pfam, E.value <0.001)
dros_focal_prop <- length(unique(focal$Gene_focal))/length(unique(dros_gp_tbl$Gene_focal))
dros_focal_pfam_raw <- length(unique(focal$Gene_focal))
dros_focal_pfam_all <- length(unique(dros_gp_tbl$Gene_focal))

#number with pfam match in focal
ortho <-left_join(dros_gp_pfam, dros_ortho_pfam, by="Gene_ortho") %>%
  filter(E.value.y <0.001)
dros_ortho_prop <- length(unique(ortho$Gene_ortho))/length(unique(dros_gp_tbl$Gene_ortho))
dros_ortho_pfam_raw <- length(unique(ortho$Gene_ortho))
dros_ortho_pfam_all <- length(unique(dros_gp_tbl$Gene_ortho))

dros_gp_pfam <- left_join(dros_gp_pfam, dros_ortho_pfam, by="Gene_ortho") %>%
  filter(E.value.x <0.001 & E.value.y < 0.001)
dros_gp_pfam_slim <- unique(dros_gp_pfam[,c("Gene_focal", "Gene_ortho", "species", "hmm_name.x", "hmm_name.y")])

dros_gp_pfam_slim$same = FALSE
dros_gp_pfam_slim[which(dros_gp_pfam_slim$hmm_name.x==dros_gp_pfam_slim$hmm_name.y), "same"] = TRUE
dros_gp_pfam_slim <- dplyr::select(dros_gp_pfam_slim, -hmm_name.x, -hmm_name.y)
dros_gp_pfam_slim <- as.data.frame(unique(dros_gp_pfam_slim))

dros_gp_pfam_slim$focal <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "len_focal"])
dros_gp_pfam_slim$ortho <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "len_ortho"])
dros_gp_pfam_slim$var <- "CDS length"

dros_main_df <- dros_gp_pfam_slim

dros_gp_pfam_slim$focal <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "disoPCT_focal"])
dros_gp_pfam_slim$ortho <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "disoPCT_ortho"])
dros_gp_pfam_slim$var <- "ISD pct"

dros_main_df <- rbind(dros_main_df, dros_gp_pfam_slim)

dros_gp_pfam_slim$focal <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "lc_focal"])
dros_gp_pfam_slim$ortho <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "lc_ortho"])
dros_gp_pfam_slim$var <- "LowComp pct"

dros_main_df <- rbind(dros_main_df, dros_gp_pfam_slim)

dros_gp_pfam_slim$focal <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "GC_focal"])
dros_gp_pfam_slim$ortho <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "GC_ortho"])
dros_gp_pfam_slim$var <- "GC pct"

dros_main_df <- rbind(dros_main_df, dros_gp_pfam_slim)

dros_gp_pfam_slim$focal <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "tm_focal"])
dros_gp_pfam_slim$ortho <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "tm_ortho"])
dros_gp_pfam_slim$var <- "TM pct"

dros_main_df <- rbind(dros_main_df, dros_gp_pfam_slim)
dros_gp_pfam_slim$focal <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "Hpct_focal"])
dros_gp_pfam_slim$ortho <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "Hpct_ortho"])
dros_gp_pfam_slim$var <- "Helix pct"

dros_main_df <- rbind(dros_main_df, dros_gp_pfam_slim)
dros_gp_pfam_slim$focal <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "Epct_focal"])
dros_gp_pfam_slim$ortho <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "Epct_ortho"])
dros_gp_pfam_slim$var <- "Strand pct"

dros_main_df <- rbind(dros_main_df, dros_gp_pfam_slim)
dros_gp_pfam_slim$focal <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "Bpct_focal"])
dros_gp_pfam_slim$ortho <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "Bpct_ortho"])
dros_gp_pfam_slim$var <- "Buried pct"


dros_main_df <- rbind(dros_main_df, dros_gp_pfam_slim)
dros_gp_pfam_slim$focal <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "Epct.1_focal"])
dros_gp_pfam_slim$ortho <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "Epct.1_ortho"])
dros_gp_pfam_slim$var <- "Exposed pct"


dros_main_df <- rbind(dros_main_df, dros_gp_pfam_slim)
dros_gp_pfam_slim$focal <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "Cpct_focal"])
dros_gp_pfam_slim$ortho <- with(dros_gp_pfam_slim, dros_gp_tbl[match(paste(Gene_focal, Gene_ortho), paste(dros_gp_tbl$Gene_focal, dros_gp_tbl$Gene_ortho)), "Cpct_ortho"])
dros_gp_pfam_slim$var <- "Coil pct"

dros_main_df <- rbind(dros_main_df, dros_gp_pfam_slim)

######

vert_focal_pfam <- read.table("Pfam_search_raw_data/human_focal_vs_pfam.txt",
                              header=TRUE,
                              stringsAsFactors = FALSE,
                              row.names=NULL,
                              fill=TRUE)
colnames(vert_focal_pfam) <- c(c("Gene_focal"), colnames(vert_focal_pfam)[2:15])
vert_ortho_pfam <- read.table("Pfam_search_raw_data/human_target_vs_pfam.txt",
                              header=TRUE,
                              stringsAsFactors = FALSE,
                              row.names=NULL,
                              fill=TRUE)
colnames(vert_ortho_pfam) <- c(c("Gene_ortho_spec"), colnames(vert_ortho_pfam)[2:15])

vert_gp_pfam <- left_join(vert_gp_tbl, vert_focal_pfam, by="Gene_focal")

#number with pfam match in focal
focal <- filter(vert_gp_pfam, E.value <0.001)
vert_focal_prop <- length(unique(focal$Gene_focal))/length(unique(vert_gp_tbl$Gene_focal))
vert_focal_pfam_raw <- length(unique(focal$Gene_focal))
vert_focal_pfam_all <- length(unique(vert_gp_tbl$Gene_focal))
#number with pfam match in focal
vert_gp_pfam$Gene_ortho_spec <- paste(vert_gp_pfam$species, vert_gp_pfam$Gene_ortho, vert_gp_pfam$len_ortho, sep="_")
ortho <-left_join(vert_gp_pfam, vert_ortho_pfam, by="Gene_ortho_spec") %>%
  filter(E.value.y <0.001)
vert_ortho_prop <- length(unique(ortho$Gene_ortho_spec))/length(unique(paste(vert_gp_tbl$Gene_ortho, vert_gp_tbl$species, vert_gp_tbl$len_ortho)))
vert_ortho_pfam_raw <- length(unique(ortho$Gene_ortho_spec))
vert_ortho_pfam_all <- length(unique(paste(vert_gp_tbl$Gene_ortho, vert_gp_tbl$species, vert_gp_tbl$len_ortho)))

vert_gp_tbl$Gene_ortho_spec <- paste(vert_gp_tbl$species, vert_gp_tbl$Gene_ortho, vert_gp_tbl$len_ortho, sep="_")


vert_gp_pfam <- left_join(vert_gp_pfam, vert_ortho_pfam, by="Gene_ortho_spec") %>%
  filter(E.value.x <0.001 & E.value.y < 0.001)
vert_gp_pfam_slim <- unique(vert_gp_pfam[,c("Gene_focal", "Gene_ortho", "species", "Gene_ortho_spec", "hmm_name.x", "hmm_name.y")])
#vert_gp_pfam_slim_spec <- unique(vert_gp_pfam[,c("Gene_focal", "Gene_ortho", "hmm_name.x", "hmm_name.y", "species")])

vert_gp_pfam_slim$same = FALSE
#vert_gp_pfam_slim_spec$same = FALSE
vert_gp_pfam_slim[which(vert_gp_pfam_slim$hmm_name.x==vert_gp_pfam_slim$hmm_name.y), "same"] = TRUE
#vert_gp_pfam_slim_spec[which(vert_gp_pfam_slim_spec$hmm_name.x==vert_gp_pfam_slim_spec$hmm_name.y), "same"] = TRUE
#vert_gp_pfam_all_dom <- unique(vert_gp_pfam_slim_spec)
vert_gp_pfam_slim <- dplyr::select(vert_gp_pfam_slim, -hmm_name.x, -hmm_name.y)
vert_gp_pfam_slim <- unique(vert_gp_pfam_slim)

vert_same_cluster <- unique(vert_gp_tbl[which(vert_gp_tbl$Same_Cluster=="Same_ortho_clusterTrue"), c("Gene_focal", "Gene_ortho", "species", "Gene_ortho_spec")])
vert_same_cluster$same <- TRUE

vert_gp_pfam_slim <- unique(bind_rows(vert_gp_pfam_slim, vert_same_cluster))

vert_gp_pfam_slim$focal <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "len_focal"])
vert_gp_pfam_slim$ortho <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "len_ortho"])
vert_gp_pfam_slim$var <- "CDS length"

vert_main_df <- vert_gp_pfam_slim

vert_gp_pfam_slim$focal <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "disoPCT_focal"])
vert_gp_pfam_slim$ortho <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "disoPCT_ortho"])
vert_gp_pfam_slim$var <- "ISD pct"

vert_main_df <- rbind(vert_main_df, vert_gp_pfam_slim)

vert_gp_pfam_slim$focal <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "lc_focal"])
vert_gp_pfam_slim$ortho <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "lc_ortho"])
vert_gp_pfam_slim$var <- "LowComp pct"

vert_main_df <- rbind(vert_main_df, vert_gp_pfam_slim)

vert_gp_pfam_slim$focal <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "GC_focal"])
vert_gp_pfam_slim$ortho <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "GC_ortho"])
vert_gp_pfam_slim$var <- "GC pct"

vert_main_df <- rbind(vert_main_df, vert_gp_pfam_slim)

vert_gp_pfam_slim$focal <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "tm_focal"])
vert_gp_pfam_slim$ortho <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "tm_ortho"])
vert_gp_pfam_slim$var <- "TM pct"

vert_main_df <- rbind(vert_main_df, vert_gp_pfam_slim)
vert_gp_pfam_slim$focal <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "Hpct_focal"])
vert_gp_pfam_slim$ortho <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "Hpct_ortho"])
vert_gp_pfam_slim$var <- "Helix pct"

vert_main_df <- rbind(vert_main_df, vert_gp_pfam_slim)
vert_gp_pfam_slim$focal <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "Epct_focal"])
vert_gp_pfam_slim$ortho <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "Epct_ortho"])
vert_gp_pfam_slim$var <- "Strand pct"

vert_main_df <- rbind(vert_main_df, vert_gp_pfam_slim)
vert_gp_pfam_slim$focal <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "Bpct_focal"])
vert_gp_pfam_slim$ortho <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "Bpct_ortho"])
vert_gp_pfam_slim$var <- "Buried pct"


vert_main_df <- rbind(vert_main_df, vert_gp_pfam_slim)
vert_gp_pfam_slim$focal <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "Epct.1_focal"])
vert_gp_pfam_slim$ortho <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "Epct.1_ortho"])
vert_gp_pfam_slim$var <- "Exposed pct"


vert_main_df <- rbind(vert_main_df, vert_gp_pfam_slim)
vert_gp_pfam_slim$focal <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "Cpct_focal"])
vert_gp_pfam_slim$ortho <- with(vert_gp_pfam_slim, vert_gp_tbl[match(paste(Gene_focal, Gene_ortho_spec), paste(vert_gp_tbl$Gene_focal, vert_gp_tbl$Gene_ortho_spec)), "Cpct_ortho"])
vert_gp_pfam_slim$var <- "Coil pct"

vert_main_df <- rbind(vert_main_df, vert_gp_pfam_slim)
vert_main_df <- select(vert_main_df, -Gene_ortho_spec)
vert_main_df$tag <- "human"
dros_main_df$tag <- "fruitfly"
cer_main_df$tag <- "yeast"

all_main_df <- rbind(cer_main_df, dros_main_df, vert_main_df)
all_main_df$focal <- round(all_main_df$focal, 3)
all_main_df$ortho <- round(all_main_df$ortho, 3)

write.csv(all_main_df, "Figure_7_source_data_1.csv", row.names = FALSE)


colnames(all_main_df) <- c("Gene_focal", "Gene_ortho", "species", "same", "focal", "ortho", "var", "tag")
all_main_df$alpha <- 0.3
all_main_df[which(all_main_df$same), "alpha"] <- 1

### get the pairs

same_pairs <- filter(all_main_df, same==TRUE, var =="TM pct")
diff_pairs <- filter(all_main_df, same==FALSE, var =="TM pct")

k=filter(all_main_df, Gene_focal=="MNE1")[c(2:20),]
t.test(k$focal, k$ortho, paired=TRUE)

#######

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

all_prop_df <- as.data.frame(rbind(c(cer_focal_prop, cer_ortho_prop),
                                   c(dros_focal_prop, dros_ortho_prop),
                                   c(vert_focal_prop, vert_ortho_prop)))

all_prop_df$dataset <- c("yeast", "fruitfly", "human")
colnames(all_prop_df) <- c("focal", "ortho", "dataset")
all_prop_df <- melt(all_prop_df)
all_prop_df$dataset <- factor(all_prop_df$dataset, levels=c("yeast", "fruitfly", "human"), ordered = TRUE)
all_prop_df$SE <- c(sqrt(cer_focal_prop*(1-cer_focal_prop)/cer_focal_pfam_all),
                    sqrt(dros_focal_prop*(1-dros_focal_prop)/dros_focal_pfam_all),
                    sqrt(vert_focal_prop*(1-vert_focal_prop)/vert_focal_pfam_all),
                    sqrt(cer_ortho_prop*(1-cer_ortho_prop)/cer_ortho_pfam_all),
                    sqrt(dros_ortho_prop*(1-dros_ortho_prop)/dros_ortho_pfam_all),
                    sqrt(vert_ortho_prop*(1-vert_ortho_prop)/vert_ortho_pfam_all))

all_raw_df <- as.data.frame(rbind(c(cer_focal_pfam_raw, cer_focal_pfam_all, cer_ortho_pfam_raw, cer_ortho_pfam_all),
                                   c(dros_focal_pfam_raw, dros_focal_pfam_all, dros_ortho_pfam_raw, dros_ortho_pfam_all),
                                   c(vert_focal_pfam_raw, vert_focal_pfam_all,vert_ortho_pfam_raw, vert_ortho_pfam_all)))
all_raw_df$dataset <- c("yeast", "fruitfly", "human")
colnames(all_raw_df) <- c("focal Pfam", "focal all", "ortho Pfam", "ortho all", "dataset")
write.csv(all_raw_df, "Figure_7-supplement_1.csv", row.names = FALSE)


chisq.test(rbind(c(cer_focal_pfam_raw, cer_focal_pfam_all-cer_focal_pfam_raw),
                 c(cer_ortho_pfam_raw, cer_ortho_pfam_all-cer_ortho_pfam_raw)))
chisq.test(rbind(c(vert_focal_pfam_raw, vert_focal_pfam_all-vert_focal_pfam_raw),
                 c(vert_ortho_pfam_raw, vert_ortho_pfam_all-vert_ortho_pfam_raw)))
chisq.test(rbind(c(dros_focal_pfam_raw, dros_focal_pfam_all-dros_focal_pfam_raw),
                 c(dros_ortho_pfam_raw, dros_ortho_pfam_all-dros_ortho_pfam_raw)))


pfam_prop_pl <- ggplot(data = all_prop_df,
                       aes(x=dataset,
                           y=value,
                           fill=variable,
                           ymin=value-SE,
                           ymax=value+SE)) +
  geom_bar(stat="identity",
           position=position_dodge(width=0.9)) +
  geom_errorbar(position=position_dodge(width=0.9),
                width=0.4) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines")) +
  xlab("") +
  ylab("Proportion with a Pfam domain") +
  scale_fill_manual(values = c("lightgray", "darkgray"))

ggsave("figures/Figure7A.pdf")


cor_plots <- ggplot() + 
  geom_point(data = all_main_df, 
             aes(x=focal, 
                 y=ortho,
                 colour=same,
                 alpha=alpha)) +
  geom_smooth(data = all_main_df, 
             aes(x=focal, 
                 y=ortho,
                 colour=same),
             method="lm") +
  facet_wrap(~var, scales="free") +
  guides(alpha="none",
         colour=guide_legend(title="Pfam domain"))

ggsave("figures/Figure7C.pdf")

sup_table_3 <- filter(all_main_df, same==TRUE) %>%
                group_by(var) %>%
                dplyr::summarize(r = round(cor.test(focal, ortho, method="spearman")$estimate, 3),
                                 p = round(cor.test(focal, ortho, method="spearman")$p.value, 5),
                                 cor_sign = 10*p<0.05)

write.csv(sup_table_3, "Figure_7-supplement_2.csv", row.names=FALSE)


#### Figure 8B ####
###################

cer_keep <- c()
for (gene in unique(cer_gp_tbl$Gene_focal))
{
  farTrip <- unique(cer_gp_oo[which(cer_gp_oo$Gene_focal==gene), "species"])
  farMis <-  unique(cer_gp_tbl[which(cer_gp_tbl$Gene_focal==gene), "species"])
  rem <- setdiff(farTrip, farMis)
  rem_div <- cer_div[rem, "V2"]
  rem_mis <- cer_div[farMis, "V2"]
  spec_with_match <- setdiff(cer_nr[which(cer_nr$V1==gene),"V2"], farMis)
  spec_with_match <- spec_with_match[which(!grepl("Saccharomyces*", spec_with_match))]
  
  if ((length(rem)==0 || min(rem_mis) >= max(rem_div)) & (length(spec_with_match) == 0) & length(unique(rem_mis)) >1)
  {
    cer_keep <- c(cer_keep, gene)
  }
}

gene_freq_cer <- filter(cer_gp_tbl, Gene_focal %in% cer_keep) %>%
  group_by(Gene_focal) %>%
  dplyr::summarize(no_sp = length(unique(div)),
                   ortho = mean(len_ortho))
gene_freq_cer$focal <- cer_gp_tbl[match(gene_freq_cer$Gene_focal, cer_gp_tbl$Gene_focal), "len_focal"]
gene_freq_cer_df <- melt(gene_freq_cer, id.vars=c("Gene_focal", "no_sp"))
gene_freq_cer_df <- gene_freq_cer_df[which(gene_freq_cer_df$no_sp>1),]
gene_freq_cer_df$dataset <- "yeast"

vert_keep <- c()
all_species_check <- c()

specOutput <- data.frame(found=character(), notFound=character())
vert_to_remove <- c("C6orf222", "CCDC195", "CCDC7", "FMR1NB", "GVQW2", "MMP24OS", "SLC51B", "SMLR1", "TAC4", "TEX29", "GAGE12B", "SUPT16H", "PROSER3")

for (gene in unique(vert_gp_tbl$Gene_focal))
{
  farTrip <- unique(vert_gp_oo[which(vert_gp_oo$Gene_focal==gene), "species"])
  farMis <-  unique(vert_gp_tbl[which(vert_gp_tbl$Gene_focal==gene), "species"])
  rem <- setdiff(farTrip, farMis)
  rem_div <- vert_div[rem, "V2"]
  rem_mis <- vert_div[farMis, "V2"]
  spec_with_match <- setdiff(vert_nr[which(vert_nr$V1==gene),"V2"], farMis)
  
  if ((length(rem)==0 || min(rem_mis) >= max(rem_div)) & length(unique(rem_mis)) >1 & (!gene %in% vert_to_remove))
  {
    all_species_check <- c(all_species_check, spec_with_match, farMis)
    vert_keep <- c(vert_keep, gene)
    specOutput <- rbind(specOutput, as.data.frame(cbind(gene, toString(spec_with_match), toString(farMis))))
  }
}


gene_freq_vert <- filter(vert_gp_tbl, Gene_focal %in% vert_keep) %>%
  group_by(Gene_focal) %>%
  dplyr::summarize(no_sp = length(unique(div)),
                   ortho = mean(len_ortho))
gene_freq_vert$focal <- vert_gp_tbl[match(gene_freq_vert$Gene_focal, vert_gp_tbl$Gene_focal), "len_focal"]
gene_freq_vert_df <- melt(gene_freq_vert, id.vars=c("Gene_focal", "no_sp"))
gene_freq_vert_df <- gene_freq_vert_df[which(gene_freq_vert_df$no_sp>1),]
gene_freq_vert_df$dataset <- "human"

dros_keep <- c()
for (gene in unique(dros_gp_tbl$Gene_focal))
{
  farTrip <- unique(dros_gp_oo[which(dros_gp_oo$Gene_focal==gene), "species"])
  farMis <-  unique(dros_gp_tbl[which(dros_gp_tbl$Gene_focal==gene), "species"])
  rem <- setdiff(farTrip, farMis)
  rem_div <- dros_div[rem, "V5"]
  rem_mis <- dros_div[farMis, "V5"]
  spec_with_match <- setdiff(dros_nr[which(dros_nr$V1==gene),"V2"], farMis)
  spec_with_match <- spec_with_match[which(!grepl("Drosophila*", spec_with_match))]
  
  if ((length(rem)==0 || min(rem_mis) >= max(rem_div)) & (length(spec_with_match) == 0) & length(unique(rem_mis)) >1)
  {
    dros_keep <- c(dros_keep, gene)
  }
}


gene_freq_dros <- filter(dros_gp_tbl, Gene_focal %in% dros_keep) %>%
  group_by(Gene_focal) %>%
  dplyr::summarize(no_sp = length(unique(div)),
                   ortho = mean(len_ortho))
gene_freq_dros$focal <- dros_gp_tbl[match(gene_freq_dros$Gene_focal, dros_gp_tbl$Gene_focal), "len_focal"]
gene_freq_dros_df <- melt(gene_freq_dros, id.vars=c("Gene_focal", "no_sp"))
gene_freq_dros_df <- gene_freq_dros_df[which(gene_freq_dros_df$no_sp>1),]
gene_freq_dros_df$dataset <- "fruitfly"



pdros <- t.test(filter(gene_freq_dros_df, variable=="ortho")$value,
                filter(gene_freq_dros_df, variable=="focal")$value, paired = TRUE)

pcer <- t.test(filter(gene_freq_cer_df, variable=="ortho")$value,
               filter(gene_freq_cer_df, variable=="focal")$value, paired = TRUE)

pvert <- t.test(filter(gene_freq_vert_df, variable=="ortho")$value,
                filter(gene_freq_vert_df, variable=="focal")$value, paired = TRUE)

df_all <- rbind(gene_freq_cer_df, gene_freq_dros_df, gene_freq_vert_df)
df_all$dataset <- factor(df_all$dataset, levels=c("yeast", "fruitfly", "human"), ordered = TRUE)

df_gathered <- spread(df_all, variable, value)
df_gathered_cp <- df_gathered
colnames(df_gathered_cp) <- c("Focal gene", "No. species with undetectable homologues", "dataset", "Mean undetectable homologue CDS length", "Focal CDS length")
write.csv(df_gathered_cp, "Figure_8-supplement_1.csv", row.names = FALSE)

dim(gene_freq_vert_df)[[1]]/2
dim(gene_freq_dros_df)[[1]]/2
dim(gene_freq_cer_df)[[1]]/2

pl_len <- ggplot(data=df_all,
                 aes(x=variable,
                     y=value)) + 
  geom_boxplot(alpha=0.2) +
  geom_point() +
  geom_segment(data=df_gathered,
               aes(x="ortho",
                   xend="focal",
                   y=ortho,
                   yend=focal),
               alpha=0.3,
               linetype=2) +
  stat_summary(fun.y = "mean",
               geom="point",
               colour="red",
               shape=8,
               size=3) +
  xlab("") + ylab("CDS length (nt)") +
  facet_grid(~dataset) +
  coord_cartesian(ylim=c(0,4000)) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14)) +
  scale_x_discrete(labels=c("target", "focal"))


ggsave("figures/Figure8B.pdf")


##### dn/ds plots

vert_all_genes_list <- read.table("full_lists_genes/all_human_genes_checked_compl.txt",
                                  stringsAsFactors = FALSE,
                                  header=TRUE)

cer_trip_genes <- unique(cer_gp_oo[which(cer_gp_oo$species %in% c("Saccharomyces_arboricola", "Saccharomyces_kudriavzevii")), "Gene_focal"])
dros_trip_genes <- unique(dros_gp_oo[which(dros_gp_oo$species=="Drosophila_simulans"), "Gene_focal"])
ver_trip_genes <- unique(vert_gp_oo[which(vert_gp_oo$species=="Mus_musculus"), "Gene_focal"])

cer_prot_gene_names <- read.table("dnds/gene_prot_names.txt")

cer_dnds <- read.table("dnds/scer_dn_ds.txt",
                       header=FALSE,
                       stringsAsFactors = FALSE)
dros_dnds <- read.table("dnds/Dmel_Dsim_analysis_results_flydivas_v1.2.txt",
                        header=TRUE,
                        stringsAsFactors = FALSE)
vert_dnds <- read.table("dnds/human_mouse_dnds.tsv",
                        header=TRUE,
                        stringsAsFactors = FALSE,
                        sep='\t')
vert_dnds <- filter(vert_dnds,
                    !is.na(dN.with.Mouse),
                    Gene.name %in% vert_all_genes_list$Gene_name) %>%
              group_by(Gene.name) %>%
              summarise(Gene.stable.ID=head(Gene.name, 1),
                        dN.with.Mouse=head(dN.with.Mouse, 1),
                        dS.with.Mouse=head(dS.with.Mouse, 1))

cer_dnds$V1 <- as.character(cer_dnds$V1)
cer_dnds$V1 <- cer_prot_gene_names[match(cer_dnds$V1, cer_prot_gene_names$V1), "V2"]

cer_dnds$inTrip <- "NO"
dros_dnds$inTrip <- "NO"
vert_dnds$inTrip <- "NO"

cer_dnds[which(cer_dnds$V1 %in% cer_trip_genes), "inTrip"] = "YES"
dros_dnds[which(dros_dnds$symbol %in% dros_trip_genes), "inTrip"] = "YES"
vert_dnds[which(vert_dnds$Gene.name %in% ver_trip_genes), "inTrip"] = "YES"

vert_dnds <- dplyr::select(vert_dnds, Gene.name, dN.with.Mouse, dS.with.Mouse, inTrip)
dros_dnds <- dplyr::select(dros_dnds, symbol, dN, dS, inTrip)
cer_dnds <- dplyr::select(cer_dnds, V1, V3, V4, inTrip)

cer_dnds$tag <- "yeast"
vert_dnds$tag <- "human"
dros_dnds$tag <- "fruitfly"

colnames(cer_dnds) <- c("Gene_name","dN", "dS", "micro-synteny", "tag")
colnames(vert_dnds) <- c("Gene_name", "dN", "dS", "micro-synteny", "tag")
colnames(dros_dnds) <- c("Gene_name", "dN", "dS", "micro-synteny", "tag")

vert_dnds$dN <- as.numeric(vert_dnds$dN)
vert_dnds$dS <- as.numeric(vert_dnds$dS)
dros_dnds$dN <- as.numeric(as.character(dros_dnds$dN))
dros_dnds$dS <- as.numeric(as.character(dros_dnds$dS))

df_all <- bind_rows(cer_dnds, vert_dnds, dros_dnds)
df_all$tag <- factor(df_all$tag, levels=c("yeast", "fruitfly", "human"), ordered = TRUE)

df_all$dNdS <- round(df_all$dN/df_all$dS,4)

df_all <- filter(df_all, dNdS != Inf)

write.csv(df_all, file="Figure_5_source_data_1.csv", row.names = FALSE)
colnames(df_all) <- c("Gene_name","dN","dS","inTrip","tag","dNdS")

dnds_stats <- df_all %>%
  group_by(tag) %>%
  summarise(dn_pval = wilcox.test(dN~inTrip)$p.value,
            ds_pval = wilcox.test(dS~inTrip)$p.value,
            dnds_pval = wilcox.test(dNdS~inTrip)$p.value,
            dn_es = abs(qnorm(dn_pval/2))/sqrt(n()),
            ds_es = abs(qnorm(ds_pval/2))/sqrt(n()),
            dnds_es = abs(qnorm(dnds_pval/2))/sqrt(n()),
            dn_Cliff = cliff.delta(dN, inTrip)$estimate,
            ds_Cliff = cliff.delta(dS, inTrip)$estimate,
            dnds_Cliff = cliff.delta(dNdS, inTrip)$estimate,
            N = n())


write.csv(dnds_stats, file="Figure5B.csv", row.names = FALSE)


df_all_d <- df_all %>%
  select(-dNdS) %>%
  filter(dN<0.5 & ((dS<0.5 & tag=="fruitfly") | (dS<2 & tag %in% c("human", "yeast")))) %>%
  gather(variable, value, dN, dS)

dens_dn_ds <- ggplot() +
  geom_density(data = df_all_d, 
               aes(x=value, 
                   fill=tag,
                   linetype=inTrip),
               alpha=0.2) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"),
        legend.position = "bottom") +
  facet_wrap(variable~tag, scales = "free")

ggsave(plot = dens_dn_ds, "figures/Figure5A.pdf")
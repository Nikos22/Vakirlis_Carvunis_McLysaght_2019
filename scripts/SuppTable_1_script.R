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

cer_df_tbl <- read.csv("scripts_submission/Fig3_cer_df.csv")
dros_df_tbl <- read.csv("scripts_submission/Fig3_dros_df.csv")
vert_df_tbl <- read.csv("scripts_submission/Fig3_vert_df.csv")

vert_div <- read.table("scripts_submission/divergence_times/vert_div.txt", row.names=1)
dros_div <- read.table("scripts_submission/divergence_times/dros_div.txt", row.names=1)
cer_div <- read.table("scripts_submission/divergence_times/cer_div.txt", row.names=1)

vert_orphan <- read.table("scripts_submission/all_genes_blast/vertebrates_forOrphan_all_-_final.txt",
                          header=TRUE,
                          stringsAsFactors = FALSE,
                          fill=TRUE)
cer_orphan <- read.table("scripts_submission/all_genes_blast/cer_all_genes_det_final.txt",
                         header=TRUE,
                         stringsAsFactors = FALSE)
dros_orphan <- read.table("scripts_submission/all_genes_blast/droso_all_genes_det_final.txt",
                          header=TRUE,
                          stringsAsFactors = FALSE)
dros_orphan <- unique(select(dros_orphan, species, Gene_focal))
cer_orphan <- unique(select(cer_orphan, species, Gene_focal))
vert_orphan <- unique(select(vert_orphan, species, Gene_focal))

### vertebrates


vert_df_tbl <- mutate(vert_df_tbl,
                      residues = ver_res[species, "V2"]) %>%
  arrange(div) # calculates new columns
vert_df_tbl$species <- factor(vert_df_tbl$species, ordered = TRUE)
vert_df_tbl$tag <- "human"

### drosophila

dros_df_tbl <- mutate(dros_df_tbl,
                      residues = dros_res[species, "V2"])
dros_df_tbl$species <- factor(dros_df_tbl$species, ordered = TRUE)
dros_df_tbl$tag <- "fruitfly"

### yeasts

cer_df_tbl <- mutate(cer_df_tbl,
                     residues = cer_res[species, "V2"])
cer_df_tbl$species <- factor(cer_df_tbl$species, ordered = TRUE)
cer_df_tbl$tag <- "yeast"


###

df <- bind_rows(cer_df_tbl, dros_df_tbl, vert_df_tbl)
df$tag <- factor(df$tag, levels=c("yeast", "fruitfly", "human"), ordered = TRUE)
df$species <- paste(substr(df$species, 1, 1), substr(unlist(strsplit(df$species, split = "_"))[2*(1:length(df$species))],1,3), sep='')
df_tbl_min_ev <- group_by(df, by=species) %>% slice(which.min(evalue))


#######


cer_fp_df <- read.table("scripts_submission/FP/yeast_FP_stats.txt", header=TRUE, stringsAsFactors = FALSE)
cer_fp_df_tbl <- tbl_df(cer_fp_df)
cer_fp_df_tbl <- mutate(cer_fp_df_tbl,
                        div = cer_div[species,"V2"],
                        residues = cer_res[species, "V2"])
cer_fp_df_tbl$tag <- "yeast"

dros_fp_df <- read.table("scripts_submission/FP/droso_FP_stats.txt", header=TRUE, stringsAsFactors = FALSE)
dros_fp_df_tbl <- tbl_df(dros_fp_df)
dros_fp_df_tbl <- mutate(dros_fp_df_tbl,
                         div = dros_div[species,"V5"],
                         residues = dros_res[species, "V2"])
dros_fp_df_tbl$tag <- "fruitfly"

ver_fp_df <- read.table("scripts_submission/FP/vert_FP_stats.txt", header=TRUE, stringsAsFactors = FALSE)
ver_fp_df_tbl <- tbl_df(ver_fp_df)
ver_fp_df_tbl <- mutate(ver_fp_df_tbl,
                        div = vert_div[species,"V2"],
                        residues = ver_res[species, "V2"])
ver_fp_df_tbl$tag <- "human"
#### these are the initial genePair files filtered through the orphan files, essentially removing
## those that have a genome-wide tblastn match

cer_tblastn_match <- read.table("scripts_submission/tblastn_genomewide/after_tblastn_cer.txt",
                                header=FALSE,
                                stringsAsFactors = FALSE,
                                row.names=NULL,
                                fill=TRUE)

vert_tblastn_match <- read.table("scripts_submission/tblastn_genomewide/after_tblastn_vert.txt",
                                 header=FALSE,
                                 stringsAsFactors = FALSE,
                                 row.names=NULL,
                                 fill=TRUE)

dros_tblastn_match <- read.table("scripts_submission/tblastn_genomewide/after_tblastn_droso.txt",
                                 header=FALSE,
                                 stringsAsFactors = FALSE,
                                 row.names=NULL,
                                 fill=TRUE)

####### F1 scores and ROC curves

cer_f1 <- mutate(cer_fp_df_tbl,
                 FPR = (found_elsewhere+found_tblastn)/total,
                 TPR = (1-((cer_df_tbl$no_match-cer_df_tbl$found_tblastn)/cer_df_tbl$total)),
                 TP = cer_df_tbl$total - (cer_df_tbl$no_match-cer_df_tbl$found_tblastn),
                 FP = found_elsewhere+found_tblastn,
                 TN = no_match-found_tblastn,
                 FN = cer_df_tbl$no_match-cer_df_tbl$found_tblastn)

dros_df_tbl <- filter(dros_df_tbl, evalue %in% dros_fp_df_tbl$evalue)
dros_f1 <- mutate(dros_fp_df_tbl,
                  FPR = (found_elsewhere+found_tblastn)/total,
                  TPR = (1-((dros_df_tbl$no_match-dros_df_tbl$found_tblastn)/dros_df_tbl$total)),
                  TP = dros_df_tbl$total - (dros_df_tbl$no_match-dros_df_tbl$found_tblastn),
                  FP = found_elsewhere+found_tblastn,
                  TN = no_match-found_tblastn,
                  FN = dros_df_tbl$no_match-dros_df_tbl$found_tblastn)

ver_f1 <- mutate(ver_fp_df_tbl,
                 FPR = (found_elsewhere+found_tblastn)/total,
                 TPR = (1-((vert_df_tbl$no_match-vert_df_tbl$found_tblastn)/vert_df_tbl$total)),
                 TP = vert_df_tbl$total - (vert_df_tbl$no_match-vert_df_tbl$found_tblastn),
                 FP = found_elsewhere+found_tblastn,
                 TN = no_match-found_tblastn,
                 FN = vert_df_tbl$no_match-vert_df_tbl$found_tblastn)

f1_score <- function(tp, fp, fn, tn)
{
  precision = tp/(tp+fp)
  recall = tp/(tp+fn)
  return((2*(precision*recall))/(precision+recall))  
}

cer_f1_final = cer_f1 %>%
  group_by(evalue, species) %>% 
  dplyr::summarize(f1score = f1_score(TP, FP, FN, TN),
                   mcc = round(mcc(TP=TP, FP=FP, FN=FN, TN=TN),2),
                   FPR=FPR) %>%
  group_by(species) %>%
  slice(which(FPR<0.05)) %>%
  group_by(species) %>%
  slice(which.max(evalue)) %>%
  mutate(div = cer_div[species, "V2"],
         residues = cer_res[species, "V2"])

dros_f1_final = dros_f1 %>%
  group_by(evalue, species) %>% 
  dplyr::summarize(f1score = f1_score(TP, FP, FN, TN),
                   mcc = round(mcc(TP=TP, FP=FP, FN=FN, TN=TN),2),
                   FPR=FPR) %>%
  group_by(species) %>%
  slice(which(FPR<0.05)) %>%
  group_by(species) %>%
  slice(which.max(evalue)) %>%
  mutate(div = dros_div[species, "V5"],
         residues = dros_res[species, "V2"])

ver_f1_final = ver_f1 %>%
  group_by(evalue, species) %>% 
  dplyr::summarize(f1score = f1_score(TP, FP, FN, TN),
                   mcc = round(mcc(TP=TP, FP=FP, FN=FN, TN=TN),2),
                   FPR=FPR) %>%
  group_by(species) %>%
  slice(which(FPR<0.05)) %>%
  group_by(species) %>%
  slice(which.max(evalue)) %>%
  mutate(div = vert_div[species, "V2"],
         residues = ver_res[species, "V2"])

#####
####

dros_df_tbl_opt <- filter(dros_df_tbl, paste(species, evalue) %in% paste(dros_f1_final$species,dros_f1_final$evalue))
toAdd <- dros_tblastn_match[match(dros_df_tbl_opt$species ,dros_tblastn_match$V2), "V1"]
toAdd[which(is.na(toAdd))] <- 0
dros_df_tbl_opt$found_tblastn_elsewhere <- (dros_df_tbl_opt$no_match - dros_df_tbl_opt$found_tblastn) - toAdd
dros_df_tbl_opt$speciesShort <- as.character(dros_df_tbl_opt$species)
dros_df_tbl_opt$speciesShort <- paste(substr(dros_df_tbl_opt$speciesShort, 1, 1), substr(unlist(strsplit(dros_df_tbl_opt$speciesShort, split = "_"))[2*(1:length(dros_df_tbl_opt$speciesShort))],1,3), sep='')
#
cer_df_tbl_opt <- filter(cer_df_tbl, paste(species, evalue) %in% paste(cer_f1_final$species,cer_f1_final$evalue))
toAdd <- cer_tblastn_match[match(cer_df_tbl_opt$species ,cer_tblastn_match$V2), "V1"]
toAdd[which(is.na(toAdd))] <- 0
cer_df_tbl_opt$found_tblastn_elsewhere <- (cer_df_tbl_opt$no_match - cer_df_tbl_opt$found_tblastn) - toAdd
cer_df_tbl_opt$speciesShort <- as.character(cer_df_tbl_opt$species)
cer_df_tbl_opt$speciesShort <- paste(substr(cer_df_tbl_opt$speciesShort, 1, 1), substr(unlist(strsplit(cer_df_tbl_opt$speciesShort, split = "_"))[2*(1:length(cer_df_tbl_opt$speciesShort))],1,3), sep='')

#
ver_df_tbl_opt <- filter(vert_df_tbl, paste(species, evalue) %in% paste(ver_f1_final$species,ver_f1_final$evalue))
toAdd <- vert_tblastn_match[match(ver_df_tbl_opt$species ,vert_tblastn_match$V2), "V1"]
toAdd[which(is.na(toAdd))] <- 0
ver_df_tbl_opt$found_tblastn_elsewhere <- (ver_df_tbl_opt$no_match - ver_df_tbl_opt$found_tblastn) - toAdd
ver_df_tbl_opt$speciesShort <- as.character(ver_df_tbl_opt$species)
ver_df_tbl_opt$speciesShort <- paste(substr(ver_df_tbl_opt$speciesShort, 1, 1), substr(unlist(strsplit(ver_df_tbl_opt$speciesShort, split = "_"))[2*(1:length(ver_df_tbl_opt$speciesShort))],1,3), sep='')

####
cer_df_tbl_opt <- arrange(cer_df_tbl_opt, div)
dros_df_tbl_opt <- arrange(dros_df_tbl_opt, div)
ver_df_tbl_opt <- arrange(ver_df_tbl_opt, div)

df <- bind_rows(cer_df_tbl_opt, dros_df_tbl_opt, ver_df_tbl_opt)
df$tag <- factor(df$tag, levels=c("yeast", "fruitfly", "human"), ordered = TRUE)

expTable <- mutate(df, foundPred = (found+found_tblastn)) %>%
  mutate(no_match = no_match - found_tblastn_elsewhere - found_tblastn) %>%
  dplyr::select(-found, -found_tblastn) %>%
  mutate(found_elsewhere = found_elsewhere + found_nn + found_tblastn_elsewhere) %>%
  select(-found_nn)

generEv <- read.csv("scripts_submission/final_scripts/file_with_general_evalues.csv")
expTable$general_evalue <- generEv[match(expTable$species, generEv$species), "evalue"]
expTable$total_checked <- 0
expTable[which(expTable$tag=="yeast"), "total_checked"] <- "5997"
expTable[which(expTable$tag=="fruitfly"), "total_checked"] <- "13929"
expTable[which(expTable$tag=="human"), "total_checked"] <- "19892"
temp_df <- rbind(as.data.frame(table(dros_orphan$species)), as.data.frame(table(cer_orphan$species)), as.data.frame(table(vert_orphan$species)))
temp_df$Var1 <- as.character(temp_df$Var1)
expTable$not_found_outside <- temp_df[match(expTable$species, temp_df$Var1), "Freq"]
expTable <- expTable[, c("tag", "species", "speciesShort", "div", "evalue", "general_evalue", "residues", "foundPred", "found_elsewhere", "no_match", "total", "not_found_outside", "total_checked")]

# remove C.elegans
expTable <- filter(expTable, species!="Caenorhabditis_elegans")
#

write.csv(expTable, "scripts_submission/final_scripts/Supplementary_Table_1.csv", row.names = FALSE)

colnames(expTable) <-c("dataset",
                       "target species",
                       "sp. abbrev.",
                       "div. time",
                       "phylostrat. E-value",
                       "general E-value",
                       "# residues",
                       "found opposite",
                       "found elsewhere",
                       "not found",
                       "total in microsynteny",
                       "not found and outside micro-synteny",
                       "total genes checked")
write.csv(expTable, "scripts_submission/final_scripts/Supplementary_Table_1_for_paper.csv", row.names = FALSE)

###


library(dplyr)
library(ggplot2)
library(mltools)
library(gridExtra)
library(reshape2)

setwd("~/Documents/research/Vakirlis_Carvunis_McLysaght_2019/")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



df_fn <- read.csv("Figure_3-source_data_1.csv")
df_fn$tag <- factor(df_fn$tag, levels=c("yeast", "fruitfly", "human"), ordered = TRUE)

#
df_fn <- filter(df_fn, species!="Caenorhabditis_elegans")
#

cer_df_tbl <- filter(df_fn, tag=="yeast")
dros_df_tbl <- filter(df_fn, tag=="fruitfly")
ver_df_tbl <- filter(df_fn, tag=="human")



df_fp <- read.csv("Figure_3-source_data_2.csv")
df_fp$tag <- factor(df_fp$tag, levels=c("yeast", "fruitfly", "human"), ordered = TRUE)

#
df_fp <- filter(df_fp, species!="Caenorhabditis_elegans")
#

cer_fp_df_tbl <- filter(df_fp, tag=="yeast")
dros_fp_df_tbl <- filter(df_fp, tag=="fruitfly")
ver_fp_df_tbl <- filter(df_fp, tag=="human")

####### F1 scores and ROC curves

cer_f1 <- mutate(cer_fp_df_tbl,
                 FPR = found/total,
                 TPR = 1-(cer_df_tbl$not_found/cer_df_tbl$total),
                 TP = cer_df_tbl$total - cer_df_tbl$not_found,
                 FP = found,
                 TN = total-found,
                 FN = cer_df_tbl$not_found)

dros_df_tbl <- filter(dros_df_tbl, evalue %in% dros_fp_df_tbl$evalue)
dros_f1 <- mutate(dros_fp_df_tbl,
                 FPR = found/total,
                 TPR = 1-(dros_df_tbl$not_found/dros_df_tbl$total),
                 TP = dros_df_tbl$total - dros_df_tbl$not_found,
                 FP = found,
                 TN = total-found,
                 FN = dros_df_tbl$not_found)

ver_f1 <- mutate(ver_fp_df_tbl,
                  FPR = found/total,
                  TPR = 1-(ver_df_tbl$not_found/ver_df_tbl$total),
                  TP = ver_df_tbl$total - ver_df_tbl$not_found,
                  FP = found,
                  TN = total-found,
                  FN = ver_df_tbl$not_found)


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
  slice(which(mcc == max(mcc))) %>%
  group_by(species) %>%
  slice(which.max(evalue))

dros_f1_final = dros_f1 %>%
  group_by(evalue, species) %>% 
  dplyr::summarize(f1score = f1_score(TP, FP, FN, TN),
                   mcc = round(mcc(TP=TP, FP=FP, FN=FN, TN=TN),2),
                   FPR=FPR) %>%
  group_by(species) %>%
  slice(which(mcc == max(mcc))) %>%
  group_by(species) %>%
  slice(which.max(evalue))

ver_f1_final = ver_f1 %>%
  group_by(evalue, species) %>% 
  dplyr::summarize(f1score = f1_score(TP, FP, FN, TN),
                   mcc = round(mcc(TP=TP, FP=FP, FN=FN, TN=TN),2),
                   FPR=FPR) %>%
  group_by(species) %>%
  slice(which(mcc == max(mcc))) %>%
  group_by(species) %>%
  slice(which.max(evalue))




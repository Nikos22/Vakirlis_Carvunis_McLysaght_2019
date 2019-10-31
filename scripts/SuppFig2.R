library(dplyr)
library(ggplot2)
library(mltools)
library(gridExtra)

setwd("~/Documents/research/ortho_buster/")

###

cer_res <- read.table("scripts_submission/no_residues/yeast_residues.txt", row.names=1)
dros_res <- read.table("scripts_submission/no_residues/droso_residues.txt", row.names=1)
ver_res <- read.table("scripts_submission/no_residues/vert_residues.txt", row.names=1)

vert_div <- read.table("scripts_submission/divergence_times/vert_div.txt", row.names=1)
dros_div <- read.table("scripts_submission/divergence_times/dros_div.txt", row.names=1)
cer_div <- read.table("scripts_submission/divergence_times/cer_div.txt", row.names=1)

df_fp <- read.csv("scripts_submission/final_scripts/fig3A2_table.csv", stringsAsFactors = FALSE)
cer_fp_df_tbl <- filter(df_fp, tag=="yeast")
dros_fp_df_tbl <- filter(df_fp, tag=="fruitfly")
vert_fp_df_tbl <- filter(df_fp, tag=="human")

df_fn <- read.csv("scripts_submission/final_scripts/fig3A1_table.csv", stringsAsFactors = FALSE)
cer_df_tbl <- filter(df_fn, tag=="yeast")
dros_df_tbl <- filter(df_fn, tag=="fruitfly")
vert_df_tbl <- filter(df_fn, tag=="human")

####### F1 scores and ROC curves

cer_f1 <- mutate(cer_fp_df_tbl,
                 FPR = (found)/total,
                 TPR = (1-((cer_df_tbl$not_found)/cer_df_tbl$total)),
                 TP = cer_df_tbl$total - cer_df_tbl$not_found,
                 FP = found,
                 TN = total-found,
                 FN = cer_df_tbl$not_found)

dros_df_tbl <- filter(dros_df_tbl, evalue %in% dros_fp_df_tbl$evalue)
dros_f1 <- mutate(dros_fp_df_tbl,
                 FPR = (found)/total,
                 TPR = (1-((dros_df_tbl$not_found)/dros_df_tbl$total)),
                 TP = dros_df_tbl$total - dros_df_tbl$not_found,
                 FP = found,
                 TN = total-found,
                 FN = dros_df_tbl$not_found)

vert_f1 <- mutate(vert_fp_df_tbl,
                 FPR = (found)/total,
                 TPR = (1-((vert_df_tbl$not_found)/vert_df_tbl$total)),
                 TP = vert_df_tbl$total - vert_df_tbl$not_found,
                 FP = found,
                 TN = total-found,
                 FN = vert_df_tbl$not_found)

f1_score <- function(tp, fp, fn, tn)
{
  precision = tp/(tp+fp)
  recall = tp/(tp+fn)
  return((2*(precision*recall))/(precision+recall))  
}

cer_f1_final = cer_f1 %>%
  group_by(evalue, species) %>% 
  dplyr::summarize(f1score = f1_score(TP, FP, FN, TN),
            mcc = round(mcc(TP=TP, FP=FP, FN=FN, TN=TN),2)) %>%
  group_by(species) %>%
  slice(which(mcc == max(mcc))) %>%
  group_by(species) %>%
  slice(which.max(evalue)) %>%
  mutate(div = cer_div[species, "V2"],
         residues = cer_res[species, "V2"])

dros_f1_final = dros_f1 %>%
  group_by(evalue, species) %>% 
  dplyr::summarize(f1score = f1_score(TP, FP, FN, TN),
            mcc = round(mcc(TP=TP, FP=FP, FN=FN, TN=TN),2)) %>%
  group_by(species) %>%
  slice(which(mcc == max(mcc))) %>%
  group_by(species) %>%
  slice(which.max(evalue)) %>%
  mutate(div = dros_div[species, "V5"],
         residues = dros_res[species, "V2"])

ver_f1_final = vert_f1 %>%
  group_by(evalue, species) %>% 
  dplyr::summarize(f1score = f1_score(TP, FP, FN, TN),
            mcc = round(mcc(TP=TP, FP=FP, FN=FN, TN=TN),2)) %>%
  group_by(species) %>%
  slice(which(mcc == max(mcc))) %>%
  group_by(species) %>%
  slice(which.max(evalue)) %>%
  mutate(div = vert_div[species, "V2"],
         residues = ver_res[species, "V2"])

########

dros_f1_final$tag <- "fruitfly"
cer_f1_final$tag <- "yeast"
ver_f1_final$tag <- "human"
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

df <- bind_rows(cer_f1_final, dros_f1_final, ver_f1_final)
df$tag <- factor(df$tag, levels=c("yeast", "fruitfly", "human"), ordered = TRUE)
write.csv(df, file="scripts_submission/final_scripts/file_with_general_evalues.csv", row.names = FALSE)

#
df <- filter(df, species!="Caenorhabditis_elegans")
#

########

cor_result_eval_noRes <- cor.test(-log(df$evalue, 10), df$residues)
# Pearson's product-moment correlation
# 
# data:  -log(df$evalue, 10) and df$residues
# t = 4.4109, df = 46, p-value = 6.151e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3089098 0.7180889
# sample estimates:
#       cor 
# 0.5451919 

cor_result_eval_div <- cor.test(-log(df$evalue, 10), df$div)
# Pearson's product-moment correlation
# 
# data:  -log(df$evalue, 10) and df$div
# t = 1.3726, df = 46, p-value = 0.1765
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.09089814  0.45675145
# sample estimates:
#       cor 
# 0.1983597 


eval_prot_pl <- ggplot() +
  geom_smooth(data = df,
              aes(x=-log(evalue,10),
                  y= residues),
              method="lm",
              colour="black",
              alpha=0.5) +
  geom_jitter(data = df,
                 aes(x= -log(evalue,10),
                     y= residues,
                     colour = tag,
                     shape = tag),
                 size=3,
                 stroke=1) +
  scale_colour_manual(values = cbPalette) +
  labs(y="Total residues in proteome",
       x="-log(E-value)") +
  theme(axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"))

###

eval_hist_pl <- ggplot() +
  geom_histogram(data = df,
                 aes(x= -log(evalue,10),
                     fill = tag),
                 position = position_dodge(),
                 binwidth=1) +
  scale_x_continuous(breaks = seq(0, 25, by=2)) +
  scale_fill_manual(values = cbPalette) +
  labs(y="Count",
       x="-log(E-value)") +
  theme(axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"))


eval_div_pl <- ggplot() +
  geom_smooth(data = df,
              aes(x=-log(evalue,10),
                  y= div),
              method="lm",
              colour="black",
              alpha=0.5) +
  geom_jitter(data = df,
              aes(x=-log(evalue,10),
                  y= div,
                  shape = tag,
                  colour=tag),
              stroke=1,
              size=3) +
  scale_colour_manual(values = cbPalette) +
  labs(x="-log(E-value)",
       y="Time since divergence (my)") +
  theme(axis.text = element_text(size=14, colour="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key.height = unit(1.4, "lines"))


pl <- grid.arrange(eval_hist_pl, eval_prot_pl, eval_div_pl,  ncol=1)

ggsave(plot = pl, "scripts_submission/final_scripts/SuppFig2.png", width=10, height = 10)


library(data.table)
library(ggplot2)
library(tidyr)
library(magrittr)
library(cowplot)
library(ggrepel)
library(ggpubr)

set.seed(4)

p_comd_coef_across_tissues <- readRDS('../plots/fig2/comd_coef_human_ndufb3.rds')
p_comd_coef_within_tissues <- readRDS('../plots/fig2/comd_coef_ndufb3_within_tissues.rds')
p_comd_coef_sdd <- readRDS('../plots/fig2/age_comd_coef.rds')


p_aranged <- ggarrange(p_comd_coef_across_tissues, NULL,NULL,p_comd_coef_within_tissues,NULL,p_comd_coef_sdd, nrow=2, ncol=3,widths = c(1, 0.08,1) )
p_aranged

save_plot('../plots/fig2/fig2_comd_coef_tpm_age_arranged.png',p_aranged, base_height=8, base_width=8)


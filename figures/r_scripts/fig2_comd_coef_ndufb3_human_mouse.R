library(data.table)
library(ggplot2)
library(tidyr)
library(magrittr)
library(cowplot)
library(ggrepel)
library(ggpubr)


data_base_path = '../figure_data/'
tpm_sdf_major_dt <- fread(paste0(data_base_path,'fig2/tpm_comd_major_tissue.csv'))
tpm_sdf_major_dt[, comd_coef_fc:=2**comd_coef]

p <- ggplot(tpm_sdf_major_dt, aes(comd_coef_fc, NDUFB3, label=tissue)) +
  geom_point(size=2, color="dodgerblue") +
  stat_cor(method = "spearman", size=3.4, cor.coef.name='rho', label.y = 2.05, label.x = 1) +
  geom_text_repel(max.overlaps =8) +
  geom_smooth(method='lm', color='blue') +
  theme_cowplot() +
  theme(panel.grid.major = element_line( size=.05, color="black" ))  +
  #scale_x_log10(breaks=c(0.7,.8, 1,1.2, 1.5, 2), limits=c(0.7,1.5)) +
  scale_y_log10() +
  labs(x='COMD coefficient [FC/codon/s]', y='NDUFB3 [TPM]') +
  annotation_logticks(sides = "l")
p

p <- ggplot(tpm_sdf_major_dt, aes(comd_coef_fc, NDUFB3, label=tissue)) +
  geom_point(size=2, color="dodgerblue") +
  geom_text_repel(max.overlaps =8) +
  geom_smooth(method='lm', color='blue') +
  theme_cowplot() +
  theme(panel.grid.major = element_line( size=.05, color="black" ))  +
  #scale_x_log10(breaks=c(0.7,.8, 1,1.2, 1.5, 2), limits=c(0.7,1.5)) +
  scale_y_log10() +
  labs(x='COMD coefficient [FC/codon/s]', y='NDUFB3 [TPM]') +
  annotation_logticks(sides = "l")
p
save_plot('../plots/fig2/comd_coef_human_ndufb3.svg',p, base_height=4.5, base_width=4.5)
save_plot('../plots/fig2/comd_coef_human_ndufb3.png',p, base_height=5, base_width=5)
saveRDS(p, '../plots/fig2/comd_coef_human_ndufb3.rds')

exp_comd_coef_mouse_dt <- fread(paste0(data_base_path,'fig2/exp_comd_coef_mouse.csv'))

p<- ggplot(exp_comd_coef_mouse_dt, aes(2**comd_coef, Ndufb3, label=tissue)) +
  geom_point(size=3, color="dodgerblue") +
  #stat_cor(method = "spearman", size=3.5, cor.coef.name='rho', label.x = 0.04) +
  geom_smooth(method='lm', color='blue') +
  geom_text_repel(size=3, max.overlaps = 4) +
  theme_cowplot() +
  theme(panel.grid.major = element_line( size=.05, color="black" ))  +
  scale_x_log10(breaks=c(.8, 1,1.2, 1.), limits=c(0.75,1.4)) +
  scale_y_log10() +
  labs(x='COMD coefficient [FC/codon/s]', y='Ndufb3 expression') 

p

save_plot('../plots/fig2/comd_coef_mouse_ndufb3.svg',p, base_height=5, base_width=5)
save_plot('../plots/fig2/comd_coef_mouse_ndufb3.png',p, base_height=5, base_width=5)


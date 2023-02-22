library(data.table)
library(ggplot2)
library(ggExtra)
library("RColorBrewer")
library(cowplot)
library(ggpubr)

data_base_path = '../figure_data/'
isch_comd_coef_dt <- fread(paste0(data_base_path,'fig3/comd_coef_isch_all_samples.csv'))


p <- ggplot(isch_comd_coef_dt[SMTS=='Esophagus'], aes(SMTSISCH, comd_coef_fc)) +
  geom_point(size=1, color="darkblue") +
  geom_smooth(method='lm', color='black') + 
  theme_cowplot() +
  theme(panel.grid.minor = element_line( size=.05, color="black" )) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~AGE, ncol=3) +
  labs(x='Ischemic time [min]', y='COMD coeffcient [FC/codon/s]', title='Esophagus') +
  stat_cor(method = "spearman", size=4, label.y=.28, cor.coef.name='rho', aes(label = ..r.label..))
p
save_plot('../plots/fig3/isch_time_comd_coef_esophagus.svg',p, base_height=5, base_width=4)
save_plot('../plots/fig3/isch_time_comd_coef_esophagus.png',p, base_height=5, base_width=4)
saveRDS(p, '../plots/fig3/isch_time_comd_coef_esophagus.rds')


library(data.table)
library(ggplot2)
library(ggExtra)
library("RColorBrewer")
library(cowplot)
library(gridExtra)
library(ggpubr)

data_base_path = '../figure_data/'
isch_comd_coef_dt <- fread(paste0(data_base_path,'fig3/comd_coef_isch_all_samples.csv'))

pdf('../plots/fig3/p_aranged_list_isch.pdf')
for(tissue in isch_comd_coef_dt[, unique(SMTSD)]){
p <- ggplot(isch_comd_coef_dt[SMTSD==tissue], aes(SMTSISCH, comd_coef_fc)) +
  geom_point(size=1, color="darkblue") +
  geom_smooth(method='lm', color='black') + 
  theme_cowplot() +
  theme(panel.grid.minor = element_line( size=.05, color="black" )) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~AGE, ncol=3) +
  labs(x='Ischemic time [min]', y='COMD coefficient (FC/codon/s)', title=tissue) +
  stat_cor(method = "spearman", size=4, label.y=.7, cor.coef.name='rho')
print(p)
}
dev.off()




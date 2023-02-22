library(data.table)
library(ggplot2)
library(ggExtra)
library("RColorBrewer")
library(cowplot)

data_base_path = '../figure_data/'
isch_comd_coef_dt <- fread(paste0(data_base_path,'fig3/corrected_comd_coef_vs_ischemic_time.csv'))


colorscale = scale_fill_gradientn(
  colors = rev(brewer.pal(9, "YlGnBu")),
  values = c(0, exp(seq(-5, 0, length.out = 100))))

p <- ggplot(isch_comd_coef_dt, aes(SMTSISCH, comd_coef_res)) +
  geom_point(color='transparent') +
  geom_hex() +
  labs(x='Ischemic time [min]', y='Adjusted COMD coefficient') +
  geom_smooth(method='lm') + 
  scale_x_log10() +
  theme_half_open() +
  background_grid() +
  colorscale
p
isch_comd_coef_dt[,cor.test(SMTSISCH, comd_coef_res, method='spearman')]
save_plot('../plots/fig3/ischemic_time_comd_coef.png',p, base_height=5, base_width=5.5)
save_plot('../plots/fig3/ischemic_time_comd_coef.svg',p, base_height=5, base_width=5.5)
saveRDS(p, '../plots/fig3/ischemic_time_comd_coef.rds')

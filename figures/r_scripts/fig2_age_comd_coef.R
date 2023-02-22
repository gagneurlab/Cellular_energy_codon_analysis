library(data.table)
library(ggplot2)
library(tidyr)
library(magrittr)
library(cowplot)
library(ggrepel)

data_base_path = '../figure_data/'
age_comd_coef_dt <- fread(paste0(data_base_path,'fig2/age_comd_coef_individual.csv'))


ggplot(age_comd_coef_dt, aes(age_value, effect_on_comd_coef)) +
  geom_point(color='darkblue', alpha=0.35) +
  geom_smooth(method='lm', color='blue') +
  stat_cor(method = "spearman", size=3.4, cor.coef.name='rho') +
  theme_cowplot() +
  theme(panel.grid.major = element_line( size=.05, color="black" ))  +
  labs(x='Age', y='Individual\'s COMD coefficient')


p <- ggplot(age_comd_coef_dt, aes(age_value, effect_on_comd_coef)) +
  geom_point(color='darkblue', alpha=0.35) +
  geom_smooth(method='lm', color='blue') +
  theme_cowplot() +
  theme(panel.grid.major = element_line( size=.05, color="black" ))  +
  labs(x='Age', y='Individual\'s COMD coefficient')
p

save_plot('../plots/fig2/age_comd_coef.svg',p, base_height=4, base_width=4)
save_plot('../plots/fig2/age_comd_coef.png',p, base_height=4, base_width=4)
saveRDS(p, '../plots/fig2/age_comd_coef.rds')

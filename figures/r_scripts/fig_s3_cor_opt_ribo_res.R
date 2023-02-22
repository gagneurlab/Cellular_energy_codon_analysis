library(data.table)
library(ggplot2)
library(ggExtra)
library("RColorBrewer")
library(cowplot)
library(ggpubr)

data_base_path = '../figure_data/'
ribo_res_dt <- fread(paste0(data_base_path,'fig3/cor_csc_dectime_ribo_res.csv'))
setnames(ribo_res_dt, c('mean A-site residency t=-5min', 'decoding rate (codon/s)'), 
         c('Mean A-site residency t=-5min', 'Decoding rate [codon/s]'))
ribo_res_melted_dt <- melt(ribo_res_dt[, .(`Mean A-site residency t=-5min`, CSC,`Decoding rate [codon/s]`,codon)], 
                           id.vars=c('codon','Mean A-site residency t=-5min'), variable.name = 'metric')

p1 <- ggplot(ribo_res_dt, aes(`Decoding rate [codon/s]`,`Mean A-site residency t=-5min`,  label=codon)) +
  geom_point(color='red') +
  geom_text() +
  geom_smooth(method='lm', color='blue') +
  theme_cowplot() +
  stat_cor(method = "spearman", cor.coef.name='rho', label.y=.6, label.x=5.5)
p1
p2 <- ggplot(ribo_res_dt, aes(CSC,`Mean A-site residency t=-5min`,  label=codon)) +
  geom_point(color='red') +
  geom_text() +
  geom_smooth(method='lm', color='blue') +
  theme_cowplot() +
  stat_cor(method = "spearman", cor.coef.name='rho', label.y=.6,label.x=-.05)
p2

p <- ggarrange(p1,p2)
p

save_plot('../plots/fig3/s3_cor_opt_ribo_res.svg',p, base_height=4, base_width=4)
save_plot('../plots/fig3/s3_cor_opt_ribo_res.png',p, base_height=4, base_width=8)
saveRDS(p, '../plots/fig3/s3_cor_opt_ribo_res.rds')

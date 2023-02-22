library(data.table)
library(ggplot2)
library(tidyr)
library(magrittr)
library(cowplot)
library(ggrepel)



data_base_path = '../figure_data/'
comd_coef_exon_df <- fread(paste0(data_base_path,'fig4/comd_coef_diff_spliced_exons_dec_rate.csv'))
comd_coef_exon_df[, comd_coef_fc:=2**comd_coef]
comd_coef_exon_df[psi_diff>=0.2, splicing:='Increased inclusion']
comd_coef_exon_df[psi_diff<= -0.2, splicing:='Increased skipping']
comd_coef_exon_df[abs(psi_diff)<=0.2, splicing:='not differentially spliced']
comd_coef_exon_df[, splicing:=factor(splicing, levels=c('Increased skipping', 'Increased inclusion', 'Not differentially spliced'), ordered = TRUE)]
comd_coef_exon_df[, avg_dec_rate:= mean(mtdr), by=.(tissue, splicing)]
comd_coef_avg_exon_df <- unique(comd_coef_exon_df[, .(avg_dec_rate, splicing, tissue,comd_coef_fc)])

p <-ggplot(comd_coef_avg_exon_df[splicing=='Increased inclusion'], aes(comd_coef_fc, avg_dec_rate)) +
  geom_point(size=2, color="darkblue") +
  geom_smooth(method='lm', color='black') + 
  theme_cowplot() +
  scale_x_log10(breaks=c(0.6, 0.7,.8, 1,1.2, 1.5, 2)) +
  labs(x='COMD coefficient [FC/codon/s]', y='average exon decoding rate') +
  stat_cor(method = "spearman", size=4, cor.coef.name='rho', alternative='greater')
p

p <-ggplot(comd_coef_avg_exon_df[splicing=='Increased skipping'], aes(comd_coef_fc, avg_dec_rate)) +
  geom_point(size=2, color="darkblue") +
  geom_smooth(method='lm', color='black') + 
  theme_cowplot() +
  scale_x_log10(breaks=c(0.7,.8, 1,1.2, 1.5, 2)) +
  labs(x='COMD coefficient [FC/codon/s]', y='average exon decoding rate [codon/s]') +
  stat_cor(method = "spearman", size=4, cor.coef.name='rho',alternative='less')
p


p <- ggplot(comd_coef_avg_exon_df[splicing=='Increased inclusion' |splicing=='Increased skipping' ], aes(comd_coef_fc, avg_dec_rate)) +
  geom_point(size=2, color="darkblue") +
  geom_smooth(method='lm', color='black') + 
  theme_cowplot() +
  scale_x_log10(breaks=c(0.7,.8, 1,1.2, 1.5, 2)) +
  facet_wrap(~splicing) +
  labs(x='COMD coefficient [FC/codon/s]', y='Mean exon decoding rate [codon/s]') 
p


save_plot('../plots/fig4/cassete_exon_comd_coef.svg',p, base_height=4, base_width=6)
save_plot('../plots/fig4/cassete_exon_comd_coef.png',p, base_height=4, base_width=6)
saveRDS(p, '../plots/fig4/cassete_exon_comd_coef.rds')

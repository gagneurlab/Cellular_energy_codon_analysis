library(data.table)
library(ggplot2)
library(tidyr)
library(magrittr)
library(cowplot)
library(ggrepel)
library("RColorBrewer")

data_base_path = '../figure_data/'
comd_coef_tpm_dt <- fread(paste0(data_base_path,'fig2/comd_coef_ndufb3_within_tissues.csv'))
comd_coef_tpm_dt[, comd_coef_fc:=2**comd_coef]
comd_coef_tpm_dt[SMTSD=='Skin - Not Sun Exposed (Suprapubic)', SMTSD:='Skin - Not Sun Exposed']
comd_coef_tpm_dt[SMTSD=='Skin - Sun Exposed (Lower leg)', SMTSD:='Skin - Sun Exposed']

tissues_of_interest=c('Brain - Cerebellar Hemisphere', 'Brain - Frontal Cortex (BA9)','Skin - Not Sun Exposed','Thyroid', 'Stomach','Adrenal Gland', 'Skin - Sun Exposed')

p<-ggplot(comd_coef_tpm_dt[SMTSD %in% tissues_of_interest], aes(comd_coef_fc, NDUFB3, color=SMTSD)) +
  geom_point(size=1, alpha=.37) +
  geom_smooth(method='lm') +
  scale_x_log10(breaks=c(0.3,.5, .8, 1, 1.5, 2, 3)) +
  scale_y_log10() +
  scale_color_manual(values =  brewer.pal(n = 7, name = 'Dark2'))+
  labs(x='COMD coefficient [FC/codon/s]', y='NDUFB3 [TPM]') +
  labs(color='') + 
  guides(colour = guide_legend(nrow = 4)) +
  theme_cowplot() +
  theme(legend.position = c(.025,.975), legend.text = element_text(size=8.5),legend.key.size = unit(.5,"line"), 
        panel.grid.major = element_line( size=.05, color="black" )) +
  annotation_logticks(sides='l')

p
save_plot('../plots/fig2/comd_coef_ndufb3_within_tissues.svg',p, base_height=4.3, base_width=4.3)
save_plot('../plots/fig2/comd_coef_ndufb3_within_tissues.png',p, base_height=4.3, base_width=4.3)
saveRDS(p, '../plots/fig2/comd_coef_ndufb3_within_tissues.rds')

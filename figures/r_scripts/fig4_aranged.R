library(data.table)
library(ggplot2)
library(tidyr)
library(magrittr)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(png)


p_skin_brain <- readRDS('../plots/fig4/skin_brain_cassette_dec_rates.rds')
p_cassette_sdf <- readRDS('../plots/fig4/cassete_exon_comd_coef.rds')


p_aranged <- ggarrange(p_skin_brain ,p_cassette_sdf, nrow=1, ncol=2, widths = c(1,1.2))
p_aranged
#dev.off()

save_plot('../plots/fig4/fig4_arranged.png',p_aranged, base_height=5, base_width=10.5)


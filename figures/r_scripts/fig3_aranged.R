library(data.table)
library(ggplot2)
library(tidyr)
library(magrittr)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(png)


p_atp_boxplots <- readRDS('../plots/fig3/atp_boxplots_codon_diff.rds')
p_atp_codon_diff <- readRDS('../plots/fig3/atp_codon_diff.rds')
p_atp_time_evol <- readRDS('../plots/fig3/atp_time_evol.rds')
p_codon_gcg <- readRDS('../plots/fig3/codon_merge_GCG_t2.rds')
p_isch_time_comd_coef_esophagus <- readRDS('../plots/fig3/isch_time_comd_coef_esophagus.rds')
p_isch_time_comd_coef <- readRDS('../plots/fig3/ischemic_time_comd_coef.rds')
isch_diagram <- readPNG("../figure_data/fig3/isch_diagram.png")

p_isch_diagram <- ggplot() + 
  background_image(isch_diagram) +
  # This ensures that the image leaves some space at the edges
  #theme(plot.margin = margin(l=1, r=1, unit = "cm")) +
  coord_fixed(.2)

p_isch_comd_coef_diagram <- ggarrange(p_isch_diagram, p_isch_time_comd_coef, nrow=2, ncol=1, heights=c(1,3))
p_isch_comd_coef_diagram_eso <- ggarrange(p_isch_comd_coef_diagram, p_isch_time_comd_coef_esophagus, widths = c(.8,1))
p_isch_comd_coef_diagram_eso


p_5pseq_cov <- ggarrange(p_atp_time_evol,p_atp_boxplots, align='h', nrow=1,
                         widths = c(.7,1))
p_5pseq_cov
p_5pseq_betas <- ggarrange(p_codon_gcg, p_atp_codon_diff, align='h', nrow=1,
                         widths = c(1, 1))
p_5pseq_betas


p_aranged <- ggarrange(p_isch_comd_coef_diagram_eso ,p_5pseq_cov, p_5pseq_betas,
                       nrow=3, ncol=1,heights=c(1.3,1,1.2))
#p_aranged
#dev.off()

save_plot('../plots/fig3/fig3_arranged_all.png',p_aranged, base_height=12.5, base_width=9.5)


library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

data_base_path = '../figure_data/'
comd_dt <- fread(paste0(data_base_path,'fig1/comd_coef_subtissues.csv'))
comd_dt[, comd_coef_fc :=2**comd_coef]

p <- ggplot(comd_dt, aes(comd_coef_fc, reorder(tissue,-comd_coef_fc))) +
   geom_point(color='darkgreen', size=3) + 
   scale_x_continuous(limits=c(.65,1.6),breaks=c(.7, 1, 1.3, 1.6)) +
   scale_y_discrete(position = "right") +
   theme_cowplot() +
   theme( axis.text = element_text(size = 12))  +
   background_grid(major = "x", minor = "none") +
   labs(x='COMD\ncoefficient', y='Tissue')

p

save_plot('../plots/fig1/comd_subtissues.png',p, base_height=9, base_width=5.5)
save_plot('../plots/fig1/comd_subtissues.svg',p, base_height=9, base_width=5.5)


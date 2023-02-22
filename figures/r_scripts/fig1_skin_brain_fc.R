library(data.table)
library(ggplot2)
library(ggExtra)
library("RColorBrewer")
library(cowplot)

data_base_path = '../figure_data/'
skin_brain_dt <- fread(paste0(data_base_path,'fig1/ei_centered_dec_rates_subtissues.csv'))
skin_brain_dt[, fc_skin_brain := 2**(`Skin - Not Sun Exposed (Suprapubic)` - `Brain - Frontal Cortex (BA9)`)]

colorscale = scale_fill_gradientn(
  colors = rev(brewer.pal(9, "YlGnBu")),
  values = c(0, exp(seq(-5, 0, length.out = 100))))


skin_brain_dt <- na.omit(skin_brain_dt, cols="fc_skin_brain")


# Add limits to plot and still show outliers

# cap values at lim depending on which upper or lower limit they exceed
lim <- 47
log2_lim <- log2(lim)
skin_brain_dt[, fc_skin_brain_capped := ifelse(abs(log2(fc_skin_brain)) > log2_lim, 2**(log2_lim*sign(log2(fc_skin_brain))), NA)] 

p <- ggplot(skin_brain_dt, aes(mtdr, fc_skin_brain)) +
  geom_point(color='transparent') +
  geom_point(aes(mtdr, fc_skin_brain_capped), size = 2L, shape=17) +
  geom_hex() +
  geom_smooth(method='lm') +
  labs(x='mRNA decoding rate\nin HEK293 [codon/s]', y='mRNA half-life ratio \n Skin Not S. E. /Brain Frontal C.') +
  scale_y_continuous(trans='log10', breaks=c(0.1,0.25,0.5,1,2,4, 10), limits=c(1/lim-0.0001,lim+0.0001)) +
  theme_half_open() +
  theme(legend.direction = "horizontal",legend.position = c(0.7, 1.1)) +
  background_grid() +
  colorscale 
p


p <- ggExtra::ggMarginal(p, type = "histogram")
p
save_plot('../plots/fig1/skin_brain_fc.svg',p, base_height=5, base_width=5.5)
save_plot('../plots/fig1/skin_brain_fc.png',p, base_height=5, base_width=5.5)

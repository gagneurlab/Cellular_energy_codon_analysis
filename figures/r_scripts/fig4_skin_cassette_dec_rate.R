library(data.table)
library(ggplot2)
library(tidyr)
library(magrittr)
library(cowplot)
library(ggrepel)
library(ggpubr)

data_base_path = '../figure_data/'
comd_coef_exon_df <- fread(paste0(data_base_path,'fig4/comd_coef_diff_spliced_exons_dec_rate.csv'))
comd_coef_exon_df[, comd_coef_fc:=2**comd_coef]
comd_coef_exon_df[psi_diff>0.2, splicing:='Increased\ninclusion']
comd_coef_exon_df[psi_diff< -0.2, splicing:='Increased\nskipping']
comd_coef_exon_df[abs(psi_diff)<=0.2, splicing:='Not differentially spliced']

comd_coef_exon_df[, splicing:=factor(splicing, levels=c('Increased\nskipping', 'Increased\ninclusion', 'Not differentially spliced'), ordered = TRUE)]

comd_coef_exon_df[, avg_dec_rate:= mean(mtdr), by=.(tissue, splicing)]
comd_coef_exon_skin_brain_df <- comd_coef_exon_df[tissue=="Sun Exposed (Lower leg) - Skin"|tissue=="Cerebellar Hemisphere - Brain"]
comd_coef_exon_skin_brain_df[tissue == "Sun Exposed (Lower leg) - Skin", tissue:="Skin - Sun Exposed"]
comd_coef_exon_skin_brain_df[tissue == "Cerebellar Hemisphere - Brain", tissue:="Brain - Cerebellar\nHemisphere"]
comd_coef_exon_skin_brain_df[, tissue:=factor(tissue, levels = c("Skin - Sun Exposed", "Brain - Cerebellar\nHemisphere"))]
comd_coef_exon_skin_brain_df <- comd_coef_exon_skin_brain_df[splicing!='Not differentially spliced']
comd_coef_exon_skin_brain_df[,  splicing_w_count:= paste0(splicing, "\n(n = ",.N,")"), by=.(splicing,tissue)]
comd_coef_exon_skin_brain_df[,  pval := wilcox.test(mtdr~splicing)$p.value, by=tissue]
light_green <- "#d6e8dd"
comd_coef_exon_skin_brain_df[, dec_rate_color:=ifelse(tissue=='Skin - Sun Exposed',
  ifelse(splicing=='Increased\ninclusion', "seagreen", light_green), 
  ifelse(splicing=='Increased\ninclusion', light_green, "seagreen"))]


p <-ggplot(comd_coef_exon_skin_brain_df, aes(splicing_w_count, mtdr)) +
  geom_boxplot(aes(fill=dec_rate_color)) +
  geom_text(data = unique(comd_coef_exon_skin_brain_df[, .(tissue, pval)]), aes(x=1.2, y=8.55,label = paste0("p = ", signif(pval,2)))) +
  scale_y_continuous(breaks=c(6, 6.5, 7.0, 7.5,8.0, 8.5))+
  theme_cowplot() +
  theme(panel.grid.major.x = element_blank(), axis.title.x=element_blank(),
        panel.grid.major.y = element_line( size=.05, color="black" ) 
  ) + 
  labs(y='Exon decoding rate [codon/s]') +
  facet_wrap(~tissue, scales = "free_x") +
  scale_fill_manual(values = c(light_green, "seagreen")) +
  guides(fill="none")
p


save_plot('../plots/fig4/skin_brain_cassette_dec_rates.svg',p, base_height=6, base_width=6.5)
save_plot('../plots/fig4/skin_brain_cassette_dec_rates.png',p, base_height=6, base_width=6.5)
saveRDS(p,'../plots/fig4/skin_brain_cassette_dec_rates.rds')


library(data.table)
library(ggplot2)
library(tidyr)
library(magrittr)
library(cowplot)
library(ggrepel)
library(dplyr)
library(ggsignif)
library(ggpubr)
library(colorspace)
library(RColorBrewer)


data_base_path = '../figure_data/'
all_betas_df <- fread(paste0(data_base_path,'fig3/all_fitted_betas.csv'))
all_betas_df[, time_point:=as.numeric(time_point)]
all_betas_df[, centered_beta_init:=mean_beta_init-mean(mean_beta_init)]
atp_dt <- data.table(time_point=c(-5,0.5,2,5,10),atp=c(3.2, 1.3, 1.2, 1.9, 2.2))

# Plotting betas relative to the mean (beta_codon_sample_i - avg_beta_codon_sample_i)

all_betas_df[, sample_centered_beta := beta_codon - mean(.SD[, beta_codon]), by=sample]
all_betas_df[, beta_relative_to_init := sample_centered_beta-centered_beta_init]


all_betas_fast_slow_df <- all_betas_df[(codon_group=='Fast codons')|(codon_group=='Slow codons')]

p <-ggplot(all_betas_fast_slow_df, aes(as.factor(time_point), exp(beta_relative_to_init),  fill=codon_group)) +
  geom_boxplot() +
  labs(fill='') +
  scale_y_log10() +
  scale_x_discrete(breaks=c(-5,0.5,2,5,10), labels = c("-5 min", "30 s", "2 min", "5 min", "10 min")) +
  scale_fill_manual(values=c("seagreen", "#d6e8dd"))+
  theme_cowplot()+
  labs(x='Time after addition of Antimycin A', y='Ribosome A-site residency\nrelative to t=-5min') 
p
p <- p + geom_pwc(
  aes(group=codon_group), tip.length = 0,
  method = "wilcox_test", label= "p.adj.signif", p.adjust.method = 'BH') 

p

save_plot('../plots/fig3/atp_boxplots_codon_diff.svg',p, base_height=5, base_width=8) 
save_plot('../plots/fig3/atp_boxplots_codon_diff.png',p, base_height=5, base_width=8) 
saveRDS(p, '../plots/fig3/atp_boxplots_codon_diff.rds')


#Plotting differences between medians:

median_diff_atp <- all_betas_fast_slow_df[, .SD[codon_group=="Slow codons", median(beta_relative_to_init)] - 
                                            .SD[codon_group=="Fast codons", median(beta_relative_to_init)], by=time_point]

setnames(median_diff_atp, 'V1', 'beta_slow_fast_diff')
median_diff_atp <- merge(median_diff_atp, atp_dt, by='time_point')

ggplot(median_diff_atp, aes(atp, beta_slow_fast_diff)) +
  geom_point()


p <- ggplot(median_diff_atp) +
   stat_cor(aes(atp, beta_slow_fast_diff), method = "pearson") +
   geom_point(aes(atp, beta_slow_fast_diff, color=factor(time_point)), size=3) +
   geom_smooth(aes(atp, beta_slow_fast_diff),method='lm') +
   theme_cowplot() +
   theme(panel.grid.major = element_line( size=.05, color="black" ))  +
   labs(x='ATP [umol/gDW]',y='Ribosome occupancy slow-fast', color='time\n[min]')
   
p

set.seed(123) 
all_betas_fast_slow_permut <- copy(all_betas_fast_slow_df[codon_group=="Slow codons"|codon_group=="Fast codons"])


get_median_diff_dt <- function(permut_dt, per_i=0){
  permut_dt <- permut_dt[, .SD[codon_group=="Slow codons", median(beta_relative_to_init)] - 
                           .SD[codon_group=="Fast codons", median(beta_relative_to_init)], by=time_point]
  setnames(permut_dt, "V1", paste("median_diff" ,per_i, sep='_'))
}

n_permuts <- 1000


for(i in 1:n_permuts){
  all_betas_fast_slow_permut <- copy(all_betas_fast_slow_df[codon_group=="Slow codons"|codon_group=="Fast codons"])
  all_betas_fast_slow_permut[, codon_group:=sample(codon_group)]
  if(i==1){
    median_diff_dt <- get_median_diff_dt(all_betas_fast_slow_permut, per_i = i)
    print(median_diff_dt)
  }
  else{
    median_diff_dt <- merge(median_diff_dt, get_median_diff_dt(all_betas_fast_slow_permut, per_i = i), by='time_point')
  }
}

median_diff_dt_melted <- melt(median_diff_dt, id.vars = 'time_point', value.name = 'median_diff', variable.name = 'per_number')
median_diff_sd <- median_diff_dt_melted[, sd(exp(median_diff)), by='time_point']
setnames(median_diff_sd, 'V1', 'median_std')
median_diff_atp <- merge(median_diff_atp, median_diff_sd, by='time_point')
median_diff_atp <-merge(median_diff_atp, data.table(time_point=c(-5,0.5,2,5,10), time_min=c('-5 min', '30 s', '2 min', '5 min', '10 min')))

median_diff_atp[, time_min:=factor(time_min, levels=c('-5 min', '30 s', '2 min', '5 min', '10 min'))]
p <- ggplot(median_diff_atp) +
  geom_point(aes(atp, exp(beta_slow_fast_diff), color=factor(time_min)), size=3) +
  geom_smooth(aes(atp, exp(beta_slow_fast_diff)), method='lm',se = FALSE) +
  geom_errorbar(aes(x=atp, ymin=exp(beta_slow_fast_diff)-median_std, ymax= exp(beta_slow_fast_diff)+median_std)) +
  theme_cowplot() +
  theme(panel.grid.major = element_line( size=.05, color="black" )) +
  labs(x = 'ATP [\U00B5mol/gDW]', y='Slow-to-fast ribosome A-site\nresidency ratio [a.u]', color='') 

p

save_plot('../plots/fig3/atp_codon_diff.svg',p, base_height=5, base_width=5) 
save_plot('../plots/fig3/atp_codon_diff.png',p, base_height=5, base_width=5) 
saveRDS(p, '../plots/fig3/atp_codon_diff.rds')



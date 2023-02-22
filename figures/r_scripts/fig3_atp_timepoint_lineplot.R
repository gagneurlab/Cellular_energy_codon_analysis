library(data.table)
library(ggplot2)
library(tidyr)
library(magrittr)
library(cowplot)
library(ggrepel)
library(dplyr)


data_base_path = '../figure_data/'
atp_dt <- data.table(time=c(-5,0,0.5,2,5,10),atp=c(3.2, 3.2, 1.3, 1.2, 1.9, 2.2))
atp_point_dt <- data.table(time=c(-5,0.5,2,5,10),atp=c(3.2, 1.3, 1.2, 1.9, 2.2))


p <- ggplot(atp_dt, aes(time, atp)) +
  geom_line(size=.8,colour='deepskyblue4') +
  geom_point(data=atp_point_dt, aes(time, atp), size=2.5, color='deepskyblue4') +
  scale_x_continuous(breaks=c(-5,0.5,2,5,10), labels = c("-5 min", "30 s", "2 min", "5 min", "10 min")) +
  theme_cowplot() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(x='Time after addition of Antimycin A', y='ATP')

p
save_plot('../plots/fig3/atp_time_evol.png',p, base_height=2.5, base_width=4) 
save_plot('../plots/fig3/atp_time_evol.svg',p, base_height=2.5, base_width=4)
saveRDS(p, '../plots/fig3/atp_time_evol.rds')


#Setup
library("tidyverse")
library(dplyr)
library("zoo") # Using index
library(readxl)
opts <- options(stringsAsFactors = F)
## data input: generates a list with samples, and keeps the frame stats, transcript descriptors and libsize for each sample. 
workdir <- "../figure_data/fig3"
fivepseq_count <- file.path(workdir,"ATP_merge")
rds <- file.path(workdir,"rds")

# Read codon_pauses_f1.txt file from 5pSeq output file
mycodon <- list()
codon_file <- function(whichframe="codon_pauses_f1.txt"){
  codon.files <- list.files(fivepseq_count, pattern = "^codon_pauses*",full.names = T, recursive = T)
  codon.files_list <- codon.files[which(basename(codon.files)==whichframe)]
  
  for (codon.f in codon.files_list) {
    codon <- read.table(codon.f, sep = "\t", header = T,  row.names=1)
    codon <- codon[,c(1:33)]
    colnames(codon) <- seq(-30,2)
    name <- basename(dirname(codon.f))
    mycodon[[name]] <- codon
  }
  return(mycodon)
}
codon_file_f1 <- codon_file(whichframe="codon_pauses.txt")

# # -----------------------------------------------------------------------
#Prepare an excel with sample name as the first column and total counts for each sample as the second row
total_codon_merge <- read_excel(file.path(workdir,"total_codon_merge.xlsx"))
total_codon_merge_list <- split(total_codon_merge$total_count,total_codon_merge$samples)
# # -----------------------------------------------------------------------
#Sum up all the count together, return with frame and count
GetFrameCount <- function(sample,counts) {
  mysum <- counts
  b <- sample %>% as_tibble() %>%
       dplyr::mutate_all(~ (. / mysum)*1000000) 
  b$codons <- rownames(sample)
  #b$codons <- unlist(sapply(rownames(sample),function(x) strsplit(as.character(x),split="_")[[1]][2]))
  b <- column_to_rownames(b,"codons")
  return(b)
}

codon_aver_f1 <- map2(codon_file_f1,total_codon_merge_list, GetFrameCount)
#saveRDS(codon_aver_f1,file.path(rds,"ATP_codon_aver_f1_indiv.rds"))
#saveRDS(codon_aver_f1,file.path(rds,"ATP_codon_aver_f1_merge.rds"))
codon_aver_f1 <- readRDS(file.path(rds,"ATP_codon_aver_f1_merge.rds"))


# # ---------Extracted samples for comparsion -----------------------------
codon_diff <- function(df,index_treat,index_contrl){
  codon_diff <- df[[index_treat]] / df[[index_contrl]] # using merged t2 - t-5
  codon_diff_mean <- codon_diff
  #codon_diff_mean <- codon_diff - rowMeans(codon_diff)
  col_pos <- as.numeric(colnames(df[[1]]))
  #col_pos_sub <- which(col_pos%%3 == 1)
  
  codon_diff_mean <- codon_diff_mean[c("ALA_GCG","ALA_GCT"),] # change the codon you want to compare
  codon_diff_mean <- as.data.frame(t(codon_diff_mean))
  
  codon_diff_mean <- rownames_to_column(codon_diff_mean,"pos")
  data_long <- gather(codon_diff_mean, codon, rpm, ALA_GCG:ALA_GCT, factor_key=TRUE)
  data_long$pos <- as.numeric(data_long$pos) 
  return(data_long)
}
df <- codon_diff(codon_aver_f1,4,1)


# # ---------Extracted the raw rpm for non-optimal codons -----------------------------
codon_rpm <- function(df,index_treat,index_contrl,codon="ALA_GCG"){
          ATP_low <- as.data.frame(t(df[[index_treat]]["ALA_GCG",]))
          ATP_high <- as.data.frame(t(df[[index_contrl]]["ALA_GCG",]))

          ATP_low$group <- "ATP_low"
          ATP_high$group <- "ATP_high"

          ATP_low <- rownames_to_column(ATP_low,"pos")
          ATP_high <- rownames_to_column(ATP_high,"pos")
          
          rownames(ATP_low) <- NULL
          rownames(ATP_high) <- NULL
          
          ATP <- rbind(ATP_low,ATP_high)
          ATP$pos <- as.numeric(ATP$pos)
          
          return(ATP)
}
ATP <- codon_rpm(codon_aver_f1,4,1,codon="ALA_GCG")

# # ----------Plot the codon raw for non-optimal codons -----------------------------------------------------
p<-ggplot(ATP, aes(x=pos, y=ALA_GCG, group=group)) +
  geom_line(aes(color=group))+
  #scale_x_continuous(breaks=c(-30, -20, -17, -10, 0)) +
  theme_cowplot() + 
  theme(legend.position = c(0.45, 1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y='5Pseq read count \nat codon GCG (rpm)', color='',x='') +
  scale_color_manual(labels = c("-5 min (ATP high)", "2 min (ATP low)"), 
                     values = c("ATP_low" = "#4D95B2FF", "ATP_high" = "#F8A16DFF")) 
  
p
save_plot('../plots/fig3/codon_merge_GCG_t2.png',p, base_height=3, base_width=5) 
save_plot('../plots/fig3/codon_merge_GCG_t2.svg',p, base_height=3, base_width=5) 
saveRDS(p, '../plots/fig3/codon_merge_GCG_t2.rds')







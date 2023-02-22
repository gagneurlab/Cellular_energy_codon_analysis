# -*- coding: utf-8 -*-
library(data.table)
library(biomaRt)
library(magrittr)
library(plyr)
library(ggplot2)
library(MASS)
library(gravity)

#Load RDS object with reads per position and sample:
data_path <- '../data/'
counts_obj <- readRDS(paste0(data_path, "2202_meta_fivepseq.count.rds"))
ensembl <- useMart("ensembl",dataset="scerevisiae_gene_ensembl")


get_counts_features <- function(condition, out_file_name, n_genes_to_sample=1000, include_all_genes=F, min_gene_counts=10, 
                                n_first_codons_to_skip=0, n_last_codons_to_skip=0){

  counts_obj_condition <- counts_obj[[condition]]
  
  coding_sequences_dt <- as.data.table(getSequence(id = counts_obj_condition$t.ass$ID, 
                                                   type="ensembl_gene_id", 
                                                   mart = ensembl, 
                                                   seqType = "coding"))
  
  #Get for every codon position (rows), the distance to start codon, the count number and the gene
  
  number_of_genes <- length(counts_obj_condition$count.f)
  gene_counts_dt_list <- vector(mode = "list", length = number_of_genes)
  
  for (gene_i in 1:number_of_genes) {
    
    # Get ID of the gene and its counts
    gene_i_counts <- counts_obj_condition$count.f[[gene_i]]
    gene_id <- counts_obj_condition$t.ass[gene_i, 'ID']
    
    # Get gene coding sequence length
    cds_gene <- coding_sequences_dt[ensembl_gene_id==gene_id,coding]
    cds_full_length <- nchar(cds_gene)
    
    # The counts correspond to the number of 
    #read 5’ endpoints mapping to each position in the transcript. The coordinates 
    #span the CDS and -100/+100 nucleotides around it. The ribosome protects 
    #the fragment on A-site -17 positions. For the start codon the tRNA methionine
    # enters in the P-site, which corresponds to the A-site of the codon after the 
    # start one. I will not consider the reads mapping to this second codon 
    # because they are affected by the initiation (Will only start considering reads -17 + 3 + 3 = -11 positions away from the start of the CDS).
    # I will also not consider the reads corresponding to the termination codon and the one immediately before it.
    
    # Get range of positions of5Pseq reads corresponding to the codons of interest
    counts_interest_start <- 101 - 11 + 3*n_first_codons_to_skip #coding sequence starts at position 101 in the vector (1-based)
    counts_interest_end <- 101 + (cds_full_length - 17 - 3) - 3  - 3*n_last_codons_to_skip
    
    if (counts_interest_end <= counts_interest_start){
      print('gene CDS is too small for specified codons to skip. Gene will be removed.')
      next 
    }
      
    # Get reads for the codons of intererest only
    codon_proteceted_positions <- seq(from=counts_interest_start, to=counts_interest_end, by=3)
    codon_pos_counts <- counts_obj_condition$count.f[[gene_i]][codon_proteceted_positions]
    
    # Get position of codon relative to the start of the codons considered in the coding sequence
    pos_rel_start <- codon_proteceted_positions - counts_interest_start + 6 + 3*n_first_codons_to_skip
    
    
    gene_counts_i_dt <- data.table(gene_id=gene_id, 
                                   position_from_start=pos_rel_start, 
                                   counts=codon_pos_counts)
    gene_counts_dt_list[gene_i] <- list(gene_counts_i_dt)
    
    if (gene_i %% 100 == 0) 
      print(gene_i)
  }
  gene_pos_counts_dt <- rbindlist(gene_counts_dt_list)
  
  gene_counts_dt <- gene_pos_counts_dt[, sum(counts), by=gene_id]
  setnames(gene_counts_dt, "V1", "total_gene_codon_counts")
  
  ggplot(gene_counts_dt, aes(total_gene_codon_counts)) +
    geom_histogram() +
    scale_x_log10() 
  
  print(paste('number of genes with more than min counts:', 
        as.character(gene_counts_dt[total_gene_codon_counts>min_gene_counts,.N])))
  
  #Discard genes with no or few counts
  genes_with_min_reads <- gene_counts_dt[total_gene_codon_counts>min_gene_counts]$gene_id
  if(include_all_genes==F){
    genes_with_min_reads <- sample(genes_with_min_reads, size=n_genes_to_sample)
  }
  gene_pos_counts_dt <- gene_pos_counts_dt[gene_id %in% genes_with_min_reads]
  
  #Only for expressed genes, get for every codon position (rows), the distance to start codon, the codon type and the gene
  coding_sequences_dt[ensembl_gene_id %in% genes_with_min_reads]
  counts_features_dt <- merge(gene_pos_counts_dt, coding_sequences_dt, by.x='gene_id', by.y='ensembl_gene_id')
  
  counts_features_dt[, a_site_codon := substr(coding, start=position_from_start+1, stop=position_from_start+3)]
  counts_features_dt[, coding:=NULL]
  write.csv(counts_features_dt, out_file_name)


}

for (condition in names(counts_obj)){
  include_all_genes<-T
  counts_features_folder <- paste0(data_path,'5p_counts_with_features/')
  #out_file_name <- paste(condition, n_genes_to_sample, 'genes', 'counts', 'features', sep='_')
  out_file_name <- paste(condition, 'counts_features', sep='_')
  get_counts_features(condition = condition, min_gene_counts = 10, include_all_genes=include_all_genes,
                      out_file_name = paste0(counts_features_folder, out_file_name,'.csv'),
                      n_first_codons_to_skip=0, n_last_codons_to_skip=0)

}



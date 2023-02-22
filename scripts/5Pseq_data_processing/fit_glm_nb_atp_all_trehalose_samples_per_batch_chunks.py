# Batchwise Fitting a model to 5Pseq reads to estimate codon dependent ribosome occupancy in 5Pseq fragments

import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
import numpy as np
import os
from plotnine import *
from scipy.stats import pearsonr, spearmanr
import seaborn as sns
import matplotlib.pyplot as plt
import time
np.random.seed(1234)


# # Trehalose time series

print('loading data...')
data_path = '../../data/5p_counts_with_features/'
all_counts_df = pd.read_csv(data_path+'trehalose_series_features_with_counts.csv', index_col=0)



samples = all_counts_df['sample'].unique()
n_samples = len(samples)



reference_gene = 'YDR077W'


all_counts_df['batch'] = all_counts_df['sample'].apply(lambda el: int(el.split('_')[-1]))

print('data was loaded')

def fit_glm_nb_sample(all_counts_df, sample_name, gene_chunk_ids, reference_gene, reference_codon='AAA', chunk_number=None):
    
    counts_df = all_counts_df[all_counts_df['sample']==sample_name]
    chunk_i = counts_df[counts_df.gene_id.isin(gene_chunk_ids)].reset_index(drop=True)
    
    #Model fitting
    glm_binom = smf.glm(formula=f"counts ~ C(a_site_codon, Treatment(reference='{reference_codon}')) + C(gene_id, Treatment(reference='{reference_gene}')) + position_from_start",
                    data=chunk_i, family=sm.families.NegativeBinomial(link=sm.families.links.log())).fit()


    # Get coefficients and process dataframe
    coefs_df = glm_binom.summary2().tables[1]
    coefs_df.index.name = 'feature'
    coefs_df.reset_index(inplace=True)
    coefs_df.feature = coefs_df.feature.apply(lambda el: el.split('T.')[1][:-1] if 'T.' in el else el)
    coefs_df.rename({'Coef.':'beta'}, axis=1, inplace=True)
    coefs_df.set_index('feature',inplace=True)
    coefs_df.columns = coefs_df.columns + f'_{sample_name}' if chunk_number is None else coefs_df.columns + f'_{sample_name}_{chunk_number}'
    
    return coefs_df


reference_codon='AAA'
coefs_df_list = []
n_chunks = 3
print('N chunks: ', n_chunks)

for batch_to_fit in range(1,4):
    
    print(f'Batch {batch_to_fit}')
        


    for sample in all_counts_df.loc[all_counts_df['batch']==batch_to_fit, 'sample'].unique():

        print(f'Sample {sample}')
        all_genes_sample = np.array(all_counts_df.loc[all_counts_df['sample']==sample, 'gene_id'].unique())
        np.random.shuffle(all_genes_sample)
        gene_chunks = np.array_split(all_genes_sample, n_chunks)
        
        for i, gene_chunk in enumerate(gene_chunks):
            print('chunk ', i)
            if reference_gene not in gene_chunk:
                  genes_to_fit = np.append(gene_chunk, reference_gene)
            else:
                genes_to_fit = gene_chunk

            # Fitting
            start = time.time()
            coefs_df_list.append(fit_glm_nb_sample(all_counts_df, sample, gene_chunk_ids=genes_to_fit,
                                                   reference_gene=reference_gene, reference_codon=reference_codon, chunk_number=i))
            end = time.time()
            print(f'sample {sample} concluded. {end-start}s elapsed')

            
series_coefs_df = pd.concat(coefs_df_list, axis=1)
series_coefs_df.to_csv(f'trehalose_series_{n_chunks}_chunks_betas_all_batches.csv')
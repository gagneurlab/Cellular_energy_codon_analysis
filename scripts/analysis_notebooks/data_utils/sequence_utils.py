import numpy as np
import pandas as pd
from itertools import product
from plotnine import *
import scipy.stats as stats
import re

ALL_CODONS_NAMES = np.array(['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA',
                             'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC',
                             'CGG',
                             'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
                             'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA',
                             'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'])

ALL_CODONS_DF = pd.DataFrame(data=np.zeros(64), index=ALL_CODONS_NAMES, columns=['codon number'])


def shift_sequence(seq, n):
    seq = seq[n:]+seq[0:n]
    return seq


def codons_count(sequence, frame=None, codon_shift=None):

    if codon_shift is not None:
        sequence = shift_sequence(sequence, codon_shift)
    if frame is not None:
        sequence = sequence[frame:]

    codon_list = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
    [existing_codons, existing_codons_count] = np.unique(codon_list, return_counts=True)
    existing_codons_df = pd.DataFrame(data=existing_codons_count, index=existing_codons, columns=['codon number'])
    existing_codons_df = ALL_CODONS_DF.add(existing_codons_df, axis='index', fill_value=0)

    return existing_codons_df


def rc_dna(seq):
    """
    Reverse complement the DNA sequence
    """
    rc_hash = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }
    return "".join([rc_hash[s] for s in reversed(seq)])


def cut_seq(seq):
    while (len(seq) % 3 != 0):
        seq = seq[:-1]
    return seq


def get_seq_gc_content(seq_str):
    return (seq_str.count('G') + seq_str.count('C')) / len(seq_str)


def get_seq_n_content(seq_str, n):
    return (seq_str.count(n) / len(seq_str))


class MotifUtils:

    def __init__(self, sequences_df):

        self.sequences_df = sequences_df

    def motif_median_fold_change(self, motif, calculate_on, mrna_section, motif_searched=False, return_medians=False):

        if not motif_searched:
            self.sequences_df['motif_in'] = self.sequences_df[mrna_section].str.contains(motif)

        hf_median_motif = np.median(self.sequences_df.loc[self.sequences_df['motif_in'] == True, calculate_on]
                                    .dropna())
        hf_median_non_motif = np.median(self.sequences_df.loc[self.sequences_df['motif_in'] == False, calculate_on]
                                        .dropna())

        median_fold_change = hf_median_motif - hf_median_non_motif
        if return_medians:
            return median_fold_change, hf_median_motif, hf_median_non_motif
        else:
            return median_fold_change


    def all_motifs_median_fold_change(self, motif_length, calculate_on, mrna_section, plot=True):
        mfc_vector = []
        motif_vector = []

        for motif_tuple in product(['A', 'C', 'T', 'G'], repeat=motif_length):
            motif_str = ''.join(motif_tuple)
            mfc_vector.append(self.motif_median_fold_change(motif_str, calculate_on=calculate_on, mrna_section
            =mrna_section))
            motif_vector.append(motif_str)

        mfc_df = pd.DataFrame({'median_fold_change': mfc_vector, 'motif': motif_vector})

        if plot:
            hist = (ggplot(mfc_df, aes('median_fold_change'))
                    + geom_histogram()
                    + labs(x='Median Fold Change', y='Count'))
            fig = hist.draw()
            fig.show()

        return mfc_df

    def motif_boxplot(self, motif, y_axis, mrna_section, motif_searched=False):

        if not motif_searched:
            self.sequences_df['motif_in'] = self.sequences_df[mrna_section].str.contains(motif)

        return (ggplot(self.sequences_df.loc[:, ['motif_in', y_axis]].dropna(),
                       aes('factor(motif_in)', y_axis))
                + geom_boxplot()
                + labs(title=motif, x='Motif number')
                + scale_x_discrete(labels=['0', '>=1'])
                )

    def wilcox_ranksum_test(self, motif, calculate_on, mrna_section, motif_searched=False):

        if not motif_searched:
            self.sequences_df['motif_in'] = self.sequences_df[mrna_section].str.contains(motif)

        motif_transcripts_hl = self.sequences_df.loc[self.sequences_df['motif_in'] == True, calculate_on]\
            .dropna().values
        non_motif_transcripts_hl = self.sequences_df.loc[self.sequences_df['motif_in'] == False, calculate_on]\
            .dropna().values

        return stats.ranksums(motif_transcripts_hl, non_motif_transcripts_hl)

    def motif_info(self, motif, calculate_on, mrna_section):

        self.sequences_df['motif_in'] = self.sequences_df[mrna_section].str.contains(motif)

        motif_mfc, motif_median, non_motif_median = self.\
            motif_median_fold_change(motif=motif, calculate_on=calculate_on, mrna_section=mrna_section,
                                     motif_searched=True, return_medians=True)

        motif_wilcox = self.wilcox_ranksum_test(motif=motif, calculate_on=calculate_on,
                                                mrna_section=mrna_section, motif_searched=True)
        boxplot = self.motif_boxplot(motif=motif, y_axis=calculate_on, mrna_section=mrna_section,
                                     motif_searched=True)
        drawn_boxplot = boxplot.draw()
        drawn_boxplot.show()

        non_na_sequences_df = self.sequences_df.loc[:, ['motif_in', calculate_on]].dropna()
        n_transcripts_w_motif = len(non_na_sequences_df.loc[non_na_sequences_df['motif_in'] == True, :])
        n_transcripts_no_motif = len(non_na_sequences_df.loc[non_na_sequences_df['motif_in'] == False, :])

        print('Median Fold Change:', motif_mfc)
        print('Wilcox ranksum test:', motif_wilcox)
        print(f'Number of transcripts with the motif:{n_transcripts_w_motif}\tmedian: {motif_median}\n '
              f'Number of transcripts without the motif:{n_transcripts_no_motif}\tmedian: {non_motif_median}')

    def find_motif_positions(self, sequence, motif):

        positions = []
        for match in re.finditer(motif, sequence):
            positions.append(match.start())

        return positions

    def get_motif_postitions(self, motif, mrna_section):

        pos_series = self.sequences_df[mrna_section].apply(self.find_motif_positions, motif=motif)
        pos_series = pos_series.explode().dropna().astype('int64')  # explode turns list elements into seperate rows,
        #  mantaining the same index. So if a transcript has 2 instances of a motif, it will appear 2 times in the
        #  pos series with 2 different position numbers
        pos_df = pos_series.to_frame()
        pos_df.rename({mrna_section: 'pos'}, inplace=True, axis='columns')
        return pos_df


class MultiTissueMotifUtils(MotifUtils):
    def __init__(self, sequences_df, tissues):
        super().__init__(sequences_df)
        self.tissues = tissues

    def all_kmers_median_fold_change(self, kmer_length, mrna_section, plot=True):
        all_kmers_df = None

        for tissue in self.tissues:
            mfc_df = super().all_motifs_median_fold_change(motif_length=kmer_length, calculate_on=tissue,
                                                           mrna_section=mrna_section, plot=False)

            mfc_df.rename(columns={'median_fold_change': tissue}, inplace=True)
            mfc_df.set_index('motif', inplace=True)

            if all_kmers_df is None:
                all_kmers_df = mfc_df
            else:
                all_kmers_df = all_kmers_df.merge(mfc_df, how='outer', left_index=True, right_index=True)

        if plot:
            self.plot_all_tissues_kmers_mfc(all_kmers_df)

        return  all_kmers_df

    def plot_all_tissues_kmers_mfc(self, all_kmers_df):

        mfc_df_melted = all_kmers_df.reset_index().melt(var_name='tissue', value_name='median_fold_change',
                                                        id_vars='motif')
        plot = (ggplot(mfc_df_melted, aes('median_fold_change'))
                     +geom_histogram(bins=50)
                     +facet_wrap('~tissue', scales='free')
                     + theme(figure_size=(20, 20))
                     + labs(x='kmer median fold change', y='Counts')
                )
        plot_drawn = plot.draw()
        plot_drawn.show()

    def rank_motifs_func(self, group_df, n):
        group_df = group_df.sort_values(by='median_fold_change', ascending=False).iloc[:n, 1:]
        # Drop previous index values of the main df
        group_df = group_df.reset_index().drop('index', axis='columns')
        group_df.index.rename('rank', inplace=True)
        return group_df

    def rank_motifs_per_tissue(self, all_kmers_df, n_rank):
        abs_mer_df = np.abs(np.log2(all_kmers_df)).T #Log in order to be able to rank the motifs by abs(log fold change)
        # this way a motif with median higher than the transcripts with non motif will be considered also
        abs_mer_df.index.rename('tissue', inplace=True)
        abs_mer_df_melted = abs_mer_df.reset_index().melt(id_vars='tissue', value_name='median_fold_change')
        ranked_mer_df = abs_mer_df_melted.groupby('tissue').apply(self.rank_motifs_func, n_rank)

        return ranked_mer_df

    def rank_metrics_func(self, group_df):
        group_df['tissues_in'] = group_df['tissue'].count()
        group_df['median_rank'] = group_df['rank'].median()
        group_df['median_median_fold_change'] = group_df['median_fold_change'].median()
        group_df = group_df.drop(['tissue', 'rank', 'median_fold_change', 'motif'], axis='columns')
        group_df = group_df.drop_duplicates()
        return group_df

    def kmer_info(self, all_kmers_df, n_rank):

        ranked_kmer_df = self.rank_motifs_per_tissue(all_kmers_df, n_rank)
        rank_info = ranked_kmer_df.reset_index().groupby('motif').apply(self.rank_metrics_func)
        rank_info = rank_info.droplevel(1).sort_values('median_rank')

        return rank_info




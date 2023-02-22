from data_utils.concise_utils import encodeSequence
from data_utils.sequence_utils import get_seq_gc_content, codons_count, get_seq_n_content
import numpy as np
import pandas as pd
import pickle


# Utilities for handling loaded data on gene sequences

class DataHandler:

    def __init__(self, loaded_data_path):

        self.data = pickle.load(open(loaded_data_path, "rb"))

        self.all_transcript_ids = [transcript['transcript_id'] for transcript in self.data]
        self.aug_start_transcripts = [transcript['transcript_id'] for transcript in self.data if
                                      transcript['non_AUG'] == 0]

    def get_sequence_onehot_df(self, mrna_section, utr_filter_size=0, stop_codon_in_3_utr=False,
                               include_non_aug_start=False, cut_n_transcripts=None, maxlen=None):
        """

        :param mrna_section: one of 3_utr, 5_utr, cds (coding sequence) or all
        :param utr_filter_size: number of minimum bases in the utr in order to include a certain transcript's utr
        :param stop_codon_in_3_utr: include the stop codon in the 3 utr or not
        :param include_non_aug_start: include transcripts which do not start with an AUG in the coding sequence
        :param cut_n_transcripts: number of transcripts you want to return, if None returns all
        :param maxlen: length of padding
        :return: dataframe with transcript ids as indexes and one hot numpy vectors in the rows for the specific
        mrna_section

        """
        sequences_df = self.get_sequence_df(mrna_section, utr_filter_size=utr_filter_size,
                                            stop_codon_in_3_utr=stop_codon_in_3_utr,
                                            include_non_aug_start=include_non_aug_start)

        if cut_n_transcripts is not None:
            sequences_df = sequences_df.iloc[:cut_n_transcripts]

        encoded_seq = encodeSequence(sequences_df[mrna_section].values, vocab=["A", "C", "G", "T"], neutral_vocab=["N"],
                                     maxlen=maxlen)
        sequences_df[mrna_section] = [[encoded_seq[i]] for i in range(encoded_seq.shape[0])]

        return sequences_df

    def get_sequence_df(self, mrna_section, utr_filter_size=0, stop_codon_in_3_utr=False, include_non_aug_start=False):
        """

        :param mrna_section: one of 3_utr, 5_utr, cds (coding sequence) or all
        :param utr_filter_size: number of minimum bases in the utr in order to include a certain transcript's utr
        :param stop_codon_in_3_utr: include the stop codon in the 3 utr or not
        :param include_non_aug_start: include transcripts which do not start with an AUG in the coding sequence
        :return: dataframe with transcript ids as indexes and string sequence in the rows for the specific mrna_section
        """

        if 'utr' in mrna_section:

            stop_codons = ['TAA', 'TAG', 'TGA']

            data_utr = [transcript for transcript in self.data if transcript['utrs'][mrna_section] is not None]

            if (mrna_section == '3_utr') & (not stop_codon_in_3_utr):
                data_utr = [transcript for transcript in data_utr if
                            len(transcript['utrs'][mrna_section]) > max(utr_filter_size + 3, 3)]
                utr_list = [transcript['utrs'][mrna_section][3:]  # Don't get stop codon
                            if (transcript['utrs'][mrna_section][:3] in stop_codons)
                            else transcript['utrs'][mrna_section] for transcript in data_utr]  # Get utr if stop codon
                # is not present
            else:
                data_utr = [transcript for transcript in data_utr if
                            len(transcript['utrs'][mrna_section]) > utr_filter_size]
                utr_list = [transcript['utrs'][mrna_section] for transcript in data_utr]

            transcript_ids = [transcript['transcript_id'] for transcript in data_utr]
            sequences_df = pd.DataFrame({mrna_section: utr_list}, index=transcript_ids)

        elif mrna_section == 'cds':

            cds_list = [transcript[mrna_section] for transcript in self.data]
            transcript_ids = self.all_transcript_ids

            sequences_df = pd.DataFrame({mrna_section: cds_list}, index=transcript_ids)

        elif mrna_section == 'all':

            utr_5_df = self.get_sequence_df(mrna_section='5_utr', utr_filter_size=utr_filter_size)
            cds_df = self.get_sequence_df(mrna_section='cds')
            utr_3_df = self.get_sequence_df(mrna_section='3_utr', utr_filter_size=utr_filter_size,
                                            stop_codon_in_3_utr=True)  # get stop codon

            merge_seq_df = utr_5_df.merge(cds_df, how='right', left_index=True, right_index=True)
            merge_seq_df = merge_seq_df.merge(utr_3_df, how='left', left_index=True, right_index=True)
            merge_seq_df['all'] = merge_seq_df['5_utr'].fillna("") + merge_seq_df['cds'] + \
                                  merge_seq_df['3_utr'].fillna("")
            merge_seq_df.drop(labels=['5_utr', 'cds', '3_utr'], axis=1, inplace=True)

            sequences_df = merge_seq_df

        if not include_non_aug_start:
            sequences_df = sequences_df.loc[sequences_df.index.intersection(pd.Index(self.aug_start_transcripts))]

        return sequences_df

    def get_n_content(self, mrna_section, n, include_non_aug_start=False):

        seq_df = self.get_sequence_df(mrna_section=mrna_section, include_non_aug_start=include_non_aug_start)
        seq_df[n+'_content_' + mrna_section] = seq_df[mrna_section].apply(get_seq_n_content, n=n)

        return seq_df.loc[:, [n+'_content_' + mrna_section]]

    def get_gc_content(self, mrna_section, include_non_aug_start=False):

        """

        :param mrna_section: one of 3_utr, 5_utr, cds (coding sequence) or all
        :param include_non_aug_start: include transcripts which do not start with an AUG in the coding sequence
        :return: dataframe with transcript ids as indexes and gc content on each row for each specific mrna_section
        """

        seq_df = self.get_sequence_df(mrna_section=mrna_section, include_non_aug_start=include_non_aug_start)
        seq_df['gc_content_' + mrna_section] = seq_df[mrna_section].apply(get_seq_gc_content)

        return seq_df.loc[:, ['gc_content_' + mrna_section]]

    def get_mrna_features(self, get_codon_frequency=True, include_non_aug_start=False, get_nucleotide_content=False,
                          compute_at_content=False):
        """

        :param get_codon_frequency: include the ratio of each codon in the coding sequence for each transcript
        :param include_non_aug_start: include transcripts which do not start with an AUG in the coding sequence
        :return: dataframe with transcript ids as indexes and transcript features as columns. Transcript
        features can be:
            - gene id
            - transcript id
            - chromosome name
            - existance of 5'UTR start codon AUG
            - transcript's cds start's with aug or not
            - 5'UTR length (natural and logged)
            - 3'UTR length (natural and logged)
            - CDS length (natural and logged)
            - Stop codon type
            - gc content 5'UTR
            - gc content 3'UTR
            - gc content CDS
            - AAA; AAC; ...; TTG; TTT codon fraction
        """

        utr_3_len_discount = 3  # account for stop codon anotated in 3 utr
        stop_codons = ['TAA', 'TAG', 'TGA']

        features_dict = {#'gene_id': [],
                         'transcript_id': [], 'chromosome': [], 'non_AUG': [], 'u_AUG': [],
                         '5_utr_length': [], '3_utr_length': [], 'cds_length': [], 'stop_codon': []}

        for transcript in self.data:
            #features_dict['gene_id'].append(transcript['gene_id'])
            features_dict['transcript_id'].append(transcript['transcript_id'])
            features_dict['chromosome'].append(transcript['chromosome'])
            features_dict['u_AUG'].append(transcript['u_aug'])
            features_dict['non_AUG'].append(transcript['non_AUG'])
            features_dict['5_utr_length'].append(len(transcript['utrs']['5_utr'])
                                                 if transcript['utrs']['5_utr'] is not None else np.nan)
            features_dict['3_utr_length'].append((max(len(transcript['utrs']['3_utr']) - utr_3_len_discount, 0))
                                                 if transcript['utrs']['3_utr'] is not None else np.nan)
            features_dict['cds_length'].append(len(transcript['cds']))
            features_dict['stop_codon'].append(np.intersect1d(transcript['utrs']['3_utr'][:3], stop_codons)
                                               if transcript['utrs']['3_utr'] is not None else np.nan)

        features_df = pd.DataFrame(data=features_dict)
        features_df = features_df.set_index('transcript_id')

        if not include_non_aug_start:
            features_df = features_df[features_df['non_AUG'] != 1]

        features_df['stop_codon'] = features_df['stop_codon'].fillna('')
        features_df['stop_codon'] = features_df['stop_codon'].apply(
            lambda y: 'stop_codon_other' if len(y) == 0 else 'stop_codon_' + y[0])
        # stop_codon_other : transcripts with no stop codon or whose anotation did not exist
        # get columns for each stop codon indicating the presence of that stop codon as a termination codon of the cds
        features_df = features_df.merge(pd.get_dummies(features_df['stop_codon']), left_index=True, right_index=True)
        features_df.drop(['stop_codon', 'stop_codon_other'], axis=1, inplace=True)

        # add one to log to account for 0 values
        features_df['log_3_utr_length'] = np.log2(1 + features_df['3_utr_length'])
        features_df['log_5_utr_length'] = np.log2(1 + features_df['5_utr_length'])
        features_df['log_cds_length'] = np.log2(1 + features_df['cds_length'])

        mrna_sections = ['5_utr', 'cds', '3_utr']
        nucleotides = ['A', 'C', 'G', 'T']

        for mrna_section in mrna_sections:

            gc_content_df = self.get_gc_content(mrna_section=mrna_section,
                                                include_non_aug_start=include_non_aug_start)
            if compute_at_content:
                gc_content_df['at_content_' + mrna_section] = 1 - gc_content_df['gc_content_'+mrna_section]
                gc_content_df.drop('gc_content_'+mrna_section, axis=1)

            gc_content_df['log_gc_content_' + mrna_section] = np.log2(gc_content_df['gc_content_'+mrna_section])
            features_df = features_df.merge(gc_content_df, right_index=True, left_index=True, how='left')

            if get_nucleotide_content:
                for nucleotide in nucleotides:
                    features_df = features_df.merge(self.get_n_content(mrna_section=mrna_section, n=nucleotide,
                                                                        include_non_aug_start=include_non_aug_start),
                                                    right_index=True, left_index=True, how='left')

        if get_codon_frequency:
            features_df = features_df.merge(self.get_codon_ratio(include_non_aug_start=include_non_aug_start),
                                            left_index=True, right_index=True)

        return features_df

    def get_computed_codon_ratio(self, codon_shift=None, transcripts=None):

        if transcripts is not None:
            sequences_df = self.get_sequence_df(mrna_section='cds')
            sequences_df = sequences_df.loc[sequences_df.index.intersection(transcripts)]
        else:
            sequences_df = self.get_sequence_df(mrna_section='cds')

        for i in range(len(sequences_df)):
            if i%1000==0:
                print(f'{i} of {len(sequences_df)} transcripts completed')
            count_df = codons_count(sequences_df.iloc[i][0], codon_shift).rename({'codon number': sequences_df.index[
                        i]}, axis=1).T
            if i == 0:
                counts_df = count_df
            else:
                counts_df = pd.concat([counts_df, count_df])

            computed_codon_ratios = counts_df.div(counts_df.sum(axis=1), axis=0)

        if codon_shift is None:
            assert computed_codon_ratios.equals(self.get_codon_ratio().loc[computed_codon_ratios.index])

        return computed_codon_ratios

    def get_codon_ratio(self, include_non_aug_start=False):

        """

        :param include_non_aug_start: include transcripts which do not start with an AUG in the coding sequence
        :return: dataframe with transcript ids as indexes and AAA; AAC; ...; TTG; TTT codon fraction on each column
        """
        # Still lacks some implemmentation to include non aug transcripts because of posible frameshift

        codons_array = np.array([transcript['codon_numbers']['codon number'].values for transcript in self.data
                                 if ((transcript['non_AUG'] == 0) & (not include_non_aug_start))])
        codons_ratio = codons_array / np.expand_dims(codons_array.sum(axis=1), axis=1)

        codons_names = self.data[0]['codon_numbers'].index.values
        transcript_id_list = self.all_transcript_ids if include_non_aug_start else self.aug_start_transcripts

        codons_df = pd.DataFrame(data=codons_ratio, columns=codons_names, index=transcript_id_list)

        return codons_df

    def split_train_test(self, data_frame, split_criterion='fraction', split_fraction=0.2,
                         split_chromosomes=None, n_split_chromosomes=4, random_seed=3):

        """

        :param data_frame: data frame you want to split. Must have transcript ids as indexes
        :param split_criterion: one of fraction or chromosomes.
        If fraction than the data frame will be splitted in
        train and test sets according to the specified fraction. Fraction must take a value between 0 and 1 and
        specifies the fraction of transcripts wanted for the test set.
        If chromsomes than the sets will be splitted by chromosome
        :param split_fraction:
        :param split_chromosomes: the specific chromosomes wanted for the test set. If none and split criterion is
        chromosomes than random chromosomes will be chosen
        :param n_split_chromosomes: number of chromsomes wanted on the test set if they are not specified
        :param random_seed: random seed number for reproducibility
        :return: 2 data frames with transcript ids as indexes - train dataframe and test dataframe
        """
        np.random.seed(random_seed)

        if split_criterion == 'fraction':
            test_df = data_frame.sample(frac=split_fraction, random_state=random_seed)
            train_indexes = data_frame.index.difference(test_df.index)
            train_df = data_frame.loc[train_indexes]

            return train_df, test_df

        if split_criterion == 'chromosome':

            transcript_id_chr = pd.DataFrame(
                {'transcript_id': self.all_transcript_ids,
                 'chromosome': [transcript['chromosome'] for transcript in self.data]}).set_index('transcript_id')

            if split_chromosomes is None:
                chromosomes = transcript_id_chr['chromosome'].unique()
                split_chromosomes = np.random.choice(chromosomes, n_split_chromosomes, replace=False)

            chromosome_df = transcript_id_chr.loc[data_frame.index]
            test_indexes = chromosome_df[chromosome_df['chromosome'].isin(split_chromosomes)].index
            test_df = data_frame.loc[test_indexes]
            train_indexes = data_frame.index.difference(test_df.index)
            train_df = data_frame.loc[train_indexes]

            print('Test set chromosomes: ', split_chromosomes)
            print('Test set split fraction: ', len(test_df) / len(data_frame))

            return train_df, test_df

        print('Split criterion is not valid. Choose between fraction or chromosome.')

    def split_cross_val(self, data_frame, split_criterion='all_transcripts', n_sets=10, n_split_chromosomes=2,
                        random_seed=2):

        """

        :param data_frame: data frame you want to split. Must have transcript ids as indexes
        :param split_criterion: one of fraction or chromosomes.
        If fraction than the data frame will be splitted in
        train and test sets according to the specified fraction. Fraction must take a value between 0 and 1 and
        specifies the fraction of transcripts wanted for the test set.
        If chromsomes than the sets will be splitted by chromosome
        :param n_sets: cross validation number of sets
        :param n_split_chromosomes:
        :param random_seed:
        :return:
        """
        np.random.seed(random_seed)

        if split_criterion == 'all_transcripts':
            data_frame = data_frame.sample(frac=1, random_state=random_seed)  # shuffle rows
            sets = np.array_split(data_frame, n_sets)

            return sets

        elif split_criterion == 'chromosome':

            transcript_id_vars = pd.DataFrame(
                {'transcript_id': self.all_transcript_ids,
                 'chromosome': [transcript['chromosome'][0] for transcript in self.data]}).set_index('transcript_id')

            chromosome_df = transcript_id_vars.loc[data_frame.index]
            chromosomes = np.array(transcript_id_vars['chromosome'].unique())
            sets = []

            while len(chromosomes) >= n_split_chromosomes:
                set_i_chromosomes = np.random.choice(chromosomes, n_split_chromosomes, replace=False)
                chromosomes = np.delete(chromosomes, np.argwhere(np.isin(chromosomes, set_i_chromosomes)))

                set_i_indexes = chromosome_df[chromosome_df['chromosome'].isin(set_i_chromosomes)].index
                sets.append(data_frame.loc[set_i_indexes])

            if len(chromosomes) > 0:  # final leftover set
                set_i_indexes = chromosome_df[chromosome_df['chromosome'].isin(chromosomes)].index
                sets.append(data_frame.loc[set_i_indexes])

            return sets

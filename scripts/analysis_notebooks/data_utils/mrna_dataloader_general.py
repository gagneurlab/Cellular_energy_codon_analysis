from pybedtools import Interval
import pyranges as pr
import pandas as pd
import numpy as np
from kipoiseq.extractors import FastaStringExtractor
from data_utils.sequence_utils import rc_dna, cut_seq, codons_count
import pickle



class GenomeCDSFetcher:

    def __init__(self, gtf, fae):

        self.fae = fae
        self.cds = (gtf
                    .query("gene_type == 'protein_coding'")
                    .query("Feature == 'CDS'")
                    .set_index('transcript_id'))
        
        self.transcripts = self.cds.index.unique()

    def __len__(self):
        return len(self.transcripts)
    
    def gtf_row2interval(self, row):
    #Note: Pyranges converts gtf coordinates to 0 based!
        return Interval(str(row.Chromosome),
                        int(row.Start),
                        int(row.End),
                        strand=str(row.Strand))

    def get_cds_exons(self, transcript_id):
        cds_exons = self.cds.loc[transcript_id]

        # get cds intervals
        if isinstance(cds_exons, pd.Series):
            # single exon
            strand = cds_exons.Strand
            intervals = [self.gtf_row2interval(cds_exons)]
        else:
            # multiple exons
            strand = cds_exons.iloc[0].Strand
            assert np.all(strand == cds_exons.Strand)

            intervals = [self.gtf_row2interval(row)
                         for i, row in cds_exons.loc[transcript_id].sort_values("Start").iterrows()]
        return intervals, strand

    # if the dna seq is not %3==0, there are unnecessary bases at the end
    # should be called only after all exons are connected!

    def get_seq(self, transcript_id):
        exons, strand = self.get_cds_exons(transcript_id)
        # merge the sequences

        seq = "".join([self.fae.extract(exon) for exon in exons])
        
        if strand == '-':
            # optionally reverse complement
            seq = rc_dna(seq)
        seq = cut_seq(seq)

        non_AUG = False
        if seq[:3] != 'ATG':
            non_AUG = True
        return seq, non_AUG

    def __getitem__(self, idx):
        return self.get_seq(self.transcripts[idx])


class MRNA_halflife_dl():

    def __init__(self, gtf_file_path: str, fasta_file_path: str, mrna_region: str = 'all',
                 get_codon_number: bool = False, u_aug_search: bool = False):

        # read gtf and fasta
        self.read_files(gtf_file_path, fasta_file_path)

        # Create CDS and UTR fetchers objects
        self.cds_fetcher = GenomeCDSFetcher(self.gtf, self.fae)
        self.utrs_fetcher = GenomeUTRFetcher(self.gtf, self.fae)

        # remove repetition of the same transcript, chromosome
        self.chromosomes_df = self.gtf.loc[:, ['Chromosome', 'transcript_id']].drop_duplicates() \
            .set_index("transcript_id")['Chromosome']
        self.gene_id_df = self.gtf.loc[:, ['gene_id', 'transcript_id']].drop_duplicates() \
            .set_index("transcript_id")['gene_id']

        self.mrna_region = mrna_region
        self.get_codon_number = get_codon_number
        self.u_aug_search = u_aug_search

        if self.mrna_region in ['cds', 'all'] :
            self.transcripts = self.cds_fetcher.transcripts
        elif self.mrna_region == 'utrs':
            self.transcripts = self.utrs_fetcher.transcripts


        all_codons_names = np.array(
            ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA',
             'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG',
             'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
             'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA',
             'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'])

        self.all_codons_df = pd.DataFrame(data=np.zeros(64), index=all_codons_names, columns=['codon number'])

    def read_files(self, gtf_file_path, fasta_file_path):
        """
        Initialization of Fasta Extractor and GTF file objects as well as reading of RNA halflife data
        :param gtf_file_path:
        :param fasta_file_path:

        """

        # read gtf and fasta data
        self.fae = FastaStringExtractor(fasta_file_path, use_strand=False, force_upper=True)
        self.gtf = pr.read_gtf(gtf_file_path, as_df=True)

    def __getitem__(self, idx):
        """

        :param idx: index of the dataloader corresponds to a different transcript
        :return: data on a particular transcript

        This features will always be present in an item of the dataloader:
        'chromosome': chromosome to which the transcript belongs to
        'transcript_id': code of the transcript

        Other features will be present depending on the region of the mRNA molecule we want to focus on:
        'all': returns sequence from the utrs and cds
        'cds': returns sequence from cds only
        'utrs': returns sequence form utrs only


        """

        transcript = self.transcripts[idx]

        data = {'chromosome': self.chromosomes_df.loc[transcript],
                'gene_id': self.gene_id_df.loc[transcript],
                'transcript_id': transcript
                }

        if self.mrna_region == 'all':
            data['cds'], data['non_AUG'] = self.cds_fetcher.get_seq(transcript)
            data['utrs'] = self.utrs_fetcher.get_seqs(transcript)

        elif self.mrna_region == 'utrs':
            data['utrs'] = self.utrs_fetcher.get_seqs(transcript)

        elif self.mrna_region == 'cds':
            data['cds'], data['non_AUG'] = self.cds_fetcher.get_seq(transcript)

        if self.get_codon_number:
            data['codon_numbers'] = self.get_codons(data['cds'])

        if self.u_aug_search:
            if data['utrs']['5_utr'] is not None:
                data['u_aug'] = int('ATG' in data['utrs']['5_utr'])
            else:
                data['u_aug'] = None

        return data

    def __len__(self):
        return len(self.transcripts)

    def get_codons(self, sequence):
        """
        :param sequence: string sequence of the transcript
        :return: data frame containing the number of each different codon on the sequence
        """
        return codons_count(sequence)

    def load_data(self, save_path=None):
        """
        Load all data into memory
        :return: list containing all the data
        """
        print('loading_data...')
        data_list = []
        for transcript_i, mrna in enumerate(self):
            data_list.append(mrna)
            if transcript_i % 1000 == 0:
                print('%f per cent completed' % (transcript_i / len(self) * 100))

        if save_path is not None:
            pickle.dump(data_list, open(save_path, "wb"))
        else:
            return data_list

class GenomeUTRFetcher:

    def __init__(self, gtf, fae):

        self.fae = fae

        self.transcript_data = (gtf
                    .query("gene_type == 'protein_coding'")
                    .query("(Feature == 'UTR') | (Feature == 'CDS')")
                    .set_index('transcript_id'))

        transcripts_utr = (gtf.query("gene_type == 'protein_coding'")
                               .query("Feature == 'UTR'")
                               .set_index('transcript_id'))

        self.transcripts = self.transcript_data.index.unique()
        # have all transcripts available have both anotated utrs and cds
        self.transcripts = self.transcripts.intersection(transcripts_utr.index).unique()

    def __len__(self):
        return len(self.transcripts)

    def get_utrs_intervals(self, transcript_id):

        utrs_plus_cds_df = self.transcript_data.loc[transcript_id]
        utrs_chunks = utrs_plus_cds_df[utrs_plus_cds_df.Feature == 'UTR'].sort_values("Start")
        cds_chunks = utrs_plus_cds_df[utrs_plus_cds_df.Feature == 'CDS'].sort_values("Start")

        if isinstance(cds_chunks, pd.Series):
            cds_min_pos = cds_chunks.Start  # gtf coordinates
            cds_max_pos = cds_chunks.End
            strand = cds_chunks.Strand
            chromosome = cds_chunks.Chromosome
        else:
            cds_min_pos = cds_chunks.iloc[0].Start  # gtf coordinates
            cds_max_pos = cds_chunks.iloc[-1].End  # gtf coordinates
            strand = cds_chunks.iloc[0].Strand
            chromosome = str(cds_chunks.iloc[0].Chromosome)

        utr_5_intervals = None
        utr_3_intervals = None

        for transcript_i, utr_chunk in utrs_chunks.iterrows():

            current_utr_chunk_start = utr_chunk.Start  # gtf coordinates
            current_utr_chunk_end = utr_chunk.End  # gtf coordinates
            # 0 based coordinates (interval from utr_start-1 until end-1 included)
            current_interval = Interval(chromosome, current_utr_chunk_start - 1,
                                        current_utr_chunk_end)

            if current_utr_chunk_start < cds_min_pos:
                # 5'UTR
                if utr_5_intervals is None:
                    utr_5_intervals = [current_interval]
                else:
                    utr_5_intervals.append(current_interval)

            elif current_utr_chunk_start > cds_max_pos:
                # 3'UTR
                if utr_3_intervals is None:
                    utr_3_intervals = [current_interval]
                    continue
                else:
                    utr_3_intervals.append(current_interval)

        if strand == '-':
            utr_5_intervals, utr_3_intervals = utr_3_intervals, utr_5_intervals  # swap variables

        return utr_5_intervals, utr_3_intervals, strand

    def get_seqs(self, transcript_id):

        if transcript_id in self.transcripts:
            utr_5_intervals, utr_3_intervals, strand = self.get_utrs_intervals(transcript_id)
        else:
            utr_5_intervals, utr_3_intervals = None, None

        # merge the sequences from 5 utr and 3 utr
        if utr_5_intervals is not None:
            seq_5_utr = "".join([self.fae.extract(interval) for interval in utr_5_intervals])
            if strand == '-':
                # reverse code if strand is negative
                seq_5_utr = rc_dna(seq_5_utr)
        else:
            seq_5_utr = None

        if utr_3_intervals is not None:
            seq_3_utr = "".join([self.fae.extract(interval) for interval in utr_3_intervals])
            if strand == '-':
                # reverse code if strand is negative
                seq_3_utr = rc_dna(seq_3_utr)
        else:
            seq_3_utr = None

        return {'5_utr': seq_5_utr,
                '3_utr': seq_3_utr}

    def __getitem__(self, idx):
        return self.get_seqs(self.transcripts[idx])


import pandas as pd
import numpy as np
from typing import List
import os
import scvelo as scv
import scanpy as sc
import anndata
import scipy
from itertools import compress, chain


class RNAHalflifeDataloader:
    def __init__(self, features_file_path: str):
        self.path = features_file_path
        self.features = pd.read_csv(features_file_path)
        self.features.loc[:, ["cds", "utr3", "utr5"]] = self.features.loc[
            :, ["cds", "utr3", "utr5"]
        ].replace(np.nan, "")

        self.features.loc[
            :, self.features.columns[self.features.columns.str.contains("gc_content")]
        ] = self.features.loc[
            :, self.features.columns[self.features.columns.str.contains("gc_content")]
        ].replace(
            np.nan, 0
        )

        self.features.set_index("transcript_id", inplace=True)

        # compute log of length
        length_cols = self.features.columns[
            self.features.columns.str.contains("length")
        ]

        for col in length_cols:
            self.features[f"log_{col}"] = np.log10(self.features[col])
            self.features[f"log_{col}"] = self.features[f"log_{col}"].replace(
                [np.nan, np.inf, -np.inf], 0
            )

    def keys(self) -> List[object]:
        return self.features.index

    def feature_names(self) -> List[object]:
        return list(self.features.columns)

    def sel(self, identifier):
        """
        extract sequence by identifier
        """
        return self.features.loc[identifier].values

    def isel(self, idx):
        """
        extract sequence by  index
        :param idx: index in the range (0, len(self))
        """
        return self.features.iloc[idx].values

    def items(self):
        for i in self.keys():
            yield i, self.sel(i)

    def __len__(self):
        return len(self.keys())

    def __getitem__(self, item):
        self.isel(item)


def convert_gene_names(data, gtf_file, copy=False):
    adata = data.copy() if copy else data
    
    import pyranges as pr
    g = pr.read_gtf(gtf_file, as_df=True, duplicate_attr=True)
    g['gene_id'] = [x.split('.')[0] for x in g.gene_id]
    g = g.loc[(g.Feature == 'gene') & (g.gene_type =='protein_coding'), ['gene_name', 'gene_id', 'gene_type']]
    
    names_to_id_dict = dict(zip(g["gene_name"].values.tolist(), g["gene_id"].values.tolist()))
    adata.var["id"] = adata.var.index.map(names_to_id_dict)

    return adata if copy else None


def strip_version_number(id):
    return id.split(".")[0]


class BaseDataLoader:
    adata: anndata.AnnData
    isoforms: pd.DataFrame
    mrnaDataloader: RNAHalflifeDataloader
    genes: dict

    def __init__(
        self,
        adata,
        isoforms_table,
        mrnaDataloader: RNAHalflifeDataloader,
        cell_type_key=None,
    ):
        """
        This base class defines dataloader attributes shared across all half life models
        """

        self.adata = adata

        self.isoforms = pd.read_csv(isoforms_table, index_col="gene_id")
        self.isoforms.index = self.isoforms.index.map(strip_version_number)
        self.filter_isoforms()

        self.mrnaDataloader = mrnaDataloader
        self.cell_type_key = cell_type_key

    def genes_to_transcript(self, genes):
        return self.isoforms.loc[genes].iloc[:, 0]

    def transcript_to_genes(self, transcripts):
        return self.isoforms.loc[self.isoforms.iloc[:, 0].isin(transcripts)].index

    def filter_isoforms(self):
        self.isoforms = self.isoforms.loc[self.gene_id()]
        return

    def gene_id(self):
        return self.adata.var["id"].values

    def names_to_id_dict(self):
        return dict(
            zip(
                self.adata.var.index.values.tolist(),
                self.adata.var["id"].values.tolist(),
            )
        )

    def id_to_names_dict(self):
        return dict(
            zip(
                self.adata.var["id"].values.tolist(),
                self.adata.var.index.values.tolist(),
            )
        )

    def cell_types(self):
        return list(self.adata.obs[self.cell_type_key].unique())

    def filter_and_normalize(self, min_shared_cells):
        scv.pp.filter_genes(
            self.adata,  min_shared_cells = min_shared_cells
        ) #, min_shared_counts=filter_genes min_cells = 38, min_cells_u=38
        if self.adata.shape[0] > 1:
            sc.pp.highly_variable_genes(
                self.adata, n_top_genes=4000, subset=False, flavor="seurat_v3"
            )
        scv.pp.normalize_per_cell(self.adata)

    def input_data(self, features=None):
        self.filter_isoforms()
        transcripts = self.isoforms.iloc[:, 0]
        if not features:
            return self.mrnaDataloader.features.loc[transcripts, :]
        return self.mrnaDataloader.features.loc[transcripts, features]

    def __len__(self):
        return len(self.gene_id())


class DataLoaderPseudoBulk(BaseDataLoader):
    def __init__(
        self,
        adata,
        isoforms_table,
        mrnaDataloader: RNAHalflifeDataloader,
        batchsize,
        celltype_key,
        filter_cells=None
    ):

        super().__init__(adata, isoforms_table, mrnaDataloader, celltype_key)

        self.batchsize = batchsize
        self.filter_cells = filter_cells
        self.subclusters = None

        # Compute pseudo bulk
        (
            self.adata,
            self.subclusters,
            self.celltypes,
            self.removed,
        ) = compute_pseudo_bulk(
            self.adata,
            cell_type_key=self.cell_type_key,
            batch=self.batchsize,
            filter_cells=self.filter_cells,
        )

        # filter genes, normalize and compute highly variable genes
        # add raw layers
        self.adata.layers["counts"] = self.adata.X.copy()
        for layer in ["unspliced", "spliced"]:
            self.adata.layers[f"{layer}_raw"] = self.adata.layers[layer].copy()
    
    def compute_half_life(
        self, 
        min_shared_cells,
        filter_highly_variable=False, 
        exclude_zeros=False,
        mask=False,
        log=True
    ):
        
        self.filter_and_normalize(
            min_shared_cells=min_shared_cells
        )

        # Filter to genes that are present in isoform table
        self.adata = self.adata[:, self.adata.var.id.isin(self.isoforms.index)]

        if filter_highly_variable:
            self.adata = self.adata[:, self.adata.var.highly_variable]

        self.filter_isoforms()

        # Compute half life
        compute_half_life(self.adata, mask)
        if log:
            sc.pp.log1p(self.adata, copy=False)
            sc.pp.log1p(self.adata, layer="unspliced", copy=False)
            sc.pp.log1p(self.adata, layer="spliced", copy=False)
            self.adata.layers["half_life"] = np.log10(self.adata.layers["half_life"])

        compute_half_life_diff(self.adata, exclude_zeros)

    def sample_label(self):
        return np.repeat(self.celltypes, [len(x) for x in self.subclusters])


def compute_pseudo_bulk(adata, cell_type_key, batch=None, filter_cells=None):
    print(filter_cells)
    celltypes = adata.obs[cell_type_key].unique()
    # TODO: Check when we should apply filter/normalization
    clusters = [
        np.where(adata.obs[cell_type_key].isin([celltype]))[0] for celltype in celltypes
    ]
    removed = None

    if filter_cells:
        suitable_clusters = np.array(
            [len(cluster) >= filter_cells for cluster in clusters]
        )
        clusters = list(compress(clusters, suitable_clusters))
        removed = celltypes[~suitable_clusters]
        celltypes = celltypes[suitable_clusters]

        print(
            f"Removed {len(removed)} celltypes as they had less cells than filter_cells {filter_cells}"
        )
        
    subclusters = [[cluster] for cluster in clusters]


    adata = aggregate_adata(
        adata, subclusters, celltypes, cell_id=cell_type_key, agg_type="sum"
    )

    label = np.repeat(celltypes, [len(x) for x in subclusters])
    ext = np.concatenate([np.arange(len(x)) for x in subclusters])
    adata.obs["bulk_label"] = [
        label.astype(str)[i] + "_" + ext.astype(str)[i] for i in range(len(label))
    ]
    return adata, subclusters, celltypes, removed


def aggregate_adata(
    adata, subclusters, celltypes, cell_id="cell_ontology_class", agg_type="sum"
):
    subclusters_flat = list(chain(*subclusters))
    print(f"Computed {len(subclusters_flat)} subcluster")

    assert agg_type in ["sum", "mean"], "agg_type can either be 'sum' or 'mean'"

    agg_func = np.sum if agg_type == "sum" else np.mean

    layer_keys = ["spliced", "unspliced"]  # list(adata.layers.keys())

    print(f"Aggregating X and layers {layer_keys}")
    X_agg = scipy.sparse.csr.csr_matrix(
        aggregate_cluster(adata.X, subclusters_flat, agg_func)
    )
    layers_agg = [
        scipy.sparse.csr.csr_matrix(
            aggregate_cluster(adata.layers[layer], subclusters_flat, agg_func)
        )
        for layer in layer_keys
    ]
    #
    labels = np.repeat(celltypes, [len(x) for x in subclusters])
    return anndata.AnnData(
        X=X_agg,
        layers=dict(zip(layer_keys, layers_agg)),
        var=adata.var,
        obs=pd.DataFrame({cell_id: labels}),
    )


def create_mask(idx, length):
    mask = np.repeat(False, length)
    if np.isfinite(idx).all():
        mask[idx] = True
    else:
        mask = mask
    return mask


def aggregate_cluster(x, subclusters, agg_func):
    x = x.A if isinstance(x, scipy.sparse.csr.csr_matrix) else x
    cluster_agg = np.asarray(
        [
            agg_func(x[create_mask(subcluster, x.shape[0])], axis=0)
            for subcluster in subclusters
        ]
    )
    return cluster_agg


def compute_half_life(adata, copy=False):

    adata = adata.copy() if copy else adata

    adata.layers["half_life"] = (adata.layers["spliced"].A + 1) / (
        adata.layers["unspliced"].A + 1
    )
    assert (adata.layers["half_life"] != 0).all(), "Half life contains zeros"

    adata.layers["half_life"] = adata.layers["half_life"]
    return adata if copy else None


def compute_half_life_diff(adata, exclude_zeros, copy=False):
    adata = adata.copy() if copy else adata

    if scipy.sparse.issparse(adata.layers["half_life"]):
        hl = adata.layers["half_life"].A
    else:
        hl = adata.layers["half_life"]

    if exclude_zeros:
        hl[hl == 0] = np.nan
        hl_means = np.nanmean(hl, axis=0)
    else:
        hl_means = np.nanmean(hl, axis=0)

    adata.varm["mean_half_life"] = hl_means
    adata.layers["half_life_diff"] = scipy.sparse.csr.csr_matrix(
        hl - hl_means[np.newaxis, :]
    )
    return adata if copy else None


def load_mouse_features(mrna_dataloader, dataloader, adata):

    features = mrna_dataloader.features
    isoform_table = dataloader.isoforms

    if adata is not None and isoform_table is not None:
        print("Filtering features table")
        transcripts = isoform_table.iloc[isoform_table.index.isin(adata.var["id"]), 0]
        features = features.loc[transcripts, :]

        def reverse_dict(dict):
            return {v: k for k, v in dict.items()}

        features["gene_id"] = features.index.map(reverse_dict(transcripts.to_dict()))
        features["gene_name"] = features["gene_id"].map(
            reverse_dict(adata.var["id"].to_dict())
        )
        features["chromosome"] = features["gene_name"].map(adata.var["Chromosome"])

        features = features.loc[
            :,
            ~(
                features.columns.str.startswith("utr3_")
                | features.columns.str.startswith("utr5_")
            ),
        ]
        return features

    
    
    
def calculate_support(transcripts):
    support =    transcripts.tag.str.contains('basic')*1 + \
                    transcripts.tag.str.contains('CCDS')*1 + \
                    transcripts.tag.str.contains('appris_principal_1')*1 + \
                    (transcripts.transcript_support_level == '1')*1 
    return support
        
        

def find_major_isoform(transcripts):
    return transcripts.sort_values(['gene_id','support'],ascending=False).groupby('gene_id').first()

def expand(major_isoform, celltypes):
    gene_id = major_isoform.index.values
    transcript_id = major_isoform.transcript_id.values
    mat = np.tile(transcript_id, (len(celltypes),1))
    return pd.DataFrame(mat.T, index=gene_id, columns=celltypes.astype('str'))
    
def create_isoform_table(adata, gtf, cds):
    celltypes = adata.obs.cell_ontology_class.unique()
    transcripts = (gtf
       .query("Feature == 'transcript'")
       .query("gene_type =='protein_coding'"))
    
    transcripts = transcripts[transcripts.transcript_id.isin(cds)]
    transcripts['support'] = calculate_support(transcripts)
    
    major_isoform = find_major_isoform(transcripts)
    # expand to celltypes
    major_isoform = expand(major_isoform, celltypes)
    return major_isoform.reset_index().rename(columns = {'index': 'gene_id'}), transcripts
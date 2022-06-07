import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd


def load_tmp(f_name, in_dir, out_dir):
    
    # load sparse matrix:
    X = io.mmread(os.path.join(in_dir, '.'.join(['counts', 'mtx'])))

    # create anndata object
    adata = anndata.AnnData(
        X=X.transpose().tocsr(),
        dtype=X.dtype
    )

    # load cell metadata:
    cell_meta = pd.read_csv(os.path.join(in_dir, '.'.join(['metadata', 'csv'])))

    # load gene names:
    with open(os.path.join(in_dir, '.'.join(['gene_names', 'csv'])), 'r') as f:
        gene_names = f.read().splitlines()

    # set anndata observations and index obs by barcodes, var by gene names
    adata.obs = cell_meta
    adata.obs.index = adata.obs['cell_id'].to_list()
    adata.var.index = gene_names

    # load dimensional reduction, set pca
    drs = glob.glob(os.path.join(in_dir,'dr_*'))

    for dr in drs:
        dr_coord = pd.read_csv(os.path.join(in_dir, '.'.join([dr, 'csv'])))
        dr_coord.index = adata.obs.index
        adata.obsm['_'.join(["X",dr])] = dr_coord.to_numpy()

    # set umap
    adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
    adata.obsm['X_umapint'] = np.vstack((adata.obs['UMAPint_1'].to_numpy(), adata.obs['UMAPint_2'].to_numpy())).T

    # save dataset as anndata format
    adata.write(os.path.join(out_dir,'.'.join([f_name,'h5ad'])))
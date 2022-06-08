#!/usr/bin/env python3

"""
Convert Seurat .h5Seurat to Anndata .h5ad.
"""
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import pandas as pd
import shutil
from pathlib import Path
import re

def seurat_to_anndata(f_name, in_dir, out_dir, purge_tmp=False):
    # # coerce type
    # in_dir = Path(in_dir)
    # out_dir = Path(out_dir)

    # load sparse matrix:
    X = io.mmread(in_dir/'counts.mtx')

    # create anndata object
    adata = anndata.AnnData(
        X=X.transpose().tocsr(),
        dtype=X.dtype
    )

    # load cell metadata:
    cell_meta = pd.read_csv(in_dir/'metadata.csv')
    # load gene names:
    with open(in_dir/'gene_names.csv', 'r') as f:
        gene_names = f.read().splitlines()

    # set anndata observations and index obs by barcodes, var by gene names
    adata.obs = cell_meta
    adata.obs.index = adata.obs['cell_id'].to_list()
    adata.var.index = gene_names

    # load dimensional reduction, set pca
    drs = in_dir.glob("dr_*")
    for dr in drs:
        dr_name = re.search('dr_(.+)', dr.stem).group(1)
        dr_coord = pd.read_csv(dr)
        dr_coord.index = adata.obs.index
        adata.obsm['_'.join(["X",dr_name])] = dr_coord.to_numpy()

    # set umap
    embs = in_dir.glob("embeddings_*")
    for emb in embs:
        emb_name = re.search('embeddings_(.+)', emb.stem).group(1)
        emb_coord = pd.read_csv(emb)
        emb_coord.index = adata.obs.index
        adata.obsm['_'.join(["X",emb_name])] = emb_coord.to_numpy()

    # save dataset as anndata format
    adata.write(out_dir/'.'.join([f_name,'h5ad']))
    
    # purge tmp files
    if purge_tmp:
        shutil.rmtree(in_dir)

if __name__ == "__main__":
    """
    Parse CML arguments
    """
    import sys
    if len(sys.argv) < 3:
        print("Insufficient argument passed. f_name=filename, in_dir=tmp dir, out_dir=output dir")
    else:
        try:
            f_name = str(sys.argv[1])
            in_dir = Path(sys.argv[2])
            out_dir = Path(sys.argv[3])
            # purge_tmp = sys.argv[4]  ## TODO: parse boolean arguments
            seurat_to_anndata(f_name, in_dir, out_dir)
        except:
            print("Error: wrong argument type.")
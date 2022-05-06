import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

f_name = "hpcs_lps_state_marked"
in_dir = "data/tmp/"
out_dir = "data/"

# load sparse matrix:
X = io.mmread(os.path.join(in_dir, '.'.join([f_name+'_counts', 'mtx'])))

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv(os.path.join(in_dir, '.'.join([f_name+'_metadata', 'csv'])))

# load gene names:
with open(os.path.join(in_dir, '.'.join([f_name+'_gene_names', 'csv'])), 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['cell_id'].to_list()
adata.var.index = gene_names

# load dimensional reduction, set pca:1
for pca in ["pca"]:
    pc_coord = pd.read_csv(os.path.join(in_dir, '.'.join([f_name + '_' + pca, 'csv'])))
    pc_coord.index = adata.obs.index
    adata.obsm['_'.join(["X",pca])] = pc_coord.to_numpy()

# set umap
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
adata.obsm['X_umapint'] = np.vstack((adata.obs['UMAPint_1'].to_numpy(), adata.obs['UMAPint_2'].to_numpy())).T

# # plot a UMAP colored by sampleID to test:
# sc.pl.umap(adata, color=['stim'], frameon=False, save=False)
# sc.pl.scatter(adata, basis='umapint', color = ['cell_type_brief'], frameon=False, save=False)

# save dataset as anndata format
adata.write(os.path.join(out_dir,'.'.join([f_name,'h5ad'])))

# # reload dataset
# adata = sc.read_h5ad('my_data.h5ad')

import celltypist
from celltypist import models
import sys
import os
import loompy
import IPython
import anndata
import scvi
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import seaborn as sns
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix

IPython.InteractiveShell.ast_node_interactivity = "all"

from IPython.display import set_matplotlib_formats
%config InlineBackend.print_figure_kwargs={'facecolor' : "w"}
%config InlineBackend.figure_format='retina'

subset_cells$barcode <- colnames(subset_cells)
subset_cells$UMAP_1 <- subset_cells@reductions$umap@cell.embeddings[,1]
subset_cells$UMAP_2 <- subset_cells@reductions$umap@cell.embeddings[,2]
write.csv(subset_cells@meta.data, file='All_cell_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(subset_cells, assay='RNA', slot='counts')
writeMM(counts_matrix, file = 'All_cell_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(subset_cells@reductions$pca@cell.embeddings, file = 'All_cell_pca.csv' , quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='All_cell_gene_names.csv',
  quote=F,row.names=F,col.names=F
)

###--------------------------------construct the anndata structure----------------------------###
# load sparse matrix:
X = io.mmread("All_cell_counts.mtx.gz")

# create anndata object
adata = anndata.AnnData(
  X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("All_cell_metadata.csv")
cell_meta["barcode"] = [x.replace("-1", ".1") for x in cell_meta["barcode"].to_list()]

# load gene names:
with open("All_cell_gene_names.csv", 'r') as f:
  gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("All_cell_pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

sc.set_figure_params(figsize=(8, 8))

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

model = models.Model.load(model = 'Immune_All_AddPIP.pkl')
predictions = celltypist.annotate(adata, model = 'Immune_All_AddPIP.pkl', majority_voting = True)

adata = predictions.to_adata()

adata.obs.to_csv("Cell_annotation.celltypist.txt")

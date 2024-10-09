import scvelo
import scanpy
import cellrank
import numpy
import pandas
import anndata
import os
import sys
import pandas as pd
import loompy
import scanpy as sc
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import scvelo as scv
import IPython
IPython.core.interactiveshell.InteractiveShell.ast_node_interactivity = "all"


###-----------------neutrophil meta.data from Seurat--------------------###
# subset_cells = subset(subset_cells, downsample = 3000)
subset_cells = subset(subset_cells, subset = Major_cluster %in% c("Neutrophil", "GMP_C28"))
subset_cells = subset(subset_cells, downsample = 3000)

subset_cells@meta.data$barcode <- colnames(subset_cells)
subset_cells@meta.data$UMAP_1 <- subset_cells@reductions$umap@cell.embeddings[,1]
subset_cells@meta.data$UMAP_2 <- subset_cells@reductions$umap@cell.embeddings[,2]
write.csv(subset_cells@meta.data, file='Neutrophil_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)

cellcycle_genes <- c(cc.genes$g2m.genes, cc.genes$g2m.genes)
selected_genes = setdiff(subset_cells %>% row.names(), cellcycle_genes)

counts_matrix <- GetAssayData(subset_cells[selected_genes,], assay = 'RNA', slot = 'counts')
writeMM(counts_matrix, file = paste0('Neutrophil_counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(subset_cells@reductions$pca@cell.embeddings, file = 'Neutrophil_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
	data.frame('gene' = rownames(counts_matrix)),file = 'Neutrophil_gene_names.csv',
	quote = F, row.names = F, col.names = F
)

###----------------------------extract neutrophis from loom file------------------------------###
adata = scvelo.read("/public/home/HeShuai/Projects/1.HCA2.0/1.Rawdata/Velocyto_looms/merged.loom", cache = True)
adata.var_names_make_unique()

Neutrophil_cell_metadata = pd.read_csv("Neutrophil_metadata.csv", sep = ",")
Neutrophil_cell_barcodes = [x.replace("-1", "") for x in Neutrophil_cell_metadata["barcode"].to_list()]

adata = scvelo.read("Neutrophil_velocity.loom")

cur_celltypes = Neutrophil_cell_barcodes

adata_subset_Neutrophil = adata[adata.obs['obs_names'].isin(cur_celltypes)]
adata_subset_Neutrophil.obs.index = adata_subset_Neutrophil.obs['obs_names'].to_list()

adata_subset_Neutrophil.write_loom("Neutrophil_velocity.loom")

ldata = scvelo.read('Neutrophil_velocity.loom', cache = True)

# load sparse matrix:
X = io.mmread("Neutrophil_counts.mtx")

# create anndata object
adata = anndata.AnnData(
	X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("Neutrophil_metadata.csv")
cell_meta["barcode"] = [x.replace("-1", "") for x in cell_meta["barcode"].to_list()]

# load gene names:
with open("Neutrophil_gene_names.csv", 'r') as f:
	gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("Neutrophil_pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color = 'seurat_clusters', frameon=False, save=True, legend_loc='on data')

# save dataset as anndata format
adata.write('Neutrophil_data.h5ad')

###----------------------------------reload dataset------------------------###
adata = sc.read_h5ad('Neutrophil_data.h5ad')
adata.obs.index = adata.obs["barcode"].to_list()

scvelo.settings.verbosity = 3
scvelo.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cellrank.settings.verbosity = 2

ldata = scvelo.read('Neutrophil_velocity.loom', cache=False)

ldata.obs.index = ldata.obs['obs_names'].to_list()
# merge matrices into the original adata object
adata = scvelo.utils.merge(adata, ldata)
# plot umap to check
sc.pl.umap(adata, color=['seurat_clusters'], frameon=False, legend_loc='on data', title='', save='_celltypes.pdf')

scv.pl.proportions(adata, groupby='Detail_Annotation_final')

# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
# scv.tl.umap(adata)

# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

# scv.tl.recover_dynamics(adata)
# scv.tl.velocity(adata, mode='dynamical')
# scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='Detail_Annotation_final', save='embedding_grid.pdf', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['Detail_Annotation_final', 'Stage'], save='embedding_stream.pdf', title='')

scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save='embedding_pseudotime.pdf')

df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
	pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

scv.get_df(adata, 'fit*', dropna=True).head()
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save='embedding_latent_time.pdf')


scv.tl.paga(adata, groups = 'Detail_Annotation_final')
# df = scv.get_df(adata, 'paga/transitions_confidence', precision = 2).T
# df.style.background_gradient(cmap = 'Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
						min_edge_width=2, node_size_scale=1.5)





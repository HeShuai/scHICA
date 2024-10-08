###---------------milopy----------------------------###
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo  ## For mouse gastrulation data 
import anndata
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
import milopy.core as milo
import milopy.plot as milopl
import bbknn
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# reload dataset
## scaled data must have the same dim with data slot or high variable feature ##
SaveH5Seurat(orig, filename = "All_HCA2.0_cells_annotation.1038993.h5Seurat", overwrite = T)
Convert("All_HCA2.0_cells_annotation.1038993.h5Seurat", dest = "h5ad", overwrite = T)

adata = sc.read_h5ad('All_HCA2.0_cells_annotation.1038993.h5ad')
# adata.obs.index = adata.obs["barcode"].to_list()

bbknn.bbknn(adata, batch_key = "Donor", n_pcs = 50)
# bbknn$bbknn(adata, batch_key = 0, n_pcs = as.integer(50), neighbors_within_batch = as.integer(1))
# sc.tl.umap(adata)
# sc.tl.leiden(adata, resolution = 1)

# sc.pl.umap(adata, color = 'leiden', save = "Myeloid_all_BBKNN.pdf")

###----------------------Construct neighbourhoods------------------###
milo.make_nhoods(adata, prop=0.1)
adata[adata.obs['nhood_ixs_refined'] != 0].obs[['nhood_ixs_refined', 'nhood_kth_distance']]

nhood_size = np.array(adata.obsm["nhoods"].sum(0)).ravel()
plt.hist(nhood_size, bins=100)

milo.count_nhoods(adata, sample_col="Donor_Name")
adata.uns["nhood_adata"]

adata.obs["Stage_continuous"] = adata.obs["Stage"].astype('category').cat.codes ## convert the object to category.
milo.DA_nhoods(adata, design="~Stage_continuous")

old_figsize = plt.rcParams["figure.figsize"]
plt.rcParams["figure.figsize"] = [10,5]
plt.subplot(1,2,1)
plt.hist(adata.uns["nhood_adata"].obs.PValue, bins=50);
plt.xlabel("P-Vals");
plt.subplot(1,2,2)
plt.plot(adata.uns["nhood_adata"].obs.logFC, -np.log10(adata.uns["nhood_adata"].obs.SpatialFDR), '.');
plt.xlabel("log-Fold Change");
plt.ylabel("- log10(Spatial FDR)");
plt.tight_layout()
plt.rcParams["figure.figsize"] = old_figsize

import milopy.utils
milopy.utils.build_nhood_graph(adata)

plt.rcParams["figure.figsize"] = [10,10]
milopl.plot_nhood_graph(adata, 
                        alpha=0.01, ## SpatialFDR level (1%) 
                        min_size=2, ## Size of smallest dot
                       save = "milopl.plot_nhood_graph_All_cells.pdf")
          
adata.obs["Global_annotation"] = adata.obs["Global_annotation"].astype(str) ## conver number to str
milopy.utils.annotate_nhoods(adata, anno_col='Global_annotation')

plt.hist(adata.uns['nhood_adata'].obs["nhood_annotation_frac"]);
plt.xlabel("Global_annotation")

adata.uns['nhood_adata'].obs.loc[adata.uns['nhood_adata'].obs["nhood_annotation_frac"] < 0.4, "nhood_annotation"] = "Mixed"

sc.pl.violin(adata.uns['nhood_adata'], "logFC", groupby="nhood_annotation", rotation=90, show=False);
plt.axhline(y=0, color='black', linestyle='--');
plt.show()

adata.uns['nhood_adata'].obs.to_csv("All_cells_milopy_cells.txt", sep = "\t")



###------------------------------plot using R----------------------------------------###
Da_results <-  read.table("All_cells_milopy_cells.txt", sep = "\t", header = T, row.names = 1)
Da_results$Nhood<- Da_results%>%row.names
png("HCA2.0_Myeloid_milo_test_dotplot_all_cells.png", height = 15, width = 25, units = "in", res = 400)
plotDAbeeswarm(Da_results, group.by = c("nhood_annotation"))
dev.off()


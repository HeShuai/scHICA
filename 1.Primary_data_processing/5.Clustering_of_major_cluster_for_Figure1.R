packages <- list("Seurat", "Matrix", "stringr", "stringi", "gridExtra", "ggplot2",
								 "plyr", "ggthemes", "cowplot", "data.table", "RColorBrewer", "ComplexHeatmap",
								 "circlize", "pheatmap", "reshape2", "scales", "rlang", "future",
								 "parallel", "ggrepel", "ggsci", "ape", "dplyr", "harmony",
								 "sigminer", "ggpubr", "rstatix", "ROGUE", "tidyverse", "reticulate", "dendextend",
								 "scRepertoire", "ggdendro", "ggtern", "viridis", "scatterpie", "SeuratData", "SeuratDisk",
								 "scibet", "scuttle", "zellkonverter", "cowplot", "gmodels") ## "descr": CrossTable

lapply(packages, library, character.only = TRUE)
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]
options(width = 240)
options(future.globals.maxSize = 400*1000 * 1024^2)
plan("multisession", workers = 1)
plan()
###----------------------------------This is an example for 109 samples (321 samples in our study)----------------------------------------------####
################################################################################################
### Merging Seurat datasets containing over approximately 1.03 million cells and 20,000+ genes can be challenging.
### As an alternative approach, consider either performing the merge using Python-based tools (e.g., Scanpy)
### or limiting the analysis in Seurat to highly variable genes instead of the full gene set.
################################################################################################

tissues_colors <- c("Ascending colon" = "#E64B35FF",
										"Sigmoid colon" = '#ff0000ff',
										"Transverse colon" = '#d62790',
										Ileum = "#91D1C2FF",
										Jejunum = '#4DBBD5FF',
										Duodenum = '#00A087FF',
										"Lesser gastric curvature" = '#F0E685FF',
										"Mesenteric lymph node" = "#FFB570FF",
										Blood = "#1B1919FF",
										PBMC = "#374E55FF",
										"Whole blood" = "#3F4041FF",
										Marrow = '#466983FF',
										Spleen = "#725663FF",
										Lung = "#6A6599FF",
										Liver = '#749B58FF',
										Kidney = '#B09C85FF',
										Skin = '#D6D6CEFF',
										Bladder = "#D49464FF",
										Esophagus = '#c5b0d5',
										Trachea = "#BA63affF",
										Heart = '#9467bd',
										Muscle = "#aec7e8",
										Pancreas = "#5050FFFF",
										Stomach = '#ADB17DFF',
										Thymus = '#8c564b'
)

####--------------------------------Step 1. enviroment seting---------------------------###
###-----------------------seting python environment--------------------------###
use_python("/public/home/HeShuai/miniconda3/bin/python")
###-----------------------loading python package-----------------##---###
anndata <- import("anndata",convert = FALSE)
bbknn <- import("bbknn", convert = FALSE)
sc <- import("scanpy",convert = FALSE)
Projectname <- "HCA2.0"
####--------------------------------Step 2. merge all the 109 samples---------------------------###
lst <- list()
i <- 0
for (singlet_data in dir("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/orig/") %>% grep(pattern = "RData", value = T)) {
	print(i)
	tmp <- get(load(paste0("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/orig/", singlet_data)))
	tmp$orig.ident <- singlet_data %>% gsub(pattern = "\\.RData", replacement = "")
	lst[[singlet_data %>% gsub(pattern = "\\.RData", replacement = "")]] <- tmp
	rm(tmp)
	print(singlet_data)
	i <- i + 1
}

subset_cells <- merge(lst[[1]], c(lst[2:109]), add.cell.ids = names(lst))
################------------------------------------------round one: batch effect correction and global clustering-----------------------------------##############################
####--------------------------------Step 3. quality control, normalization and scale data---------------------------###
mito.features <- grep(pattern = "^MT-", x = rownames(x = subset_cells), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = subset_cells, slot = 'counts')[mito.features, ])/Matrix::colSums(x = GetAssayData(object = subset_cells, slot = 'counts'))
subset_cells[["percent.mito"]] <- percent.mito
subset_cells <- subset(x = subset_cells, subset = nFeature_RNA >= 500 & nFeature_RNA < 7000 & percent.mito <= 0.1 & nCount_RNA >= 1000 & nCount_RNA <= 50000)

###---------------------------normalization and scale------------------###
subset_cells <- NormalizeData(subset_cells, verbose = TRUE)
subset_cells <- FindVariableFeatures(subset_cells, nfeatures = 2000, selection.method = 'vst', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf))
subset_cells <- ScaleData(subset_cells, verbose = TRUE,  vars.to.regress = c("nCount_RNA") )
subset_cells <- RunPCA(subset_cells, features = setdiff(subset_cells@assays$RNA@var.features, row.names(subset_cells) %>% grep(pattern = "^MT-|^RPL|^RPS|^IGK|^IGV|^IGL|^IGH", v = T)) , npcs = 50, verbose = TRUE)

###---------------------------loading metadata-------------------------------------###
metadat <- read.table("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/Sample_information_of_HCA2.0.csv", sep = ",", header = T, row.names = 1, stringsAsFactors = F)
subset_cells$Name <- mapvalues(subset_cells@meta.data$orig.ident, from = metadat %>% row.names, to = metadat$Name)
subset_cells$Donor <- mapvalues(subset_cells@meta.data$orig.ident, from = metadat %>% row.names, to = metadat$Donor)
subset_cells$Stage <- mapvalues(subset_cells@meta.data$orig.ident, from = metadat %>% row.names, to = metadat$Stage)
####--------------------------------Step 4. batch effect correction by BBKNN---------------------------###
PCA <- subset_cells@reductions$pca@cell.embeddings
batch <- subset_cells@meta.data$orig.ident

adata <- anndata$AnnData(X = PCA, obs = batch)
sc$tl$pca(adata, n_comps = as.integer(50), svd_solver = 'auto')
adata$obsm$X_pca <- PCA
bbknn$bbknn(adata, batch_key = 0, n_pcs = as.integer(50))
sc$tl$umap(adata)

sc$tl$leiden(adata, resolution = 1)
sc$pl$umap(adata, color = 'leiden', save = paste0("BBKNN_EBV_GC_", Projectname, "_data.pdf"))

umap <-  py_to_r(adata$obsm[["X_umap"]])
clusters <-  py_to_r(adata$obs)
subset_cells$seurat_clusters <- clusters$leiden

subset_cells@reductions$umap <- subset_cells@reductions$pca
subset_cells@reductions$umap@cell.embeddings <- umap
row.names(subset_cells@reductions$umap@cell.embeddings) <- row.names(subset_cells@reductions$pca@cell.embeddings)
colnames(subset_cells@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
subset_cells@reductions$umap@key <- "UMAP_"
subset_cells@reductions$umap@feature.loadings <- matrix(0)
##-----------change the pca embeding---------####
pca <-  py_to_r(adata$obsm[["X_pca"]])
subset_cells@reductions$pca@cell.embeddings <- pca
row.names(subset_cells@reductions$pca@cell.embeddings) <- row.names(subset_cells@reductions$umap@cell.embeddings)
colnames(subset_cells@reductions$pca@cell.embeddings) <- paste0("PC_", 1:50)
subset_cells@reductions$pca@key <- "PC_"
subset_cells@reductions$pca@feature.loadings <- matrix(0)

Idents(subset_cells) <- subset_cells$seurat_clusters
plan("multisession", workers = 1)
subset_cell <- subset_cells[, WhichCells(subset_cells, downsample = 500)]

result <- mclapply(as.numeric(levels(subset_cell@active.ident)),
									 FUN =  function(x) {FindMarkers(subset_cell, ident.1 = x, ident.2 = NULL, max.cells.per.ident = 500)},
									 mc.cores = 24)
RESULT <- result

roundN <- 1
while(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
	if(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
		recalculate_clusters <- which(mapply(length, result, SIMPLIFY = TRUE)!=5)-1
		print(recalculate_clusters)
		result1 <- mclapply(recalculate_clusters,
												FUN =  function(x) {FindMarkers(subset_cell, ident.1 = x, ident.2 = NULL, max.cells.per.ident = 500)},
												mc.cores = 35)
	}
	print(roundN + 1)
	for(i in 1:length(recalculate_clusters)){
		result[c(recalculate_clusters+1)[i]] <- result1[i]
	}
}

all_markers <- do.call(rbind, result)
all_markers$gene <- unlist(mapply(rownames, result))
all_markers$cluster <- rep(levels(subset_cell@active.ident), times = mapply(dim, result, SIMPLIFY = TRUE)[1,])
subset_cells.markers <- all_markers

subset_cells.markers %>% TOP_N(50, pct.1 = 0.2) -> top50
subset_cells.markers <- subset_cells.markers %>% TOP_N(5000)

write.table(top50, file = paste0(Projectname, "_top50_DEGs.csv"), sep = ",", row.names = T, quote = F)
write.table(subset_cells.markers, file = paste0(Projectname, "_culster_all_DEGs.csv"), sep = ",", row.names = T, quote = F)

###-----------------------------plot---------------------------------###

Dimpolt_fixed <- function(obj = subset_cells, file_names = "file_names", group.by = "ident", pt.size = 0.1, label.size = 5,
													raster = F, downsample = 100000, label = TRUE, cols = color_used, width = 15, height = 15, res = 300, pdf = FALSE, size = 5, ncol = 1, drawlegend = TRUE){
	txt <- "p <- DimPlot(object = obj[, WhichCells(obj, downsample = downsample)], reduction = 'umap', label = label, cols = cols, group.by = group.by,
								 pt.size = pt.size, label.size = label.size, raster = raster, repel = T) + guides(colour = guide_legend(ncol = ncol, override.aes = list(size = size)))
		        print(p)"
	txt1 <- "p <- DimPlot(object = obj[, WhichCells(obj, downsample = downsample)], reduction = 'umap', label = label, cols = cols, group.by = group.by,
								 pt.size = pt.size, label.size = label.size, raster = raster, repel = T) + NoLegend()
						print(p)"
	if(pdf == FALSE & drawlegend == TRUE){
		png(paste0(file_names, ".png"), width = width, height = height, units = "in", res = res)
		eval(parse(text = txt))
		dev.off()
	} else if (pdf == TRUE & drawlegend == TRUE) {
		pdf(paste0(file_names, ".pdf"), width = width, height = height)
		eval(parse(text = txt))
		dev.off()
	} else if (pdf == FALSE & drawlegend == FALSE) {
		png(paste0(file_names, ".png"), width = width, height = height, units = "in", res = res)
		eval(parse(text = txt1))
		dev.off()
	}
}

pt.size <- 0.1
Dimpolt_fixed(obj = subset_cells, file_names = paste0(Projectname, "_celltype_ident"), group.by = "ident", pt.size = pt.size)
Dimpolt_fixed(obj = subset_cells, file_names = paste0(Projectname, "_patient"), group.by = "Patients", pt.size = pt.size)
Dimpolt_fixed(obj = subset_cells, file_names = paste0(Projectname, "_porig.ident"), group.by = "orig.ident", pt.size = pt.size)
Dimpolt_fixed(obj = subset_cells, file_names = paste0(Projectname, "_Detail_Annotation_final"), group.by = "Detail_Annotation_final", pt.size = pt.size, label.size = 7)



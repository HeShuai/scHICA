library('SoupX')
library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)
library(parallel)
library(stringi)
library(data.table)
library(RColorBrewer)
library(DoubletFinder)
library(future)
library(ggsci)

color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jco()(10), pal_nejm()(8))[-8]

###-----------------------Step 1. romove the ambient RNA using Soupx-----------------------####
ags <- commandArgs(trailingOnly = T)
setwd("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/HCA2.0/All_sample_contamination_removed")

samplename <- ags[1] %>% str_split(pattern = "/", simplify = T) %>% `[`(10) %>% gsub(pattern = "_cDNA", replacement = "")

sc <-  load10X(ags[1])
sc <-  autoEstCont(sc)
out <-  adjustCounts(sc) %>% round()

###-----------------------Step 2. quality control-----------------------####
Seurat_object <- CreateSeuratObject(out, min.cells = 1, min.features = 0)

mito.features <- grep(pattern = "^MT-", x = rownames(x = Seurat_object), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = Seurat_object, slot = 'counts')[mito.features, ])/Matrix::colSums(x = GetAssayData(object = Seurat_object, slot = 'counts'))
TenXdat <- Seurat_object
TenXdat[["percent.mito"]] <- percent.mito
TenXdat <- subset(x = TenXdat, subset = nFeature_RNA >= 500 & nFeature_RNA < 7000 & percent.mito <= 0.1 & nCount_RNA >= 1000 & nCount_RNA <= 50000)

TenXdat <- NormalizeData(object = TenXdat, normalization.method = "LogNormalize", scale.factor = 1e4)
TenXdat <- FindVariableFeatures(object = TenXdat, nfeatures = 2000)
TenXdat <- ScaleData(object = TenXdat, features = VariableFeatures(object = TenXdat), vars.to.regress = c("nCount_RNA"))###only regress counts and variable featrues was used.
TenXdat <- RunPCA(object = TenXdat, features = VariableFeatures(object = TenXdat), verbose = FALSE)

subset_cells <- TenXdat

dim.use <- 30
res.use <- 1

subset_cells <- FindNeighbors(object = subset_cells, dims = 1:dim.use)
subset_cells <- FindClusters(object = subset_cells, resolution = res.use)

##------------------------5). Run the uMAP
subset_cells <- RunUMAP(object = subset_cells, dims = 1:dim.use, umap.method = "uwot")

##------------------------## plot the first run tSNE
pdf(paste0(samplename, "_", dim.use, "_", res.use, "_First_Run_tSNE.pdf"))
DimPlot(object = subset_cells, reduction = 'umap', label = TRUE, cols = color_used, pt.size = 1)
dev.off()

write.table(data.frame(Tissue = samplename, genes = dim(subset_cells@assays$RNA@data)[1], cells = dim(subset_cells@assays$RNA@data)[2]),
            file = paste0(samplename, "_before_dobuletfinder.txt"),
            sep = "\t", row.names = F, quote = F)

##------------------------Step 3. doublets removing using doubletfinder---------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(subset_cells, PCs = subset_cells@commands$FindNeighbors.RNA.pca$dims)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- bcmvn[which(max(bcmvn$BCmetric) == bcmvn$BCmetric),2] %>% as.character %>% as.numeric

##------------------------1). Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(subset_cells@active.ident)
nExp_poi <- round(0.07*length(colnames(subset_cells)))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

##------------------------2). Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
subset_cells <- doubletFinder_v3(subset_cells, PCs = subset_cells@commands$FindNeighbors.RNA.pca$dims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE)
subset_cells <- doubletFinder_v3(subset_cells, PCs = subset_cells@commands$FindNeighbors.RNA.pca$dims, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = grep("pANN", names(subset_cells@meta.data), value = T))

##------------------------3). Plot results --------------------------------------------------------------------------------------------------------------
high_of_low <- (subset_cells@meta.data[, grep("^DF\\.classifications", names(subset_cells@meta.data), value = T)] == "Singlet")+0
subset_cells@meta.data[high_of_low[, 1] + high_of_low[, 2] == 2, "DF_hi.lo"] <- "Singlet"
subset_cells@meta.data[high_of_low[, 1] + high_of_low[, 2] == 1, "DF_hi.lo"] <- "Doublet_lo"
subset_cells@meta.data[high_of_low[, 1] + high_of_low[, 2] == 0, "DF_hi.lo"] <- "Doublet_hi"

##------------------------4).remove the doublets and recluster-----------------#######
# subset_cells <- subset_cells[, subset_cells@meta.data[grepl(subset_cells$DF_hi.lo, pattern = "Singlet"), ] %>% row.names()]
subset_cells <- TenXdat[, subset_cells@meta.data[grepl(subset_cells$DF_hi.lo, pattern = "Singlet"), ] %>% row.names()]
##------------------------. recluster-----------------#######

write.table(data.frame(Tissue = samplename, genes = dim(subset_cells@assays$RNA@data)[1], cells = dim(subset_cells@assays$RNA@data)[2]),
            file = paste0(samplename, "_after_dobuletfinder.txt"),
            sep = "\t", row.names = F, quote = F)

subset_cells[["tissue"]] <- samplename
assign(samplename, subset_cells)

###it is so import that "list" args was used in save function.
save(list = samplename, file = paste0(samplename, ".RData"))



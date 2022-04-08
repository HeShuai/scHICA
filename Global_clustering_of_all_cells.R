###-------------------------------The script performed contaminated gene reomove, batch effect correction, and clustering

library(Seurat)
library(Matrix)
library(stringr)
library(gridExtra)
library(ggplot2)
library(stringi)
library(plyr)
library(ggthemes)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(reshape2)
library(scales)
library(rlang)
library(future)
library(parallel)
library(ggrepel)
library(ggsci)
library(ape)
library(dplyr)
library(harmony)
library(sigminer)
library(ggpubr)
library(rstatix)
library(ROGUE)
library(tidyverse)
library(ggsci)
library(reticulate)

####--------------------------------Step 1. enviroment seting---------------------------###

###-----------------------seting python environment--------------------------###
use_python("/public/home/HeShuai/anaconda3/envs/R4/bin/python")

###-----------------------loading python package-----------------##---###
anndata <- import("anndata",convert = FALSE)
bbknn <- import("bbknn", convert = FALSE)
sc <- import("scanpy",convert = FALSE)

options(future.globals.maxSize = 400*1000 * 1024^2)
plan("multiprocess", workers = 30)
plan()

color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10))[-8]

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
subset_cells <- FindVariableFeatures(merged_data, nfeatures = 2000, selection.method = 'vst', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf))
subset_cells <- ScaleData(subset_cells, verbose = TRUE,  vars.to.regress = c("nCount_RNA") )
subset_cells <- RunPCA(subset_cells, features = setdiff(subset_cells@assays$RNA@var.features, row.names(subset_cells) %>% grep(pattern = "^MT-|^RPL|^RPS|^IGK|^IGV|^IGL|^IGH", v = T)) , npcs = 50, verbose = TRUE)

###---------------------------loading metadata-------------------------------------###
metadat <- read.table("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/Sample_information_of_HCA2.0.csv", sep = ",", header = T, row.names = 1, stringsAsFactors = F)
subset_cells$Name <- mapvalues(subset_cells@meta.data$orig.ident, from = metadat %>% row.names, to = metadat$Name)
subset_cells$Donor <- mapvalues(subset_cells@meta.data$orig.ident, from = metadat %>% row.names, to = metadat$Donor)
subset_cells$Stage <- mapvalues(subset_cells@meta.data$orig.ident, from = metadat %>% row.names, to = metadat$Stage)

####--------------------------------Step 4. batch effect correction by BBKNN---------------------------###

###----------------------------batch effect correction by BBKNN------------------------------###
PCA <- subset_cells@reductions$pca@cell.embeddings
Batch <- subset_cells@meta.data$Donor

adata <- anndata$AnnData(X = PCA, obs = Batch)
sc$tl$pca(adata)
adata$obsm$X_pca <- PCA
bbknn$bbknn(adata, batch_key = 0)
sc$tl$umap(adata)

sc$tl$leiden(adata, resolution = 2)
sc$pl$umap(adata, color = 'leiden', save = "BBKNN.pdf")

###----------------------------mapping the BBKNN objects to corresponding slots in Seurat------------------------------###
umap <-  py_to_r(adata$obsm[["X_umap"]])
clusters <-  py_to_r(adata$obs)

subset_cells$seurat_clusters <- clusters$leiden

subset_cells@reductions$umap <- subset_cells@reductions$pca
subset_cells@reductions$umap@cell.embeddings <- umap
row.names(subset_cells@reductions$umap@cell.embeddings) <- row.names(subset_cells@reductions$pca@cell.embeddings)
colnames(subset_cells@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
subset_cells@reductions$umap@key <- "UMAP_"
subset_cells@reductions$umap@feature.loadings <- matrix(0)

Idents(subset_cells) <- subset_cells$seurat_clusters

png("HCA2.0_6samples_colored_by_cluster.png",
    width = 15, height = 15, units = "in", res = 300)
p5 <- DimPlot(object = subset_cells[, WhichCells(subset_cells, downsample = 3000)], reduction = 'umap', label = TRUE, cols = colors_dxs, pt.size = 0.0001, label.size = 5, raster = F) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
print(p5)
dev.off()

png("HCA2.0_6samples_colored_by_tissue_namer.png",
    width = 15, height = 15, units = "in", res = 300)
p5 <- DimPlot(object = subset_cells[, WhichCells(subset_cells, downsample = 30000)], reduction = 'umap', label = TRUE, cols = color_used, group.by = "Name", pt.size = 0.0001, label.size = 5, raster = F) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
print(p5)
dev.off()

png("HCA2.0_6samples_colored_by_adult_or_fetal.png",
    width = 15, height = 15, units = "in", res = 300)
p5 <- DimPlot(object = subset_cells[, WhichCells(subset_cells, downsample = 30000)], reduction = 'umap', label = TRUE, cols = color_used, group.by = "Stage", pt.size = 0.0001, label.size = 5, raster = F) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
print(p5)
dev.off()

png("HCA2.0_6samples_colored_by_donor.png",
    width = 15, height = 15, units = "in", res = 300)
p5 <- DimPlot(object = subset_cells[, WhichCells(subset_cells, downsample = 30000)], reduction = 'umap', label = TRUE, cols = color_used, group.by = "Donor", pt.size = 0.0001, label.size = 5, raster = F) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
print(p5)
dev.off()

####--------------------------------Step 5. DEGs calculation--------------------------###

###-------------------------------------calculate the DEGs-------------------------------###
plan("multiprocess", workers = 1)

subset_cell <- subset_cells[, WhichCells(subset_cells, downsample = 500)]###downsampling before calculating DEGs

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

write.table(top50,
            file = "HCA2.0_subset_cells_culster_top50_DEGs_allcell.csv",
            sep = ",",
            row.names = T,
            quote = F)

write.table(subset_cells.markers,
            file = "HCA2.0_allcell_subset_cells_culster_all_DEGs_allcell.csv",
            sep = ",",
            row.names = T,
            quote = F)

png(paste0("featureplot_all_cells.png"), height = 20, width = 25, units = "in", res = 600)
FeaturePlot(subset_cells[, WhichCells(subset_cells, downsample = 3000)],
            features =  c("PTPRC", "CD3D" , "CD3E", "MS4A1", "CD79A", "CD14", "C1QC", "MKI67", 
                          "CLEC10A", "CLEC9A", "PPBP", "CSF3R"), ## active and repression marker genes
            # features = c("MKI67", "PCNA", "TYMS"),
            pt.size = 0.0001,
            ncol = 4,
            # cols = c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12))[],
            reduction = "umap",
            raster = F)  
dev.off()


################------------------------------------------round two:remove the potential contamination, and then re-run round one for final clustering------------------------------------------##############################

removed_gene_list <- c('ACTA1',
                       'ACTN2',
                       'ATP2A1',
                       'CA3',
                       'CKM',
                       'COX6A2',
                       'COX7A1',
                       'CRYAB',
                       'DES',
                       'EEF1A2',
                       'ENO3',
                       'FXYD1',
                       'HSPB6',
                       'KLHL41',
                       'MYBPC1',
                       'MYL1',
                       'MYL2',
                       'MYLPF',
                       'MYOZ1',
                       'NEB',
                       'PGAM2',
                       'RYR1',
                       'SLN',
                       'TCAP',
                       'TCEA3',
                       'TNNC1',
                       'TNNC2',
                       'TNNI1',
                       'TNNI2',
                       'TNNT3',
                       'TPM1',
                       'TPM2',
                       'ACTC1',
                       'MYL5',
                       'MYOZ2',
                       'MYO7A')

Fetal_muscle_cells <- subset_cells@meta.data %>% subset(Name == "Muscle" & Stage == "Fetal") %>% row.names()
subset_cells@assays$RNA@counts[removed_gene_list, Fetal_muscle_cells] <- 0

###-------------------then reomeve the doublets cluster and run the "Round one" step for the final clustering.

###-----------------------------------------------------------------------------------------------------------------------------###
colors_dxs <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', 
                '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', 
                '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', 
                '#C1E6F3', '#6778AE', '#91D0BE', 
                '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175',
                '#FEE453', '#FE9955', '#FC6FB3', '#72F717', '#D0EAA4', '#C8FF8A', 
                '#53ECFE','#FBD1B9', '#DAB8F9', '#A849FD', '#54A4DE', '#43E3B2', 
                '#F3B3DB', '#F49089', '#FF8787', 
                '#B3CDFE', '#E4C755', '#6565FE', '#66B2CD','#F1BB72', '#E63863', '#5EF196', 
                '#DF57B1', '#6778AE', '#3750A8', '#B53E2B', '#CA6542',
                '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175',
                '#FEE453', '#FE9955', '#53ECFE', '#DAB8F9', '#A849FD', '#54A4DE', '#43E3B2', 
                '#F3B3DB', '#F49089', '#FF8787', 
                '#B3CDFE', '#E4C755', '#6565FE', '#66B2CD','#F1BB72', '#E63863', '#5EF196', 
                '#DF57B1', '#6778AE', '#3750A8', '#B53E2B', '#CA6542',
                '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')


TOP_N <- function(x, n, pct.1 = 0.1, first = "avg_log2FC", second = "p_val_adj", sig.padj = NULL, fc.threshold = c(-0.25, 0.25)){
  if(table(names(x) %in% c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster"))["TRUE"] == 7){
    reordered_daf <- data.frame()
    if (!is.null(pct.1) & (pct.1 > 0 & pct.1 < 1)) {
      x <- x[x$pct.1 >= pct.1, ]
    }
    if (!is.null(sig.padj) & is.numeric(sig.padj)) {
      x <- x[x$p_val_adj <= sig.padj, ]
    }
    if (length(fc.threshold) == 2) { 
      fc.threshold <- sort(fc.threshold)
      x <- x[(x$avg_log2FC < fc.threshold[1] | x$avg_log2FC > fc.threshold[2]), ]
    } else if (length(fc.threshold) == 1) {
      x <- x[x$avg_log2FC >= fc.threshold, ]
    }
    
    for (i in unique(x$cluster)) {
      if(n > dim(x[x$cluster == i, ])[1]){
        message("FBI warning, n < ", n)
        tmp <- dplyr::arrange(.data = x[x$cluster == i, ], desc(get(first)), desc(get(second)))[1:(dim(x[x$cluster == i, ])[1]), ]
        reordered_daf <- rbind(reordered_daf, tmp)
      } else {
        tmp <- dplyr::arrange(.data = x[x$cluster == i, ], desc(get(first)), desc(get(second)))[1:n, ]
        reordered_daf <- rbind(reordered_daf, tmp)
      }
    }
  }
  else {
    message("be careful !")
  }
  return(reordered_daf)
}


sessionInfo()
R version 4.1.2 (2021-11-01)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /public/home/HeShuai/anaconda3/envs/R4/lib/libopenblasp-r0.3.18.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  grid      stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] reticulate_1.24       forcats_0.5.1         purrr_0.3.4          
 [4] readr_2.1.2           tidyr_1.2.0           tibble_3.1.6         
 [7] tidyverse_1.3.1       ROGUE_1.0             rstatix_0.7.0        
[10] ggpubr_0.4.0          sigminer_2.1.2        Biobase_2.54.0       
[13] BiocGenerics_0.40.0   harmony_0.1.0         Rcpp_1.0.8.3         
[16] dplyr_1.0.8           ape_5.6-1             ggsci_2.9            
[19] ggrepel_0.9.1         future_1.24.0         rlang_1.0.2          
[22] scales_1.1.1          reshape2_1.4.4        pheatmap_1.0.12      
[25] circlize_0.4.14       ComplexHeatmap_2.10.0 RColorBrewer_1.1-2   
[28] data.table_1.14.2     cowplot_1.1.1         ggthemes_4.2.4       
[31] plyr_1.8.6            stringi_1.7.6         ggplot2_3.3.5        
[34] gridExtra_2.3         stringr_1.4.0         Matrix_1.4-0         
[37] SeuratObject_4.0.4    Seurat_4.1.0         

loaded via a namespace (and not attached):
  [1] readxl_1.3.1          backports_1.4.1       NMF_0.23.0           
  [4] igraph_1.2.11         lazyeval_0.2.2        splines_4.1.2        
  [7] listenv_0.8.0         scattermore_0.8       gridBase_0.4-7       
 [10] digest_0.6.29         foreach_1.5.2         htmltools_0.5.2      
 [13] fansi_1.0.2           magrittr_2.0.2        tensor_1.5           
 [16] cluster_2.1.2         doParallel_1.0.17     ROCR_1.0-11          
 [19] tzdb_0.2.0            globals_0.14.0        modelr_0.1.8         
 [22] matrixStats_0.61.0    spatstat.sparse_2.1-0 colorspace_2.0-3     
 [25] rvest_1.0.2           haven_2.4.3           crayon_1.5.0         
 [28] jsonlite_1.8.0        spatstat.data_2.1-2   survival_3.3-1       
 [31] zoo_1.8-9             iterators_1.0.14      glue_1.6.2           
 [34] polyclip_1.10-0       registry_0.5-1        gtable_0.3.0         
 [37] leiden_0.3.9          GetoptLong_1.0.5      car_3.0-12           
 [40] future.apply_1.8.1    shape_1.4.6           abind_1.4-5          
 [43] DBI_1.1.2             rngtools_1.5.2        spatstat.random_2.1-0
 [46] miniUI_0.1.1.1        viridisLite_0.4.0     xtable_1.8-4         
 [49] clue_0.3-60           spatstat.core_2.4-0   stats4_4.1.2         
 [52] htmlwidgets_1.5.4     httr_1.4.2            ellipsis_0.3.2       
 [55] ica_1.0-2             pkgconfig_2.0.3       dbplyr_2.1.1         
 [58] uwot_0.1.11           deldir_1.0-6          here_1.0.1           
 [61] utf8_1.2.2            tidyselect_1.1.2      later_1.3.0          
 [64] cellranger_1.1.0      munsell_0.5.0         tools_4.1.2          
 [67] cli_3.2.0             generics_0.1.2        broom_0.7.12         
 [70] ggridges_0.5.3        fastmap_1.1.0         goftest_1.2-3        
 [73] fs_1.5.2              fitdistrplus_1.1-8    RANN_2.6.1           
 [76] pbapply_1.5-0         nlme_3.1-155          mime_0.12            
 [79] xml2_1.3.3            rstudioapi_0.13       compiler_4.1.2       
 [82] plotly_4.10.0         png_0.1-7             ggsignif_0.6.3       
 [85] spatstat.utils_2.3-0  reprex_2.0.1          lattice_0.20-45      
 [88] vctrs_0.3.8           pillar_1.7.0          lifecycle_1.0.1      
 [91] furrr_0.2.3           spatstat.geom_2.3-2   lmtest_0.9-40        
 [94] GlobalOptions_0.1.2   RcppAnnoy_0.0.19      irlba_2.3.5          
 [97] httpuv_1.6.5          patchwork_1.1.1       R6_2.5.1             
[100] promises_1.2.0.1      KernSmooth_2.23-20    IRanges_2.28.0       
[103] parallelly_1.30.0     codetools_0.2-18      MASS_7.3-55          
[106] assertthat_0.2.1      rprojroot_2.0.2       pkgmaker_0.32.2      
[109] rjson_0.2.21          withr_2.5.0           sctransform_0.3.3    
[112] S4Vectors_0.32.3      hms_1.1.1             mgcv_1.8-39          
[115] rpart_4.1.16          carData_3.0-5         Rtsne_0.15           
[118] lubridate_1.8.0       shiny_1.7.1          

packages <- list("Seurat", "Matrix", "stringr", "stringi", "gridExtra", "ggplot2",
								 "plyr", "ggthemes", "cowplot", "data.table", "RColorBrewer", "ComplexHeatmap",
								 "circlize", "pheatmap", "reshape2", "scales", "rlang", "future",
								 "parallel", "ggrepel", "ggsci", "ape", "dplyr", "harmony",
								 "sigminer", "ggpubr", "rstatix", "ROGUE", "tidyverse", "reticulate", "dendextend",
								 "scRepertoire", "ggdendro", "ggtern", "viridis", "scatterpie", "SeuratData", "SeuratDisk",
								 "scibet", "scuttle", "zellkonverter")
lapply(packages, library, character.only = TRUE)
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]
options(width = 240)
source("/public/home/HeShuai/Projects/Script/1.TOP_N.R")
options(future.globals.maxSize = 400*1000 * 1024^2)
plan("multisession", workers = 20)
plan()
###-----------------------seting python environment--------------------------###
use_python("/public/home/HeShuai/miniconda3/bin/python")
###-----------------------loading python package-----------------##---###
anndata <- import("anndata",convert = FALSE)
bbknn <- import("bbknn", convert = FALSE)
sc <- import("scanpy",convert = FALSE)

####--------------------------------------------------------------Construct the signature gene set for naive CD4/8 and Treg using published dataset--------------------------------------------------------###
###-----------------------loading python package done-----------------##---###
ad <- readH5AD("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/SCIENCE_Adult_and_Fetal_2022/PAN.A01.v01.raw_count.20210429.NKT.embedding.h5ad")
adata_Seurat <- as.Seurat(ad, counts = "X", data = NULL)
subset_cells <- subset(adata_Seurat, subset = celltype_annotation %in% c("CD4+T", "CD8+T", "DP(P)_T", "DP(Q)_T", "TREG", "ABT(ENTRY)", "DN(P)_T", "DN(Q)_T", "DN(early)_T", "TREG"))

tmp <- CreateSeuratObject(subset_cells@assays$originalexp@counts, project = "Fetal_all_T", min.cells = 1, min.features = 0)
tmp@meta.data <- subset_cells@meta.data
subset_cells <- tmp
###---------------------optional step end---------------------------##
subset_cells <- NormalizeData(subset_cells, verbose = T)
subset_cells <- FindVariableFeatures(subset_cells, nfeatures = 2000, selection.method = 'vst', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf))
subset_cells <- ScaleData(subset_cells)#, vars.to.regress = c("nCount_originalexp"), verbose = TRUE)#"nCount_RNA"
subset_cells <- RunPCA(subset_cells, features = setdiff(subset_cells@assays$RNA@var.features, row.names(subset_cells) %>% grep(pattern = "^MT-|^RPL|^RPS|^IGK|^IGV|^IGL|^IGH", value = T)) , npcs = 50, verbose = TRUE)
###----------------------------batch effect correction by BBKNN------------------------------###
subset_cells$annotation <- subset_cells$celltype_annotation
Idents(subset_cells) <- subset_cells$annotation
subset_cells <- RenameIdents(subset_cells, "TREG" = "Treg",
														 "DN(P)_T" = "Precursor", "DN(Q)_T" = "Precursor", "DN(early)_T" = "Precursor", "DP(P)_T" = "Precursor",
														 "ABT(ENTRY)" = "Precursor", "DP(Q)_T" = "Precursor",
														 "CD4+T" = "CD4", "CD8+T" = "CD8")
subset_cells$annotation <- Idents(subset_cells)
selected_cells <- setdiff(subset_cells %>% colnames, subset_cells@meta.data %>% filter(annotation == "Treg", organ == "TH") %>% rownames())
subset_cells <- subset_cells[, selected_cells]

Idents(subset_cells) <- subset_cells$seurat_clusters
plan("multisession", workers = 1)

subset_cell <- subset_cells[, WhichCells(subset_cells, downsample = 1000)]

result <- mclapply(levels(subset_cell@active.ident),
									 FUN =  function(x) {FindMarkers(subset_cell, ident.1 = x, ident.2 = NULL, max.cells.per.ident = 1000)},
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

write.table(top50, file = "Science_top50_DEGs_publication_atlas_T_precursor.csv", sep = ",", row.names = T, quote = F)
write.table(subset_cells.markers, file = "Science_subset_cells_culster_all_DEGs_publication_atlas_T_precursor.csv", sep = ",", row.names = T, quote = F)

###--------------------------save the data----------------------------###
subset_cells@assays$RNA@data <- as.matrix(0)
subset_cells@assays$RNA@scale.data <- as.matrix(0)
save(subset_cells, file = "Science_all_Fetal_all_T.RData")


####-----------------------------------------------------------------calculate the module score--------------------------------------------------------##
###----------------------------------Add modulescore-------------------------------###
load("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/Figure.1/2024.08.04.Selected_T_cells.58476.RData")
DefaultAssay(subset_cells) <- "RNA"
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features = VariableFeatures(subset_cells)[1])
subset_cells <- subset_cells %>% subset(Global_annotation %in% c("CD4 T TCF7", "Precursor T ITM2A"))
subset_cells$Annotation <- subset_cells$Global_annotation %>% as.character()
subset_cells@meta.data$Annotation[subset_cells$Name != "Thymus" &  subset_cells$Annotation == "Precursor T ITM2A"] <- "Other_tissue"

subset_cells <- AddModuleScore(subset_cells, features = list(CD4 = subset_cells.markers %>% filter(cluster == "CD4") %>% filter(avg_log2FC > 0.5) %>% select(gene) %>% unlist %>% unname,
Precursor = subset_cells.markers %>% filter(cluster == "Precursor") %>% filter(avg_log2FC > 0.5) %>% select(gene) %>% unlist %>% unname,
Treg = subset_cells.markers %>% filter(cluster == "Treg") %>% filter(avg_log2FC > 0.5) %>% select(gene) %>% unlist %>% unname), name = c("CD4_score", "Precursor_score", "Treg"))

subset_cells@meta.data %>% write.table(file = "Modulescore.precursor.mature.T.txt", sep = "\t", col.names=T, row.names =T, quote = F)

png('Modulescore.png', res = 300, units = "in", height = 8, width = 6)
VlnPlot(subset_cells, features = c("CD4_score1", "Precursor_score2", "Treg3"), group.by = "Annotation", pt.size = 0)
dev.off()

Modulescore_data <- read.table("Modulescore.precursor.mature.T.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
selected_data <- Modulescore_data %>% select(Annotation, CD4_score1, Precursor_score2) %>% melt()
colnames(selected_data)[2:3] <- c("Cell_type", "Modulescore")
selected_data$Annotation <- factor(selected_data$Annotation, levels = c("Precursor T ITM2A", "Other_tissue", "CD4 T TCF7"))

p <- ggviolin(selected_data, x = "Annotation", y = "Modulescore",
					fill = "Annotation", add = "boxplot",
				 facet.by = "Cell_type", add.params = list(fill = "white"),
				 trim = F, size = 0) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
	scale_y_continuous(expand = c(0.01, 0))

pdf("ModuleScore_T_precursor.pdf", height = 7, width = 6)
p +  geom_pwc(
	aes(group = Annotation), tip.length = 0,
	method = "t_test", label = "p.adj.format",
	p.adjust.method = "fdr",
	bracket.nudge.y = -0.08,
	p.adjust.by = "group",
	hide.ns = FALSE) +
	scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) #+
# coord_flip() 
dev.off()


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
#####------------------------------define a function to returen data ordered by considering avg firstly---------------###### 
load("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/HCA2.0_alltissues_final_annotation_1038993.RData")
DefaultAssay(subset_cells) <- "RNA"
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features = VariableFeatures(subset_cells)[1])

selected_cells <- subset_cells %>% subset(Global_annotation %in% c("Precursor T ITM2A", "Treg", "CD4 T TCF7") & Stage == "Fetal")
selected_cells$Name <- selected_cells$Name %>% gsub(pattern = "^Colon$", replacement = "Jejunum")
selected_cells$Name <- selected_cells$Name %>% gsub(pattern = "^Stomach Protein$", replacement = "Stomach")
subset_cells <- selected_cells

subset_cells@assays$RNA@data <- as.matrix(0)
subset_cells@assays$RNA@scale.data <- as.matrix(0)
save(subset_cells, file = "2024.08.04.Selected_T_cells.58476.RData")

####------------------------------------------only show the developmental trajectory from TP to ET cells--------------------------------------------##
###-------------------------------------------------------------------Monocle3------------------------------------------------------------------######
load("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/Figure.1/2024.08.04.Selected_T_cells.58476.RData")
DefaultAssay(subset_cells) <- "RNA"
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features = VariableFeatures(subset_cells)[1])

subset_cells <- subset_cells %>% subset(Global_annotation == "Precursor T ITM2A")
expression_matrix <- subset_cells@assays$RNA@counts %>% as.matrix()
cell_metadata <- subset_cells@meta.data
cell_metadata$orig <- cell_metadata$Name
cell_metadata$orig[cell_metadata$orig != "Thymus"] <- "Other_tissue"

gene_annotation <- data.frame(gene_short_name = row.names(subset_cells@assays$RNA@counts), row.names = row.names(subset_cells@assays$RNA@counts), stringsAsFactors = F)

cds <- new_cell_data_set(expression_matrix,
												 cell_metadata = cell_metadata,
												 gene_metadata = gene_annotation)

cells <- cell_metadata %>% filter(Global_annotation == "Precursor T ITM2A") %>% rownames()
cds <- cds[, cells]

cds <- preprocess_cds(cds, num_dim = 30) 
cds <- align_cds(cds, alignment_group = "tissue", umap.min_dist = 0.5, umap.n_neighbors = 20L)
cds <- reduce_dimension(cds)

reducedDims(cds)[["UMAP"]] <- Embeddings(object = subset_cells[["umap"]])[cells, ]
colnames(reducedDims(cds)[["UMAP"]]) <- NULL

cds@reduce_dim_aux$UMAP$model$umap_model$embedding <- Embeddings(object = subset_cells[["umap"]])[cells, ]
colnames(cds@reduce_dim_aux$UMAP$model$umap_model$embedding) <- NULL

cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = T)

cds <- cds[, cds@reduce_dim_aux$UMAP@listData$model$umap_model$embedding[, 1] < 4.5 & cds@reduce_dim_aux$UMAP@listData$model$umap_model$embedding[, 2] > 4.5]

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin = "Thymus"){
	cell_ids <- which(colData(cds)[, "orig"] == time_bin)
	closest_vertex <-
		cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
	closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
	root_pr_nodes <-
		igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
																															(which.max(table(closest_vertex[cell_ids,]))))]
	root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))

cds$Global_annotation <- factor(cds$Global_annotation, levels = c("Thymus", "Other_tissue"))

pdf("Trajectory_monocle3.cluster.pdf", height = 10, width = 10)
# png("Trajectory_monocle3.cluster.png", height = 10, width = 10, res = 400, units = "in")
plot_cells(cds,
					 label_groups_by_cluster = T,
					 color_cells_by = "orig",
					 label_cell_groups = F,
					 label_leaves = T,
					 label_branch_points = T,
					 cell_size = 0.1,
					 trajectory_graph_segment_size = 0.3,
					 graph_label_size = 2) +
	scale_color_manual(values = color_used)
dev.off()

pdf("monocyte3_Precursor_psedutime_test_1.pdf", height = 10, width = 10)
# png("monocyte3_Precursor_psedutime_test_1.png", height = 10, width = 10, res = 300, units = "in")
plot_cells(cds,
					 color_cells_by = "pseudotime",
					 label_cell_groups = FALSE,
					 label_leaves = FALSE,
					 label_branch_points = FALSE,
					 graph_label_size = 1.5,
					 cell_size = 0.1,
					 trajectory_graph_segment_size = 0.3)
dev.off()

png("monocyte3_Precursor_trajectory_annotation_test.png", height = 10, width = 10, res = 300, units = "in")
# pdf("monocyte3_Precursor_annotation.pdf", height = 10, width = 10)
plot_cells(cds,
					 color_cells_by = "orig",
					 label_cell_groups = F,
					 label_leaves = T,
					 label_branch_points = T,
					 cell_size = 0.3,
					 graph_label_size = 2
) + scale_color_manual(values = color_used)
# scale_color_manual(values = eval(parse(text = color_used)))
dev.off()

# png("monocyte3_Precursor_trajectory_annotation_tissue.png", height = 10, width = 10, res = 300, units = "in")
pdf("monocyte3_Precursor_trajectory_annotation_tissue.pdf", height = 10, width = 10)
plot_cells(cds,
					 color_cells_by = "Name",
					 label_cell_groups = F,
					 label_leaves = F,
					 label_branch_points = F,
					 graph_label_size = 2,
					 cell_size = 0.1,
					 trajectory_graph_segment_size = 0.3
) + scale_color_manual(values = color_used)
# scale_color_manual(values = eval(parse(text = color_used)))
dev.off()

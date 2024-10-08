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
###------------------------HSC---------------------------###
# load("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/HCA2.0_alltissues_final_annotation_1038993.RData")
DefaultAssay(subset_cells) <- "RNA"
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features = VariableFeatures(subset_cells)[1])
Idents(subset_cells) <- subset_cells$Global_annotation
selected_cells <- subset_cells %>% subset(downsample  = 1000)

###----------------------RORE calculation-----------------------------###
library(gmodels)
meta.data <- read.table("All_cells.meta.data.txt", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
selected <- meta.data %>% select(Global_annotation, Name, Stage) %>% as.tibble()
selected[selected$Name %>% grepl(pattern = "^Colon$"), "Name"] <- "Jejunum"
selected[selected$Name %>% grepl(pattern = "Stomach Protein"), "Name"] <- "Stomach"
selected$Samples <- paste0(selected$Stage, "_", selected$Name)

cell_number <- table(selected$Samples, selected$Global_annotation) %>% unclass() %>% as.data.frame()
RO_RE_cell_number <- CrossTable(cell_number %>% as.matrix(), expected = T)
RO_RE <- (RO_RE_cell_number$chisq$observed/RO_RE_cell_number$chisq$expected)
RO_RE[(RO_RE_cell_number$chisq$observed <= 10)] <- 0
# RO_RE[(RO_RE_cell_number$CST$observed <= 10) + (RO_RE_cell_number$CST$expected) != 0 ] <- 0

RO_RE[is.na(RO_RE)] <- 0
RO_RE[is.infinite(RO_RE %>% as.matrix())] <- 0
# RO_RE <- RO_RE[, Order_Myeloid_RORE]
RO_RE[RO_RE > 10] <- 10

RO_REnum <- RO_RE %>% round(2)
RO_REnum[RO_REnum == 0] <- ""

annotation_row <- data.frame(Stage = rep(c("Adult", "Fetus"), times = c(22, 16)))
row.names(annotation_row) <- row.names(RO_RE)
	
annotation_col <- data.frame(Celltype = rep(x = c("B/plasma", "T/NK/ILC", "T/NK/ILC", "DC",  "T/NK/ILC", "Other Precursor",
																							 "T/NK/ILC", "DC", "Monocyte/macrophage", "Mast cell", "Monocyte/macrophage",
																							 "B/plasma", "Neutrophil", "T/NK/ILC", "DC", "B/plasma", "B/plasma",
																							 "Other Precursor", "Neutrophil", "T/NK/ILC", "B/plasma", "T/NK/ILC"),
																				 times = c(1, 4, 5, 3, 2, 1, 1, 1, 7, 1, 2, 3, 3, 2, 1, 3, 2, 1, 1, 2, 1, 4)))
row.names(annotation_col) <- colnames(RO_RE)

annotation_color <- list(Celltype = c("Monocyte/macrophage" = color_used[43] %>% gsub(pattern = "FF", replacement = ""),
														 "T/NK/ILC" = color_used[44] %>% gsub(pattern = "FF", replacement = ""),
														 "Neutrophil" = color_used[45] %>% gsub(pattern = "FF", replacement = ""),
														 "DC" = color_used[46] %>% gsub(pattern = "FF", replacement = ""),
														 "Mast cell" = color_used[47] %>% gsub(pattern = "FF", replacement = ""),
														 "B/plasma" = color_used[49] %>% gsub(pattern = "FF", replacement = ""),
														 "Other Precursor" = color_used[48] %>% gsub(pattern = "FF", replacement = "")),
														 Stage = c(Adult = color_used[c(6)] %>% gsub(pattern = "FF", replacement = ""),
														 					Fetus = c(color_used[c(7)] %>% gsub(pattern = "FF", replacement = ""))))

row.names(annotation_col) <- colnames(RO_RE)

pdf("RORE_all_cells.pdf", width = 10.7, height = 6)
pheatmap(RO_RE,
				 # color = colorRampPalette(c("#181495", "blue", "white", "red", "#930107"))(100),
				 color = colorRampPalette(c( "navy", "white", "firebrick3"))(50),
				 # color = colorRampPalette(c("dodgerblue4", "peachpuff", "deeppink4"))(100),
				 breaks = seq(-5, 5, length.out = 50),
				 scale = "none",
				 cluster_rows = T,
				 cluster_cols = T,
				 fontsize_row = 8,
				 fontsize_col = 8,
				 clustering_distance_cols = "manhattan",
				 clustering_distance_rows = "euclidean",
				 # display_numbers = RO_REnum,
				 annotation_col = annotation_col,
				 annotation_row = annotation_row,
				 annotation_colors = annotation_color,
				 treeheight_row = 18,
				 treeheight_col = 18,
				 # gaps_col = c(6, 16, 22, 29, 31, 32, 34),
				 # gaps_row = c(16, 21),
				 # cutree_col = 2,
				 # cutree_rows = 3,
				 border_color = "white")
dev.off()

packages <- list("Seurat", "Matrix", "stringr", "stringi", "gridExtra", "ggplot2",
								 "plyr", "ggthemes", "cowplot", "data.table", "RColorBrewer", "ComplexHeatmap",
								 "circlize", "pheatmap", "reshape2", "scales", "rlang", "future",
								 "parallel", "ggrepel", "ggsci", "ape", "dplyr", "harmony",
								 "sigminer", "ggpubr", "rstatix", "tidyverse", "reticulate", "dendextend",
								 "scRepertoire", "ggdendro", "ggtern", "viridis", "scatterpie")

lapply(packages, library, character.only = TRUE)
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]
options(width = 240)

tissues_colors <- c("Ascending colon" = "#E64B35FF", "Sigmoid colon" = '#ff0000ff', "Transverse colon" = '#d62790',
										Ileum = "#91D1C2FF", Jejunum = '#4DBBD5FF', Duodenum = '#00A087FF',
										"Lesser gastric curvature" = '#F0E685FF', "Mesenteric lymph node" = "#FFB570FF",
										Blood = "#1B1919FF", PBMC = "#374E55FF", "Whole blood" = "#3F4041FF",	Marrow = '#466983FF', Spleen = "#725663FF",
										Lung = "#6A6599FF", Liver = '#749B58FF', Kidney = '#B09C85FF', Skin = '#D6D6CEFF', Bladder = "#D49464FF",
										Esophagus = '#c5b0d5', Trachea = "#BA63affF", Heart = '#9467bd', Muscle = "#aec7e8", Pancreas = "#5050FFFF",
										Stomach = '#ADB17DFF', Thymus = '#8c564b')

color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]
options(width = 220)



###-------------------------------------The following script is an demo for analysing TCR sharing between organs for CD4 and CD8 T cells in adults-------------------------------------###
###------Note: A similar strategy was applied to analyze TCR sharing between organs for CD4/CD8 TEFF/EM and TRM cells in adults---------------------------------###


#####-----------------------------------------TCR tracking analysis tissues--------------------------------###
TCR_meta.data <- read.table("Merged_TCR.information_at_least_one_paired_VDJ_cellranger6.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
TCR_meta.data$renamed_Samlename <- TCR_meta.data$Name
TCR_meta.data$renamed_Samlename[TCR_meta.data$renamed_Samlename %>% grepl(pattern = "^Colon$")] <- "Jejunum"
TCR_meta.data$renamed_Samlename[TCR_meta.data$renamed_Samlename %>% grepl(pattern = "Stomach Protein")] <- "Stomach"
TCR_meta.data$renamed_Samlename[TCR_meta.data$renamed_Samlename %>% grepl(pattern = ".*blood|.*PBMC.*")] <- "Blood"

TCR_meta.data <- TCR_meta.data %>% filter(! Detail_Annotation_final %in% c("DG T1", "DG T2", "DG T3", "ILCs", "MAIT",
																																					 "NK CD16 KIR2DL3", "NK CD16 MYOM2", "NK CD56 CCL3", "NK CD56 KIR2DL4", "Pro ILC", "Pro NK"))
###--------------------TCR composition of each tissue--------------------###
T_cell_colon_unique_barcode <- TCR_meta.data %>% select(barcode, customer_clone, sample, renamed_Samlename, Detail_Annotation_final, Stage) %>% unique
CD4 <- T_cell_colon_unique_barcode %>% filter(grepl(Detail_Annotation_final, pattern = "CD4|Treg")) %>% filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8")) 
CD8 <- T_cell_colon_unique_barcode %>% filter(grepl(Detail_Annotation_final, pattern = "CD8")) %>% filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8"))
###--------------------TCR expansion of each tissue--------------------###
TCR <- CD8 %>%
	dplyr::select(c(barcode, customer_clone, renamed_Samlename, sample, Detail_Annotation_final, Stage)) %>%
	# filter(Detail_Annotation_final %>% grepl(., pattern = "TRM")) %>% filter(! (renamed_Samlename %>% grepl(., pattern = "Blood|Kidney|Marrow"))) %>%
	# filter(Detail_Annotation_final %>% grepl(., pattern = "TEM|TEFF")) %>% filter(Stage == "Adult") %>%  
	setnames(old = c("barcode", "customer_clone", "renamed_Samlename", "sample", "Detail_Annotation_final"), new = c("Cell_Name", "clone.id", "majorCluster", "patient", "loc")) %>% unique()

###-----------------------------------Adult CD8-----------------------------------###
TCR$loc <- 1 ### set loc equal to 1, and using transition index between organs instead of migration
obj <- Startrac.run(TCR %>% filter(Stage == "Adult"), proj = "HCA2.0", verbose = T)
Startrac::plot(obj, index.type = "cluster.all", byPatient = T)

dat_for_plot_CD8_ad_migration_ave <- obj@pIndex.tran %>% filter(aid == "HCA2.0")

tmp_CD8_Ad_AVE <- dat_for_plot_CD8_ad_migration_ave %>% `[`(, -c(1))
tmp_CD8_Ad_AVE <- column_to_rownames(tmp_CD8_Ad_AVE, var = "majorCluster")
tmp_CD8_Ad_AVE[is.na(tmp_CD8_Ad_AVE)] <- 0

overlapped_names <- intersect(tissues_colors %>% names(), row.names(tmp_CD8_Ad_AVE))
tmp_CD8_Ad_AVE <- tmp_CD8_Ad_AVE[overlapped_names, overlapped_names]

trans_score <- obj@cluster.data %>% filter(aid == "HCA2.0") %>% select(tran)
trans_score <- trans_score[overlapped_names,, drop = F]

migration_across_tissue <- tmp_CD8_Ad_AVE 
row.names(migration_across_tissue) <- colnames(migration_across_tissue)
pheatmap::pheatmap(migration_across_tissue)

uppertri <- migration_across_tissue
uppertri[!upper.tri(uppertri)] <- 10
uppertri <- as.matrix(uppertri)
tmp <- as.vector(t(uppertri))
weight <- tmp[!((tmp == 10)|(tmp == 0)) ]

###-----------------------------------Adult CD4-----------------------------------###
TCR$loc <- 1 ### set loc equal to 1, and using transition index instead of migration
obj <- Startrac.run(TCR %>% filter(Stage == "Adult"), proj = "HCA2.0", verbose = T)
Startrac::plot(obj, index.type = "cluster.all", byPatient = T)

dat_for_plot_CD4_ad_migration_ave <- obj@pIndex.tran %>% filter(aid == "HCA2.0")

tmp_CD4_Ad_AVE <- dat_for_plot_CD4_ad_migration_ave %>% `[`(, -c(1))
tmp_CD4_Ad_AVE <- column_to_rownames(tmp_CD4_Ad_AVE, var = "majorCluster")
tmp_CD4_Ad_AVE[is.na(tmp_CD4_Ad_AVE)] <- 0

overlapped_names <- intersect(tissues_colors %>% names(), row.names(tmp_CD4_Ad_AVE))
tmp_CD4_Ad_AVE <- tmp_CD4_Ad_AVE[overlapped_names, overlapped_names]

trans_score <- obj@cluster.data %>% filter(aid == "HCA2.0") %>% select(tran)
trans_score <- trans_score[overlapped_names,, drop = F]

migration_across_tissue <- tmp_CD4_Ad_AVE 
row.names(migration_across_tissue) <- colnames(migration_across_tissue)
pheatmap::pheatmap(migration_across_tissue)

uppertri <- migration_across_tissue
uppertri[!upper.tri(uppertri)] <- 10
uppertri <- as.matrix(uppertri)
tmp <- as.vector(t(uppertri))
weight <- tmp[!((tmp == 10)|(tmp == 0)) ]

###---------------------------------igraph plot： plot the CD4 and CD8 individually------------------------------###
shared_matrix <- unclass(with(TCR %>% filter(Stage == "Adult"), table(clone.id, majorCluster)))
shared_matrix <- shared_matrix %>% unclass() %>% as.data.frame()
rownames(shared_matrix) <- 1:dim(shared_matrix)[1]
shared_matrix[shared_matrix > 200] <- 200
# pheatmap(shared_matrix)
shared_matrix <- shared_matrix[, overlapped_names]

# write.table(shared_matrix, "clone_type_numbers_in_each_tissue_of_CD4.txt", sep = "\t", col.names = T, row.names = T, quote = F)
shared_matrix[ shared_matrix > 0 ] <- 1
shared_m <- t(shared_matrix) %*% (shared_matrix %>% as.matrix())
cluster.gr <- igraph::graph_from_adjacency_matrix(shared_m/sum(shared_m), 
																									mode = "undirected", weighted = TRUE, diag = FALSE)

# E(cluster.gr)$weight <- weight*5
E(cluster.gr)$weight <- weight*10
E(cluster.gr)$width <- E(cluster.gr)$weight

# V(cluster.gr)$size <- trans_score$tran*5 #Fetal*400
V(cluster.gr)$size <- trans_score$tran*8
V(cluster.gr)$frame.color <- tissues_colors[trans_score %>% row.names()]
V(cluster.gr)$color <- tissues_colors[trans_score %>% row.names()]

# png("Tissue_TCR_shared_network_CD4.png", height = 10, width = 10, res = 400, units = "in")
pdf("Tissue_TCR_shared_network_CD4_adult.pdf", height = 10, width = 10)
# pdf("Cluster_TCR_shared_network_CD4_Fetal.pdf", height = 10, width = 10)
set.seed(2)
pt <- base::plot(cluster.gr,
								 # edge.width = igraph::E(cluster.gr)$weight, #Fetal*100
								 # vertex.color = unname(tissues_colors)[match(x = names(V(cluster.gr)), table = names(tissues_colors))],
								 vertex.color = tissues_colors[trans_score %>% row.names()],
								 # vertex.color = color_used[1:14],
								 vertex.label.dist = 3,
								 # vertex.label.color = unname(tissues_colors)[match(x = names(V(cluster.gr)), table = names(tissues_colors))],
								 vertex.label.color = tissues_colors[trans_score %>% row.names()],
								 vertex.label.cex	= 1,
								 layout = layout_in_circle(cluster.gr))
legend( title = "TCR_sharing*10", list(x = -0.2, y = -1.1), legend = c(0.0005, 0.1, 0.214),
				col = "grey", lty = 1, lwd = c(0.005, 1, 2.144))
legend( title = "Migration_ability*100", list(x = 0.2, y = -1.1),
				legend = c(0.19, 0.56, 0.92),# V(g)$degree。
				col = "black", pch = 1, pt.lwd = 1, pt.cex = c(0.75, 2.23, 3.72) 
				# 设置相同数值，igraph与legend生成的节点大小不一样，
				# 图例节点差不多是图节点的2倍大，这里通过设置图节点大小为degree*2解决，需要手动调整。
)
dev.off()












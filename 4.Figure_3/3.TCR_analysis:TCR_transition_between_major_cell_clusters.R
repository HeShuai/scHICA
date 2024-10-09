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
###--------------------TCR transition between different cell clusters--------------------### Analyzing CD8 T cells as an example
TCR <- CD8 %>%
	dplyr::select(c(barcode, customer_clone, renamed_Samlename, sample, Detail_Annotation_final, Stage)) %>%
	# filter(Detail_Annotation_final %>% grepl(., pattern = "TRM")) %>% filter(! (renamed_Samlename %>% grepl(., pattern = "Blood|Kidney|Marrow"))) %>%
	# filter(Detail_Annotation_final %>% grepl(., pattern = "TEM|TEFF")) %>% filter(Stage == "Adult") %>%  
	setnames(old = c("barcode", "customer_clone", "renamed_Samlename", "sample", "Detail_Annotation_final"), new = c("Cell_Name", "clone.id", "loc", "patient", "majorCluster")) %>% unique()


obj <- Startrac.run(TCR, proj = "HCA2.0", verbose = T)
Startrac::plot(obj, index.type = "cluster.all", byPatient = T)

dat_for_plot_CD8 <- obj@cluster.data %>% filter(aid != "HCA2.0")
dat_for_plot_CD8$Stage <- dat_for_plot_CD8$aid %>% gsub(pattern = "[1-3]", replacement = "")
dat_for_plot_CD8$CD4_CD8 <- "CD8"

dat_for_plot_CD8 <- obj@pIndex.tran %>% filter(aid != "HCA2.0")

tmp1 <- dat_for_plot_CD8 %>% filter(aid == "Adult3") %>% `[`(, -c(1))
tmp1 <- column_to_rownames(tmp1, var = "majorCluster")
tmp1[is.na(tmp1)] <- 0

pdf("pheatmap_transition_cluster_Adult3.pdf", height = 6, width = 6)
pheatmap(tmp1,
				 color = colorRampPalette(c(col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4'))(1000),
				 treeheight_col = 10, treeheight_row = 10, border_color = "white")
dev.off()

packages <- list("Seurat", "Matrix", "stringr", "stringi", "gridExtra", "ggplot2",
								 "plyr", "ggthemes", "cowplot", "data.table", "RColorBrewer", "ComplexHeatmap",
								 "circlize", "pheatmap", "reshape2", "scales", "rlang", "future",
								 "parallel", "ggrepel", "ggsci", "ape", "dplyr", "harmony",
								 "sigminer", "ggpubr", "rstatix", "ROGUE", "tidyverse", "reticulate", "dendextend",
								 "scRepertoire", "ggdendro", "ggtern", "viridis", "scatterpie", "SeuratData", "SeuratDisk",
								 "scibet", "scuttle", "zellkonverter", "cowplot")

lapply(packages, library, character.only = TRUE)
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]
options(width = 240)

options(future.globals.maxSize = 400*1000 * 1024^2)
plan("multisession", workers = 1)
plan()

tissues_colors <- c("Ascending colon" = "#E64B35FF", "Sigmoid colon" = '#ff0000ff', "Transverse colon" = '#d62790',
										Ileum = "#91D1C2FF", Jejunum = '#4DBBD5FF', Duodenum = '#00A087FF',
										"Lesser gastric curvature" = '#F0E685FF', "Mesenteric lymph node" = "#FFB570FF",
										Blood = "#1B1919FF", PBMC = "#374E55FF", "Whole blood" = "#3F4041FF",	Marrow = '#466983FF', Spleen = "#725663FF",
										Lung = "#6A6599FF", Liver = '#749B58FF', Kidney = '#B09C85FF', Skin = '#D6D6CEFF', Bladder = "#D49464FF",
										Esophagus = '#c5b0d5', Trachea = "#BA63affF", Heart = '#9467bd', Muscle = "#aec7e8", Pancreas = "#5050FFFF",
										Stomach = '#ADB17DFF', Thymus = '#8c564b')

#####-----------------------------------------TCR tracking analysis tissues--------------------------------###
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]
options(width = 180)
###---------------------------------------TCR clone size (bar plot) in each cell cluster, Fetal--------------------------###
UMAP <- read.table("HCA2.0_TNK_UMAP.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
TCR_meta.data <- read.table("20230404_TNK_meta.data_with_TCR_information_.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
TCR_meta.data$renamed_Samlename <- TCR_meta.data$Name
TCR_meta.data$renamed_Samlename[TCR_meta.data$renamed_Samlename %>% grepl(pattern = "^Colon$")] <- "Jejunum"
TCR_meta.data$renamed_Samlename[TCR_meta.data$renamed_Samlename %>% grepl(pattern = "Stomach Protein")] <- "Stomach"
TCR_meta.data$renamed_Samlename[TCR_meta.data$renamed_Samlename %>% grepl(pattern = ".*blood|.*PBMC.*")] <- "Blood"

TCR_meta.data <- TCR_meta.data %>% filter(! Detail_Annotation_final %in% c("DG T1", "DG T2", "DG T3", "ILCs",
																																					 "MAIT", "NK CD16 KIR2DL3","NK CD16 MYOM2",
																																					 "NK CD56 CCL3", "NK CD56 KIR2DL4", "Pro ILC", "Pro NK"))
###------------choose CD4 and CD8 cell separately----------###
T_cell_colon_unique_barcode <- TCR_meta.data %>% select(row.names, customer_clone, Name:Stage, renamed_Samlename, Detail_Annotation_final) %>% unique
CD4 <- T_cell_colon_unique_barcode %>% filter(grepl(Detail_Annotation_final, pattern = "CD4|Treg")) %>% filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8")) 
CD8 <- T_cell_colon_unique_barcode %>% filter(grepl(Detail_Annotation_final, pattern = "CD8")) %>% filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8"))
Others <- T_cell_colon_unique_barcode %>% filter(Detail_Annotation_final %in% c("DP RAG1", "DP TOX2", "High pro DP", "CD4CD8 low"))

###------------CD4 CD8 clone state calculation separately------------------###
clone_ID_CD4 <- CD4$customer_clone %>% table %>% as.data.frame() %>% filter(`.` != "Not")
CD4$Clone_Size <- mapvalues(CD4$customer_clone, from = clone_ID_CD4$".", to = clone_ID_CD4$Freq, warn_missing = F) %>% as.numeric()
clone_ID_CD8 <- CD8$customer_clone %>% table %>% as.data.frame() %>% filter(`.` != "Not")
CD8$Clone_Size <- mapvalues(CD8$customer_clone, from = clone_ID_CD8$".", to = clone_ID_CD8$Freq, warn_missing = F) %>% as.numeric()
clone_ID_Others <- Others$customer_clone %>% table %>% as.data.frame() %>% filter(`.` != "Not")
Others$Clone_Size <- mapvalues(Others$customer_clone, from = clone_ID_Others$".", to = clone_ID_Others$Freq, warn_missing = F) %>% as.numeric()

TCR_meta.data <- rbind(CD4, CD8, Others) %>% as.data.frame()
TCR_meta.data <- TCR_meta.data %>% filter(!is.na(Clone_Size))

TCR_meta.data$UMAP_1 <- mapvalues(TCR_meta.data$row.names, from = UMAP %>% row.names(), to = UMAP$UMAP_1, warn_missing = F) %>% as.numeric()
TCR_meta.data$UMAP_2 <- mapvalues(TCR_meta.data$row.names, from = UMAP %>% row.names(), to = UMAP$UMAP_2, warn_missing = F) %>% as.numeric()
TCR_meta.data$Clone_strategy <- with(TCR_meta.data, ifelse(Clone_Size == 1, "Single",
																													 ifelse(Clone_Size > 1 & Clone_Size <= 5, "2-5",
																													 			 ifelse(Clone_Size > 5 & Clone_Size <= 10, "6-10",
																													 			 			 ifelse(Clone_Size > 10 & Clone_Size <= 50, "11-50",
																													 			 			 			 ifelse(Clone_Size > 50 & Clone_Size <= 100, "51-100",
																													 			 			 			 			 ifelse(Clone_Size > 100 & Clone_Size <= 250, "101-250",
																													 			 			 			 			 			 ifelse(Clone_Size > 250 & Clone_Size <= 500, "251-500",
																													 			 			 			 			 			 			 ifelse(Clone_Size > 500 & Clone_Size <= 1000, "501-1000",
																													 			 			 			 			 			 			 			 ifelse(Clone_Size > 1000, "1001+", "NON"))))))))))

Clone_state <- c("Single", "2-5", "6-10", "11-50", "51-100", "101-250", "251-500", "501-1000", "1001+")
TCR_meta.data$Clone_strategy <- factor(TCR_meta.data$Clone_strategy %>% as.character(), levels = rev(Clone_state))

TCR_color <- c("#019875", "#D2C49F", "#FEC6A2", "#FE9C7F", "#FE715C", "#FE4739", "#F70206", "#CC112B", "#96255B")

###-------------------------------CD4 CD8 clone state in fetal subclusters----------------------------------###
Fetal_TCR <- TCR_meta.data %>% filter(Stage == "Fetal")

dats <- table(Fetal_TCR$Clone_strategy, Fetal_TCR$Detail_Annotation_final) %>% t %>% unclass() %>% as.data.frame() %>% rownames_to_column(var = "Cell_type") 
order_Fetal_celltype <- dats %>% mutate(multiple = c(`1001+` + `501-1000` + `251-500` + `101-250` + `51-100` + `11-50` + `6-10` + `2-5`),
																				multiple_proportion = multiple/Single, All = multiple + Single) %>% filter(All >= 50) %>%
	arrange(desc(multiple_proportion)) %>% select(Cell_type) %>% unlist %>% unname()
dats <- melt(dats) %>% filter(variable != "NON") %>% filter(value != 0)

dats$Cell_type <- factor(dats$Cell_type, levels = c(order_Fetal_celltype %>% grep(pattern = "CD4|Treg", value = T),
																										order_Fetal_celltype %>% grep(pattern = "CD8", value = T),
																										order_Fetal_celltype %>% grep(pattern = "DP RAG1|DP TOX2|High pro DP", value = T)))

pdf("Fetal_TCR_clone_structrue_each_cluster_CD4_CD8_separate.pdf", width = 6, height = 4)
ggplot(dats, aes(x = Cell_type, y = value, fill = variable, color = variable)) + 
	geom_bar(stat = "identity", position = "fill") +
	scale_fill_manual(values = rev(TCR_color[1:4]), name = "Clone size") +
	scale_color_manual(values = rep("white", 4)) +
	scale_y_continuous(expand = c(0, 0.01)) +
	theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
dev.off()
###------------Fetal_TCR_clone_structrue in each cluster RE/EO analysis--------------------###
Fetal_TCR_meta <- table(Fetal_TCR$Clone_strategy, Fetal_TCR$Detail_Annotation_final) %>% t %>% unclass %>% as.data.frame()
Fetal_TCR_meta <- Fetal_TCR_meta %>% mutate(multiple = `1001+` + `501-1000` + `251-500` + `101-250` + `51-100` + `11-50` + `6-10` + `2-5`) %>% select(Single, multiple) %>% filter(Single >= 50)
Fetal_TCR_meta_CD4 <- Fetal_TCR_meta[Fetal_TCR_meta %>% row.names() %>% grep(pattern = "CD4|Treg"), ]
Fetal_TCR_meta_CD8 <- Fetal_TCR_meta[Fetal_TCR_meta %>% row.names() %>% grep(pattern = "CD8"), ]
Fetal_TCR_meta_DP <- Fetal_TCR_meta[Fetal_TCR_meta %>% row.names() %>% grep(pattern = "DP"), ]

R_O_E_CD4 <- CrossTable(Fetal_TCR_meta_CD4 %>% as.matrix(), expected = T)
Multiple_CD4 <- (R_O_E_CD4$CST$observed/R_O_E_CD4$CST$expected)[, 2]
R_O_E_CD8 <- CrossTable(Fetal_TCR_meta_CD8 %>% as.matrix(), expected = T)
Multiple_CD8 <- (R_O_E_CD8$CST$observed/R_O_E_CD8$CST$expected)[, 2]
rbind_RERO_TCR <- c(Multiple_CD4, Multiple_CD8) %>% as.matrix() %>% t
R_O_E_DP <- CrossTable(Fetal_TCR_meta_DP %>% as.matrix(), expected = T)
Multiple_DP <- (R_O_E_DP$CST$observed/R_O_E_DP$CST$expected)[, 2]

rbind_RERO_TCR <- rbind_RERO_TCR %>% as.data.frame() %>% t 
rbind_RERO_TCR <- rbind_RERO_TCR[c(order_Fetal_celltype %>% grep(pattern = "CD4|Treg", value = T),
																	 order_Fetal_celltype %>% grep(pattern = "CD8", value = T)), ]
numbers <- round(rbind_RERO_TCR %>% unname(), 2)

pheatmap::pheatmap(rbind_RERO_TCR %>% t,
									 cluster_rows = F, cluster_cols = F,
									 color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
									 display_numbers = round(rbind_RERO_TCR %>% unname(), 2) %>% as.matrix() %>% t)
###---------------------------------------TCR clone size (bar plot) in each cell cluster, Adult--------------------------###
Adult_TCR <- TCR_meta.data %>% filter(Stage == "Adult")
dats <- table(Adult_TCR$Clone_strategy, Adult_TCR$Detail_Annotation_final) %>% t %>% unclass() %>% as.data.frame() %>% rownames_to_column(var = "Cell_type") 
order_Adult_celltype <- dats %>% mutate(multiple = c(`1001+` + `501-1000` + `251-500` + `101-250` + `51-100` + `11-50` + `6-10` + `2-5`),
																				multiple_proportion = multiple/Single, All = multiple + Single) %>% filter(All >= 80) %>%
	arrange(desc(multiple_proportion)) %>% select(Cell_type) %>% unlist %>% unname()
dats <- melt(dats) %>% filter(variable != "NON") %>% filter(value != 0)
dats$Cell_type <- factor(dats$Cell_type, levels = c(order_Adult_celltype %>% grep(pattern = "CD4|Treg", value = T),
																										order_Adult_celltype %>% grep(pattern = "CD8", value = T),
																										order_Adult_celltype %>% grep(pattern = "DP RAG1|DP TOX2|High pro DP", value = T)))
dats <- dats %>% filter(!is.na(Cell_type))

pdf("Adult_TCR_clone_structrue_each_cluster_CD4_CD8_separate.pdf", width = 6, height = 4)
ggplot(dats, aes(x = Cell_type, y = value, fill = variable, color = variable)) + 
	geom_bar(stat = "identity", position = "fill") +
	scale_fill_manual(values = rev(TCR_color), name = "Clone size") +
	scale_color_manual(values = rep("white", 9)) +
	scale_y_continuous(expand = c(0, 0.01)) +
	theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
dev.off()

Adult_TCR_meta <- table(Adult_TCR$Clone_strategy, Adult_TCR$Detail_Annotation_final) %>% t %>% unclass %>% as.data.frame()
Adult_TCR_meta <- Adult_TCR_meta %>% mutate(multiple = `1001+` + `501-1000` + `251-500` + `101-250` + `51-100` + `11-50` + `6-10` + `2-5`, All = multiple + Single) %>% filter(All >= 80) %>% select(Single, multiple) 
Adult_TCR_meta_CD4 <- Adult_TCR_meta[Adult_TCR_meta %>% row.names() %>% grep(pattern = "CD4|Treg"), ]
Adult_TCR_meta_CD8 <- Adult_TCR_meta[Adult_TCR_meta %>% row.names() %>% grep(pattern = "CD8"), ]
Adult_TCR_meta_DP <- Adult_TCR_meta[Adult_TCR_meta %>% row.names() %>% grep(pattern = "DP"), ]

R_O_E_CD4 <- CrossTable(Adult_TCR_meta_CD4 %>% as.matrix(), expected = T)
Multiple_CD4 <- (R_O_E_CD4$CST$observed/R_O_E_CD4$CST$expected)[, 2]
R_O_E_CD8 <- CrossTable(Adult_TCR_meta_CD8 %>% as.matrix(), expected = T)
Multiple_CD8 <- (R_O_E_CD8$CST$observed/R_O_E_CD8$CST$expected)[, 2]
rbind_RERO_TCR <- c(Multiple_CD4, Multiple_CD8) %>% as.matrix() %>% t
R_O_E_DP <- CrossTable(Adult_TCR_meta_DP %>% as.matrix(), expected = T)
Multiple_DP <- (R_O_E_DP$CST$observed/R_O_E_DP$CST$expected)[, 2]

rbind_RERO_TCR <- rbind_RERO_TCR %>% as.data.frame() %>% t 
rbind_RERO_TCR <- rbind_RERO_TCR[c(order_Adult_celltype %>% grep(pattern = "CD4|Treg", value = T),
																	 order_Adult_celltype %>% grep(pattern = "CD8", value = T)), ]
numbers <- round(rbind_RERO_TCR %>% unname(), 2)

pheatmap::pheatmap(rbind_RERO_TCR %>% t,
									 cluster_rows = F, cluster_cols = F,
									 color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
									 display_numbers = round(rbind_RERO_TCR %>% unname(), 2) %>% as.matrix() %>% t)

###------------------pairwise comparison---------------###
ovlerlap_celltype <- intersect(Fetal_TCR_meta %>% row.names(), Adult_TCR_meta %>% row.names())
ovlerlap_celltype_Fetal <- Fetal_TCR_meta[ovlerlap_celltype, ]
ovlerlap_celltype_Fetal <- rownames_to_column(ovlerlap_celltype_Fetal, var = "Cell_type")
ovlerlap_celltype_Fetal$Stage <- "Fetal"
ovlerlap_celltype_Fetal <- melt(ovlerlap_celltype_Fetal)

ovlerlap_celltype_Adult <- Adult_TCR_meta[ovlerlap_celltype, ]
ovlerlap_celltype_Adult <- rownames_to_column(ovlerlap_celltype_Adult, var = "Cell_type")
ovlerlap_celltype_Adult$Stage <- "Adult"
ovlerlap_celltype_Adult <- melt(ovlerlap_celltype_Adult)

rbind_overlap <- rbind(ovlerlap_celltype_Fetal, ovlerlap_celltype_Adult)
rbind_overlap$Stage <- factor(rbind_overlap$Stage %>% as.character(), levels = c("Fetal", "Adult"))

dat_pvalue <- c()
for(i in rbind_overlap$Cell_type %>% unique()){
	tmp <- rbind_overlap %>% filter(Cell_type %in% i) %>% select(2:4) %>% dcast(Stage ~  variable) %>% select(2:3) %>% chisq.test()
	pvalue <- tmp$p.value
	names(pvalue) <- i
	dat_pvalue <- c(dat_pvalue, pvalue)
}
dat_pvalue <- dat_pvalue %>% as.data.frame() %>% arrange(desc(`.`))

rbind_overlap$Cell_type <- factor(rbind_overlap$Cell_type %>% as.character(), levels = dat_pvalue %>% row.names())

pdf("pairwised_comparison_TCR_clone_CD4_CD8_separate.pdf", width = 6, height = 4)
ggplot(data = rbind_overlap, aes(x = Stage, y = value, group = Stage, fill = variable)) +
	geom_bar(stat = "identity", position = "fill") +
	scale_y_continuous(expand = c(0, 0.01)) +
	scale_fill_manual(values = alpha(c("#019875", "firebrick3"), 0.6)) +
	facet_grid( ~ Cell_type) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
				panel.spacing.x = unit(0, "lines"), 
				panel.border = element_rect(color = "black", fill = NA, size = 1),
				strip.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0))
dev.off()
###---------------------cell cycle analysis using addmodule score-------------------------####
load("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/TNK/Cells/2023.4.4.TNK_cells_final.465018_annotation.RData")

DefaultAssay(subset_cells) <- "RNA"
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features = VariableFeatures(subset_cells)[1])

gmt_file <- read.gmt("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/TNK/Cells/TCR/HALLMARK_G2M_CHECKPOINT.v2023.1.Hs.gmt")
subset_cells <- AddModuleScore(subset_cells, features = list(gmt_file$gene))
write.table(subset_cells@meta.data %>% select(Name:Stage, Detail_Annotation_final, Cluster1), file = "proliferation.TNK_meta.data.txt", col.names = T, row.names = T, sep = "\t", quote = F )

proliferation <- read.table("proliferation.TNK_meta.data.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
CD4_TCR <- TCR_meta.data %>% filter(Detail_Annotation_final %>% grepl(pattern = "CD4|Treg")) %>% filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8"))
CD4_TCR$Cycle <- mapvalues(CD4_TCR$row.names, from = proliferation %>% row.names(), to = proliferation$Cluster1, warn_missing = F)
CD8_TCR <- TCR_meta.data %>% filter(Detail_Annotation_final %>% grepl(pattern = "CD8")) %>% filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8"))
CD8_TCR$Cycle <- mapvalues(CD8_TCR$row.names, from = proliferation %>% row.names(), to = proliferation$Cluster1, warn_missing = F)

CD4_cycle <- list(Single = CD4_TCR %>% filter(Stage == "Fetal") %>% filter(Clone_strategy == "Single") %>% select(Cycle) %>% unlist() %>% as.numeric(),
									Multiple = CD4_TCR %>% filter(Stage == "Fetal")  %>% filter(Clone_strategy != "Single") %>% select(Cycle) %>% unlist() %>% as.numeric())

CD8_cycle <- list(Single = CD8_TCR %>% filter(Stage == "Fetal")  %>% filter(Clone_strategy == "Single") %>% select(Cycle) %>% unlist() %>% as.numeric(),
									Multiple = CD8_TCR %>% filter(Stage == "Fetal")  %>% filter(Clone_strategy != "Single") %>% select(Cycle) %>% unlist() %>% as.numeric())

t.test(CD4_cycle$Single, CD4_cycle$Multiple)
t.test(CD8_cycle$Single, CD8_cycle$Multiple)

rbind_cycle <- rbind(CD4_TCR %>% select(Clone_strategy, Cycle, Stage) %>%
										 	mutate(Cell_type = "CD4", Multiple = ifelse(Clone_strategy != "Single", "Multiple", "Single"), Cluster = paste0(Stage, "_", Cell_type)),
										 CD8_TCR %>% select(Clone_strategy, Cycle, Stage) %>%
										 	mutate(Cell_type = "CD8", Multiple = ifelse(Clone_strategy != "Single", "Multiple", "Single"), Cluster = paste0(Stage, "_", Cell_type)))
rbind_cycle$Multiple <- factor(rbind_cycle$Multiple %>% as.character(), levels = c("Single", "Multiple"))
rbind_cycle$Cluster <- factor(rbind_cycle$Cluster %>% as.character(), levels = c("Fetal_CD4", "Fetal_CD8", "Adult_CD4", "Adult_CD8"))

pdf("cycle_module_score_CD4_CD8_single_multiple_clones.pdf", width = 6, height = 4)
ggplot(data = rbind_cycle, aes(x = Multiple, y = Cycle %>% as.numeric())) +
	geom_violin(aes(fill = factor(Multiple))) + 
	facet_grid( ~ Cluster) +
	scale_fill_manual(values = alpha(c("#019875", "firebrick3"), 0.6)) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
				panel.spacing.x = unit(0, "lines"), 
				panel.border = element_rect(color = "black", fill = NA, size = 1),
				strip.text.x.top = element_text(angle = 0, hjust = 0.5, vjust = 0)) +
	labs(x = "Module score of cycle")
dev.off()



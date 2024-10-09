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

TCR_meta.data <- read.table("Merged_TCR.information_at_least_one_paired_VDJ_cellranger6.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
TCR_meta.data$renamed_Samlename <- TCR_meta.data$Name
TCR_meta.data$renamed_Samlename[TCR_meta.data$renamed_Samlename %>% grepl(pattern = "^Colon$")] <- "Jejunum"
TCR_meta.data$renamed_Samlename[TCR_meta.data$renamed_Samlename %>% grepl(pattern = "Stomach Protein")] <- "Stomach"
TCR_meta.data$renamed_Samlename[TCR_meta.data$renamed_Samlename %>% grepl(pattern = ".*blood|.*PBMC.*")] <- "Blood"


TCR_meta.data <- TCR_meta.data %>% filter(! Detail_Annotation_final %in% c("DG T1", "DG T2", "DG T3", "ILCs",
									"MAIT", "NK CD16 KIR2DL3","NK CD16 MYOM2", "NK CD56 CCL3", "NK CD56 KIR2DL4", "Pro ILC", "Pro NK", "DP RAG1", "DP TOX2", "High pro DP", "CD4CD8 low"))

###--------------------TCR composition of each tissue--------------------###
T_cell_colon_unique_barcode <- TCR_meta.data %>% select(barcode, customer_clone, sample, renamed_Samlename, Detail_Annotation_final) %>% unique
CD4 <- T_cell_colon_unique_barcode %>% filter(grepl(Detail_Annotation_final, pattern = "CD4|Treg")) %>% filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8")) 
CD8 <- T_cell_colon_unique_barcode %>% filter(grepl(Detail_Annotation_final, pattern = "CD8")) %>% filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8"))

dat <- data.frame()
for(i in CD4$sample %>% unique()){
	tmp <- CD4 %>% filter(sample == i)
	clone_cellnumber <- tmp$customer_clone %>% table() %>% sort(decreasing = T) %>% as.data.frame()
	clone_cellnumber <- rownames_to_column(clone_cellnumber, var = "clone_rank")
	names(clone_cellnumber)[2] <- "Clone_type"
	dat <- rbind(dat, clone_cellnumber)
}
CD4$clone_rank <- mapvalues(CD4$customer_clone, from = dat$Clone_type, to = dat$clone_rank, warn_missing = F)
CD4$clone_rank <- CD4$clone_rank %>% as.numeric()

dat1 <- data.frame()
for(i in CD8$sample %>% unique()){
	tmp <- CD8 %>% filter(sample == i)
	clone_cellnumber <- tmp$customer_clone %>% table() %>% sort(decreasing = T) %>% as.data.frame()
	clone_cellnumber <- rownames_to_column(clone_cellnumber, var = "clone_rank")
	names(clone_cellnumber)[2] <- "Clone_type"
	dat1 <- rbind(dat1, clone_cellnumber)
}
CD8$clone_rank <- mapvalues(CD8$customer_clone, from = dat1$Clone_type, to = dat1$clone_rank, warn_missing = F)
CD8$clone_rank <- CD8$clone_rank %>% as.numeric()

dat2 <- data.frame()
for(i in CD8$sample %>% unique()){
	tmp <- CD8 %>% filter(sample == i)
	tmp1 <- tmp %>% group_by(renamed_Samlename) %>% summarise_at(.vars = "clone_rank", .funs = function(x)(return(cut(x, c(0, 10, 100, 1000, 100000)) %>% as.data.frame())))
	tmp1$clone_rank <- tmp1$clone_rank$. %>% as.character()
	tmp2 <- table(tmp1) %>% unclass() %>% prop.table(margin = 1) %>% unclass() %>% as.data.frame()
	colnames(tmp2) <- c("one_ten", "Ten_onehundred", "oneh_thousand", "thousand_more")
	tmp2 <- rownames_to_column(tmp2, "tissues")
	tmp2 <- melt(tmp2)
	tmp2$Donor <- i
	dat2 <- rbind(dat2, tmp2)
}

dat3 <- data.frame()
for(i in CD4$sample %>% unique()){
	tmp <- CD4 %>% filter(sample == i)
	tmp1 <- tmp %>% group_by(renamed_Samlename) %>% summarise_at(.vars = "clone_rank", .funs = function(x)(return(cut(x, c(0, 10, 100, 1000, 100000)) %>% as.data.frame())))
	tmp1$clone_rank <- tmp1$clone_rank$. %>% as.character()
	tmp2 <- table(tmp1) %>% unclass() %>% prop.table(margin = 1) %>% unclass() %>% as.data.frame()
	colnames(tmp2) <- c("one_ten", "Ten_onehundred", "oneh_thousand", "thousand_more")
	tmp2 <- rownames_to_column(tmp2, "tissues")
	tmp2 <- melt(tmp2)
	tmp2$Donor <- i
	dat3 <- rbind(dat3, tmp2)
}

pdf("pie_char_TCR_CD8.pdf", height = 8, width = 18)
ggplot(dat2, aes(x = "", y = value, fill = variable)) +
	geom_bar(width = 1, stat = "identity", color = "white") +  coord_polar("y") +
	facet_grid( Donor ~ tissues) +
	scale_fill_manual(values = color_used) +
	theme(panel.spacing = unit(0, "lines")) +
	theme_classic2()
dev.off()

pdf("pie_char_TCR_CD4.pdf", height = 8, width = 18)
ggplot(dat3, aes(x = "", y = value, fill = variable)) +
	geom_bar(width = 1, stat = "identity", color = "white") +  coord_polar("y") +
	facet_grid( Donor ~ tissues) +
	scale_fill_manual(values = color_used) +
	theme(panel.spacing = unit(0, "lines")) +
	theme_classic2()
dev.off()

###--------------------TCR expansion of each tissue: bar plot--------------------###
TCR <- CD8 %>%
	dplyr::select(c(barcode, customer_clone, renamed_Samlename, sample, Detail_Annotation_final)) %>%
	setnames(old = c("barcode", "customer_clone", "renamed_Samlename", "sample", "Detail_Annotation_final"), new = c("Cell_Name", "clone.id", "majorCluster", "patient", "loc")) %>% unique()

dat_CD8 <- data.frame()
for(i in TCR$patient %>% unique()){
	tmp <- TCR %>% filter(patient == i)
	obj <- new("Startrac", tmp, aid = "HCA2.0")
	obj <- calIndex(obj)
	tmp1 <- obj@cluster.data[2:5] %>% melt()
	tmp1$Donor <- i
	dat_CD8 <- rbind(dat_CD8, tmp1)
}

TCR <- CD4 %>%
	dplyr::select(c(barcode, customer_clone, renamed_Samlename, sample, Detail_Annotation_final)) %>%
	setnames(old = c("barcode", "customer_clone", "renamed_Samlename", "sample", "Detail_Annotation_final"), new = c("Cell_Name", "clone.id", "majorCluster", "patient", "loc")) %>% unique()

dat_CD4 <- data.frame()
for(i in TCR$patient %>% unique()){
	tmp <- TCR %>% filter(patient == i)
	obj <- new("Startrac", tmp, aid = "HCA2.0")
	obj <- calIndex(obj)
	tmp1 <- obj@cluster.data[2:5] %>% melt()
	tmp1$Donor <- i
	dat_CD4 <- rbind(dat_CD4, tmp1)
}

dat_CD4$CD4_CD8 <- "CD4"
dat_CD8$CD4_CD8 <- "CD8"

CD4_CD8 <- rbind(dat_CD4, dat_CD8)
CD4_CD8$Stage <- CD4_CD8$Donor %>% gsub(pattern = "[1-3]", replacement = "")
CD4_CD8$Stage <- factor(CD4_CD8$Stage %>% as.character(), levels = c("Fetal", "Adult"))

pdf("T_cell_expansion.pdf", height = 6, width = 15)
ggplot(CD4_CD8 %>% filter(variable == "expa"), aes(x = Stage, y = value, fill = Stage)) + 
	geom_boxplot() +
	geom_point() +
	scale_fill_manual(values = alpha(color_used[c(6, 7)], 0.9)) +
	scale_y_continuous(expand = c(0, 0.01)) +
	facet_grid(CD4_CD8 ~ majorCluster, space = "free_y", scales = "free_y") +
	theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
				panel.spacing.x = unit(0, "lines"), 
				panel.border = element_rect(color = "black", fill = NA, size = 1),
				strip.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0))
dev.off()

###---------------------------------------single or multiple TCR clone in different tissues: bar plot--------------------------###
T_cell_colon_unique_barcode <- TCR_meta.data %>% select(barcode, customer_clone, sample, renamed_Samlename, Detail_Annotation_final, Stage) %>% unique
CD4 <- T_cell_colon_unique_barcode %>% filter(grepl(Detail_Annotation_final, pattern = "CD4|Treg")) %>% filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8")) 
CD8 <- T_cell_colon_unique_barcode %>% filter(grepl(Detail_Annotation_final, pattern = "CD8")) %>% filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8"))

result <- function(Fetal_Adult, CD4_CD8){res <- lapply(CD4_CD8$renamed_Samlename %>% unique(), FUN = function(x){
	tmp <- CD4_CD8 %>% filter(Stage == Fetal_Adult) %>% filter(renamed_Samlename == x)
	tmp1 <- tmp$customer_clone %>% table %>% unclass %>% as.numeric()
	tmp2 <- data.frame(mutiple_clone_proportion = 1 - sum(tmp1 == 1)/sum(tmp1),
										 Single_clone_proportion = sum(tmp1 == 1)/sum(tmp1),
										 Tissue = x)
	return(tmp2)
})
return(res)}

CD4_clone_adult <- do.call(rbind, result(Fetal_Adult = "Adult", CD4_CD8 = CD4))
CD4_clone_adult$Stage <- "CD4_clone_Adult"
CD8_clone_adult <- do.call(rbind, result(Fetal_Adult = "Adult", CD4_CD8 = CD8))
CD8_clone_adult$Stage <- "CD8_clone_Adult"
CD4_clone_Fetal <- do.call(rbind, result(Fetal_Adult = "Fetal", CD4_CD8 = CD4))
CD4_clone_Fetal$Stage <- "CD4_clone_Fetal"
CD8_clone_Fetal <- do.call(rbind, result(Fetal_Adult = "Fetal", CD4_CD8 = CD8))
CD8_clone_Fetal$Stage <- "CD8_clone_Fetal"

rbind_TCR_clone_state <- reduce(list(CD4_clone_adult, CD8_clone_adult, CD4_clone_Fetal, CD8_clone_Fetal), rbind)
rbind_TCR_clone_state <- rbind_TCR_clone_state %>% mutate(CD4_CD8 = str_split(Stage, pattern = "_", simplify = T) %>% `[`(, 1),
				Stage = str_split(Stage, pattern = "_", simplify = T) %>% `[`(, 3))
rbind_TCR_clone_state <- rbind_TCR_clone_state %>% melt()
rbind_TCR_clone_state$Stage <- factor(rbind_TCR_clone_state$Stage %>% as.character(), levels = c("Fetal", "Adult"))

pdf("T_cell_expansion_tissues.pdf", height = 6, width = 15)
ggplot(rbind_TCR_clone_state, aes(x = Stage, y = value, fill = variable)) + 
	geom_bar(stat = "identity", position = "fill") +
	scale_fill_manual(values = rev(alpha(c("#019875", "#FEC6A2"), 0.8))) +
	scale_y_continuous(expand = c(0, 0.01)) +
	facet_grid(CD4_CD8 ~ Tissue, space = "free_y", scales = "free_y") +
	theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
				panel.spacing.x = unit(0, "lines"), 
				panel.border = element_rect(color = "black", fill = NA, size = 1),
				strip.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0))
dev.off()
###---------------------------------------TCR clone size plot:pie chart-----------------------------------------###
UMAP <- read.table("HCA2.0_TNK_UMAP.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
TCR_meta.data <- read.table("20230404_TNK_meta.data_with_TCR_information_.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
TCR_meta.data$renamed_Samlename <- TCR_meta.data$Name
TCR_meta.data$renamed_Samlename[TCR_meta.data$renamed_Samlename %>% grepl(pattern = "^Colon$")] <- "Jejunum"
TCR_meta.data$renamed_Samlename[TCR_meta.data$renamed_Samlename %>% grepl(pattern = "Stomach Protein")] <- "Stomach"
TCR_meta.data$renamed_Samlename[TCR_meta.data$renamed_Samlename %>% grepl(pattern = ".*blood|.*PBMC.*")] <- "Blood"

TCR_meta.data <- TCR_meta.data %>% filter(! Detail_Annotation_final %in% c("DG T1", "DG T2", "DG T3", "ILCs", "CD4CD8 low",
									"MAIT", "NK CD16 KIR2DL3","NK CD16 MYOM2",
									"NK CD56 CCL3", "NK CD56 KIR2DL4", "Pro ILC",
									"Pro NK"))

T_cell_colon_unique_barcode <- TCR_meta.data %>% select(row.names, customer_clone, Name:Stage, renamed_Samlename, Detail_Annotation_final) %>% unique
CD4 <- T_cell_colon_unique_barcode %>% filter(grepl(Detail_Annotation_final, pattern = "CD4|Treg")) %>% filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8")) 
CD8 <- T_cell_colon_unique_barcode %>% filter(grepl(Detail_Annotation_final, pattern = "CD8")) %>% filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8"))
Others <- T_cell_colon_unique_barcode %>% filter(Detail_Annotation_final %in% c("DP RAG1", "DP TOX2", "High pro DP", "CD4CD8 low"))

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
# subset_cells$Clone_state <- factor(TCR_meta.data$Clone_strategy %>% as.character(), levels = Clone_state)

TCR_color <- c("#019875", "#D2C49F", "#FEC6A2", "#FE9C7F", "#FE715C", "#FE4739", "#F70206", "#CC112B", "#96255B") # c("#EFEFEFFF"),

dat_for_plot <- with(TCR_meta.data , table(Clone_strategy %>% as.character(), Donor %>% as.character())) %>%
	unclass() %>% prop.table(margin = 2) %>% as.data.frame()

dat_for_plot <- rownames_to_column(dat_for_plot, var = "Clone_state") %>% melt()
setnames(dat_for_plot, old = c("variable", "value"), new = c("Donor", "Proportion"))
dat_for_plot$Clone_state <- factor(dat_for_plot$Clone_state %>% as.character(), levels = rev(Clone_state))
dat_for_plot$Donor <- factor(dat_for_plot$Donor %>% as.character(), levels = c("Fetal1", "Fetal2", "Fetal3", "Adult1", "Adult2", "Adult3"))

###-----------------CD4 and CD8 cell separated: barplot for TCR clone structure--------------###
dat_for_plot_CD4CD8 <- TCR_meta.data %>% filter(Clone_strategy != "NON") %>%
	filter(! Detail_Annotation_final %in% c("CD4CD8 low", "DG T1", "High pro DP", "DP RAG1", "DP TOX2", "DG T2",
		"DG T3", "ILCs", "NK CD16 KIR2DL3", "NK CD16 MYOM2", "NK CD56 CCL3", "NK CD56 KIR2DL4", "Pro ILC", "Pro NK", "MAIT"))

CD4_TCR <- dat_for_plot_CD4CD8 %>% filter(grepl(Detail_Annotation_final, pattern = "CD4|Treg")) %>%
	filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8")) %>% filter(as.character(Clone_strategy) != "NON")
CD8_TCR <- dat_for_plot_CD4CD8 %>% filter(grepl(Detail_Annotation_final, pattern = "CD8")) %>%
	filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8")) %>% filter(as.character(Clone_strategy) != "NON")

dat_CD4_plot <- table(CD4_TCR$Clone_strategy %>% as.character(), CD4_TCR$Donor) %>% unclass() %>% prop.table(margin = 2) %>% as.data.frame()
dat_CD4_plot <- rownames_to_column(dat_CD4_plot, var = "Clone_state") %>% melt()
setnames(dat_CD4_plot, old = c("variable", "value"), new = c("Donor", "Proportion"))
dat_CD4_plot$Clone_state <- factor(dat_CD4_plot$Clone_state %>% as.character(), levels = rev(Clone_state))
dat_CD4_plot$Donor <- factor(dat_CD4_plot$Donor %>% as.character(), levels = c("Fetal1", "Fetal2", "Fetal3", "Adult1", "Adult2", "Adult3"))
dat_CD4_plot$CD4_CD8 <- "CD4"

dat_CD8_plot <- table(CD8_TCR$Clone_strategy %>% as.character(), CD8_TCR$Donor) %>% unclass() %>% prop.table(margin = 2) %>% as.data.frame()
dat_CD8_plot <- rownames_to_column(dat_CD8_plot, var = "Clone_state") %>% melt()
setnames(dat_CD8_plot, old = c("variable", "value"), new = c("Donor", "Proportion"))
dat_CD8_plot$Clone_state <- factor(dat_CD8_plot$Clone_state %>% as.character(), levels = rev(Clone_state))
dat_CD8_plot$Donor <- factor(dat_CD8_plot$Donor %>% as.character(), levels = c("Fetal1", "Fetal2", "Fetal3", "Adult1", "Adult2", "Adult3"))
dat_CD8_plot$CD4_CD8 <- "CD8"

CD4_CD8_levels <- c(paste0("CD4_Fetal", c(1:3)), paste0("CD4_Adult", c(1:3)), paste0("CD8_Fetal", c(1:3)), paste0("CD8_Adult", c(1:3)))
dat_CD8_CD8_plot <- rbind(dat_CD4_plot, dat_CD8_plot) %>% as.data.frame()
dat_CD8_CD8_plot$Cluster <- factor(paste0(dat_CD8_CD8_plot$CD4_CD8, "_", dat_CD8_CD8_plot$Donor),
				levels = CD4_CD8_levels)

pdf("TCR_clone_structrue_each_donor_CD4_CD8_separated.pdf", width = 5, height = 6)
ggplot(data = dat_CD8_CD8_plot, aes(x = Cluster, y = Proportion, fill = Clone_state, color = Clone_state)) +
	geom_bar(stat = "identity", position = "fill") +
	# coord_polar(theta = "y") +
	scale_fill_manual(values = rev(TCR_color),
										name = "Clone size", labels = rev(Clone_state)) +
	scale_color_manual(values = rep("white", 9)) +
	scale_y_continuous(expand = c(0, 0.01)) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

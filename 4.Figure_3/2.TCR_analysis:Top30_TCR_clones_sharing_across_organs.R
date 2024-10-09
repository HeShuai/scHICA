packages <- list("Seurat", "Matrix", "stringr", "stringi", "gridExtra", "ggplot2",
								 "plyr", "ggthemes", "cowplot", "data.table", "RColorBrewer", "ComplexHeatmap",
								 "circlize", "pheatmap", "reshape2", "scales", "rlang", "future",
								 "parallel", "ggrepel", "ggsci", "ape", "dplyr", "harmony",
								 "sigminer", "ggpubr", "rstatix", "tidyverse", "dendextend",
								 "scRepertoire", "ggdendro", "ggtern", "viridis", "scatterpie", "zellkonverter")

lapply(packages, library, character.only = TRUE)
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]
options(width = 230)

tissues_colors <- c("Ascending colon" = "#E64B35FF", "Sigmoid colon" = '#ff0000ff', "Transverse colon" = '#d62790',
										Ileum = "#91D1C2FF", Jejunum = '#4DBBD5FF', Duodenum = '#00A087FF',
										"Lesser gastric curvature" = '#F0E685FF', "Mesenteric lymph node" = "#FFB570FF",
										Blood = "#1B1919FF", PBMC = "#374E55FF", "Whole blood" = "#3F4041FF",	Marrow = '#466983FF', Spleen = "#725663FF",
										Lung = "#6A6599FF", Liver = '#749B58FF', Kidney = '#B09C85FF', Skin = '#D6D6CEFF', Bladder = "#D49464FF",
										Esophagus = '#c5b0d5', Trachea = "#BA63affF", Heart = '#9467bd', Muscle = "#aec7e8", Pancreas = "#5050FFFF",
										Stomach = '#ADB17DFF', Thymus = '#8c564b')
#####-----------------------------------------TCR tracking analysis subclusters--------------------------------###
###---------------------------------------This is a demo for ploting TRM TCR sharing across organs in adult1. A similar strategy can be applied for other donor accordingly.---------------------------####
###---------------------------------------TCR clone size plot--------------------------###
UMAP <- read.table("HCA2.0_TNK_UMAP.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
TCR_meta.data_1 <- read.table("20230404_TNK_meta.data_with_TCR_information_.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
TCR_meta.data_1$renamed_Samlename <- TCR_meta.data_1$Name
TCR_meta.data_1$renamed_Samlename[TCR_meta.data_1$renamed_Samlename %>% grepl(pattern = "^Colon$")] <- "Jejunum"
TCR_meta.data_1$renamed_Samlename[TCR_meta.data_1$renamed_Samlename %>% grepl(pattern = "Stomach Protein")] <- "Stomach"
TCR_meta.data_1$renamed_Samlename[TCR_meta.data_1$renamed_Samlename %>% grepl(pattern = ".*blood|.*PBMC.*")] <- "Blood"

TCR_meta.data <- TCR_meta.data_1 %>% filter(! Detail_Annotation_final %in% c("DG T1", "DG T2", "DG T3", "ILCs", "CD4CD8 low",
																																					 "MAIT", "NK CD16 KIR2DL3","NK CD16 MYOM2",
																																					 "NK CD56 CCL3", "NK CD56 KIR2DL4", "Pro ILC",
																																					 "Pro NK"))

T_cell_colon_unique_barcode <- TCR_meta.data %>% select(row.names, customer_clone, Name:Stage, renamed_Samlename, Detail_Annotation_final) %>% unique
CD4 <- T_cell_colon_unique_barcode %>% filter(grepl(Detail_Annotation_final, pattern = "CD4|Treg")) %>%
	filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8")) %>% filter(Donor == "Adult1")
CD8 <- T_cell_colon_unique_barcode %>% filter(grepl(Detail_Annotation_final, pattern = "CD8")) %>%
	filter(! grepl(Detail_Annotation_final, pattern = "DP|CD4CD8")) %>% filter(Donor == "Adult1")

clone_ID_CD4 <- CD4$customer_clone %>% table %>% as.data.frame() %>% filter(`.` != "Not")
clone_ID_CD4_top30 <- clone_ID_CD4 %>% arrange(desc(Freq)) %>% head(30) %>% `[`(".") %>% unlist() %>% unname() %>% as.character()
CD4$Clone_Size <- mapvalues(CD4$customer_clone, from = clone_ID_CD4$".", to = clone_ID_CD4$Freq, warn_missing = F) %>% as.numeric()

clone_ID_CD8 <- CD8$customer_clone %>% table %>% as.data.frame() %>% filter(`.` != "Not")
clone_ID_CD8_top30 <- clone_ID_CD8 %>% arrange(desc(Freq)) %>% head(30) %>% `[`(".") %>% unlist() %>% unname() %>% as.character()
CD8$Clone_Size <- mapvalues(CD8$customer_clone, from = clone_ID_CD8$".", to = clone_ID_CD8$Freq, warn_missing = F) %>% as.numeric()

CD4_TCR <- CD4 %>% filter(grepl(Detail_Annotation_final, pattern = "CD4|Treg")) %>% filter(customer_clone %in% clone_ID_CD4_top30)
CD8_TCR <- CD8 %>% filter(grepl(Detail_Annotation_final, pattern = "CD8")) %>% filter(customer_clone %in% clone_ID_CD8_top30)
CD8_TCR$num <- 1
CD4_TCR$num <- 1

###------------replace CD8_TCR with CD4_TCR to plot TCR sharing across organs for CD4 T cells--------####
top30_CD8_TCR <- CD8_TCR %>% select(customer_clone, Detail_Annotation_final, renamed_Samlename, num)

sum_CD8_TCR <- top30_CD8_TCR %>% group_by(customer_clone, Detail_Annotation_final, renamed_Samlename) %>% summarise_at(.vars = "num", .funs = sum)
sum_CD8_TCR$customer_clone <- sum_CD8_TCR$customer_clone %>% as.factor() %>% as.numeric()
sum_CD8_TCR$customer_clone <- mapvalues(sum_CD8_TCR$customer_clone %>% as.numeric(), from = 1:30, to = seq(1, 30, 1))#seq(1, 45, 1.5))

sum_CD8_TCR_1 <- dcast(sum_CD8_TCR, customer_clone + renamed_Samlename ~ Detail_Annotation_final, fill = 0)
sum_CD8_TCR_1$size <- sum_CD8_TCR_1[, -c(1:2)] %>% rowSums()
sum_CD8_TCR_1$Organ <- sum_CD8_TCR_1$renamed_Samlename
sum_CD8_TCR_1$renamed_Samlename <- sum_CD8_TCR_1$renamed_Samlename %>% as.factor() %>% as.numeric()

overlapped_names <- intersect(tissues_colors %>% names(), sum_CD8_TCR_1$Organ)
sum_CD8_TCR_1$Organ <- factor(sum_CD8_TCR_1$Organ, levels = overlapped_names)
sum_CD8_TCR_1$renamed_Samlename <- mapvalues(sum_CD8_TCR_1$Organ %>% as.character(),
																						 from = sum_CD8_TCR_1$Organ %>% levels(), to = seq(1, 1*length(sum_CD8_TCR_1$Organ %>% levels()), 19/19))
sum_CD8_TCR_1$renamed_Samlename <- sum_CD8_TCR_1$renamed_Samlename %>% as.numeric()
sum_CD8_TCR_1$size[sum_CD8_TCR_1$size > 600] <- 600
sum_CD8_TCR_1$size[sum_CD8_TCR_1$size < 5] <- 0
# sum_CD8_TCR_1 <- sum_CD8_TCR_1 %>% select(customer_clone:renamed_Samlename, Organ, size, `CD8 TN/CM`,
# 																					`CD8 TEM`, `CD8 TEM GZMK`, `CD8 TEFF`, `CD8 TRM GPR15`, `CD8 TRM ZNF683`, `Pro CD8 T`) #, `CD8 TN/CM` `Pro CD8 T`
sum_CD8_TCR_1 <- sum_CD8_TCR_1 %>% select(customer_clone:renamed_Samlename, Organ, size, `CD4 TN/CM CCR7 low`, `CD4 Tfh CXCR5`, `CD4 TEM`, `CD4 TEM IL2`,
													`CD4 TEFF`, `CD4 Th17`, `CD4 TRM CCR6`, `CD4 TRM EMP1`, `Treg naive`, `Pro Treg`) #, `CD8 TN/CM` `Pro CD8 T`
# `CD4 TN/CM CCR7 low`, `CD4 TEM`, `CD4 TEM IL2`,
# `CD4 Tfh CXCR5`, `CD4 TEFF`, `CD4 Th2`, `CD4 Th17`, `CD4 TRM CCR6`, `CD4 TRM EMP1`, `Treg naive`, `Pro Treg`, `Pro CD4 T`
####----------------------order the TCR, from universal to specific-----------------------##
order_TCR <- sum_CD8_TCR_1 %>% filter(size > 0) %>% select(customer_clone) %>% table %>% sort(decreasing = T) %>% names() %>% as.numeric()
sum_CD8_TCR_1 <- sum_CD8_TCR_1 %>% filter(customer_clone %in% order_TCR)
sum_CD8_TCR_1$customer_clone <- mapvalues(sum_CD8_TCR_1$customer_clone, from = order_TCR %>% as.numeric(), to = seq(1, 30, 30/30)) #seq(1, 45, 1.5))

pdf("CD4_TCR_sharing_accross_organs_Adult1.pdf", height = 10, width = 8)
ggplot() + geom_scatterpie(aes(x = customer_clone, y = renamed_Samlename, r = 0.15*sqrt(size/2/pi)), 
								data = sum_CD8_TCR_1, cols = c(sum_CD8_TCR_1 %>% colnames())[5:14], long_format = F, color = NA) + #coord_equal() +
	scale_x_continuous(breaks = unique(sum_CD8_TCR_1$customer_clone), labels = unique(sum_CD8_TCR_1$customer_clone)) +
	scale_y_reverse(breaks = unique(sum_CD8_TCR_1$renamed_Samlename), labels = unique(sum_CD8_TCR_1$Organ)) +
	coord_fixed() +
	theme(panel.grid.major = element_line(linewidth = 0.2, color = "grey70", linetype = 1),
				panel.background = element_blank(),
				axis.ticks = element_blank()) +
	scale_fill_manual(values = color_used) +
	geom_scatterpie_legend(0.15*sqrt(sum_CD8_TCR_1$size/2/pi),  x = 30, y = 1.5)# x = 35, y = 1)
dev.off()



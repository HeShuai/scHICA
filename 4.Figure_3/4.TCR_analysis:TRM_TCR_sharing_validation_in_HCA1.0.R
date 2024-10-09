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

options(future.globals.maxSize = 400*1000 * 1024^2)
plan("multisession", workers = 20)
plan()
####-------------------------HCA1.0 validation----------------------#####
library(readxl)
TCR_CD8_HCA1.0 <- read_excel("D:/20200715_GB_RV/12. GB.RV2/Figures_and_Supplementary_materials/GB_RV3_V1/Final_Version/Additional file 4_for_HCA2.0_validation.xlsx",
														 sheet = "Table S27 CD8")
TCR_CD8_HCA1.0 <- read_excel("D:/20200715_GB_RV/12. GB.RV2/Figures_and_Supplementary_materials/GB_RV3_V1/Final_Version/Additional file 4_for_HCA2.0_validation.xlsx",
														 sheet = "Table S26 CD4")

selected_cells <- TCR_CD8_HCA1.0 %>% select(cell_barcode, tissue, `cell type`, raw_clonotype_id) %>% unique() %>% setnames(old = "cell type", new = "Cell_type")
CD8 <- selected_cells
setnames(CD8, old = c("tissue", "Cell_type", "raw_clonotype_id"), new = c("renamed_Samlename", "Detail_Annotation_final", "customer_clone"))
###----------------for CD8 T cell-----------------------###
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*TEM.*", replacement = "CD8 TEM")
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*TEFF.*", replacement = "CD8 TEFF")
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*TRM.*", replacement = "CD8 TRM")
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*IEL.*", replacement = "CD8 TRM")
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*MAIT.*", replacement = "MAIT")
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*TN/CM.*", replacement = "CD8 TN/CM")
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*TN.*", replacement = "CD8 TN/CM")

###----------------for CD4 T cell-----------------------###
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*CTLA4_Treg.*", replacement = "CD4 Treg")
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*TN/CM.*", replacement = "CD4 TN/CM")
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*TN.*", replacement = "CD4 TN/CM")
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*TCM.*", replacement = "CD4 TN/CM")
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*Th1.*", replacement = "CD4 TEFF")
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*TRM.*", replacement = "CD4 TRM")
CD8$Detail_Annotation_final <- CD8$Detail_Annotation_final %>% gsub(pattern = ".*TEM.*", replacement = "CD4 TEM")

CD8 <- CD8 %>% filter(Detail_Annotation_final != "MAIT")

clone_ID_CD8 <- CD8$customer_clone %>% table %>% as.data.frame() %>% filter(`.` != "Not")
clone_ID_CD8_top30 <- clone_ID_CD8 %>% arrange(desc(Freq)) %>% head(30) %>% `[`(".") %>% unlist() %>% unname() %>% as.character()
CD8$Clone_Size <- mapvalues(CD8$customer_clone, from = clone_ID_CD8$".", to = clone_ID_CD8$Freq, warn_missing = F) %>% as.numeric()

CD8_TCR <- CD8 %>% filter(customer_clone %in% clone_ID_CD8_top30)
CD8_TCR$num <- 1
top30_CD8_TCR <- CD8_TCR %>% select(customer_clone, Detail_Annotation_final, renamed_Samlename, num)

sum_CD8_TCR <- top30_CD8_TCR %>% group_by(customer_clone, Detail_Annotation_final, renamed_Samlename) %>% summarise_at(.vars = "num", .funs = sum)
sum_CD8_TCR$customer_clone <- sum_CD8_TCR$customer_clone %>% as.factor() %>% as.numeric()
sum_CD8_TCR$customer_clone <- mapvalues(sum_CD8_TCR$customer_clone %>% as.numeric(), from = 1:30, to = seq(1, 30, 1))#seq(1, 45, 1.5))

sum_CD8_TCR_1 <- dcast(sum_CD8_TCR, customer_clone + renamed_Samlename ~ Detail_Annotation_final, fill = 0)
sum_CD8_TCR_1$size <- sum_CD8_TCR_1[, -c(1:2)] %>% rowSums()
sum_CD8_TCR_1$Organ <- sum_CD8_TCR_1$renamed_Samlename %>% as.factor()

sum_CD8_TCR_1$renamed_Samlename <- mapvalues(sum_CD8_TCR_1$Organ %>% as.character(),
																						 from = sum_CD8_TCR_1$Organ %>% levels(), to = seq(1, 1*length(sum_CD8_TCR_1$Organ %>% levels), 14/14))
sum_CD8_TCR_1$renamed_Samlename <- sum_CD8_TCR_1$renamed_Samlename %>% as.numeric()
sum_CD8_TCR_1$size[sum_CD8_TCR_1$size > 600] <- 600
sum_CD8_TCR_1$size[sum_CD8_TCR_1$size < 5] <- 0
# sum_CD8_TCR_1 <- sum_CD8_TCR_1 %>% select(customer_clone:renamed_Samlename, Organ, size,
# 																					`CD8 TEM`, `CD8 TEFF`, `CD8 TRM`) #, `CD8 TN/CM` `Pro CD8 T`
sum_CD8_TCR_1 <- sum_CD8_TCR_1 %>% select(customer_clone:renamed_Samlename, Organ, size, `CD4 TN/CM`, `CD4 TEM`,
																					`CD4 TEFF`, `CD4 Treg`, `CD4 TRM`) #, `CD8 TN/CM` `Pro CD8 T`
# `CD4 TN/CM CCR7 low`, `CD4 TEM`, `CD4 TEM IL2`,
# `CD4 Tfh CXCR5`, `CD4 TEFF`, `CD4 Th2`, `CD4 Th17`, `CD4 TRM CCR6`, `CD4 TRM EMP1`, `Treg naive`, `Pro Treg`, `Pro CD4 T`
####----------------------order the TCR, from universal to specific-----------------------##
order_TCR <- sum_CD8_TCR_1 %>% filter(size > 0) %>% select(customer_clone) %>% table %>% sort(decreasing = T) %>% names() %>% as.numeric()
sum_CD8_TCR_1 <- sum_CD8_TCR_1 %>% filter(customer_clone %in% order_TCR)
sum_CD8_TCR_1$customer_clone <- mapvalues(sum_CD8_TCR_1$customer_clone, from = order_TCR %>% as.numeric(), to = seq(1, 6, 6/6)) #seq(1, 45, 1.5))

# pdf("CD4_TCR_sharing_accross_organs_Adult3.pdf", height = 10, width = 8)
pdf("CD4_TCR_sharing_accross_organs_HCA1.0.pdf", height = 10, width = 8)
ggplot() + geom_scatterpie(aes(x = customer_clone, y = renamed_Samlename, r = 0.3*sqrt(size/2/pi)), 
													 data = sum_CD8_TCR_1, cols = c(sum_CD8_TCR_1 %>% colnames())[5:9], long_format = F, color = NA) + #coord_equal() +
	scale_x_continuous(breaks = unique(sum_CD8_TCR_1$customer_clone), labels = unique(sum_CD8_TCR_1$customer_clone)) +
	scale_y_reverse(breaks = unique(sum_CD8_TCR_1$renamed_Samlename), labels = unique(sum_CD8_TCR_1$Organ)) +
	coord_fixed() +
	theme(panel.grid.major = element_line(linewidth = 0.2, color = "grey70", linetype = 1),
				panel.background = element_blank(),
				axis.ticks = element_blank()) +
	scale_fill_manual(values = color_used) +
	geom_scatterpie_legend(0.3*sqrt(sum_CD8_TCR_1$size/2/pi),  x = 6, y = 1.5)# x = 35, y = 1)
dev.off()

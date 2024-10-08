packages <- list("Seurat", "Matrix", "stringr", "stringi", "gridExtra", "ggplot2",
								 "plyr", "ggthemes", "cowplot", "data.table", "RColorBrewer", "ComplexHeatmap",
								 "circlize", "pheatmap", "reshape2", "scales", "rlang", "future",
								 "parallel", "ggrepel", "ggsci", "ape", "dplyr", "harmony",
								 "sigminer", "ggpubr", "rstatix", "ROGUE", "tidyverse", "reticulate", "dendextend",
								 "scRepertoire", "ggdendro", "ggtern", "viridis", "scatterpie", "SeuratData", "SeuratDisk",
								 "scibet", "scuttle", "zellkonverter", "cowplot", "gmodels", "ggalt", "ggalluvial", "xtable") ## "descr": CrossTable

lapply(packages, library, character.only = TRUE)
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]
options(width = 240)
options(future.globals.maxSize = 400*1000 * 1024^2)
plan("multisession", workers = 1)
plan()

###------------------------------------------------###
load("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/HCA2.0_alltissues_final_1038993.RData")

DefaultAssay(subset_cells) <- "RNA"
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features = VariableFeatures(subset_cells)[1])

subset_cells$renamed_Samlename <- subset_cells$Name
subset_cells$renamed_Samlename[subset_cells$renamed_Samlename %>% grepl(pattern = "^Colon$")] <- "Jejunum"
subset_cells$renamed_Samlename[subset_cells$renamed_Samlename %>% grepl(pattern = "Stomach Protein")] <- "Stomach"
# subset_cells$renamed_Samlename[subset_cells$renamed_Samlename %>% grepl(pattern = ".*blood|.*PBMC.*")] <- "Blood"
subset_cells$blood_split <- subset_cells$Name
subset_cells$blood_split[subset_cells$blood_split %>% grepl(pattern = "^Colon$")] <- "Jejunum"
subset_cells$blood_split[subset_cells$blood_split %>% grepl(pattern = "Stomach Protein")] <- "Stomach"

subset_cells$Donor_Name <- paste0(subset_cells$Donor, "_", subset_cells$blood_split)

###------------------average expression PCA analysis--------------------------###

Idents(subset_cells) <- subset_cells$Donor_Name
ave_subset_cells <- AverageExpression(subset_cells, return.seurat = T)
ave_subset_cells$Donor_Name <- row.names(ave_subset_cells@meta.data)

ave_subset_cells <- FindVariableFeatures(ave_subset_cells, nfeatures = 3000, selection.method = 'vst', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf))
ave_subset_cells <- ScaleData(ave_subset_cells, verbose = TRUE)

ave_subset_cells <- RunPCA(ave_subset_cells, features = setdiff(ave_subset_cells@assays$RNA@var.features, row.names(ave_subset_cells) %>% grep(pattern = "^MT-|^RPL|^RPS|^IGK|^IGV|^IGL|^IGH", v = T)) , npcs = 50, verbose = TRUE)

pdf(paste0("HCA2.0_all_cells", "_PCA", ".pdf"),
		width = 20, height = 15)
DimPlot(object = ave_subset_cells, reduction = 'pca', label = F, cols = c(color_used, colors_dxs), pt.size = 5, group.by = "orig.ident", label.size = 5, raster = F) +
	guides(colour = guide_legend(ncol = 4, override.aes = list(size = 5))) 
dev.off()

ave_subset_cells$Donor_Name <- ave_subset_cells$Donor_Name %>% gsub(pattern = "1_|2_|3_", replacement = "_")
###------------------------PCA plot add the circle-------------------------###
PCA <- Embeddings(object = ave_subset_cells, reduction = "pca") %>% as.data.frame()
PCA$Donor_Name <- ave_subset_cells$Donor_Name
PCA$orig.ident <- ave_subset_cells$orig.ident

pdf("PCA_map_Donor_Name.pdf", height = 10, width = 13)
ggplot(PCA, aes(x = PC_1, y = PC_2, color = Donor_Name, group = Donor_Name)) +
	labs(x = "PC 1",
			 y = "PC 2",
			 title = "PCA") +
	geom_point(size = 5) + 
	# geom_encircle(aes(fill = Donor_Name), alpha = 0.1, show.legend = F) +
	theme_classic() + coord_fixed(1) +
	scale_color_manual(values = color_used)
dev.off()

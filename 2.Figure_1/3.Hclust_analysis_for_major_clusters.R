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

###------------------------------------------clustering ----------------------------------
### as organs
cluster.averages <- AverageExpression(subset_cells, return.seurat = T)

DEGs <- read.table("HCA2.0_allcell_subset_cells_culster_all_DEGs_allcell.csv", sep = ",", header = T, row.names = 1, stringsAsFactors = F) %>%
  # mutate(Cluster = plyr::mapvalues(cluster, from = top50$cluster, to = top50$Cluster)) %>% filter(grepl(Cluster, pattern = "")) %>% 
  TOP_N(500, pct.1 = 0.2, sig.padj = 0.05, fc.threshold = 0.25)

DEGs_top50 <- read.table("HCA2.0_subset_cells_culster_top50_DEGs_allcell.csv", sep = ",", header = T, row.names = 1, stringsAsFactors = F) %>%
  # mutate(Cluster = plyr::mapvalues(cluster, from = top50$cluster, to = top50$Cluster)) %>% filter(grepl(Cluster, pattern = "")) %>% 
  TOP_N(50, pct.1 = 0.2, sig.padj = 0.05, fc.threshold = 0.2)

spearman_cor <- cor(cluster.averages@assays$RNA@data %>% as.matrix() %>% `[`(DEGs$gene %>% unique(), ), method = "pearson")
row.names(spearman_cor) <- mapvalues(row.names(spearman_cor), from = DEGs_top50$cluster %>% as.character, to = DEGs_top50$Cell_type %>% as.character)
colnames(spearman_cor) <- mapvalues(colnames(spearman_cor), from = DEGs_top50$cluster %>% as.character, to = DEGs_top50$Cell_type %>% as.character)

# row.names(spearman_cor) <- mapvalues(row.names(spearman_cor), from = tissue_to_number %>% names, to = tissue_to_number)
# colnames(spearman_cor) <- mapvalues(colnames(spearman_cor), from = tissue_to_number %>% names, to = tissue_to_number)

distance_of_cluster <- as.dist(1-spearman_cor)

hc <- hclust(distance_of_cluster, method = "complete")
hc <- as.dendrogram(hc)
# dend_data <- dendro_data(hc, type = "rectangle")
asggden <- as.ggdend(hc)

pdf("HCA_Hcluster_all_tissue.pdf", width = 15, height = 15)
ggplot(asggden, horiz = TRUE, theme = NULL, offset_labels = -0.01) +
  ylim(0.8, -0.5) +
  # theme_base()+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 

dev.off()

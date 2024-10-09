packages <- list("Seurat", "Matrix", "stringr", "stringi", "gridExtra", "ggplot2",
								 "plyr", "ggthemes", "cowplot", "data.table", "RColorBrewer", "ComplexHeatmap",
								 "circlize", "pheatmap", "reshape2", "scales", "rlang", "future",
								 "parallel", "ggrepel", "ggsci", "ape", "dplyr", "harmony",
								 "sigminer", "ggpubr", "rstatix", "ROGUE", "tidyverse", "reticulate", "dendextend",
								 "scRepertoire", "ggdendro", "ggtern", "viridis", "scatterpie", "SeuratData", "SeuratDisk",
								 "scibet", "scuttle", "zellkonverter", "cowplot", "CellChat")

lapply(packages, library, character.only = TRUE)
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7), pal_lancet()(9), pal_frontiers()(9))[-8]
options(width = 240)
source("/public/home/HeShuai/Projects/Script/1.TOP_N.R")
options(future.globals.maxSize = 400*1000 * 1024^2)
plan("multisession", workers = 1)
plan()

load('2024.06.08.adult_nonimmune_cells_1500_cells_HCA1.0.RData')
load('2024.06.08.HCA2.0_annotation_adult_T_cells.RData')

load('2024.06.08.fetal_nonimmune_cells_1500_cells_SCIENCE.RData')
load('2024.06.08.HCA2.0_annotation_fetal_T_cells.RData')

DefaultAssay(Fetal_T) <- "RNA"
Fetal_T@assays$RNA@data <- Fetal_T@assays$RNA@counts
Fetal_T <- NormalizeData(Fetal_T)
Fetal_T <- ScaleData(Fetal_T, features = VariableFeatures(Fetal_T)[1])

DefaultAssay(Adult_T) <- "RNA"
Adult_T@assays$RNA@data <- Adult_T@assays$RNA@counts
Adult_T <- NormalizeData(Adult_T)
Adult_T <- ScaleData(Adult_T, features = VariableFeatures(Adult_T)[1])

DefaultAssay(Adult_nonimmune) <- "RNA"
Adult_nonimmune@assays$RNA@data <- Adult_nonimmune@assays$RNA@counts
Adult_nonimmune <- NormalizeData(Adult_nonimmune)
Adult_nonimmune <- ScaleData(Adult_nonimmune, features = VariableFeatures(Adult_nonimmune)[1])

DefaultAssay(Fetal_nonimmune) <- "RNA"
Fetal_nonimmune@assays$RNA@data <- Fetal_nonimmune@assays$RNA@counts
Fetal_nonimmune <- NormalizeData(Fetal_nonimmune)
Fetal_nonimmune <- ScaleData(Fetal_nonimmune, features = VariableFeatures(Fetal_nonimmune)[1])


Adult_T@meta.data[, 1:5] <- NULL
Adult_T@meta.data[, 4:13] <- NULL
colnames(Adult_nonimmune@meta.data)[1] <- "Name"
Adult_nonimmune@meta.data[, c(2:6, 8)] <- NULL
colnames(Adult_nonimmune@meta.data)[2] <- "Major_clusters"
Adult_nonimmune@meta.data$Stage <- "Adult"
Adult_nonimmune@meta.data$Donor <- NULL
Fetal_nonimmune@meta.data[, 1:5] <- NULL
Fetal_nonimmune@meta.data[, c(3, 4, 6, 7, 8)] <- NULL
Fetal_nonimmune@meta.data[, c(2)] <- NULL
Fetal_nonimmune$Stage <- "Fetal"
colnames(Fetal_nonimmune@meta.data)[2] <- "Major_clusters"
colnames(Fetal_nonimmune@meta.data)[1] <- "Name"
Fetal_T@meta.data[, 1:5] <- NULL
Fetal_T@meta.data[, c(2, 4:13)] <- NULL

Adult <- merge(Adult_T, Adult_nonimmune)
Fetus <- merge(Fetal_T, Fetal_nonimmune)

data.input = Fetus@assays$RNA@data  # normalized data matrix
meta = Fetus@meta.data # a dataframe with rownames containing cell mata data

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Major_clusters")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 40) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

pdf("Cell_chat_fetal_immune_vs_stroma.pdf", width = 15, height = 7.5)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
								 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
								 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

save(cellchat, file = "Fetus.immune_vs_stroma.cellchat.RData")
save(cellchat, file = "Adult.immune_vs_stroma.cellchat.RData")

adults_chat <- get(load("Adult.immune_vs_stroma.cellchat.RData"))
fetus_chat <- get(load("Fetus.immune_vs_stroma.cellchat.RData"))

cco.list <- list(adult = adults_chat, fetus = fetus_chat)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("Overview_number_strength_immune_stroma.pdf", p, width = 6, height = 4)

pdf("Cell_chat_Fetus_adult_Ni_immune_stroma.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
	netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge = F, 
									 edge.weight.max = weight.max[2], edge.width.max = 5, 
									 title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}
# save as Counts_Compare_net.pdf
dev.off()

pdf("Cell_chat_Fetus_adult_W_immune_stroma.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","weight"))
for (i in 1:length(cco.list)) {
	netVisual_circle(cco.list[[i]]@net$weight, weight.scale = T, label.edge = F, 
									 edge.weight.max = weight.max[2], edge.width.max = 5, 
									 title.name = paste0("Strength of interactions - ", names(cco.list)[i]))
}
# save as Counts_Compare_net.pdf
dev.off()

## comparison of signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("Compare_pathway_strengh_immune_stroma.pdf", p, width = 10, height = 6)

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
#netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
p <- rankSimilarity(cellchat, type = "structural") + ggtitle("Structural similarity of pathway")
ggsave("Pathway_Similarity.pdf", p, width = 8, height = 5)

pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, 
																				title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
																				title = names(cco.list)[2], width = 8, height = 10)

p <- ht1 + ht2
ggsave("Compare_pathway_overall_strengh.pdf", p, width = 10, height = 6)

# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_all.pdf  10*6
pairLRuse <- data.frame(interaction_name = "PGE2-PTGES3_PTGER4")

levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(3:5), targets.use = c(1:2),  comparison = c(1, 2), angle.x = 45)
ggsave("Compare_LR_bubble_Fetus_adult_immune_stroma_CD4_CD8.pdf", p, width = 5, height = 10)

p <- netVisual_bubble(cellchat, sources.use = c(3:5), targets.use = c(1),  comparison = c(1, 2), angle.x = 45, pairLR.use = pairLRuse)
ggsave("Compare_LR_bubble_CD4_selected.pdf", p, width = 5, height = 3)

p <- netVisual_bubble(cellchat, sources.use = c(3:5), targets.use = c(2),  comparison = c(1, 2), angle.x = 45, pairLR.use = pairLRuse)
ggsave("Compare_LR_bubble_CD8_selected.pdf", p, width = 5, height = 3)

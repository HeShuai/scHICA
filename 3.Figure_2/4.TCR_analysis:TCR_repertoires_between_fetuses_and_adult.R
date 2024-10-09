packages <- list("Seurat", "Matrix", "stringr", "stringi", "gridExtra", "ggplot2",
								 "plyr", "ggthemes", "cowplot", "data.table", "RColorBrewer", "ComplexHeatmap",
								 "circlize", "pheatmap", "reshape2", "scales", "rlang", "future",
								 "parallel", "ggrepel", "ggsci", "ape", "dplyr", "harmony",
								 "sigminer", "ggpubr", "rstatix", "ROGUE", "tidyverse", "reticulate", "dendextend",
								 "scRepertoire", "ggdendro", "ggtern", "viridis", "scatterpie", "SeuratData", "SeuratDisk",
								 "scibet", "scuttle", "zellkonverter", "pbapply", "ggalt")
lapply(packages, library, character.only = TRUE)
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]
options(width = 240)

options(future.globals.maxSize = 400*1000 * 1024^2)
plan("multisession", workers = 1)
plan()
###-----------------------seting python environment--------------------------### 
###---------loading the files-------------###
meta.data <- read.table("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/All_cells.meta.data.txt", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
meta.data$Name <- meta.data$Name %>% gsub(pattern = "^Colon$", replacement = "Jejunum")
meta.data$Name <- meta.data$Name %>% gsub(pattern = "^Stomach Protein$", replacement = "Stomach")
TCRinform <- read.table("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/TNK/Cells/TCR_rawdata_cellranger6/Merged_TCR.information_at_least_one_paired_VDJ_cellranger6.txt",
												sep = "\t", header = T, row.names = 1, stringsAsFactors = F)									 		

lst <- list()
for (i in dir() %>% grep(pattern = "TCR_filtered_contig_annotations.csv", value = T)) {
	sample_name <- i %>% gsub(pattern = "_TCR_filtered_contig_annotations.csv", replacement = "")
	tmp <- read.table(i, header = T, row.names = NULL, sep = ",", stringsAsFactors = F)
	tmp$samples <- mapvalues(sample_name, from = meta.data$orig.ident, to = meta.data$Donor, warn_missing = F)
	tmp <- tmp %>% filter(barcode %in% TCRinform$barcode)
	lst[[sample_name]] <- tmp
	print(i)}

combined.TCR <- combineTCR(lst, samples = mapvalues(names(lst), from = meta.data$orig.ident, to = meta.data$Donor, warn_missing = F), removeNA = TRUE,  removeMulti = FALSE,  filterMulti = FALSE)

percentVJ_fixed <- function(input.data, chain = "TRB", group.by = NULL, order.by = NULL, 
					exportTable = TRUE, palette = "inferno") 
{
	input.data <- scRepertoire:::.data.wrangle(input.data, group.by, "CTgene", 
															chain)
	input.data <- scRepertoire:::.groupList(input.data, group.by)
	
	if (chain %in% c("TRA", "TRG", "IGL")) {
		positions <- c(1, 2)
	}
	else {
		positions <- c(1, 3)
	}
	gene_counts <- lapply(input.data, function(x) {
		tmp <- unlist(str_split(x[, "CTgene"], ";"))
		tmp <- str_split(tmp, "[.]", simplify = TRUE)
		strings <- paste0(tmp[, positions[1]], ";", tmp[, positions[2]])
		strings <- strings[which(strings != "NA;")]
	})
	gene.dictionary <- unique(unlist(gene_counts))
	coordinates <- strsplit(gene.dictionary, ";")
	V_values <- unique(unlist(lapply(coordinates, function(coord) coord[1])))
	J_values <- unique(unlist(lapply(coordinates, function(coord) coord[2])))
	summary <- lapply(gene_counts, function(x) {
		percentages <- unlist(prop.table(table(x)))
		genes.to.add <- setdiff(gene.dictionary, names(percentages))
		if (length(genes.to.add) >= 1) {
			percentages.to.add <- rep(0, length(genes.to.add))
			names(percentages.to.add) <- genes.to.add
			percentages <- c(percentages, percentages.to.add)
		}
		percentages <- percentages[match(str_sort(names(percentages), 
																							numeric = TRUE), names(percentages))]
	})
	if (exportTable == TRUE) {
		summary.matrix <- do.call(rbind, summary)
		return(summary.matrix)
	}
}

pdf("TCR.information.AA.TRA.pdf", height = 10, width = 10)
percentAA_fixed(combined.TCR, 
					chain = "TRA", 
					aa.length = 20,
					palette = color_used)
dev.off()

###-----------------------------------VJ pair usage------------------------------###
TCRinform <- read.table("E:/8.Everyday_working/R/202200221_HCA2.0_analysis/TCR/First_run/Merged_TCR.information_at_least_one_paired_VDJ_cellranger6.txt",
												sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
TCRinform$Name <- TCRinform$Name %>% gsub(pattern = "^Colon$", replacement = "Jejunum")
TCRinform$Name <- TCRinform$Name %>% gsub(pattern = "^Stomach Protein$", replacement = "Stomach")
VJ_pools <- TCRinform %>% dplyr::select(v_gene, j_gene) %>% mutate(VJ = paste0(v_gene, ",", j_gene)) %>% select(VJ) %>% unlist() %>% unique

counts_VJ <- lapply(TCRinform$Donor %>% unique, function(x) {tmp <- TCRinform %>% filter(Donor == x) %>%
	dplyr::select(v_gene, j_gene) %>% mutate(VJ = paste0(v_gene, ",", j_gene)) %>% select(VJ) %>% unlist() %>% unname %>% table %>% as.data.frame
colnames(tmp)[1] <- "VJ"
tmp$Donor <- x
return(tmp)
})
names(counts_VJ) <- TCRinform$Donor %>% unique
merged_data <- do.call(rbind, counts_VJ)
merged_data <- merged_data[, c(1, 3, 2)] %>% as.data.frame()
merged_data$VJ <- merged_data$VJ %>% as.character()
merged_data <- dcast(merged_data, Donor ~ VJ) %>% as.data.frame() %>% column_to_rownames(var = "Donor") %>% as.matrix()
merged_data[is.na(merged_data)] <- 0
merged_data <- merged_data %>% as.data.frame()

pc <- prcomp(merged_data)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(pc$x[,1:2]), c(rep(c("Adult", "Fetus"), each = 3))))
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

pdf("PCA_of_VJ_adult_fetus.pdf", height = 10, width = 10)
ggplot(df, aes(x = PC1, y = PC2, Group = V4)) + 
	geom_point(aes(fill = df[,3]), shape = 21, size = 5) + 
	guides(fill = guide_legend(title = "Samples")) + 
	scale_fill_manual(values = color_used) +
	geom_encircle(aes(fill = V4), alpha = 0.1, show.legend = F) +
	coord_fixed(1) +
	theme_classic()
dev.off()

###-------------------------antigen recognization-----------------------------###
database_TRA <- read.table("E:/8.Everyday_working/R/202200221_HCA2.0_analysis/TCR/Trex.database.TRA.txt", sep = "\t", header = T, row.names = NULL, stringsAsFactors = F, quote = "")	
database_TRB <- read.table("E:/8.Everyday_working/R/202200221_HCA2.0_analysis/TCR/Trex.database.TRB.txt", sep = "\t", header = T, row.names = NULL, stringsAsFactors = F, quote = "")									 		
HCA_TRA <- TCRinform %>% filter(chain == "TRA") %>% dplyr::select(cdr3, Donor, Name)
HCA_TRB <- TCRinform %>% filter(chain == "TRB") %>% dplyr::select(cdr3, Donor, Name)

HCA_TRA$Antigen <- mapvalues(HCA_TRA$cdr3, from = database_TRA$CDR3, to = database_TRA$Epitope.species, warn_missing = F)
HCA_TRA$Antigen[!HCA_TRA$Antigen %in% database_TRA$Epitope.species] <- "No"
HCA_TRA_type <- HCA_TRA$Antigen %>% str_split(pattern = ";")
HCA_TRA[lapply(HCA_TRA_type, length) %>% do.call(rbind, .) %>% c != 1, "Antigen"] <- "Multple"

HCA_TRB$Antigen <- mapvalues(HCA_TRB$cdr3, from = database_TRB$CDR3, to = database_TRB$Epitope.species, warn_missing = F)
HCA_TRB$Antigen[!HCA_TRB$Antigen %in% database_TRB$Epitope.species] <- "No"
HCA_TRB_type <- HCA_TRB$Antigen %>% str_split(pattern = ";")
HCA_TRB[lapply(HCA_TRB_type, length) %>% do.call(rbind, .) %>% c != 1, "Antigen"] <- "Multple"


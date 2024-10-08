packages <- list("Seurat", "Matrix", "stringr", "stringi", "gridExtra", "ggplot2",
								 "plyr", "ggthemes", "cowplot", "data.table", "RColorBrewer", "ComplexHeatmap",
								 "circlize", "pheatmap", "reshape2", "scales", "rlang", "future",
								 "parallel", "ggrepel", "ggsci", "ape", "dplyr", "harmony",
								 "sigminer", "ggpubr", "rstatix", "ROGUE", "tidyverse", "reticulate", "dendextend",
								 "scRepertoire", "ggdendro", "ggtern", "viridis", "scatterpie", "SeuratData", "SeuratDisk",
								 "scibet", "scuttle", "zellkonverter")
s
lapply(packages, library, character.only = TRUE)
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]
options(width = 240)

options(future.globals.maxSize = 400*1000 * 1024^2)
plan("multisession", workers = 20)
plan()

###---------change the cell barcodes-------------###
for (i in dir() %>% grep(pattern = "_TCR_filtered_contig_annotations.csv", value = T)) {
	sample_name <- i %>% gsub(pattern = "_TCR_filtered_contig_annotations.csv", replacement = "")
	tmp <- read.table(i, header = T, row.names = NULL, sep = ",", stringsAsFactors = F)
	tmp$barcode <- paste0(sample_name, "_", tmp$barcode)
	tmp$contig_id <- paste0(sample_name, "_", tmp$contig_id)
	write.table(tmp, file = paste0(sample_name, "_TCR_filtered_contig_annotations.csv"), row.names = F, col.names = T, sep = ",", quote = F)
}

###cat * > all_TCR_filtered_contig_annotations_cellranger6.tsv
###sed -n '1p' all_TCR_filtered_contig_annotations_cellranger6.tsv > header
###grep -Ev "barcode" all_TCR_filtered_contig_annotations_cellranger6.tsv > tmp
###cat header tmp >all_TCR_filtered_contig_annotations_cellranger6.tsv

#####--------------- read T cell metadat and TCR clone types information-----------------###########
T_cell_clone_raw <- read.table("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/TNK/Cells/TCR_rawdata_cellranger6/all_TCR_filtered_contig_annotations_cellranger6.tsv",
															 sep = ",", header = T, row.names = NULL, stringsAsFactors = F)

T_cells_meta.data <- read.table("/public/home/HeShuai/Projects/1.HCA2.0/2.Analysis/All_cells/All_sample_contamination_removed/TNK/Cells/20230404_TNK_meta.data.txt",
																header = T, row.names = NULL, stringsAsFactors = F, sep = "\t", comment.char = "", encoding = "UTF-8")

overlapped_T_cell_clone <- T_cell_clone_raw %>% subset(barcode %in% T_cells_meta.data$row.names ) %>% select(barcode, chain, v_gene, j_gene, c_gene, cdr3)
##-------------------remove the cells with only one chain---------#################
T_cell_clone <- overlapped_T_cell_clone
multiple_chains <- T_cell_clone[duplicated(T_cell_clone$barcode), ] %>% `[`(, "barcode") %>% unique()
T_cell_clone <- T_cell_clone[T_cell_clone$barcode %in% multiple_chains, ]

barcodes_cells <- unique(T_cell_clone$barcode)

library(parallel)
cl <- makeCluster(40)
clusterExport(cl = cl, varlist = c("barcodes_cells", "T_cell_clone"))
clusterEvalQ(cl, library(dplyr))
qualited_T_cells <- pblapply(cl = cl, X = T_cell_clone, FUN = function(x){
	tmp <- T_cell_clone[T_cell_clone$barcode %in% x, "chain"] %>% unlist
	if(all("TRB" %in% tmp, "TRA" %in% tmp)){
		return(x)
	}
}) %>% do.call(rbind, .)
stopCluster(cl)

T_cell_clone_uniq <- T_cell_clone[T_cell_clone$barcode %in% qualited_T_cells, ] ##
T_cell_clone_uniq <- subset(T_cell_clone_uniq, c_gene != "")

###-------------------remove the cells with only one chain---------#################
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, "barcode") %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] ##262625 qualited T cells
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))

# select(T_cell_clone_uniq, c(raw_barcode, chain, unique_clone)) %>% write.table("TCR_info.txt", sep = "\t", col.names = T, row.names = F) ## save the TCR informations.
TCR_info <- select(T_cell_clone_uniq, c(barcode, chain, unique_clone))
# TCR_info <- read.table("TCR_info.txt", header = T, row.names = NULL, sep = "\t", stringsAsFactors = F)

library(parallel)
cl <- makeCluster(30) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))

customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
	tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
	clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
	# dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
	# dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### ??һ????????һ??????һs????
	dat <- c(x, clonetype_customer)
}) %>% do.call(rbind, .)
stopCluster(cl)

colnames(customer_clone) <- c("cellbarcode", "clonetype_customer")
customer_clone <- customer_clone %>% as.data.frame()
# customer_clone %>% write.table("customer_clone.txt", sep = "\t", col.names = T, row.names = T)
# customer_clone <- read.table("customer_clone.txt", sep = "\t", header = T, row.names = NULL, stringsAsFactors = F)

T_cell_clone_uniq$customer_clone <- mapvalues(T_cell_clone_uniq$barcode, from = customer_clone$cellbarcode, to = customer_clone$clonetype_customer, warn_missing = F)
meta.dat <- read.table("Sample_information_of_HCA2.0_H_to_Donor.csv", header = T, row.names = NULL, sep = ",", stringsAsFactors = F)
T_cell_clone_uniq$sample <- mapvalues(T_cell_clone_uniq$barcode %>% str_split(pattern = "_", simplify = T) %>% `[`(, 1), from = meta.dat$Tissue, to = meta.dat$Donor, warn_missing = F)
T_cell_clone_uniq$customer_clone <- paste0(T_cell_clone_uniq$customer_clone, "_", T_cell_clone_uniq$sample)

unique_TCR <- T_cell_clone_uniq %>% select(barcode, customer_clone) %>% unique
multiple_clones <- (c((unique_TCR$customer_clone %>% table) > 1) %>% names)[c(unique_TCR$customer_clone %>% table) > 1]
T_cell_clone_uniq$multiple_clones <- mapvalues(T_cell_clone_uniq$customer_clone, from = multiple_clones, to = rep("Multiple", time = length(multiple_clones)), warn_missing = FALSE)
T_cell_clone_uniq$multiple_clones <- gsub(T_cell_clone_uniq$multiple_clones, pattern = "^TR.*", replacement = "Single")

unique_TCR_pri <- T_cell_clone_uniq %>% select(barcode, customer_clone) %>% unique
unique_TCR_pri$customer_clone <- paste0(unique_TCR_pri$customer_clone) ###customer_clone: TRA+TRB+(both CDR3)
unique_TCR_pri <- unique_TCR_pri %>% unique
unique_TCR_pri$Progect <- "HCA2.0"

dats <- data.frame()
x <- 1
for(variable in unique_TCR_pri$Progect %>% unique){
	unique_TCR <- subset(unique_TCR_pri, Progect == "HCA2.0")
	Two_Ten <- c(table(unique_TCR$customer_clone) %>% names)[table(unique_TCR$customer_clone) >= 2 & table(unique_TCR$customer_clone) < 10]
	Ten_fifty <- c(table(unique_TCR$customer_clone) %>% names)[table(unique_TCR$customer_clone) >= 10 & table(unique_TCR$customer_clone) < 50]
	Fifty_100 <- c(table(unique_TCR$customer_clone) %>% names)[table(unique_TCR$customer_clone) >= 50 & table(unique_TCR$customer_clone) < 100]
	one100_200 <- c(table(unique_TCR$customer_clone) %>% names)[table(unique_TCR$customer_clone) >= 100 & table(unique_TCR$customer_clone) < 200]
	Two200_500 <- c(table(unique_TCR$customer_clone) %>% names)[table(unique_TCR$customer_clone) >= 200 & table(unique_TCR$customer_clone) < 500]
	Two500_1000 <- c(table(unique_TCR$customer_clone) %>% names)[table(unique_TCR$customer_clone) >= 500 & table(unique_TCR$customer_clone) < 1000]
	Two1000 <- c(table(unique_TCR$customer_clone) %>% names)[table(unique_TCR$customer_clone) >= 1000]
	
	df1 <- data.frame(customer_clones = c(Two_Ten, Ten_fifty, Fifty_100, one100_200, Two200_500, Two500_1000, Two1000),
										state = rep(c("2-10", "10-50", "50-100", "100-200", "200-500", "500-1000", ">1000"), times = c(sapply(list(Two_Ten, Ten_fifty, Fifty_100, one100_200, Two200_500, Two500_1000, Two1000), length))))
	unique_TCR$clone_state <- mapvalues(unique_TCR$customer_clone, from = df1$customer_clones, to = df1$state, warn_missing = F)
	unique_TCR$clone_state <- gsub(unique_TCR$clone_state, pattern = "^TR.*", replacement = "Single")
	dats <- rbind(unique_TCR, dats)
	print(variable)
	print(x)
	x <- x + 1
	rm(df1, unique_TCR, Two_Ten, Ten_fifty, Fifty_100, one100_200, Two200_500, Two500_1000, Two1000)
}

T_cell_clone_uniq$clone_state <- mapvalues(T_cell_clone_uniq$customer_clone, from = dats$customer_clone, to = dats$clone_state, warn_missing = F)
# T_cell_clone_uniq$clone_state <- gsub(T_cell_clone_uniq$clone_state, pattern = "2-10|10-50|50-100|100-200|>200", replacement = "Multiple")
T_cell_clone_uniq$orig.ident <- mapvalues(T_cell_clone_uniq$barcode, from = T_cells_meta.data$row.names, to = T_cells_meta.data$orig.ident, warn_missing = F)
T_cell_clone_uniq$nCount_RNA <- mapvalues(T_cell_clone_uniq$barcode, from = T_cells_meta.data$row.names, to = T_cells_meta.data$nCount_RNA, warn_missing = F)
T_cell_clone_uniq$Name <- mapvalues(T_cell_clone_uniq$barcode, from = T_cells_meta.data$row.names, to = T_cells_meta.data$Name, warn_missing = F)
T_cell_clone_uniq$Donor <- mapvalues(T_cell_clone_uniq$barcode, from = T_cells_meta.data$row.names, to = T_cells_meta.data$Donor, warn_missing = F)
T_cell_clone_uniq$Stage <- mapvalues(T_cell_clone_uniq$barcode, from = T_cells_meta.data$row.names, to = T_cells_meta.data$Stage, warn_missing = F)

T_cell_clone_uniq$seurat_clusters <- mapvalues(T_cell_clone_uniq$barcode, from = T_cells_meta.data$row.names, to = T_cells_meta.data$seurat_clusters, warn_missing = F)
T_cell_clone_uniq$Global_annotation <- mapvalues(T_cell_clone_uniq$barcode, from = T_cells_meta.data$row.names, to = T_cells_meta.data$Global_annotation, warn_missing = F)
T_cell_clone_uniq$Major_cluster <- mapvalues(T_cell_clone_uniq$barcode, from = T_cells_meta.data$row.names, to = T_cells_meta.data$Major_cluster, warn_missing = F)
T_cell_clone_uniq$Samples <- mapvalues(T_cell_clone_uniq$barcode, from = T_cells_meta.data$row.names, to = T_cells_meta.data$Samples, warn_missing = F)
T_cell_clone_uniq$Tmp_annotation <- mapvalues(T_cell_clone_uniq$barcode, from = T_cells_meta.data$row.names, to = T_cells_meta.data$Tmp_annotation, warn_missing = F)
T_cell_clone_uniq$Detail_Annotation_final <- mapvalues(T_cell_clone_uniq$barcode, from = T_cells_meta.data$row.names, to = T_cells_meta.data$Detail_Annotation_final, warn_missing = F)

T_cell_clone_uniq %>% write.table(file = "Merged_TCR.information_at_least_one_paired_VDJ_cellranger6.txt", sep = "\t", col.names = T, row.names = T)

T_cells_meta.data$cloned_Tcells <- plyr::mapvalues(T_cells_meta.data$row.names, from = T_cell_clone_uniq$barcode, to = T_cell_clone_uniq$clone_state, warn_missing = F)
T_cells_meta.data$multiple_clones <- mapvalues(T_cells_meta.data$row.names, from = T_cell_clone_uniq$barcode, to = T_cell_clone_uniq$multiple_clones, warn_missing = F)

T_cells_meta.data$multiple_clones <- gsub(T_cells_meta.data$multiple_clones, pattern = ".*-.*", replacement = "Not")
T_cells_meta.data$cloned_Tcells <- gsub(T_cells_meta.data$cloned_Tcells, pattern = ".*_.*", replacement = "Not")

T_cells_meta.data$customer_clone <- mapvalues(T_cells_meta.data$row.names, from = T_cell_clone_uniq$barcode, to = T_cell_clone_uniq$customer_clone, warn_missing = F)
T_cells_meta.data$customer_clone <- gsub(T_cells_meta.data$customer_clone, pattern = ".*-1$", replacement = "Not")

T_cells_meta.data %>% write.table("20230404_TNK_meta.data_with_TCR_information_.txt", col.names = T, row.names = T, sep = "\t")


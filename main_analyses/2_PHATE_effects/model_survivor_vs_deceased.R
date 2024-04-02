library(limma)
library(edgeR)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(data.table)

current = getwd()
folder = "2_PHATE_effects"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

## inputs
residuals_dir <- "inputs/"
permutations_dir <- "inputs/permutations/"

## outputs
results_dir <- paste0("outputs/",folder,"/")

## source function
setwd("common_functions")
source("permFDR.R")
setwd(current)

## read in meta data from individuals
md <- fread(paste0("inputs/meta_data.csv"), data.table = FALSE)
md <- subset(md, select = c(Assigned_Name, Center_LocalCode, DECEASED_less60, Center, Age, Sex, BMI, FlowCell_ID, PHATE_cluster, DSO11, bulk_RNASeq, CBC_PC1, CBC_PC2))
md <- subset(md, DSO11 == 1 & bulk_RNASeq == 1)
md$Assigned_Name <- factor(md$Assigned_Name)
md$Center <- factor(md$Center)
md$DSO11 <- factor(md$DSO11)
md$bulk_RNASeq <- factor(md$bulk_RNASeq)
md$Sex <- factor(md$Sex)
md$FlowCell_ID <- factor(md$FlowCell_ID)
md$PHATE_cluster <- factor(md$PHATE_cluster, levels = c("1","2","3","4"))
md$DECEASED_less60 <- factor(md$DECEASED_less60, levels = c("1","0"))
rownames(md) <- md$Assigned_Name

## fill in missing BMI with average and scale numeric covariates
unique_indivs <- subset(md, select = c(Center_LocalCode, BMI))
unique_indivs <- unique_indivs[!duplicated(unique_indivs$Center_LocalCode), ]
md$BMI <- ifelse(is.na(md$BMI), mean(unique_indivs$BMI, na.rm = TRUE), md$BMI)
md$BMI_scale <- scale(md$BMI)
md$Age_scale <- scale(md$Age)

## subset on ONLY CLUSTER 1 individuals
md <- subset(md, PHATE_cluster == 1)

## contrasts 
contrasts <- c("survivor_vs_deceased")
clusters1 <- c("1")
clusters2 <- c("0")

for(z in 1:length(contrasts)){

	contrast_z <- contrasts[z]
	cluster1 <- clusters1[z]
	cluster2 <- clusters2[z]

	## subset metadata on this
	meta_data_i <- md
	meta_data_i <- subset(meta_data_i, DECEASED_less60 %in% c(cluster1, cluster2))
	meta_data_i$DECEASED_less60 <- factor(meta_data_i$DECEASED_less60, levels = c(cluster1, cluster2))

	# read in permutations
	permutations <- read.table(paste0(permutations_dir,"permutations_PHATE_",contrast_z,".txt"), header = TRUE, sep = ",")

	## read in corrected expression
	residuals <- read.table(paste0(residuals_dir,"corrected_expression_Center.txt"), header = TRUE, sep = ",")
	residuals <- residuals[colnames(residuals) %in% meta_data_i$Assigned_Name]

	## read in weights
	weights <- read.table(paste0(residuals_dir,"weights.txt"), header = TRUE, sep = ",")
	weights <- weights[colnames(weights) %in% meta_data_i$Assigned_Name]

	reorder_names <- rownames(meta_data_i)
	if(length(reorder_names) == dim(residuals)[2]){
		residuals <- residuals[reorder_names]
		weights <- weights[reorder_names]
	}else{
		correct_names <- colnames(residuals)
		meta_data_i <- meta_data_i[rownames(meta_data_i) %in% correct_names,]
		weights <- weights[,colnames(weights) %in% correct_names]
		reorder_names <- rownames(meta_data_i)
		residuals <- residuals[reorder_names]
		weights <- weights[reorder_names]
	}
	length(which(colnames(residuals)!=rownames(meta_data_i)))

	## model infection differential expression
	design = model.matrix(~ Age_scale + BMI_scale + Sex + FlowCell_ID + CBC_PC1 + CBC_PC2 + DECEASED_less60, data = meta_data_i)
	design <- design[, colSums(design != 0) > 0]

	vfit <- lmFit(residuals, weights = weights, design)
	vfit <- eBayes(vfit)

	## collect outputs
	betas = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% paste0("DECEASED_less60",cluster2))]); colnames(betas)[1] <- "betas"
	p_values = as.data.frame(vfit$p.value[, which(colnames(vfit$coefficients) %in% paste0("DECEASED_less60",cluster2))]); colnames(p_values)[1] <- "pvalues"
	fdrs = as.data.frame(p.adjust(p_values[,1], method = "BH")); colnames(fdrs)[1] <- "fdrs"
	t_stats = as.data.frame(cbind(rownames(topTable(vfit, coef = paste0("DECEASED_less60",cluster2), number = Inf)), topTable(vfit, coef = paste0("DECEASED_less60",cluster2), number = Inf)$t))
	colnames(t_stats) <- c("genes","t_stat")

	results <- cbind(betas, p_values, fdrs)
	results$genes <- rownames(results)
	results <- join(results, t_stats, by = "genes")

	## match row names for permutations
	perms <- permutations[match(results$genes, permutations$genes),]
	length(which(results$genes!=perms$genes))
	rownames(results) <- results$genes; results$genes <- NULL
	rownames(perms) <- perms$genes; perms$genes <- NULL
		
	## calculate qvalues using a permutation-based null
	perm_fdrs = permFDR(full_data = results, full_column_id = "pvalues", perm_data = perms, perm_column_ids = "all", output_name = paste0(results_dir))
	write.table(perm_fdrs$fdrs, file = paste0(results_dir,"results_with_qvalues_PHATE_",contrast_z,".txt"), quote = FALSE, sep = ",")

	corr_fdrs <- subset(perm_fdrs$fdrs)

  	numGenes <- dim(results)[1]
  	numDEgenes_10 <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.10))
 	numDEgenes_05 <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.05))
  	numDEgenes_01 <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.01))

	numDEgenes_per_CT <- as.data.frame(cbind(numDEgenes_10, numDEgenes_05, numDEgenes_01, numGenes))
	colnames(numDEgenes_per_CT) <- c("numDEgenes_fdr10","numDEgenes_fdr05","numDEgenes_fdr01","num_genes_tested")
	write.table(numDEgenes_per_CT, paste0(results_dir,"DEG_qvalues_PHATE_",contrast_z,".txt"), quote = FALSE)
}


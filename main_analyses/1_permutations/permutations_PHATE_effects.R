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
folder = "1_permutations"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

## inputs
residuals_dir <- paste0("inputs/")
num_permutations <- c(1:100)
set.seed(2022)

## outputs
results_dir <- paste0("outputs/",folder,"/")

## read in meta data from individuals
md <- fread(paste0("inputs/meta_data.csv"), data.table = FALSE)
md <- subset(md, select = c(Assigned_Name, Center_LocalCode, Center, Age, Sex, BMI, FlowCell_ID, PHATE_cluster, DSO11, bulk_RNASeq, CBC_PC1, CBC_PC2))
md <- subset(md, DSO11 == 1 & bulk_RNASeq == 1)
md$Assigned_Name <- factor(md$Assigned_Name)
md$Center <- factor(md$Center)
md$DSO11 <- factor(md$DSO11)
md$bulk_RNASeq <- factor(md$bulk_RNASeq)
md$Sex <- factor(md$Sex)
md$FlowCell_ID <- factor(md$FlowCell_ID)
md$PHATE_cluster <- factor(md$PHATE_cluster, levels = c("1","2","3","4"))
rownames(md) <- md$Assigned_Name

## fill in missing BMI with average and scale numeric covariates
unique_indivs <- subset(md, select = c(Center_LocalCode, BMI))
unique_indivs <- unique_indivs[!duplicated(unique_indivs$Center_LocalCode), ]
md$BMI <- ifelse(is.na(md$BMI), mean(unique_indivs$BMI, na.rm = TRUE), md$BMI)
md$BMI_scale <- scale(md$BMI)
md$Age_scale <- scale(md$Age)

## contrasts 
contrasts <- c("1_vs_2","3_vs_4","1_vs_3","1_vs_4","2_vs_3","2_vs_4")
clusters1 <- c("2","4","3","4","3","4")
clusters2 <- c("1","3","1","1","2","2")

for(z in 1:length(contrasts)){

	contrast_z <- contrasts[z]
	cluster1 <- clusters1[z]
	cluster2 <- clusters2[z]

	## subset metadata on this
	meta_data_i <- md
	meta_data_i <- subset(meta_data_i, PHATE_cluster %in% c(cluster1, cluster2))
	meta_data_i$PHATE_cluster <- factor(meta_data_i$PHATE_cluster, levels = c(cluster1, cluster2))

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

	for (i in 1:length(num_permutations)){

		perm_number <- num_permutations[i]

		## PERMUTE SEVERITY ACROSS INDIVIDUALS
		meta_data_PERM <- meta_data_i
		meta_data_PERM$perm_PHATE <- sample(meta_data_PERM$PHATE_cluster)

		reorder_names <- rownames(meta_data_PERM)
		residuals <- residuals[reorder_names]
		weights <- weights[reorder_names]

		length(which(colnames(residuals)!=rownames(meta_data_PERM)))
		length(which(colnames(weights)!=rownames(meta_data_PERM)))
		length(which(colnames(residuals)!=colnames(weights)))

		## differential expression
		design = model.matrix(~ Age_scale + BMI_scale + Sex + FlowCell_ID + CBC_PC1 + CBC_PC2 + perm_PHATE, data = meta_data_PERM)
		## remove columns that are all 0s
		design <- design[, colSums(design != 0) > 0]

		vfit <- lmFit(residuals, weights = weights, design)
		vfit <- eBayes(vfit)

		pvalues = as.data.frame(vfit$p.value[, which(colnames(vfit$coefficients) %in% c(paste0("perm_PHATE",cluster2)))]); colnames(pvalues)[1] <- paste0("pvalues_perm",perm_number)
		pvalues$genes <- rownames(pvalues)

		if(i == 1){
			results <- pvalues
		}else{
			results <- merge(results, pvalues, by = "genes")
		}

		if(i %% 10 == 0){print(i)}
	}

	write.table(results, paste0(results_dir,"permutations_PHATE_",contrast_z,".txt"), quote = FALSE, sep = ",", row.names = FALSE)
	print(contrast_z)
}


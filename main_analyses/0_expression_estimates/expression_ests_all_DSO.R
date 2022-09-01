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
folder = "0_expression_estimates"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

## outputs
results_dir <- paste0("outputs/",folder,"/")

## read in meta data from individuals
md <- fread(paste0("inputs/meta_data.csv"), data.table = FALSE)
md <- subset(md, select = c(Assigned_Name, Center_LocalCode, Center, Age, Sex, BMI, FlowCell_ID, PHATE_cluster, DSO11, bulk_RNASeq, CBC_PC1, CBC_PC2))
md <- subset(md, bulk_RNASeq == 1)
md$Assigned_Name <- factor(md$Assigned_Name)
md$Center <- factor(md$Center)
md$DSO11 <- factor(md$DSO11)
md$bulk_RNASeq <- factor(md$bulk_RNASeq)
md$Sex <- factor(md$Sex)
md$FlowCell_ID <- factor(md$FlowCell_ID)
md$PHATE_cluster <- factor(md$PHATE_cluster, levels = c("1","2","3","4"))

## fill in missing BMI with average and scale numeric covariates
unique_indivs <- subset(md, select = c(Center_LocalCode, BMI))
unique_indivs <- unique_indivs[!duplicated(unique_indivs$Center_LocalCode), ]
md$BMI <- ifelse(is.na(md$BMI), mean(unique_indivs$BMI, na.rm = TRUE), md$BMI)
md$BMI_scale <- scale(md$BMI)
md$Age_scale <- scale(md$Age)

## read in expression data
reads <- read.table(paste0("inputs/raw_counts.txt"), header = TRUE)
reads <- reads[!duplicated(reads$genes), ]
rownames(reads) <- reads$genes
reads$genes <- NULL

## reorder
rownames(md) <- md$Assigned_Name
reorder_names <- rownames(md)
if(length(reorder_names) == dim(reads)[2]){
	reads <- reads[reorder_names]
}else{
	correct_names <- colnames(reads)
	md <- md[rownames(md) %in% correct_names,]
	reorder_names <- rownames(md)
	reads <- reads[,colnames(reads) %in% reorder_names]
	reads <- reads[reorder_names]
}
length(which(colnames(reads)!=rownames(md)))

## subset on protein-coding genes only
coding_ids <- read.table("inputs/coding_genes.txt")
reads <- reads[which(rownames(reads) %in% coding_ids$V2), ]

## remove lowly-expressed genes
dge <- DGEList(counts = reads)
dge <- calcNormFactors(dge)
design = model.matrix(~ Age_scale + BMI_scale + Sex + FlowCell_ID + Center + CBC_PC1 + CBC_PC2, data = md)
## remove columns that are all 0s
design <- design[, colSums(design != 0) > 0]
v <- voom(dge, design, plot = TRUE)

tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
tab = tab[order(-tab$medians), ]
tab$order_by_median = 1:nrow(tab)
tab = tab[order(tab$order), ]
reads <- reads[which(tab$medians > 1), ]

## voom after removal of lowly expressed genes
dge <- DGEList(counts = reads)
dge <- calcNormFactors(dge)
design <- model.matrix(~ Age_scale + BMI_scale + Sex + FlowCell_ID + Center + CBC_PC1 + CBC_PC2, data = md)
v <- voom(dge, design, plot = TRUE)

## collect
expression <- v$E
weights <- v$weights
colnames(weights) <- colnames(expression)
rownames(weights) <- rownames(expression)

## write uncorrected expression
write.table(expression, paste0(results_dir,"uncorrected_expression.txt"), quote = FALSE, sep = ",")

## write weights
write.table(weights, paste0(results_dir,"weights.txt"), quote = FALSE, sep = ",")

## correct for covariates to get GE matrix for PCA, plotting
design <- model.matrix(~ Age_scale + BMI_scale + Sex + FlowCell_ID + Center + CBC_PC1 + CBC_PC2, data = md)
v <- voom(dge, design, plot = FALSE)
vfit <-lmFit(v, design)
vfit <- eBayes(vfit)
corrected_expression <- v$E - vfit$coefficients[,"Age_scale"]%*%t(design[,"Age_scale"]) - vfit$coefficients[,"BMI_scale"]%*%t(design[,"BMI_scale"]) - vfit$coefficients[,"SexM"]%*%t(design[,"SexM"]) - vfit$coefficients[,"FlowCell_ID1557"]%*%t(design[,"FlowCell_ID1557"]) - vfit$coefficients[,"FlowCell_ID1576"]%*%t(design[,"FlowCell_ID1576"]) - vfit$coefficients[,"FlowCell_ID1708"]%*%t(design[,"FlowCell_ID1708"]) - vfit$coefficients[,"FlowCell_ID1773"]%*%t(design[,"FlowCell_ID1773"]) - vfit$coefficients[,"FlowCell_ID1802"]%*%t(design[,"FlowCell_ID1802"]) - vfit$coefficients[,"FlowCell_ID1851"]%*%t(design[,"FlowCell_ID1851"]) - vfit$coefficients[,"CenterJGH"]%*%t(design[,"CenterJGH"]) - vfit$coefficients[,"CBC_PC1"]%*%t(design[,"CBC_PC1"]) - vfit$coefficients[,"CBC_PC2"]%*%t(design[,"CBC_PC2"])
write.table(corrected_expression, paste0(results_dir,"corrected_expression_Age_BMI_Sex_Flowcell_Center_CBC.txt"), quote = FALSE, sep = ",")

corrected_expression_Center <- v$E - vfit$coefficients[,"CenterJGH"]%*%t(design[,"CenterJGH"])
write.table(corrected_expression_Center, paste0(results_dir,"corrected_expression_Center.txt"), quote = FALSE, sep = ",")


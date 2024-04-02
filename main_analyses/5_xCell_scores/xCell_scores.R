library(dplyr)
library(plyr)
library(reshape2)
library(tidyr)
library(data.table)
library(gridExtra)
library(xCell)
library(data.table)

current = getwd()
folder = "5_xCell_scores"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

## inputs
residuals_dir <- paste0("inputs/")

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
rownames(md) <- md$Assigned_Name

## read in exp
exp <- as.data.frame(read.table(paste0(residuals_dir,"/corrected_expression_Age_BMI_Sex_Flowcell_Center_CBC.txt"), header = TRUE, sep = ","))

correct_names <- colnames(exp)
md <- md[rownames(md) %in% correct_names,]
reorder_names <- rownames(md)
exp <- exp[reorder_names]
length(which(colnames(exp)!=as.character(md$Assigned_Name)))

## function
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

## calculate xcell scores
xcell <- xCellAnalysis(exp, cell.types.use = c("B-cells", "Basophils", "CD4+ memory T-cells", "CD4+ naive T-cells", "CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem", "CD8+ naive T-cells", "CD8+ T-cells", "CD8+ Tcm", "CD8+ Tem", "Class-switched memory B-cells", "pDC", "Eosinophils", "Mast cells", "Megakaryocytes", "Memory B-cells", "Monocytes", "naive B-cells", "Neutrophils", "NK cells", "NKT", "Plasma cells", "Platelets", "Th1 cells", "Th2 cells", "Tregs"), parallel.sz = 1)

## QN scores
xcell_QN <- as.data.frame(quantile_normalisation(t(xcell)))
xcell_QN$Assigned_Name <- rownames(xcell_QN)
xcell_QN <- xcell_QN %>% select(Assigned_Name, everything())
write.table(xcell_QN, paste0(results_dir,"xCell_scores_QN.txt"), quote = FALSE, sep = ",", row.names = FALSE)


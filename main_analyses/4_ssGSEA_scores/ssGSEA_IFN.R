library(ggplot2)
library(dplyr)
library(plyr)
library(reshape2)
library(tidyr)
library(data.table)
library(gridExtra)
library(cowplot)
library(mgsub)
library(GSVA)
library(qusage)

current = getwd()
folder = "4_ssGSEA_scores"
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

## IFN genes
hallmark_IFN_gamma <- as.character(unlist(read.gmt(paste0("inputs/genelists/hallmark_IFN_gamma.gmt"))))
hallmark_IFN_alpha <- as.character(unlist(read.gmt(paste0("inputs/genelists/hallmark_IFN_alpha.gmt"))))
IFN <- unique(c(hallmark_IFN_gamma, hallmark_IFN_alpha))
IFN_gene_list <- list(IFN)
names(IFN_gene_list) <- "IFN"

## read in
exp <- as.data.frame(read.table(paste0(residuals_dir,"/corrected_expression_Age_BMI_Sex_Flowcell_Center_CBC.txt"), header = TRUE, sep = ",", check.names = FALSE))

correct_names <- colnames(exp)
md <- md[rownames(md) %in% correct_names,]
reorder_names <- rownames(md)
exp <- exp[reorder_names]
length(which(colnames(exp)!=as.character(md$Assigned_Name)))

## subset expression matrices 
exp_gvsa <- exp

scores <- gsva(as.matrix(exp_gvsa),
				gset.idx.list = IFN_gene_list,
				method = "ssgsea",
				verbose = FALSE,
				parallel.sz = 1)

scores <- as.data.frame(t(scores))
colnames(scores) <- paste0("IFN_score")
scores$Assigned_Name <- rownames(scores)
scores <- scores[,c(2,1)]

## write combined scores
write.table(scores, paste0(results_dir,"GSVA_IFN_score.txt"), sep = ",", row.names = FALSE, quote = FALSE)

## process for merge
md <- subset(md, select = c(Assigned_Name, PHATE_cluster, DSO11))
md <- subset(md, DSO11 == 1)

## plot
df <- join(md, scores, by = "Assigned_Name", type = "inner")
df$PHATE_cluster <- factor(df$PHATE_cluster, levels = c("1","2","3","4"))

plot <- ggplot(df, aes(x = PHATE_cluster, y = IFN_score, color = PHATE_cluster, fill = PHATE_cluster)) +
			geom_boxplot(alpha = 0.3, width = 0.5, outlier.shape = NA) +
			geom_jitter(width = 0.2, alpha = 0.7) +
			theme_bw() +
			scale_color_manual(values = c("#DC267F","#FFB000","#7395EC","#4C2DDC")) +
			scale_fill_manual(values = c("#DC267F","#FFB000","#7395EC","#4C2DDC")) +
			theme(panel.border = element_rect(colour = "black", fill = NA, size = .75), 
				  plot.title = element_text(size = 10), 
				  panel.grid.major = element_blank(), 
				  panel.grid.minor = element_blank(),
				  axis.title.x = element_text(size = 10),
				  axis.title.y = element_text(size = 10),
				  legend.position = "none") +
			xlab("PHATE cluster") +
			ylab("ssGSEA IFN score") +
			ggtitle(paste0("bulk whole blood, DSO 11"))

pdf(paste0(results_dir,"by_cluster_IFN_score.pdf"), width = 4, height = 3.5)
print(plot)
dev.off()

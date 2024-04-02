library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(data.table)
library(fgsea)
library(qusage)
library(RColorBrewer)
library(ggplot2)

current = getwd()
folder = "3_GSEA"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

## inputs
set.seed(2022)
residuals_dir <- paste0("inputs/")
DE_dir <- paste0("inputs/DE_results/")
genesets_dir <- paste0("inputs/genelists/")
celltypes <- c("bulk_whole_blood")

## outputs
results_dir <- paste0("outputs/",folder,"/")

## input genesets
hallmark <- read.gmt(paste0(genesets_dir,"/h.all.v7.1.symbols.gmt"))

plot_GSEA <- function(geneset, geneset_name, contrast){

	## directories
	OUTPUTS_dir <- paste0(results_dir,"/enrichments/")
	system(paste0("mkdir -p ", OUTPUTS_dir,"/fgsea_heatmaps/"))
	system(paste0("mkdir -p ", OUTPUTS_dir,"/fgsea_tables/"))
	system(paste0("mkdir -p ", OUTPUTS_dir,"/fgsea_tables/",geneset_name,"/"))
	
	## read in tstats
	readIn <- function(celltype){
		df <- read.table(paste0(DE_dir,"results_with_qvalues_PHATE_",contrast,".txt"), header = TRUE, sep = ",", check.names = FALSE)
		df$genes <- rownames(df)
		df <- subset(df, select = c(genes, t_stat))
		colnames(df)[2] <- paste0(celltype)
		return(df)
	}

	blood <- readIn("blood")

	## combine
	t_statistics <- blood
	rownames(t_statistics) <- t_statistics$genes; t_statistics$genes <- NULL

	## function to calculate gsea
	for(i in 1:ncol(t_statistics)){

		name <- colnames(t_statistics)[i]

		## take the t stats of condition i and rank them
		t_i <- subset(t_statistics, select = name); colnames(t_i)[1] <- "t_stat"
		t_i$gene <- rownames(t_i)
		t_i <- data.table(t_i[,c(2,1)])
		t_i_rank <- t_i[order(t_stat), ]

		input <- setNames(t_i_rank$t_stat, t_i_rank$gene)

		fgseaRes <- fgsea(pathways = geneset, 
		                  stats = input,
		                  minSize=15,
		                  maxSize=500,
		                  nperm=100000)

		fgseaRes_csv <- data.frame(lapply(fgseaRes, as.character), stringsAsFactors = FALSE)
		write.csv(fgseaRes_csv, file = paste0(OUTPUTS_dir,"/fgsea_tables/",geneset_name,"/",geneset_name,"_",name,"_",contrast,".csv"))

		## format for downstream
		colnames(fgseaRes) <- paste0(colnames(fgseaRes),"_",name)
		fgseaRes <- fgseaRes[,c(1,3,5)]
		colnames(fgseaRes)[1] <- "pathway"
		assign(paste0(name), fgseaRes)
	}

	## combine dfs 
	merged <- blood

	## subset pathways by threshold adjusted pval
	if(geneset_name == "hallmark"){
		thresh = 1
	}
	
	merged_subset <- as.data.frame(merged[which(-log10(merged$padj_blood) > thresh),])
	merged_subset_m <- merged_subset

	merged_subset_m$pathway <- gsub("HALLMARK_", "", merged_subset_m$pathway)
	merged_subset_m$pathway <- gsub("_"," ",merged_subset_m$pathway)
	merged_subset_m$pathway <- tolower(merged_subset_m$pathway)

	## plot
	base_size = 4.5
	breaks_pval = c(signif(min(-log10(merged_subset_m$padj_blood)), 2), signif(max(-log10(merged_subset_m$padj_blood)), 2))
	limits_pval = c(signif(min(-log10(merged_subset_m$padj_blood)), 2), signif(max(-log10(merged_subset_m$padj_blood)), 2))
	labels_pval = c(signif(min(-log10(merged_subset_m$padj_blood)), 2), signif(max(-log10(merged_subset_m$padj_blood)), 2))

	breaks_NES = c(signif(min(merged_subset_m$NES_blood), 1)-0.05, signif(max(merged_subset_m$NES_blood), 1)+0.05)
	limits_NES = c(signif(min(merged_subset_m$NES_blood), 1)-0.05, signif(max(merged_subset_m$NES_blood), 1)+0.05)
	labels_NES = c(signif(min(merged_subset_m$NES_blood), 1)-0.05, signif(max(merged_subset_m$NES_blood), 1)+0.05)

	## color dots that are not significant gray
	merged_subset_m$value_NES_withNA <- ifelse(merged_subset_m$padj_blood < 0.10, merged_subset_m$NES_blood, merged_subset_m$NES_blood == NA)
	merged_subset_m$CT <- paste0("blood")

	plot <- ggplot(merged_subset_m, aes(x = CT, y = pathway, fill = value_NES_withNA, size = -log10(padj_blood), color = CT)) + 
			geom_point(shape = 21, stroke = 0.8) +
			scale_color_manual(values = rep("black", 7), name = NULL, guide = FALSE) +
	   		scale_size(name = "-log10(padj)") +
	   		scale_fill_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title = "enrichment score", frame.linewidth = 0.5, ticks = FALSE), breaks = breaks_NES, limits = limits_NES, labels = labels_NES, na.value = "grey75") +
			guides(colour = FALSE) + 
			theme_bw(base_size = base_size)+
		  	coord_equal(ratio=0.8) +
		  	theme(legend.position = "bottom",
		  		legend.box = "vertical",
		        legend.text = element_text(size=base_size*1.5),
		        legend.title = element_text(size=base_size*1.5),
		        axis.ticks.x = element_blank(),
		        axis.ticks.y = element_blank(),
		        axis.text.x = element_text(size = base_size*1.5, angle = 45, hjust = 1, vjust = 1, colour = "black"),
		        axis.text.y = element_text(size = base_size*1.5),
		        panel.border = element_blank(), 
		        panel.grid.minor = element_blank(),
		        panel.grid.major = element_blank(),
		        legend.key.height = unit(0.5, "cm"),
		        legend.key = element_rect(fill = "white"),
		        legend.key.width = unit(0.6, "cm"),
		        plot.title = element_text(size = 10)) +
		  	ggtitle(paste0(geneset_name)) + 
		  	xlab("") +
		  	ylab("")

	if(geneset_name == "hallmark"){
			pdf(paste0(OUTPUTS_dir,"/fgsea_heatmaps/",geneset_name,"_PHATE_",contrast,".pdf"), width = 6, height = 5, useDingbats=FALSE)
		}
	print(plot)
	dev.off()
}

plot_GSEA(hallmark, "hallmark", "3_vs_4")
plot_GSEA(hallmark, "hallmark", "1_vs_2")
plot_GSEA(hallmark, "hallmark", "survivor_vs_deceased")


#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
inputrds <- args[1]
outputdir<-args[2]
group<-args[3]
#inputrds<-"/cluster/home/liuji/T1D/multiome/analysis/Cyt49_S4D5/scRNAseq/Clusters_Combined.QC_filtered.normalized.rds"
#inputrds<-"/cluster/home/liuji/T1D/multiome/analysis/Cyt49_S4D5/clustering/Cyt49_S4D5.sce.rd"
#group<-"Clusters_Combined"
library(pathfindR)
library(pheatmap)
library(scater)
library(scran)
library(SingleCellExperiment)
library(scuttle) # For logNormCounts
library(dplyr)
library(glmGamPoi)
library(reshape2)


sce<-readRDS(inputrds)
#is.mito <- any(seqnames(location)=="MT")
#rownames(sce)[is.mito]
#colData(sce)$orig.ident
sampleids<-unique(colData(sce)[[group]])
###### marker detection
# AUC represents the probability that a randomly chosen observation from our cluster of interest is greater than a randomly chosen observation from the other cluster.
# AUC: 1- upregulation, 0.5 - no change, 0 - downregulation
# Cohen’s d: number of standard deviations that separate the means of the two groups
# Cohen’s d: positive - upregulated, negative - downregulation, zero - little difference
# logFC.detected: log-fold change in the proportion of cells with detected expression between clusters

#function to identify pathways in the gene lists
pathfind<-function(prefix,selected.vgs){
  #group<-"S5D3_R01"
  #selected.vgs<-data.frame(top.ranked.rowdata)
  selected.vgs<-data.frame(selected.vgs)
  output_selected.vgs <- run_pathfindR(selected.vgs,gene_sets="Reactome",min_gset_size=3)
  #enrichment_chart(result_df = output_selected.vgs,top_terms = 15)
  #head(output_selected.vgs,3)
  input_processed <- input_processing(selected.vgs,p_val_threshold = 0.05)
  #gg_list <- visualize_terms(result_df = output_topranked,input_processed = input_processed)  # this function returns a list of ggraph objects (named by Term ID)
  output_topranked_clustered <- cluster_enriched_terms(output_selected.vgs, plot_dend = FALSE, plot_clusters_graph = FALSE)
  #head(output_topranked_clustered)
  #output_topranked_clustered[output_topranked_clustered$ID=="hsa04512",]
  selected_clusters <- output_topranked_clustered[output_topranked_clustered$Status=="Representative",]
  selected_pathway<-file.path(outputdir,paste(prefix,"mean.logFC.pathway","png",sep="."))
  png(selected_pathway, width=10, height=12, units="in", res=500)
  p<-enrichment_chart(selected_clusters, plot_by_cluster = TRUE)
  print(p)
  dev.off()
  outfile<-file.path(outputdir, paste(prefix, "mean.logFC.pathway","csv",sep="."))
  write.csv(selected_clusters, outfile,row.names = FALSE)
  
  # gg_list <- visualize_terms(result_df = output_selected.vgs,input_processed = input_processed)  # this function returns a list of ggraph objects (named by Term ID)
  
}

marker.info <- scoreMarkers(sce, sce[[group]])
## plot the markers 
#outputdir<-"/Users/liuji/Downloads/scmultiome/SEM01_S4D5/analysis/marker"
for(cluster in unique(sce[[group]])){
  #print(group)
  #group<-"C4"
  chosen <- marker.info[[cluster]]
  #ordered <- chosen[order(chosen$median.logFC.cohen,decreasing=TRUE),]
  ordered <- chosen[order(chosen$mean.AUC,decreasing=TRUE),]
  #dim(ordered[ordered$mean.AUC>0.5,])
  #head(ordered[,1:4]) # showing basic stats only, for brevity.
  top.ranked <- head(ordered,70)
  
  hvgs_expression_fig1<-file.path(outputdir,paste(group,cluster,"mean.AUC.heatmap","png",sep="."))
  png(hvgs_expression_fig1, width=6, height=12, units="in", res=500)
  fig1<-plotGroupedHeatmap(sce, features=rownames(top.ranked), group=group, center=TRUE, zlim=c(-3, 3))
  print(fig1)
  dev.off()
  
  hvgs_expression_fig2<-file.path(outputdir,paste(group,cluster,"mean.AUC.violin","png",sep="."))
  png(hvgs_expression_fig2, width=12, height=10, units="in", res=500)
  genes<-head(rownames(ordered),24)
  fig2<-plotExpression(sce, features=genes, ncol = 3,
                       x=group, colour_by=group,point_size=0.01)+  
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) 
  print(fig2)
  dev.off()
  
}

## pathway analysis
for(cluster in unique(sce[[group]])){
  #print(cluster)
  #cluster<-"C1"
  chosen <- marker.info[[cluster]]
  #ordered <- chosen[order(chosen$median.logFC.cohen,decreasing=TRUE),]
  ordered <- chosen[order(chosen$mean.AUC,decreasing=TRUE),]
  #head(ordered,10)
  top.ranked.rowdata<-ordered[ordered$mean.logFC.detected>0.1 & ordered$mean.AUC>0.5 & ordered$mean.logFC.cohen >0,
                              c("mean.logFC.detected","rank.logFC.detected")]
  top.ranked.rowdata<-head(top.ranked.rowdata,200)
  colnames(top.ranked.rowdata) <- c("FC", "rank.logFC.detected")
  top.ranked.rowdata$ADJ.PVAL <- 0.01
  top.ranked.rowdata[] <- lapply(top.ranked.rowdata, as.numeric)
  top.ranked.rowdata$Gene.symbol<-rownames(top.ranked.rowdata)
  top.ranked.rowdata<-data.frame(top.ranked.rowdata[,c("Gene.symbol","FC","ADJ.PVAL")])
  pathfind(paste(group,cluster,sep="."),top.ranked.rowdata)
}


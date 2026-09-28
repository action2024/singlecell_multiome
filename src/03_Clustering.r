#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
inputrds<-args[1]
outputdir<-args[2]
prefix<-args[3]

#source ~/miniforge3/etc/profile.d/conda.sh
#conda activate multiome

library(ArchR)
library(Signac)
library(Seurat)
#library(EnsDb.Hsapiens.v86)
library(stringr)
library(scater)
library(devtools)
library(ggplot2) # for plotting
library(reshape2)
library(msigdbr) # for gathering gene sets  #- NOT INSTALLED
#library(SeuratData)
# library(future) # for parallel computing
# library(future.apply) # for parallel computing
library(dplyr)
library(patchwork)
library("Matrix")
library("readr")
library(scuttle)
library(robustbase)
library(scran)
library(scRNAseq)
library(cluster)
library(dendextend)
library(AUCell)
library(clustree)
library(PCAtools)
library(celldex)
library(SC3)
library("ggalluvial")
library(future)  # for parallel computing
library(future.apply)  # for parallel computing
#install.packages("anticlust", version = "0.6.1")
library(SingleCellExperiment)
library(DropletUtils) #- NOT INSTALLED
library(scMerge)
#library(pathfindR)
library(pheatmap)
library(bluster)

source("/home/l128405/multiome/src/R/variables/colors.R")
source("/home/l128405/multiome/src/R/variables/devgenes.R")
source("/home/l128405/multiome/src/R/variables/celltypegenes.R")
source("/home/l128405/multiome/src/R/functions/violin_plot_qc.R")
source("/home/l128405/multiome/src/R/functions/LSI_dimplot.R")
source("/home/l128405/multiome/src/R/functions/stacked_barplot.R")
source("/home/l128405/multiome/src/R/functions/markerdetect.R")
source("/home/l128405/multiome/src/R/functions/reducedimplot.R")


# archRdirpath <- c("/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen/S39")
# outputdir<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen/S39"
# group<-"Clusters_Combined"
# prefix<-"S39"
# inputrds<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen/S39/S39.normalized.sce.rds"
# 
# archRdirpath <- c("/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen/S40")
# outputdir<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen_v2_05152026/S40_clust"
# group<-"Clusters_Combined"
# prefix<-"S40"
# inputrds<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen_v2_05152026/S40_qc/S40.normalized.sce.rds"

dir.create(outputdir)

setwd(outputdir)
#projMulti_filtered <- loadArchRProject(path = archRdirpath)
#cellIDs <- BiocGenerics::which(!(projMulti_filtered$Sample %in% c("S38A","S39B")))
#projSubset <- subsetArchRProject(ArchRProj = projMulti_filtered,cells = projMulti_filtered$cellNames[cellIDs],outputDirectory = outputdir_filtered,dropCells = TRUE,force = TRUE)

# remove mitocondria prior to clustering
sce<-readRDS(inputrds)
#sce<-sce[,which(!(sce$Sample %in% c("S38A","S39B")))]
sce<-sce[which(!grepl("^mt-", rownames(sce))),]

#top variable genes
dec <- modelGeneVar(sce)
hvg.var <- getTopHVGs(dec, n=2500)
AAVgenes<-rownames(sce)[grepl("transcript_type", rownames(sce))]
genelist<-unique(genelist[genelist%in% rownames(sce)])
hvg.var<-unique(c(hvg.var,AAVgenes,genelist))
#str(hvg.var)
sce <- runPCA(sce, subset_row=hvg.var)
#dim(reducedDim(sce, "PCA"))
percent.var <- attr(reducedDim(sce), "percentVar")
#####plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
#nn.clusters <- clusterCells(sce, use.dimred="PCA")
#table(nn.clusters)
#nn.clust.info <- clusterCells(sce, use.dimred="PCA", full=TRUE)
#nn.clust.info$objects$graph
#reducedDim(sce, "force") <- igraph::layout_with_fr(nn.clust.info$objects$graph)

#out <- RunHarmony(sce, group.by.vars = "")
#stopifnot(all.equal(colnames(sce), colnames(out)))
#reducedDim(sce, "harmony") <- reducedDim(out, "HARMONY")


nn.clusters <- clusterCells(sce, use.dimred="PCA", 
                            BLUSPARAM=SNNGraphParam(k=10, type="rank", cluster.fun="walktrap"))
colLabels(sce) <- factor(nn.clusters)
#table(nn.clusters)
sce <- runTSNE(sce, dimred="PCA")
sce <- runUMAP(sce, dimred="PCA")

colData(sce)$MergedSample <- str_sub(as.data.frame(colData(sce))$Sample, end = -2)
## append col data from sce to archR
samplesheet_file<-"/home/l128405/multiome/src/bash/samplesheet.txt"
samplesheet<-read.csv(samplesheet_file, sep = "\t",header=TRUE)
#colData(sce)<-merge(as.data.frame(colData(sce)),samplesheet,by="MergedSample",all= TRUE)
matched_rows <- match(sce$MergedSample, samplesheet$MergedSample)
# 2. Map and add the target column from dataframe2 into your experiment object
# Replace 'target_col' with the column name you want to bring over
colData(sce)$treatment <- samplesheet$treatment[matched_rows]
colData(sce)$time <- samplesheet$time[matched_rows]
colData(sce)$batch <- samplesheet$batch[matched_rows]

library(harmony)
library(BiocSingular)
set.seed(230616)
batchlist<-unique(colData(sce)$batch)
if(length(batchlist)>1){
  out <- RunHarmony(sce, group.by.vars = "batch")
  stopifnot(all.equal(colnames(sce), colnames(out)))
  reducedDim(sce, "harmony") <- reducedDim(out, "HARMONY")
  
  nn.clusters <- clusterCells(sce, use.dimred="harmony", 
                              BLUSPARAM=SNNGraphParam(k=10, type="rank", cluster.fun="walktrap"))
  sce$label_harmony <- factor(nn.clusters)
  
  sce <- runTSNE(sce, dimred="harmony", name = "TSNE_harmony")
  sce <- runUMAP(sce, dimred="harmony", name = "UMAP_harmony")}

#save clustering sce
clust.rds<-file.path(outputdir,paste(prefix,"clust","rds",sep="."))
saveRDS(sce, clust.rds)

#hclust.dyn <- clusterCells(sce, use.dimred="PCA",BLUSPARAM=HclustParam(method="ward.D2", cut.dynamic=TRUE, 
#                                                                           cut.params=list(minClusterSize=10, deepSplit=1)),full=TRUE)
#table(hclust.dyn)
#colLabels(sce) <- factor(hclust.dyn$clusters)
#sce$hclust <- factor(hclust.dyn$clusters)
#sce$hclust <- factor(hclust.dyn$clusters)
#quantify the proportion of clusters
df_clusters<-as.data.frame(table(nn.clusters))
colnames(df_clusters)<-c("cluster","cellnum")
df_clusters$cellprop<-round(df_clusters$cellnum/sum(df_clusters$cellnum)*100,1)
outfile<-file.path(outputdir, paste(prefix, "clust.num","csv",sep="."))
write.csv(df_clusters, outfile,row.names = FALSE,quote=FALSE)

#table(colData(sce)$Clusters_RNA)
if(length(batchlist)>1){group<-"label_harmony"}else{group<-"label"}
cellclusters<-as.data.frame(colData(sce)[,c("Clusters_Combined",group)])
cellclusters$Clusters_Combined <- sub("^C", "", cellclusters$Clusters_Combined)
#cellclusters$Clusters_Combined <- as.numeric(as.character(cellclusters$Clusters_Combined))
cellclusterscount<-cellclusters %>% group_by_all() %>% summarise(COUNT = n())

p1<-stacked_barplot(cellclusterscount,group,"COUNT","Clusters_Combined",group)
p2<-stacked_barplot(cellclusterscount,"Clusters_Combined","COUNT",group,"Clusters_Combined")

plotPDF(p1, p2,name = "scRNA_scMultiome_RNA-Clusters_Combined-stackedbar", addDOC = FALSE)

#table(colData(sce)$Clusters_RNA)
cellclusters<-as.data.frame(colData(sce)[,c("Clusters_RNA",group)])
cellclusters$Clusters_RNA <- sub("^C", "", cellclusters$Clusters_RNA)
#cellclusters$Clusters_Combined <- as.numeric(as.character(cellclusters$Clusters_Combined))
cellclusterscount<-cellclusters %>% group_by_all() %>% summarise(COUNT = n())

p1<-stacked_barplot(cellclusterscount,group,"COUNT","Clusters_RNA",group)
p2<-stacked_barplot(cellclusterscount,"Clusters_RNA","COUNT",group,"Clusters_RNA")

plotPDF(p1, p2,name = "scRNA_scMultiome_RNA-Clusters_RNA-stackedbar", addDOC = FALSE)

cellclusters<-as.data.frame(colData(sce)[,c("MergedSample",group)])
cellclusterscount<-cellclusters %>% group_by_all() %>% summarise(COUNT = n())
p1<-stacked_barplot(cellclusterscount,group,"COUNT","MergedSample",group)
p2<-stacked_barplot(cellclusterscount,"MergedSample","COUNT",group,"MergedSample")

plotPDF(p1, p2,name = "scRNA_label-MergedSample-stackedbar", addDOC = FALSE)


colData(sce)$Sample <- as.data.frame(colData(sce))$Sample

pca_plot<-reducedimplot(sce,"PCA",group,group,3)
pca_ident_plot<-reducedimplot(sce,"PCA","Sample","Sample",3)
#harmony_plot<-reducedimplot(sce,"harmony","label","label",3)
#harmony_ident_plot<-reducedimplot(sce,"harmony","Sample","Sample",3)
tsne_plot<-reducedimplot(sce,"TSNE",group,group,2)
tsne_ident_plot<-reducedimplot(sce,"TSNE","Sample","Sample",2)
umap_plot<-reducedimplot(sce,"UMAP",group,group,2)
umap_ident_plot<-reducedimplot(sce,"UMAP","Sample","Sample",2)
plotPDF(pca_plot,pca_ident_plot,tsne_plot,tsne_ident_plot, umap_plot,umap_ident_plot, name = "scRNA-PCA_TSNE_UMAP", addDOC = FALSE)
#plotPDF(pca_plot, tsne_plot,umap_plot,name = "scRNA-PCA_TSNE_UMAP", addDOC = FALSE)

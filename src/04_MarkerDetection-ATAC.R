#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
inputrds<-args[1]
outputdir<-args[2]
prefix<-args[3]
### plot expression level of selected gene sets and identify markers of all clusters### 
#source ~/miniforge3/etc/profile.d/conda.sh
#conda activate multiome
#outputdir<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen_v3/S13/clust"
#inputrds<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen_v3/S13/clust/S13.clust.rds"
#prefix<-"S13"

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
source("/home/l128405/multiome/src/R/functions/sctypeanno.R")


# outputdir<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen/S39"
# group<-"Clusters_Combined"
# prefix<-"S39"
# inputrds<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen/S39/S39.normalized.sce.rds"
# 
# outputdir<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen/S40"
# group<-"Clusters_Combined"
# prefix<-"S26-S27"
# inputrds<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen_v2_05152026/S26-S27_clust/S26-S27.label.anno.rds"

dir.create(outputdir)
dir.create(file.path(outputdir,"Plots"))

setwd(outputdir)
#cellIDs <- BiocGenerics::which(!(projMulti_filtered$Sample %in% c("S38A","S39B")))
#projSubset <- subsetArchRProject(ArchRProj = projMulti_filtered,cells = projMulti_filtered$cellNames[cellIDs],outputDirectory = outputdir_filtered,dropCells = TRUE,force = TRUE)

# remove mitocondria prior to clustering
sce<-readRDS(inputrds)
library(gridExtra)
genelist<-unique(genelist[genelist%in% rownames(sce)])
# plot heatmap plot
if(any(str_detect(colnames(colData(sce)), "label"))){
  genelist_heatmap<-plotGroupedHeatmap(sce, features=genelist, group="label", 
                                       center=TRUE, zlim=c(-3, 3)) 
  
  heatmap_fig<-file.path(outputdir,"Plots",paste(prefix,"label","scRNA.GE.heatmap","pdf",sep="."))
  pdf(heatmap_fig, width=9, height=12,0)
  print(genelist_heatmap)
  dev.off()
  
  devgenelist<-unique(devgenelist[devgenelist%in% rownames(sce)])
  genelist_heatmap<-plotGroupedHeatmap(sce, features=devgenelist, group="label", 
                                       center=TRUE, zlim=c(-3, 3)) 
  
  heatmap_fig<-file.path(outputdir,"Plots",paste(prefix,"label","scRNA.devGE.heatmap","pdf",sep="."))
  pdf(heatmap_fig, width=9, height=12,0)
  print(genelist_heatmap)
  dev.off()
}
if(any(str_detect(colnames(colData(sce)), "label_harmony"))){
  genelist_heatmap<-plotGroupedHeatmap(sce, features=genelist, group="label_harmony", 
                                       center=TRUE, zlim=c(-3, 3)) 
  
  heatmap_fig<-file.path(outputdir,"Plots",paste(prefix,"label_harmony","scRNA.GE.heatmap","pdf",sep="."))
  pdf(heatmap_fig, width=9, height=12,0)
  print(genelist_heatmap)
  dev.off()
  
  devgenelist<-unique(devgenelist[devgenelist%in% rownames(sce)])
  genelist_heatmap<-plotGroupedHeatmap(sce, features=devgenelist, group="label_harmony", 
                                       center=TRUE, zlim=c(-3, 3)) 
  
  heatmap_fig<-file.path(outputdir,"Plots",paste(prefix,"label_harmony","scRNA.devGE.heatmap","pdf",sep="."))
  pdf(heatmap_fig, width=9, height=12,0)
  print(genelist_heatmap)
  dev.off()
}


# plot heatmap for multiome: Clusters_Combined
genelist_heatmap<-plotGroupedHeatmap(sce, features=genelist, group="Clusters_Combined", 
                                     center=TRUE, zlim=c(-3, 3)) 

heatmap_fig<-file.path(outputdir,"Plots",paste(prefix,"multiome.Clusters_Combined.GE.heatmap","pdf",sep="."))
pdf(heatmap_fig, width=9, height=12,0)
print(genelist_heatmap)
dev.off()

# plot heatmap for multiome: Clusters_RNA
genelist_heatmap<-plotGroupedHeatmap(sce, features=genelist, group="Clusters_RNA", 
                                     center=TRUE, zlim=c(-3, 3)) 

heatmap_fig<-file.path(outputdir,"Plots",paste(prefix,"multiome.Clusters_RNA.GE.heatmap","pdf",sep="."))
pdf(heatmap_fig, width=9, height=12,0)
print(genelist_heatmap)
dev.off()

pca_plot<-reducedimplot(sce,"PCA","label","label",2)
tsne_plot<-reducedimplot(sce,"TSNE","label","label",2)
samplenum<-length(genelist)

if(any(str_detect(reducedDimNames(sce), "harmony"))){
  harmony_plot<-reducedimplot(sce,"harmony","label_harmony","label_harmony",2)
  tsne_harmony_plot<-reducedimplot(sce,"TSNE_harmony","label_harmony","label_harmony",2)
  sample_plots<-list(tsne_harmony_plot)
  
  for(gene in genelist){
    #gene<-'DLK1'
    #print(gene)
    p<-plotReducedDim(sce, "TSNE_harmony", 
                      colour_by=gene,point_size=0.002) 
    sample_plots<-c(sample_plots,list(p))
  }
  geneexpression_fig<-file.path(outputdir,"Plots",paste(prefix,"GE.TSNE_harmony","pdf",sep="."))
  pdf(geneexpression_fig, width=3*round(sqrt(samplenum),0)+9, height=3*round(sqrt(samplenum),0))
  grid.arrange(grobs = sample_plots, 3,0)## display plot
  dev.off()
}

sample_plots<-list(pca_plot)
for(gene in genelist){
  #gene<-'DLK1'
  #print(gene)
  p<-plotReducedDim(sce, "PCA", 
                    colour_by=gene,point_size=0.002) 
  sample_plots<-c(sample_plots,list(p))
}

geneexpression_fig<-file.path(outputdir,"Plots",paste(prefix,"GE.PCA","pdf",sep="."))
pdf(geneexpression_fig, width=3*round(sqrt(samplenum),0)+9, height=3*round(sqrt(samplenum),0))
grid.arrange(grobs = sample_plots, 3,0)## display plot
dev.off()

sample_plots<-list(tsne_plot)
for(gene in genelist){
  #gene<-'DLK1'
  #print(gene)
  p<-plotReducedDim(sce, "TSNE", 
                    colour_by=gene,point_size=0.002) 
  sample_plots<-c(sample_plots,list(p))
}

geneexpression_fig<-file.path(outputdir,"Plots",paste(prefix,"GE.TSNE","pdf",sep="."))
pdf(geneexpression_fig, width=3*round(sqrt(samplenum),0)+9, height=3*round(sqrt(samplenum),0))
grid.arrange(grobs = sample_plots, 3,0)## display plot
dev.off()







if(any(str_detect(colnames(colData(sce)), "label"))){
  violin<-plotExpression(sce, features=genelist,
                         x=I(colData(sce)$label),color_by =I(colData(sce)$label),  ncol = 3,point_size=0.01)  + 
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    facet_wrap(~Feature, scales = "free_y") + 
    scale_color_manual(values=customcolor)
  geneexpression_fig<-file.path(outputdir,"Plots",paste(prefix,"label.GE.violin","pdf",sep="."))
  pdf(geneexpression_fig, width=3*round(sqrt(samplenum),0)+9, height=3*round(sqrt(samplenum),0))
  print(violin)
  dev.off()}

if(any(str_detect(colnames(colData(sce)), "label_harmony"))){
  violin<-plotExpression(sce, features=devgenelist,
                         x=I(colData(sce)$label_harmony),color_by =I(colData(sce)$label_harmony),  ncol = 3,point_size=0.01)  + 
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    facet_wrap(~Feature, scales = "free_y") + 
    scale_color_manual(values=customcolor)
  geneexpression_fig<-file.path(outputdir,"Plots",paste(prefix,"label_harmony.GE.violin","pdf",sep="."))
  pdf(geneexpression_fig, width=3*round(sqrt(samplenum),0)+9, height=3*round(sqrt(samplenum),0))
  print(violin)
  dev.off()
}
###### marker detection
# AUC represents the probability that a randomly chosen observation from our cluster of interest is greater than a randomly chosen observation from the other cluster.
# AUC: 1- upregulation, 0.5 - no change, 0 - downregulation
# Cohen’s d: number of standard deviations that separate the means of the two groups
# Cohen’s d: positive - upregulated, negative - downregulation, zero - little difference
# logFC.detected: log-fold change in the proportion of cells with detected expression between clusters
#inputrds<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen/S40/S40.clust.rds"
#sce<-readRDS(inputrds)

if(any(str_detect(colnames(colData(sce)), "label"))){markerdetect(sce,"label")}
if(any(str_detect(colnames(colData(sce)), "label_harmony"))){markerdetect(sce,"label_harmony")}

#markerdetect(sce,"Clusters_Combined")

counts_matrix <- assays(sce)$logcounts
genelist_GE <- t(as.data.frame(counts_matrix[genelist, ]))
genelist_GE <- cbind( t(as.data.frame(counts_matrix[genelist, ])),as.data.frame(colData(sce)))

grouplist<-c("label","label_harmony")
for(group in grouplist){
  ## count AAV by group
  if(any(str_detect(colnames(colData(sce)),group))){
    #genelist_GE[genelist_GE['3mir-hATOH1-hPOU4F3-3Flag-3mir&#8221;; transcript_type ']>0,]
    aavcount1<-as.data.frame(table(genelist_GE[genelist_GE["3mir-hATOH1-hPOU4F3-3Flag-3mir&#8221;; transcript_type "]>0,][[group]]))
    colnames(aavcount1) <-c("cluster","hATOH1-hPOU4F3")
    aav_gfi_counts<-genelist_GE[genelist_GE["hGFI1-HA&#8221;; transcript_type "]>0 | genelist_GE["hGFI1-HA-mScarlet&#8221;; transcript_type "]>0 ,]
    aavcount2<-as.data.frame(table(aav_gfi_counts[[group]]))
    colnames(aavcount2) <-c("cluster","hGFI1")
    aavcount<- merge(aavcount1, aavcount2, by = "cluster", all= TRUE)
    #genelist_GE[genelist_GE['hGFI1-HA&#8221;; transcript_type ']>0 | genelist_GE['3mir-hATOH1-hPOU4F3-3Flag-3mir&#8221;; transcript_type ']>0,]
    aav1_IDs<-rownames(genelist_GE[genelist_GE["3mir-hATOH1-hPOU4F3-3Flag-3mir&#8221;; transcript_type "]>0,])
    aav2_IDs1<-rownames(genelist_GE[genelist_GE["hGFI1-HA&#8221;; transcript_type "]>0,])
    aav2_IDs2<-rownames(genelist_GE[genelist_GE["hGFI1-HA-mScarlet&#8221;; transcript_type "]>0,])
    aav2_IDs<-unique(c(aav2_IDs1,aav2_IDs2))
    aav_count_bygroup_file<-file.path(outputdir,paste(prefix,group,"scRNA.aav_count","csv",sep="."))
    aavcount %>%
      write.table(file = aav_count_bygroup_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)
    aav1_cellIDs_file<-file.path(outputdir,paste(prefix,group,"hATOH1-hPOU4F3.cellIDs","txt",sep="."))
    aav1_IDs %>% write.table(aav1_cellIDs_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=FALSE)
    aav2_cellIDs_file<-file.path(outputdir,paste(prefix,group,"hGfi1.cellIDs","txt",sep="."))
    aav2_IDs %>% write.table(aav1_cellIDs_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=FALSE)
    
    ### annotation of clusters
    
    sctype_scores<-sctypeanno(sce,group,db_,"P20 sc-Cochlea")
    
    sctype_anno_file<-file.path(outputdir,paste(prefix,group,"clusters.anno","csv",sep="."))
    sctype_scores %>% write.csv(sctype_anno_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)
    SC_cluster<-as.data.frame(sctype_scores[sctype_scores$type %in% c("CC-IS-OS","Pillar cells"),]$cluster)
    HC_cluster<-as.data.frame(sctype_scores[sctype_scores$type %in% c("Outer hair cells","Inner hair cells"),]$cluster)
    colnames(SC_cluster)<-"cluster"
    colnames(HC_cluster)<-"cluster"
    HCSC_cluster<-rbind(SC_cluster,HC_cluster)
    SC_cluster_file<-file.path(outputdir,paste(prefix,group,"SC.clusters","txt",sep="."))
    SC_cluster %>% write.csv(SC_cluster_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)
    HC_cluster_file<-file.path(outputdir,paste(prefix,group,"HC.clusters","txt",sep="."))
    HC_cluster %>% write.csv(HC_cluster_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)
    HCSC_cluster_file<-file.path(outputdir,paste(prefix,group,"HCSC.clusters","txt",sep="."))
    HCSC_cluster %>% write.csv(HCSC_cluster_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)
    
    SC_cellIDs<-rownames(colData(sce)[colData(sce)[[group]] %in% SC_cluster$cluster,])
    HC_cellIDs<-rownames(colData(sce)[colData(sce)[[group]] %in% HC_cluster$cluster,])
    HCSC_cellIDs<-rownames(colData(sce)[colData(sce)[[group]] %in% HCSC_cluster$cluster,])
    
    SC_cellIDs_file<-file.path(outputdir,paste(prefix,group,"SC.cellIDs","txt",sep="."))
    SC_cellIDs %>% write.table(SC_cellIDs_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=FALSE)
    HC_cellIDs_file<-file.path(outputdir,paste(prefix,group,"HC.cellIDs","txt",sep="."))
    HC_cellIDs %>% write.table(HC_cellIDs_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=FALSE)
    HCSC_cellIDs_file<-file.path(outputdir,paste(prefix,group,"HCSC.cellIDs","txt",sep="."))
    HCSC_cellIDs %>% write.table(HCSC_cellIDs_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=FALSE)
    
    if(group=="label"){
      # Transfer label from sctype_scores to sce
      match_idx <- match(sce$label, sctype_scores$cluster)
      colData(sce)$celltype <- sctype_scores$type[match_idx]
      #Save clustering sce
      clust.rds<-file.path(outputdir,paste(prefix,"label.anno","rds",sep="."))
      saveRDS(sce, clust.rds)
    }
    
    if(group=="label_harmony"){
      # Transfer label from sctype_scores to sce
      match_idx <- match(sce$label_harmony, sctype_scores$cluster)
      colData(sce)$celltype <- sctype_scores$type[match_idx]
      #Save clustering sce
      clust.rds<-file.path(outputdir,paste(prefix,"label_harmony.anno","rds",sep="."))
      saveRDS(sce, clust.rds)
    }
    sctype_scores_aavcount<- merge(as.data.frame(sctype_scores), aavcount, by = "cluster", all= TRUE) %>%
      arrange(type)
    
    aav_count_byanno_file<-file.path(outputdir,paste(prefix,group,"clusters.aav_count.anno","csv",sep="."))
    sctype_scores_aavcount %>% write.table(aav_count_byanno_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)
    
    sctype_scores<-sctypeanno(sce,group,customdb_,"Cochlear")
    sctype_anno_file<-file.path(outputdir,paste(prefix,group,"clusters.custom_anno","csv",sep="."))
    sctype_scores %>% write.csv(sctype_anno_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)
    
    
  }}

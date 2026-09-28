#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
inputrds<-args[1]
selectedgroup<-args[2]
dimreduc_method<-args[3]
selected_clust_file<-args[4]
outputdir<-args[5]
prefix<-args[6]

set.seed(230616)

#library(ArchR)
#library(Signac)
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

set.seed(230616)

#outputdir<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen_v3/S12/clust/S12_HCSC"
#inputrds<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen_v3/S12/clust/S12.label.anno.rds"
#prefix<-"S12"
#selected_clust_file<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen_v3/S12/clust/S12.label.HCSC.clusters.txt"
#selectedgroup<-"label"
#dimreduc_method<-"PCA"
dir.create(outputdir)
setwd(outputdir)

sce<-readRDS(inputrds)
library(scDblFinder)

selected_clust<-read.csv(selected_clust_file)$cluster
sce_scRNA_selected_cells<-colnames(sce[,sce[[selectedgroup]] %in% selected_clust])
sce_selected<-sce[,sce[[selectedgroup]] %in% selected_clust]

cellIDfile<-file.path(outputdir,paste(prefix,selectedgroup,"cellIDs","txt",sep="."))
colnames(sce) %>%  write.table(file = cellIDfile, sep = "\t", quote=FALSE,row.names = FALSE, col.names=FALSE)

sce<-sce_selected
#top variable genes
dec <- modelGeneVar(sce)
hvg.var <- getTopHVGs(dec, n=2000)
#str(hvg.var)
sce <- runPCA(sce, subset_row=hvg.var)
percent.var <- attr(reducedDim(sce), "percentVar")

nn.clusters <- clusterCells(sce, use.dimred=dimreduc_method, 
                            BLUSPARAM=SNNGraphParam(k=5, type="rank", cluster.fun="walktrap"))
#table(nn.clusters)
sce$reclust <- factor(nn.clusters)
sce <- runTSNE(sce, dimred=dimreduc_method, name=paste("TSNE",dimreduc_method,sep="_"))
sce <- runUMAP(sce, dimred=dimreduc_method, name=paste("UMAP",dimreduc_method,sep="_"))
hclust.dyn <- clusterCells(sce, use.dimred=dimreduc_method,BLUSPARAM=HclustParam(method="ward.D2", cut.dynamic=TRUE, 
                                                                                 cut.params=list(minClusterSize=5, deepSplit=1)),full=TRUE)
#table(hclust.dyn)
sce$hclust <- factor(hclust.dyn$clusters)

## append col data from sce to archR
samplesheet_file<-"/home/l128405/polyA/src/bash/REGEN-AK059.samplesheet.txt"
samplesheet<-read.csv(samplesheet_file, sep = "\t",header=TRUE)
#colData(sce)<-merge(as.data.frame(colData(sce)),samplesheet,by="MergedSample",all= TRUE)
matched_rows <- match(sce$MergedSample, samplesheet$MergedSample)
# 2. Map and add the target column from dataframe2 into your experiment object
# Replace 'target_col' with the column name you want to bring over
colData(sce)$treatment <- samplesheet$treatment[matched_rows]
colData(sce)$time <- samplesheet$time[matched_rows]
colData(sce)$batch <- samplesheet$batch[matched_rows]

library("openxlsx")
library(dplyr)
library(Seurat)
library(HGNChelper)
library(scater)
library(scran)
library(ArchR)

sctype_scores<-sctypeanno(sce,"reclust",db_,"P20 sc-Cochlea")
sctype_anno_file<-file.path(outputdir,paste(prefix,"reclust","clusters.anno","csv",sep="."))
sctype_scores %>% write.csv(sctype_anno_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)
# Transfer label from sctype_scores to sce
match_idx <- match(sce$reclust, sctype_scores$cluster)
colData(sce)$celltype <- sctype_scores$type[match_idx]
#save clustering sce
clust.rds<-file.path(outputdir,paste(prefix,"reclust","rds",sep="."))
saveRDS(sce, clust.rds)

# count AAV in each group
counts_matrix <- assays(sce)$logcounts
genelist<-unique(genelist[genelist%in% rownames(sce)])
genelist_GE <- t(as.data.frame(counts_matrix[genelist, ]))
genelist_GE <- merge( t(as.data.frame(counts_matrix[genelist, ])),as.data.frame(colData(sce)),by=0)
df_hATOH1_hPOU4F3<-genelist_GE[genelist_GE['3mir-hATOH1-hPOU4F3-3Flag-3mir&#8221;; transcript_type ']>0,]
df_hGFI1<-genelist_GE[genelist_GE['hGFI1-HA&#8221;; transcript_type ']>0 | genelist_GE["hGFI1-HA-mScarlet&#8221;; transcript_type "]>0,]
colData(sce)$AAV<-"NA"
overlap_hATOH1_hPOU4F3_hGFI1<-intersect(df_hATOH1_hPOU4F3$Row.names, df_hGFI1$Row.names)
unique_hATOH1_hPOU4F3<-setdiff(df_hATOH1_hPOU4F3$Row.names, df_hGFI1$Row.names)
unique_hGFI1<-setdiff(df_hGFI1$Row.names,df_hATOH1_hPOU4F3$Row.names)
if(length(overlap_hATOH1_hPOU4F3_hGFI1)>0){colData(sce)[colnames(sce) %in% overlap_hATOH1_hPOU4F3_hGFI1,]$AAV<-"hATOH1_hPOU4F3_hGFI1"}
if(length(unique_hATOH1_hPOU4F3)>0){colData(sce)[colnames(sce) %in% unique_hATOH1_hPOU4F3,]$AAV<-"hATOH1_hPOU4F3"}
if(length(unique_hGFI1)>0){colData(sce)[colnames(sce) %in% unique_hGFI1,]$AAV<-"hGFI1"}

aavcount1<-as.data.frame(table(df_hATOH1_hPOU4F3[["reclust"]]))
colnames(aavcount1) <-c("cluster","hATOH1-hPOU4F3")
aavcount2<-as.data.frame(table(df_hGFI1[["reclust"]]))
colnames(aavcount2) <-c("cluster","hGFI1")
aavcount<- merge(aavcount1, aavcount2, by = "cluster", all= TRUE)
#genelist_GE[genelist_GE['hGFI1-HA&#8221;; transcript_type ']>0 | genelist_GE['3mir-hATOH1-hPOU4F3-3Flag-3mir&#8221;; transcript_type ']>0,]
sctype_scores_aavcount<- merge(as.data.frame(sctype_scores), aavcount, by = "cluster", all= TRUE) %>%
  arrange(type)

aav_count_byanno_file<-file.path(outputdir,paste(prefix,"reclust","scRNA.aav_count.anno","csv",sep="."))
sctype_scores_aavcount %>% write.table(aav_count_byanno_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)

SC_cluster<-as.data.frame(sctype_scores[sctype_scores$type %in% c("CC-IS-OS","Pillar cells","IBC-IPhC-HenSC","Deiter's cells"),]$cluster)
HC_cluster<-as.data.frame(sctype_scores[sctype_scores$type %in% c("Outer hair cells","Inner hair cells"),]$cluster)
colnames(SC_cluster)<-"cluster"
colnames(HC_cluster)<-"cluster"
HCSC_cluster<-rbind(SC_cluster,HC_cluster)
group<-"reclust"
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



## Plot distribution of clusters in stacked bar plot format
cellclusters<-as.data.frame(colData(sce)[,c("reclust","MergedSample")])
#cellclusters$Clusters_Combined <- as.numeric(as.character(cellclusters$Clusters_Combined))
cellclusterscount<-cellclusters %>% group_by_all() %>% dplyr::summarise(COUNT = n())

cellgroup1count<-cellclusters %>% group_by(reclust) %>% dplyr::summarise(reclust_Total = n())
cellgroup2count<-cellclusters %>% group_by(MergedSample) %>% dplyr::summarise(MergedSample_Total = n())
#normalize by total number of cell per sample (remove the effect of uneual cell num between treated vs untreated)
cellcount_merged<-merge(cellclusterscount,cellgroup1count,by="reclust") %>% merge(cellgroup2count,by="MergedSample")
cellcount_merged$MergedSample_ratio <- cellcount_merged$COUNT/cellcount_merged$MergedSample_Total*100 
#normalize by cluster size
cellgroup3count<-cellcount_merged %>% group_by(reclust) %>% dplyr::summarise(totalratio = sum(MergedSample_ratio))
cellcount_merged<-merge(cellcount_merged,cellgroup3count,by="reclust")
cellcount_merged$norm <- cellcount_merged$MergedSample_ratio/cellcount_merged$totalratio*100 

p1<-stacked_barplot(cellclusterscount,"reclust","COUNT","MergedSample","reclust")
p2<-stacked_barplot(cellclusterscount,"MergedSample","COUNT","reclust","MergedSample")
p1_norm<-stacked_barplot(cellcount_merged,"MergedSample","norm","reclust","NormalizedbyTotalCell")

cellclusterscount_file<-file.path(outputdir,paste(prefix,"scRNA_reclust-MergedSample.count","csv",sep="."))
cellclusterscount %>%
  write.table(file = cellclusterscount_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)

## Plot distribution of clusters in stacked bar plot format
cellclusters<-as.data.frame(colData(sce)[,c("reclust","treatment")])
#cellclusters$Clusters_Combined <- as.numeric(as.character(cellclusters$Clusters_Combined))
cellclusterscount<-cellclusters %>% group_by_all() %>% dplyr::summarise(COUNT = n())

cellgroup1count<-cellclusters %>% group_by(reclust) %>% dplyr::summarise(reclust_Total = n())
cellgroup2count<-cellclusters %>% group_by(treatment) %>% dplyr::summarise(MergedSample_Total = n())
#normalize by total number of cell per sample (remove the effect of uneual cell num between treated vs untreated)
cellcount_merged<-merge(cellclusterscount,cellgroup1count,by="reclust") %>% merge(cellgroup2count,by="treatment")
cellcount_merged$MergedSample_ratio <- cellcount_merged$COUNT/cellcount_merged$MergedSample_Total*100 
#normalize by cluster size
cellgroup3count<-cellcount_merged %>% group_by(reclust) %>% dplyr::summarise(totalratio = sum(MergedSample_ratio))
cellcount_merged<-merge(cellcount_merged,cellgroup3count,by="reclust")
cellcount_merged$norm <- cellcount_merged$MergedSample_ratio/cellcount_merged$totalratio*100 

p3<-stacked_barplot(cellclusterscount,"treatment","COUNT","reclust","treatment")
p4<-stacked_barplot(cellclusterscount,"treatment","COUNT","reclust","treatment")
p2_norm<-stacked_barplot(cellcount_merged,"treatment","norm","reclust","NormalizedbyTotalCell")

cellclusterscount_file<-file.path(outputdir,paste(prefix,"scRNA_reclust-treatment.count","csv",sep="."))
cellclusterscount %>%
  write.table(file = cellclusterscount_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)

## Plot distribution of clusters in stacked bar plot format
cellclusters<-as.data.frame(colData(sce)[,c("reclust","time")])
#cellclusters$Clusters_Combined <- as.numeric(as.character(cellclusters$Clusters_Combined))
cellclusterscount<-cellclusters %>% group_by_all() %>% dplyr::summarise(COUNT = n())

cellgroup1count<-cellclusters %>% group_by(reclust) %>% dplyr::summarise(reclust_Total = n())
cellgroup2count<-cellclusters %>% group_by(time) %>% dplyr::summarise(MergedSample_Total = n())
#normalize by total number of cell per sample (remove the effect of uneual cell num between treated vs untreated)
cellcount_merged<-merge(cellclusterscount,cellgroup1count,by="reclust") %>% merge(cellgroup2count,by="time")
cellcount_merged$MergedSample_ratio <- cellcount_merged$COUNT/cellcount_merged$MergedSample_Total*100 
#normalize by cluster size
cellgroup3count<-cellcount_merged %>% group_by(reclust) %>% dplyr::summarise(totalratio = sum(MergedSample_ratio))
cellcount_merged<-merge(cellcount_merged,cellgroup3count,by="reclust")
cellcount_merged$norm <- cellcount_merged$MergedSample_ratio/cellcount_merged$totalratio*100 
p5<-stacked_barplot(cellclusterscount,"time","COUNT","reclust","time")
p6<-stacked_barplot(cellclusterscount,"time","COUNT","reclust","time")
p3_norm<-stacked_barplot(cellcount_merged,"time","norm","reclust","NormalizedbyTotalCell")

plotPDF(p1, p2,p1_norm,p3, p4,p2_norm,p5, p6,p3_norm,name = "scRNA_reclust-MergedSample_stackedbar", addDOC = FALSE)

cellclusterscount_file<-file.path(outputdir,paste(prefix,"scRNA_reclust-time.count","csv",sep="."))
cellclusterscount %>%
  write.table(file = cellclusterscount_file, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)

#quantify the proportion of clusters
hclust_cluster<-as.data.frame(table(hclust.dyn$clusters))
colnames(hclust_cluster)<-c("cluster","cellnum")
hclust_cluster$cellprop<-round(hclust_cluster$cellnum/sum(hclust_cluster$cellnum)*100,1)
outfile<-file.path(outputdir, paste(prefix, "hclust.num","csv",sep="."))
write.csv(hclust_cluster, outfile,row.names = FALSE,quote=FALSE)

##plot dendrogram
tree <- hclust.dyn$objects$hclust
tree$labels<-seq_along(tree$labels)
dend <- as.dendrogram(tree)
# reorder the dendrogram labels to match the cluster color and label
set_col <- customcolor[as.numeric(hclust.dyn$clusters)]
set_col <- set_col[order.dendrogram(dend)]
set_col <- factor(set_col, unique(set_col))
dend <- as.dendrogram(tree) %>%
  color_branches(clusters = as.numeric(set_col), col = levels(set_col), 
                 groupLabels=TRUE)# %>%
#set("branches_lwd", 1) 
#dend_data <-dendro_data(dend, type = "rectangle")
ggd1 <- as.ggdend(dend)
hclustdend_plot<-ggplot(ggd1, labels=FALSE) 


pca_plot<-reducedimplot(sce,"PCA",selectedgroup,selectedgroup,3)
pca_2dim_plot<-reducedimplot(sce,"PCA",selectedgroup,selectedgroup,2)
pca_2dim_plot_2<-reducedimplot(sce,"PCA","reclust","reclust",2)
dim_plot<-reducedimplot(sce,dimreduc_method,selectedgroup,selectedgroup,3)
dim_recluster_plot<-reducedimplot(sce,dimreduc_method,"reclust","reclust",3)
dim_ident_plot<-reducedimplot(sce,paste("TSNE",dimreduc_method,sep="_"),"MergedSample","MergedSample",2)
dim_ident_plot2<-reducedimplot(sce,paste("TSNE",dimreduc_method,sep="_"),"treatment","treatment",2)
dim_ident_plot3<-reducedimplot(sce,paste("TSNE",dimreduc_method,sep="_"),"time","time",2)

pca_ident_plot<-reducedimplot(sce,"PCA","MergedSample","MergedSample",2)
recluster_tsne_plot<-reducedimplot(sce,paste("TSNE",dimreduc_method,sep="_"),"reclust","reclust",2)
recluster_umap_plot<-reducedimplot(sce,paste("UMAP",dimreduc_method,sep="_"),"reclust","reclust",2)
hclust_plot<-reducedimplot(sce,dimreduc_method,"hclust","hclust",2)
hclust_tsne_plot<-reducedimplot(sce,paste("TSNE",dimreduc_method,sep="_"),"hclust","hclust",2)
hclust_umap_plot<-reducedimplot(sce,paste("UMAP",dimreduc_method,sep="_"),"hclust","hclust",2)


## plot AAV distribution 
color_values <- rep("#8A969E", length(unique(colData(sce)$AAV)))
names(color_values) <- unique(colData(sce)$AAV)
# Set only your group of interest to a specific color
color_values["hATOH1_hPOU4F3_hGFI1"] <- "#511207"
color_values["hGFI1"] <- "#E1251B"
color_values["hATOH1_hPOU4F3"] <- "#F58E7D"

aav_pca_plot<-plotReducedDim(sce, dimred="PCA", colour_by="AAV",point_size=1) +
  guides(colour = guide_legend(override.aes = list(size=4)))  + 
  theme(legend.title = element_blank()) + scale_colour_manual(values = color_values)
aav_tsne_plot<-plotReducedDim(sce, dimred=paste("TSNE",dimreduc_method,sep="_"), colour_by="AAV",point_size=1) +
  guides(colour = guide_legend(override.aes = list(size=4)))  + 
  theme(legend.title = element_blank()) + scale_colour_manual(values = color_values)
aav_umap_plot<-plotReducedDim(sce, dimred=paste("UMAP",dimreduc_method,sep="_"), colour_by="AAV",point_size=1) +
  guides(colour = guide_legend(override.aes = list(size=4)))  + 
  theme(legend.title = element_blank()) + scale_colour_manual(values = color_values)

quality_pca_plot<-plotReducedDim(sce, dimred=paste("UMAP",dimreduc_method,sep="_"), colour_by="Gex_nUMI",point_size=1) +
  theme(legend.title = element_blank())+ scale_color_continuous(trans = "log10", type = "viridis")
quality_tsne_plot<-plotReducedDim(sce, dimred=paste("TSNE",dimreduc_method,sep="_"), colour_by="Gex_nUMI",point_size=1) +
  theme(legend.title = element_blank())+ scale_color_continuous(trans = "log10", type = "viridis")
quality_umap_plot<-plotReducedDim(sce, dimred=paste("UMAP",dimreduc_method,sep="_"), colour_by="Gex_nUMI",point_size=1) +
  theme(legend.title = element_blank())+ scale_color_continuous(trans = "log10", type = "viridis")

## plot sample distribution (each sample highlighted by red per figure)
mergedsamples<-unique(colData(sce)$MergedSample)
sample_plots<-list()
for(samples in mergedsamples){
  color_values <- rep("#8A969E", length(mergedsamples))
  names(color_values) <- mergedsamples
  # Set only your group of interest to a specific color
  color_values[samples] <- "#E1251B"
  p<-plotReducedDim(sce, dimred=paste("TSNE",dimreduc_method,sep="_"), colour_by="MergedSample",point_size=1) +
    guides(colour = guide_legend(override.aes = list(size=4)))  + 
    theme(legend.title = element_blank()) + scale_colour_manual(values = color_values)
  sample_plots<-c(sample_plots,list(p))
}
plotPDF(pca_plot, dim_plot, pca_2dim_plot,pca_2dim_plot_2,dim_recluster_plot,dim_ident_plot,dim_ident_plot2,dim_ident_plot3,pca_ident_plot,recluster_tsne_plot,recluster_umap_plot,
        hclustdend_plot,hclust_plot,hclust_tsne_plot,hclust_umap_plot,sample_plots,
        aav_pca_plot,aav_tsne_plot,aav_umap_plot,quality_pca_plot,quality_tsne_plot,quality_umap_plot,name = paste("reclust",dimreduc_method,sep="."), addDOC = FALSE)

# plot heatmap plot
genelist<-unique(genelist[genelist%in% rownames(sce)])
genelist_heatmap<-plotGroupedHeatmap(sce, features=genelist, group="reclust", 
                                     center=TRUE, zlim=c(-3, 3)) 

heatmap_fig<-file.path(outputdir,"Plots",paste(prefix,"scRNA.reclust.GE.heatmap","pdf",sep="."))
pdf(heatmap_fig, width=9, height=12,0)
print(genelist_heatmap)
dev.off()

genelist<-unique(genelist[genelist%in% rownames(sce)])
genelist_heatmap<-plotGroupedHeatmap(sce, features=genelist, group="hclust", 
                                     center=TRUE, zlim=c(-3, 3)) 

heatmap_fig<-file.path(outputdir,"Plots",paste(prefix,"scRNA.hclust.GE.heatmap","pdf",sep="."))
pdf(heatmap_fig, width=9, height=12,0)
print(genelist_heatmap)
dev.off()

genelist<-unique(genelist[genelist%in% rownames(sce)])
counts_matrix <- assays(sce)$logcounts
genelist_GE <- t(as.data.frame(as.matrix(counts_matrix[genelist, ])))
genelist_GE <- merge( t(as.data.frame(as.matrix(counts_matrix[genelist, ]))),as.data.frame(colData(sce)),by=0)
gene_cluster<-genelist_GE[,c("reclust",genelist)] %>% 
  group_by(reclust)  %>%  
  melt(id.vars = c("reclust")) %>%  
  group_by(reclust,variable) %>% 
  summarise_all(c(count = mean, cell_exp_ct = ~ sum(. > 0),cell_ct=~ n()))
gene_cluster<-gene_cluster %>%
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
  filter(count > 0, `% Expressing` > 1) 
outfile<-file.path(outputdir, paste(prefix,"reclust.genequant","csv",sep="."))
write.csv(gene_cluster, outfile,row.names = FALSE,quote=FALSE)

dotplot<-ggplot(gene_cluster,aes(x=variable, y = reclust, color = count, size = `% Expressing`)) + 
  geom_point() + 
  scale_size(limits = c(0,max(10,max(gene_cluster$`% Expressing`)))) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 40, hjust=1),    axis.text.y = element_text(margin = margin(r = 0.1))) +
  ylab('') + xlab("") +
  scale_color_gradientn(colours = colorRampPalette(c("#0F3A85","#99BFE5","#E4EBF1","#8A969E", "#F58E7D","#E1251B", "#511207"))(10), 
                        limits = c(0,min(1,max(gene_cluster$count))), 
                        labels = ~sprintf("%g", .x),
                        oob = scales::squish, name = 'log2 (count + 1)') +
  theme(legend.position = "bottom")

dotplot_file<-file.path(outputdir,"Plots",paste(prefix,"reclust.GE.dotplot",dimreduc_method,"pdf",sep="."))
pdf(dotplot_file, width=14, height=14)
print(dotplot)
dev.off()

sample_plots<-list(recluster_tsne_plot,pca_2dim_plot_2)
for(gene in genelist){
  #gene<-'DLK1'
  #print(gene)
  p<-plotReducedDim(sce, paste("TSNE",dimreduc_method,sep="_"), 
                    colour_by=gene,point_size=1) +
    theme(legend.position="top")
  sample_plots<-c(sample_plots,list(p))
}

samplenum<-length(genelist)
geneexpression_fig<-file.path(outputdir,"Plots",paste(prefix,"reclust.GE",paste("TSNE",dimreduc_method,sep="_"),"pdf",sep="."))
pdf(geneexpression_fig, width=1.5*round(sqrt(samplenum),0)+9, height=1.5*round(sqrt(samplenum)+6,0))
grid.arrange(grobs = sample_plots, 3,0)## display plot
dev.off()

sample_plots<-list(recluster_tsne_plot,pca_2dim_plot_2)
for(gene in genelist){
  #gene<-'DLK1'
  #print(gene)
  p<-plotReducedDim(sce, "PCA", 
                    colour_by=gene,point_size=1) +
    theme(legend.position="top")
  sample_plots<-c(sample_plots,list(p))
}

samplenum<-length(genelist)
geneexpression_fig<-file.path(outputdir,"Plots",paste(prefix,"reclust.GE","PCA","pdf",sep="."))
pdf(geneexpression_fig, width=1.5*round(sqrt(samplenum),0)+9, height=1.5*round(sqrt(samplenum)+6,0))
grid.arrange(grobs = sample_plots, 3,0)## display plot
dev.off()

# plot heatmap plot
devgenelist<-unique(devgenelist[devgenelist%in% rownames(sce)])
genelist_heatmap<-plotGroupedHeatmap(sce, features=devgenelist, group="reclust", 
                                     center=TRUE, zlim=c(-3, 3)) 

heatmap_fig<-file.path(outputdir,"Plots",paste(prefix,"scRNA.reclust.devGE.heatmap","pdf",sep="."))
pdf(heatmap_fig, width=9, height=12,0)
print(genelist_heatmap)
dev.off()

genelist_heatmap<-plotGroupedHeatmap(sce, features=devgenelist, group="hclust", 
                                     center=TRUE, zlim=c(-3, 3)) 

heatmap_fig<-file.path(outputdir,"Plots",paste(prefix,"scRNA.hclust.devGE.heatmap","pdf",sep="."))
pdf(heatmap_fig, width=9, height=12,0)
print(genelist_heatmap)
dev.off()
# plot violin plot
violin<-plotExpression(sce, features=devgenelist,
                       x=I(colData(sce)$reclust),color_by =I(colData(sce)$reclust),  ncol = 3,point_size=1)  + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  facet_wrap(~Feature, scales = "free_y") + 
  scale_color_manual(values=customcolor)
geneexpression_fig<-file.path(outputdir,"Plots",paste(prefix,"reclust.devGE.violin","pdf",sep="."))
pdf(geneexpression_fig, width=1*round(sqrt(samplenum),0)+3, height=1*round(sqrt(samplenum),0))
print(violin)
dev.off()

violin<-plotExpression(sce, features=genelist,
                       x=I(colData(sce)$reclust),color_by =I(colData(sce)$reclust),  ncol = 3,point_size=1)  + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  facet_wrap(~Feature, scales = "free_y") + 
  scale_color_manual(values=customcolor)
geneexpression_fig<-file.path(outputdir,"Plots",paste(prefix,"reclust.GE.violin","pdf",sep="."))
pdf(geneexpression_fig, width=1*round(sqrt(samplenum),0)+9, height=1*round(sqrt(samplenum),0))
print(violin)
dev.off()


counts_matrix <- assays(sce)$logcounts
genelist_GE <- t(as.data.frame(as.matrix(counts_matrix[devgenelist, ])))
genelist_GE <- merge( t(as.data.frame(as.matrix(counts_matrix[devgenelist, ]))),as.data.frame(colData(sce)),by=0)

gene_cluster<-genelist_GE[,c("reclust",devgenelist)] %>% 
  group_by(reclust)  %>%  
  melt(id.vars = c("reclust")) %>%  
  group_by(reclust,variable) %>% 
  summarise_all(c(count = mean, cell_exp_ct = ~ sum(. > 0),cell_ct=~ n()))
gene_cluster<-gene_cluster %>%
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
  filter(count > 0, `% Expressing` > 1) 
outfile<-file.path(outputdir, paste(prefix,"reclust.devgenequant","csv",sep="."))
write.csv(gene_cluster, outfile,row.names = FALSE,quote=FALSE)

dotplot<-ggplot(gene_cluster,aes(x=variable, y = reclust, color = count, size = `% Expressing`)) + 
  geom_point() + 
  scale_size(limits = c(0,max(10,max(gene_cluster$`% Expressing`)))) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 40, hjust=1),    axis.text.y = element_text(margin = margin(r = 0.1))) +
  ylab('') + xlab("") +
  scale_color_gradientn(colours = colorRampPalette(c("#0F3A85","#99BFE5","#E4EBF1","#8A969E", "#F58E7D","#E1251B", "#511207"))(10), 
                        limits = c(0,min(1,max(gene_cluster$count))), 
                        labels = ~sprintf("%g", .x),
                        oob = scales::squish, name = 'log2 (count + 1)') +
  theme(legend.position = "bottom")

dotplot_file<-file.path(outputdir,"Plots",paste(prefix,"reclust.devGE.dotplot",dimreduc_method,"pdf",sep="."))
pdf(dotplot_file, width=14, height=14)
print(dotplot)
dev.off()


#sample_plots<-list(pca_plot,pca_recluster_plot)
sample_plots<-list(recluster_tsne_plot,pca_2dim_plot_2)

for(gene in devgenelist){
  #gene<-'DLK1'
  #print(gene)
  p<-plotReducedDim(sce, paste("TSNE",dimreduc_method,sep="_"), 
                    colour_by=gene,point_size=1) +
    theme(legend.position="top")
  sample_plots<-c(sample_plots,list(p))
}

samplenum<-length(genelist)
geneexpression_fig<-file.path(outputdir,"Plots",paste(prefix,"reclust.devGE",paste("TSNE",dimreduc_method,sep="_"),"pdf",sep="."))
pdf(geneexpression_fig, width=1*round(sqrt(samplenum),0)+9, height=1*round(sqrt(samplenum)+6,0))
grid.arrange(grobs = sample_plots, 3,0)## display plot
dev.off()

sample_plots<-list(recluster_tsne_plot,pca_2dim_plot_2)

for(gene in devgenelist){
  #gene<-'DLK1'
  #print(gene)
  p<-plotReducedDim(sce, "PCA", 
                    colour_by=gene,point_size=1) +
    theme(legend.position="top")
  sample_plots<-c(sample_plots,list(p))
}

samplenum<-length(genelist)
geneexpression_fig<-file.path(outputdir,"Plots",paste(prefix,"reclust.devGE","PCA","pdf",sep="."))
pdf(geneexpression_fig, width=1*round(sqrt(samplenum),0)+9, height=1*round(sqrt(samplenum)+6,0))
grid.arrange(grobs = sample_plots, 3,0)## display plot
dev.off()

markerdetect(sce,"reclust")

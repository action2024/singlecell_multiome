#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
archRdirpath <- args[1]
outputdir<-args[2]
group<-args[3]

library(ArchR)
library(Signac)
library(Seurat)
#library(EnsDb.Hsapiens.v86)
library(stringr)
library(scater)

#archRdirpath = "$PATH/$SAMPLE/clustering"
#outputdir = "$PATH/$SAMPLE/MarkersRNA"
#group<-"Clusters_Combined"
projMulti_filtered <- loadArchRProject(path = archRdirpath)
# Identify peaks are unique to an individual cluster or a small group of clusters
se <- getMarkerFeatures(ArchRProj = projMulti_filtered,
                        groupBy = group,
                        useMatrix = "GeneScoreMatrix", 
                        bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"))

heatmap_gex <- plotMarkerHeatmap(
  seMarker = se, 
  log2Norm = TRUE,scaleRows = TRUE,
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
  nLabel = 6,
  transpose = TRUE
)

setwd(outputdir)
plotPDF(heatmap_gex, name = "GeneMarker-Heatmap", ArchRProj = projMulti_filtered)

se_markerList <- getMarkers(se, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
for(cluster in names(se_markerList)){
  se_markerList_cluster<-se_markerList[[cluster]]
  se_markerList_cluster<-se_markerList_cluster[order(-se_markerList_cluster$Log2FC), ]
  write.csv(se_markerList_cluster, file.path(outputdir,  paste0(cluster,".genescoremarkers.csv")), row.names = FALSE)
  
  }

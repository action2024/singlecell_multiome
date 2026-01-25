#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
archRdirpath <- args[1]
outputdir<-args[2]
group<-args[3]
library(ArchR)
library(Signac)
library(Seurat)
#BiocManager::install("EnsDb.Hsapiens.v86")
#library(EnsDb.Hsapiens.v86)
library(stringr)
library(scater)
library(BSgenome.Hsapiens.UCSC.hg38)

#archRdirpath<-"/cluster/home/liuji/T1D/multiome/analysis/Cyt49_S4D5/peakcalling/Clusters_Combined"
#outputdir<-"/cluster/home/liuji/T1D/multiome/analysis/Cyt49_S4D5/PeakMarker"
#group<-"Clusters_Combined"
projMulti_filtered <- loadArchRProject(path = archRdirpath)
#marker detection based on peak calling
markersPeaks <- getMarkerFeatures(
  ArchRProj = projMulti_filtered, 
  useMatrix = "PeakMatrix", 
  groupBy = group,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

se_markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
for(cluster in names(se_markerList)){
  se_markerList_cluster<-se_markerList[[cluster]]
  se_markerList_cluster<-se_markerList_cluster[order(-se_markerList_cluster$Log2FC), ]
  write.csv(se_markerList_cluster, file.path(outputdir,  paste0(cluster,".peakmarkers.csv")), row.names = FALSE)
  
}

# Plot marker peaks
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
setwd(outputdir)
plotPDF(heatmapPeaks, name = paste0(group,"-PeakMarker-Heatmap"), width = 8, height = 6, ArchRProj = projMulti_filtered, addDOC = FALSE)
# Motif Enrichment in Marker Peaks
# getCellColData(projMulti_filtered)
# add motif annotations
if("Motif" %ni% names(projMulti_filtered@peakAnnotation)){
  projMulti_filtered <- addMotifAnnotations(ArchRProj = projMulti_filtered, motifSet = "cisbp", name = "Motif")
}
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projMulti_filtered,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

if(ncol(colData(enrichMotifs)) > 0){
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
  #ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projMulti_filtered, addDOC = FALSE)
  df <- data.frame(TF = rownames(enrichMotifs), mlog10Padj = assay(enrichMotifs)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))
  write.csv(df, file.path(outputdir,  paste0("Motifs.PeakMarker.csv")), row.names = FALSE)
  
}

if("EncodeTFBS" %ni% names(projMulti_filtered@peakAnnotation)){
  projMulti_filtered <- addArchRAnnotations(ArchRProj = projMulti_filtered, collection = "EncodeTFBS")}
  enrichEncode <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projMulti_filtered,
  peakAnnotation = "EncodeTFBS",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
if(ncol(colData(enrichEncode)) > 0){
  heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
  #ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projMulti_filtered, addDOC = FALSE)
  df <- data.frame(TF = rownames(enrichEncode), mlog10Padj = assay(enrichEncode)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))
  write.csv(df, file.path(outputdir,  paste0("EncodeTFBS.PeakMarker.csv")), row.names = FALSE)
}

# ChromVAR Deviatons Enrichment 
# add a set of background peaks which are used in computing deviations
projMulti_filtered <- addBgdPeaks(projMulti_filtered)
projMulti_filtered <- addDeviationsMatrix(
  ArchRProj = projMulti_filtered, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(projMulti_filtered, name = "MotifMatrix", plot = FALSE)
write.csv(data.frame(plotVarDev), file.path(outputdir,  paste0("chromVAR.csv")), row.names = FALSE)
plotVarDev <- getVarDeviations(projMulti_filtered, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "ChromVAR", width = 8, height = 6, ArchRProj = projMulti_filtered, addDOC = FALSE)




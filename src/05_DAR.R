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
#library(ArchRtoSignac) # {Link: GitHub https://github.com/swaruplabUCI/ArchRtoSignac}

#archRdirpath<-"/cluster/home/liuji/T1D/multiome/analysis/Cyt49_S4D5/peakcalling/Clusters_Combined"
# group<-"Clusters_Combined"
# outputdir<-"/cluster/home/liuji/T1D/multiome/analysis/Cyt49_S4D5/PeakMarker"

projMulti_filtered <- loadArchRProject(path = archRdirpath)
if("Motif" %ni% names(projMulti_filtered@peakAnnotation)){
projMulti_filtered <- addMotifAnnotations(ArchRProj = projMulti_filtered, motifSet = "cisbp", name = "Motif")}

allclusters<-unique(getCellColData(projMulti_filtered)[,group])
for(i in 1:length(allclusters)){
  test<-allclusters[i]
  control<-allclusters[!allclusters %in% test]
  markerTest <- getMarkerFeatures(
    ArchRProj = projMulti_filtered, 
    useMatrix = "PeakMatrix",
    groupBy = group,
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = test,
    bgdGroups = control
  )
  pma <- markerPlot(seMarker = markerTest, name = test, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")

  motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projMulti_filtered,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  updf <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
  updf <- updf[order(updf$mlog10Padj, decreasing = TRUE),]
  updf$rank <- seq_len(nrow(updf))
  write.csv(updf, file.path(outputdir,  paste0(test,".DAR.upregulated.csv")), row.names = FALSE)
  
  ggUp <- ggplot(updf, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = updf[rev(seq_len(30)), ], 
      aes(x = rank, y = mlog10Padj, label = TF), 
      size = 3,
      nudge_x = 2,
      color = "black"
    ) + theme_ArchR() + 
    ylab("-log10(P-adj) Motif Enrichment") + 
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  plotPDF(pma, ggUp, name = paste0(test,".DAR.upregulatedpeaks",".pdf"), addDOC = FALSE)

  motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projMulti_filtered,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )
  downdf <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
  downdf <- downdf[order(downdf$mlog10Padj, decreasing = TRUE),]
  downdf$rank <- seq_len(nrow(downdf))
  write.csv(downdf, file.path(outputdir,  paste0(test,".DAR.downregulated.csv")), row.names = FALSE)
  
  ggDo <- ggplot(downdf, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = downdf[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
      size = 1.5,
      nudge_x = 2,
      color = "black"
    ) + theme_ArchR() + 
    ylab("-log10(FDR) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  plotPDF(pma, ggDo, name = paste0(test,".DAR.downregulatedpeaks",".pdf"), addDOC = FALSE)
  
}

## plot gene expression level with selected markers

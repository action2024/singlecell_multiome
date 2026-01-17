#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
archRdirpath <- args[1]
outputdir<-args[2]

library(ArchR)
library(Signac)
library(Seurat)
#BiocManager::install("EnsDb.Hsapiens.v86")
#library(EnsDb.Hsapiens.v86)
library(stringr)
library(scater)
library(BSgenome.Hsapiens.UCSC.hg38)
projMulti_filtered <- loadArchRProject(path = archRdirpath)

# Example: Keeping only cells from a specific list of cell IDs
# Optional: Get a subset of cell names from existing project metadata
# getCellColData(projMulti_filtered)
# cluster<-"C1"
for(cluster in unique(projMulti_filtered$Clusters_Combined)){
idxPass <- which(projMulti_filtered$Clusters_Combined == cluster)
table(projMulti_filtered$Clusters_Combined)
cellsPass <- projMulti_filtered$cellNames[idxPass]

# Create the new, subsetted ArchR project
#outputdir_sub<-file.path(outputdir,paste(cluster,"subset",sep="_"))
# projMulti_filtered_subset <- subsetArchRProject(ArchRProj = projMulti_filtered, 
#                                  cells = cellsPass,force = TRUE)
projMulti_filtered_subset <- projMulti_filtered[cellsPass]

}


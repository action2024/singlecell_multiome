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
projMulti_filtered <- loadArchRProject(path = archRdirpath)
macs2path <- findMacs2()
# macs2path<-"/usr/local/Caskroom/miniforge/base/envs/macs2_env/bin/macs2"
# group<-"Clusters_Combined"
#archRdirpath<-"$path/$sample/peakcalling/Clusters_Combined"
#Calling Peaks w/ Macs2
# 1. create pseudobulk replicates 
projMulti_filtered <- addGroupCoverages(ArchRProj = projMulti_filtered, groupBy = group, verbose = FALSE)
# 2. peak calling with MACS2
#  peaks are called across the replicates for each group indicated by the groupBy argument
# score column: normalized -log10(pval), replicateScoreQuantile: quantile rank of that score for the individual pseudobulk replicate 
projMulti_filtered <- addReproduciblePeakSet(ArchRProj = projMulti_filtered, groupBy = group, pathToMacs2 = macs2path)
#retrieve the peak set stored in your ArchRProject as a GRanges object
getPeakSet(projMulti_filtered)
# add peak matrix: compute counts for each peak per cell 
projMulti_filtered <- addPeakMatrix(ArchRProj = projMulti_filtered)
# save archR project
projMulti_filtered <- saveArchRProject(ArchRProj = projMulti_filtered, outputDirectory = outputdir, load = TRUE)


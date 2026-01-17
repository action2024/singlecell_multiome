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

#group<-"Clusters_Combined"
## Identification of Positive TF-Regulators
# Step 1. Identify Deviant TF Motifs
seGroupMotif <- getGroupSE(ArchRProj = projMulti_filtered, 
                           useMatrix = "MotifMatrix", 
                           groupBy = group)

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

# Step 2. Identify Correlated TF Motifs and TF Gene Score/Expression
allclusters<-unique(getCellColData(projMulti_filtered)[,group])

for(i in 1:length(allclusters)){
  cluster<-allclusters[i]
  #cluster="C11"
  idxPass <- which(getCellColData(projMulti_filtered)[,group] == cluster)
  #table(projMulti_filtered$Clusters_Combined)
  cellsPass <- projMulti_filtered$cellNames[idxPass]
  if(length(cellsPass)>=100){
  # Create the new, subsetted ArchR project
  #outputdir_sub<-file.path(outputdir,paste(cluster,"subset",sep="_"))
  # projMulti_filtered_subset <- subsetArchRProject(ArchRProj = projMulti_filtered, 
  #                                  cells = cellsPass,force = TRUE)
  projMulti_filtered_subset <- projMulti_filtered[cellsPass]
  
  corGSM_MM <- correlateMatrices(
    ArchRProj = projMulti_filtered_subset,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "LSI_Combined"
  )
  
  # Step 3. Add Maximum Delta Deviation to the Correlation Data Frame
  corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
  # Step 4. Identify Positive TF Regulators
  corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
  corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
  corGSM_MM$TFRegulator <- "NO"
  corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
  df<-as.data.frame(sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",]))
  write.csv(df, file.path(outputdir,  paste0(cluster,".positiveTFs.csv")), row.names = FALSE)
  }
}
  
  

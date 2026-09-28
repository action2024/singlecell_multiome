#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
archRdirpath <- args[1]
inputrds<-args[2]
outputdir<-args[3]
prefix<-args[4]

source("/home/l128405/multiome/src/R/variables/colors.R")
source("/home/l128405/multiome/src/R/functions/violin_plot_qc.R")
source("/home/l128405/multiome/src/R/functions/LSI_dimplot.R")
source("/home/l128405/multiome/src/R/functions/stacked_barplot.R")


library(ArchR)
library(Signac)
library(Seurat)
#BiocManager::install("EnsDb.Hsapiens.v86")
#library(EnsDb.Hsapiens.v86)
library(stringr)
library(scater)
library(SingleCellExperiment)
library(dplyr)
library(scran)
addArchRThreads(16)
## Setting default number of Parallel threads to 8.
addArchRLocking(locking = TRUE)

#archRdirpath<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen/S41A/"
#inputrds<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen/S41A/S41A.unnormalized.sce.rds"

dir.create(outputdir)
setwd(outputdir)
projMulti_filtered <- loadArchRProject(path = archRdirpath)
sce<-readRDS(inputrds)

## filter out low quality cells with RNA transcripts <500, gene number <200
mito_cutoff<-20
genenum_cutoff<-200
transcriptnum_cutoff<-500

gene_expr_matrix <- getMatrixFromProject(projMulti_filtered, useMatrix = "GeneExpressionMatrix")
sampleQC_df_RNA<-colData(gene_expr_matrix)
# remove low-quality cells
lowQcellIDs_mito<-rownames(sampleQC_df_RNA[sampleQC_df_RNA$subsets_Mito_percent>mito_cutoff,])
lowQcellIDs_genes<-rownames(sampleQC_df_RNA[sampleQC_df_RNA$Gex_nGenes<genenum_cutoff,])
lowQcellIDs_transcripts<-rownames(sampleQC_df_RNA[sampleQC_df_RNA$Gex_nUMI<transcriptnum_cutoff,])
lowQcellIDs<-unique(c(lowQcellIDs_mito,lowQcellIDs_genes,lowQcellIDs_transcripts))
cellsToKeep <- getCellNames(projMulti_filtered)[which(!getCellNames(projMulti_filtered) %in% lowQcellIDs)]

projMulti_filtered <- projMulti_filtered[cellsToKeep,]
sce<-sce[,cellsToKeep]

#table(assay(sce,"counts")['3mir-hATOH1-hPOU4F3-3Flag-3mir&#8221;; transcript_type ',])
#table(assay(sce,"counts")['hGFI1-HA&#8221;; transcript_type ',])

# dimensionality reduction with IterativeLSI
# LSI was performed on scATAC-seq data via the TileMatrix and on the scRNA-seq data via the GeneExpressionMatrix
# scATAC-seq
projMulti_filtered <- addIterativeLSI(
  ArchRProj = projMulti_filtered, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC",
  force = TRUE
)

# scRNA-seq
projMulti_filtered <- addIterativeLSI(
  ArchRProj = projMulti_filtered, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  selectionMethod= "var",
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA",
  force = TRUE
)

#find and remove cells that do not overlap between LSI_ATAC and LSI_RNA
nonoverlap_cellIDs<-setdiff(rownames(projMulti_filtered@reducedDims$LSI_ATAC$matSVD), rownames(projMulti_filtered@reducedDims$LSI_RNA$matSVD)) 
cellsToKeep <- getCellNames(projMulti_filtered)[which(!getCellNames(projMulti_filtered) %in% nonoverlap_cellIDs)]
projMulti_filtered <- subsetCells(projMulti_filtered, cellNames = cellsToKeep)
sce<-sce[,cellsToKeep]

cellIDfile<-file.path(outputdir,paste(prefix,"cellIDs","txt",sep="."))
cellsToKeep %>%  write.table(file = cellIDfile, sep = "\t", quote=FALSE,row.names = FALSE, col.names=FALSE)
#getGenes(projMulti_filtered)$symbol
#projMulti_filtered@reducedDims$LSI_ATAC$matSVD %>% dim()
#projMulti_filtered@reducedDims$LSI_RNA$matSVD %>% dim()
#getGenes(projMulti_filtered)$symbol[grepl("^mt-", getGenes(projMulti_filtered)$symbol)]
# create a dimensionality reduction that uses information from both the scATAC-seq and scRNA-seq data
projMulti_filtered <- addCombinedDims(projMulti_filtered, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
# create UMAP embeddings for each of these dimensionality reductions
projMulti_filtered <- addUMAP(projMulti_filtered, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
projMulti_filtered <- addUMAP(projMulti_filtered, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
projMulti_filtered <- addUMAP(projMulti_filtered, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
# call clusters for each
projMulti_filtered <- addClusters(projMulti_filtered, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.4, force = TRUE)
projMulti_filtered <- addClusters(projMulti_filtered, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 0.4, force = TRUE)
projMulti_filtered <- addClusters(projMulti_filtered, reducedDims = "LSI_Combined", name = "Clusters_Combined", resolution = 0.4, force = TRUE)
#save project to outputdir
projMulti_filtered <- saveArchRProject(ArchRProj = projMulti_filtered, outputDirectory = outputdir, load = TRUE)


# save as single-cell experiment for single-cell-only processing
#getAvailableMatrices(projMulti_filtered)
#normalization
lib.sce <- librarySizeFactors(sce)
clust <- quickCluster(sce) 
sce <- computeSumFactors(sce, cluster=clust, min.mean=0.1)
sce <- logNormCounts(sce)
#add cluster info to sce
gene_expr_matrix <- getMatrixFromProject(projMulti_filtered, useMatrix = "GeneExpressionMatrix")
colData(sce) <- colData(gene_expr_matrix)

# add dimention reduction from archR project
reducedDim(sce, "LSI_Combined") <- getReducedDims(projMulti_filtered,"LSI_Combined")
reducedDim(sce, "LSI_RNA") <- getReducedDims(projMulti_filtered,"LSI_RNA")
reducedDim(sce, "LSI_ATAC") <- getReducedDims(projMulti_filtered,"LSI_ATAC")
saveRDS(sce, file = file.path(outputdir, paste0(prefix,".normalized.sce.rds")))
#rownames(sce)[which(grepl("^mt-", rownames(sce)))]

#cell proportion for each cluster
rna_clusters<-as.data.frame(table(projMulti_filtered$Clusters_RNA))
names(rna_clusters)<-c("clusters","cellnum")
rna_clusters$prop<-round(rna_clusters$cellnum/sum(rna_clusters$cellnum)*100,1)

atac_clusters<-as.data.frame(table(projMulti_filtered$Clusters_ATAC))
names(atac_clusters)<-c("clusters","cellnum")
atac_clusters$prop<-round(atac_clusters$cellnum/sum(atac_clusters$cellnum)*100,1)

combined_clusters<-as.data.frame(table(projMulti_filtered$Clusters_Combined))
names(combined_clusters)<-c("clusters","cellnum")
combined_clusters$prop<-round(combined_clusters$cellnum/sum(combined_clusters$cellnum)*100,1)

setwd(outputdir)
plotPDF(grid.arrange(top="ATAC", tableGrob(atac_clusters)),
        grid.arrange(top="RNA", tableGrob(rna_clusters)),
        grid.arrange(top="Combined", tableGrob(combined_clusters)),name = "Table-scATAC-scRNA-Combined", addDOC = FALSE)


p1_ATAC <- plotEmbedding(projMulti_filtered, name = "Clusters_ATAC", embedding = "UMAP_ATAC", size = 0.1)
p1 <- plotEmbedding(projMulti_filtered, name = "Clusters_Combined", embedding = "UMAP_ATAC", size = 0.1)
p2_RNA <- plotEmbedding(projMulti_filtered, name = "Clusters_RNA", embedding = "UMAP_RNA", size = 0.1)
p2 <- plotEmbedding(projMulti_filtered, name = "Clusters_Combined", embedding = "UMAP_RNA", size = 0.1)
p3 <- plotEmbedding(projMulti_filtered, name = "Clusters_Combined", embedding = "UMAP_Combined", size = 0.1)
p4<- plotEmbedding(projMulti_filtered, name = "Sample", embedding = "UMAP_Combined", size = 0.1)

dim1<-LSI_dimplot(projMulti_filtered,"LSI_Combined","Clusters_Combined","LSI1","LSI2")
dim2<-LSI_dimplot(projMulti_filtered,"LSI_Combined","Clusters_Combined","LSI1","LSI3")
dim3<-LSI_dimplot(projMulti_filtered,"LSI_Combined","Clusters_Combined","LSI2","LSI3")

dim1_RNA<-LSI_dimplot(projMulti_filtered,"LSI_RNA","Clusters_RNA","LSI1","LSI2")
dim2_RNA<-LSI_dimplot(projMulti_filtered,"LSI_RNA","Clusters_RNA","LSI1","LSI3")
dim3_RNA<-LSI_dimplot(projMulti_filtered,"LSI_RNA","Clusters_RNA","LSI2","LSI3")

setwd(outputdir)
plotPDF(p1_ATAC,p1, p2_RNA,p2, p3,p4,dim1,dim2,dim3,dim1_RNA,dim2_RNA,dim3_RNA, name = "UMAP_LSI-scATAC-scRNA-Combined", addDOC = FALSE)


# visualize differences in cluster residence of cells between scATAC-seq, scRNA-seq and combined
cM_atac_rna <- confusionMatrix(paste0(projMulti_filtered$Clusters_ATAC), paste0(projMulti_filtered$Clusters_RNA))
cM_atac_rna <- cM_atac_rna / Matrix::rowSums(cM_atac_rna)
library(pheatmap)
p_atac_rna <- pheatmap::pheatmap(
  mat = as.matrix(cM_atac_rna), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)


cM_combined_rna <- confusionMatrix(paste0(projMulti_filtered$Clusters_Combined), paste0(projMulti_filtered$Clusters_RNA))
cM_combined_rna <- cM_combined_rna / Matrix::rowSums(cM_combined_rna)
p_combined_rna <- pheatmap::pheatmap(
  mat = as.matrix(cM_combined_rna), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

cM_combined_atac <- confusionMatrix(paste0(projMulti_filtered$Clusters_Combined), paste0(projMulti_filtered$Clusters_ATAC))
cM_combined_atac <- cM_combined_atac / Matrix::rowSums(cM_combined_atac)
p_combined_atac <- pheatmap::pheatmap(
  mat = as.matrix(cM_combined_atac), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
plotPDF(p_atac_rna$gtable, p_combined_rna$gtable, p_combined_atac$gtable, name = "Heatmap-scATAC-scRNA-Combined", addDOC = FALSE)



cellclustersref<-as.data.frame(getCellColData(projMulti_filtered)[c("Clusters_ATAC","Clusters_Combined")])
cellclusterscount<-cellclustersref %>% group_by_all() %>% summarise(COUNT = n())
p1<-stacked_barplot(cellclusterscount,"Clusters_ATAC","COUNT","Clusters_Combined","Combined-ATAC")
p2<-stacked_barplot(cellclusterscount,"Clusters_Combined","COUNT","Clusters_ATAC","ATAC-Combined")

cellclustersref<-as.data.frame(getCellColData(projMulti_filtered)[c("Clusters_RNA","Clusters_Combined")])
cellclusterscount<-cellclustersref %>% group_by_all() %>% summarise(COUNT = n())
p3<-stacked_barplot(cellclusterscount,"Clusters_RNA","COUNT","Clusters_Combined","Combined-RNA")
p4<-stacked_barplot(cellclusterscount,"Clusters_Combined","COUNT","Clusters_RNA","RNA-Combined")

plotPDF(p1,p2,p3,p4,name = "UMAP-scATAC-scRNA-Combined-stackedbar", addDOC = FALSE)



sampleQC_df_multiome_filtered<-as.data.frame(getCellColData(projMulti_filtered))

scmultiome_qc_sum_filtered<-file.path(outputdir,paste(prefix,"scmultiome.filtered.sample.qc.sum","csv",sep="."))
sampleQC_df_multiome_filtered %>%
  group_by(Sample) %>% 
  summarise(cells=n(),
            median_transcripts = round(median(Gex_nUMI),0),
            median_genes = round(median(Gex_nGenes),0),
            median_frag = round(median(nFrags),0),
            median_TSSenrichment = round(median(TSSEnrichment),0),
            quantile_10_transcripts = quantile(Gex_nUMI, probs = 0.1),
            quantile_90_transcripts = quantile(Gex_nUMI, probs = 0.9),
            quantile_10_genes = quantile(Gex_nGenes, probs = 0.1),
            quantile_90_genes = quantile(Gex_nGenes, probs = 0.9),
            quantile_10_frag = quantile(nFrags, probs = 0.1),
            quantile_90_frag = quantile(nFrags, probs = 0.9),
            quantile_10_TSSenrichment = quantile(TSSEnrichment, probs = 0.1),
            quantile_90_TSSenrichment = quantile(TSSEnrichment, probs = 0.9)) %>% 
  write.table(file = scmultiome_qc_sum_filtered, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)


#sum_threshold<-20000
#TSSEnrichment_threshold<-max(50,quantile(sampleQC_df_scRNA_filtered$TSSEnrichment, probs = 0.5)*2,quantile(sampleQC_df_scRNA_filtered$TSSEnrichment, probs = 0.9))
#detected_threshold<-10000
#nFrags_threshold<-max(20000,quantile(sampleQC_df_scRNA_filtered$nFrags, probs = 0.5)*2,quantile(sampleQC_df_scRNA_filtered$nFrags, probs = 0.9))
#nGenes_threshold<-max(5000,quantile(sampleQC_df_scRNA_filtered$Gex_nGenes, probs = 0.5)*2,quantile(sampleQC_df_scRNA_filtered$Gex_nGenes, probs = 0.9))
#nTranscripts_threshold<-max(10000,quantile(sampleQC_df_scRNA_filtered$Gex_nUMI, probs = 0.5)*2,quantile(sampleQC_df_scRNA_filtered$Gex_nUMI, probs = 0.9))
df<-sampleQC_df_multiome_filtered[,c("TSSEnrichment","nFrags","Gex_nGenes","Gex_nUMI","Sample")]
df_cutoff<-data.frame(variable = c("TSSEnrichment", "nFrags","Gex_nGenes","Gex_nUMI"), cutoff = c(4,1000,200,500))
#df<-sampleQC_df_RNA[sampleQC_df_RNA$TSSEnrichment<TSSEnrichment_threshold & sampleQC_df_RNA$nFrags<nFrags_threshold & sampleQC_df_RNA$Gex_nGenes<nGenes_threshold & sampleQC_df_RNA$Gex_nUMI<nTranscripts_threshold,c("Sample","TSSEnrichment","nFrags","Gex_nUMI","Gex_nGenes")]
df_melt<-melt(df, id.vars = c("Sample"))
df <- merge(df_melt,df_cutoff,by="variable")
#df<-sampleQC_df_scRNA_filtered[sampleQC_df_scRNA_filtered$TSSEnrichment<TSSEnrichment_threshold & sampleQC_df_scRNA_filtered$nFrags<nFrags_threshold & sampleQC_df_scRNA_filtered$Gex_nGenes<nGenes_threshold & sampleQC_df_scRNA_filtered$Gex_nUMI<nTranscripts_threshold, c("Sample","TSSEnrichment","nFrags","Gex_nUMI","Gex_nGenes")]
df$value <- as.numeric(df$value)
sampleids<-unique(df$Sample)
atac_qc_sum_violin<-file.path(outputdir,"Plots",paste(prefix,"scmultiome.sample.stats.sum","pdf",sep="."))
pdf(atac_qc_sum_violin, height=20,width = 6+1*length(sampleids))
print(violin_plot_qcstats(df))
dev.off()

is.mito <- grep("^mt-", rowData(sce)$name, ignore.case=TRUE)
sce <- addPerCellQC(sce, subsets=list(Mito=is.mito))
sampleQC_df_scRNA_filtered<-as.data.frame(colData(sce))

df<-sampleQC_df_scRNA_filtered[,c("sum","detected","subsets_Mito_percent","Sample")]
#sum_threshold<-max(10000,quantile(df$sum, probs = 0.5)*2,quantile(df$sum, probs = 0.9))
#detected_threshold<-max(5000,quantile(df$detected, probs = 0.5)*2,quantile(df$detected, probs = 0.9))
#subsets_Mito_percent_threshold<-max(10,quantile(df$subsets_Mito_percent, probs = 0.5)*2,quantile(df$subsets_Mito_percent, probs = 0.9))
df_cutoff<-data.frame(variable = c("sum", "detected","subsets_Mito_percent"), cutoff = c(500, 200,20))

df<-melt(df, id = c("Sample"))
df <- merge(df,df_cutoff,by="variable")

#df<-df[df$sum<sum_threshold & df$detected<detected_threshold & df$subsets_Mito_percent<subsets_Mito_percent_threshold,]
levels(df$variable) <- c("transcripts", "genes", "mitocondria(%)")
dir.create(file.path(outputdir,"Plots"))
scRNA_qc_sum_violin<-file.path(outputdir,"Plots",paste(prefix,"scRNA.filtered.sample.stats.sum","pdf",sep="."))
sampleids<-unique(df$Sample)
pdf(scRNA_qc_sum_violin, height=12,width = 6+1*length(sampleids))
print(violin_plot_qcstats(df))
dev.off()

sampleQC_df_scRNA_filtered_sum<-file.path(outputdir,paste(prefix,"scRNA.filtered.sample.qc.sum","csv",sep="."))
sampleQC_df_scRNA_filtered %>%
  group_by(Sample) %>% 
  summarise(cells=n(),
            median_transcripts = round(median(sum),0),
            median_genes = round(median(detected),0),
            median_mito = round(median(subsets_Mito_percent),2),
            quantile_10_transcripts = quantile(sum, probs = 0.1),
            quantile_90_transcripts = quantile(sum, probs = 0.9),
            quantile_10_genes = quantile(detected, probs = 0.1),
            quantile_90_genes = quantile(detected, probs = 0.9),
            quantile_10_mito = quantile(subsets_Mito_percent, probs = 0.1),
            quantile_90_mito = quantile(subsets_Mito_percent, probs = 0.9)) %>% 
  write.table(file = sampleQC_df_scRNA_filtered_sum, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)


aavcount1<-as.data.frame(table(assay(sce,"counts")["3mir-hATOH1-hPOU4F3-3Flag-3mir&#8221;; transcript_type ",]))
colnames(aavcount1) <-c("Freq","hATOH1-hPOU4F3")
aav_gfi_counts<-assay(sce,"counts")[c("hGFI1-HA&#8221;; transcript_type ","hGFI1-HA-mScarlet&#8221;; transcript_type "),]
aavcount2<-as.data.frame(table(colSums(aav_gfi_counts) ))
colnames(aavcount2) <-c("Freq","hGFI1")

aavcount<- as.data.frame(merge(aavcount1, aavcount2, by = "Freq", all= TRUE))  %>%  arrange(as.numeric(Freq))
aavcount[is.na(aavcount)] <- 0
aavcount<-aavcount %>% mutate(across(everything(), as.character))%>% mutate(across(everything(), as.numeric))
aav_countfile<-file.path(outputdir,paste(prefix,"AAV.rawcount","csv",sep="."))
aavcount %>% write.table(file = aav_countfile, sep = "\t", quote=FALSE,row.names = FALSE,col.names = TRUE)

print(paste("Total num of cells that express hGFI1 is:",as.character(sum(aavcount[aavcount['Freq'] >0,]$hGFI1))))
print(paste("Total num of cells that express hATOH1-hPOU4F3 is:",as.character(sum(aavcount[aavcount['Freq'] >0,]$`hATOH1-hPOU4F3`))))


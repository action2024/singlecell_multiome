#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
inputpath <- args[1]
outputdir<-args[2]

library(ArchR)
library(Signac)
library(Seurat)
#BiocManager::install("EnsDb.Hsapiens.v86")
#library(EnsDb.Hsapiens.v86)
library(stringr)
library(scater)
library(EnsDb.Hsapiens.v86)
library(SingleCellExperiment)

#library(ArchRtoSignac) # {Link: GitHub https://github.com/swaruplabUCI/ArchRtoSignac}
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.22")
# BiocManager::install("biovizBase")
# devtools::install_github("swaruplabUCI/ArchRtoSignac", dependencies = TRUE)
addArchRGenome("hg38")
## Setting default genome to Hg38.
addArchRThreads(16)
## Setting default number of Parallel threads to 8.
addArchRLocking(locking = TRUE)
## Setting ArchRLocking to TRUE.
set.seed(1)

# inputpath <- c("/Users/liuji/Downloads/scmultiome/SEM01_S4D5")
#atacFiles <- inputFiles[grep(pattern = "\\.fragments.tsv.gz$", x = inputFiles)]
#rnaFiles <- inputFiles[grep(pattern = "\\.filtered_feature_bc_matrix.h5$", x = inputFiles)]
atacFiles <- list.files(path=inputpath,pattern = "atac_fragments.tsv.gz$", full.names = TRUE)
rnaFiles <- list.files(path=inputpath,pattern = "filtered_feature_bc_matrix.h5$", full.names = TRUE)
names(atacFiles)<-basename(dirname(inputpath))
names(rnaFiles)<-basename(dirname(inputpath))
# inputpath <- c("/Users/liuji/Downloads/scmultiome/Cyt49_S0C")
# rnaFiles <- c("SEM01_S4D5" = list.files(path = inputpath,pattern = "filtered_feature_bc_matrix.h5$",full.names =TRUE))
# atacFiles <- c("SEM01_S4D5" = list.files(path = inputpath,pattern = "atac_fragments.tsv.gz$",full.names =TRUE))

#create ArrowFiles from the scATAC-seq fragment files
Multiome_ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

#create an ArchRProject object from those ArrowFiles
projMulti <- ArchRProject(ArrowFiles = Multiome_ArrowFiles)

# scRNA-seq data load 
seRNA <- import10xFeatureMatrix(
  input = rnaFiles,
  names = names(rnaFiles),
  strictMatch = TRUE
)

#rescue mitcondria genes
seRNA <- import10xFeatureMatrix(
  input = rnaFiles,
  names = names(rnaFiles),
  strictMatch = TRUE,
  features = genes(EnsDb.Hsapiens.v86)
)


#  add this scRNA-seq data to our ArchRProject via the addGeneExpressionMatrix() function
#length(which(getCellNames(projMulti) %ni% colnames(seRNA)))
# length(colnames(seRNA))
cellsToKeep <- which(getCellNames(projMulti) %in% colnames(seRNA))
# length(cellsToKeep)
#keep cells that pass scRNA-seq quality control and scATAC-seq quality control 
projMulti_filtered <- subsetArchRProject(ArchRProj = projMulti, cells = getCellNames(projMulti)[cellsToKeep], 
                                         outputDirectory = outputdir, force = TRUE)
#add the gene expression data to our project
projMulti_filtered <- addGeneExpressionMatrix(input = projMulti_filtered, seRNA = seRNA,scaleTo = 10000, strictMatch = TRUE, force = TRUE)
#filter out any doublets
#projMulti_filtered <- addDoubletScores(projMulti_filtered)
projMulti_filtered <- addDoubletScores(
  input = projMulti_filtered,force = FALSE,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
projMulti_filtered <- filterDoublets(projMulti_filtered)

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
  name = "LSI_ATAC"
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
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)

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
# save GeneExpressionMatrix to sce format 
gene_expr_matrix <- getMatrixFromProject(projMulti_filtered, useMatrix = "GeneExpressionMatrix")
sce <- as(gene_expr_matrix, "SingleCellExperiment")
# add gene names to sce rowname
rownames(sce) <- rowData(sce)$name
# note: GeneExpressionMatrix is scaled and raw counts is preferred for following analysis,we proceed by adding raw counts from the seRNA loading dataset
#counts <- as.matrix(assay(sce, "GeneExpressionMatrix"))
#length(cellsToKeepIDs)
#subset raw counts from seRNA input by GeneExpressionMatrix geneIDs and cellIDs
cellsToKeepIDs<-colnames(sce)
genesToKeepIDs<-rownames(sce)
seRNA_filtered<-seRNA[genesToKeepIDs,cellsToKeepIDs]
#sce2 <- SingleCellExperiment(assays = list(counts = counts))
# add raw counts from seRNA to sce derived from GeneExpressionMatrix 1) sort by genename and cellIDs 2) add as raw counts
counts<- assay(seRNA_filtered_sorted,"data")
counts <- counts[ rownames(sce),  colnames(sce)]
assay(sce, "counts") <- counts
# add log counts
libSizes <- colSums(counts)
sizeFactors <- libSizes/mean(libSizes)
assays(sce)$logcounts <- log2(t(t(counts)/sizeFactors) + 1)
#colData(sce) <- getCellColData(projMulti_filtered)
# add dimention reduction from archR project
reducedDim(sce, "LSI_Combined") <- getReducedDims(projMulti_filtered,"LSI_Combined")
reducedDim(sce, "LSI_RNA") <- getReducedDims(projMulti_filtered,"LSI_RNA")
reducedDim(sce, "LSI_ATAC") <- getReducedDims(projMulti_filtered,"LSI_ATAC")
saveRDS(sce, file = file.path(outputdir, paste0(basename(dirname(inputpath)),".sce.rds")))


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

# plot dimensionality reductioins
p1_ATAC <- plotEmbedding(projMulti_filtered, name = "Clusters_ATAC", embedding = "UMAP_ATAC", size = 0.1)
p1 <- plotEmbedding(projMulti_filtered, name = "Clusters_Combined", embedding = "UMAP_ATAC", size = 0.1)
p2_RNA <- plotEmbedding(projMulti_filtered, name = "Clusters_RNA", embedding = "UMAP_RNA", size = 0.1)
p2 <- plotEmbedding(projMulti_filtered, name = "Clusters_Combined", embedding = "UMAP_RNA", size = 0.1)
p3 <- plotEmbedding(projMulti_filtered, name = "Clusters_Combined", embedding = "UMAP_Combined", size = 0.1)

setwd(outputdir)
plotPDF(p1_ATAC,p1, p2_RNA,p2, p3, name = "UMAP-scATAC-scRNA-Combined", addDOC = FALSE)

# p <- lapply(list(p1,p2,p3), function(x){
#   x + guides() + 
#     theme_ArchR(baseSize = 6.5) +
#     theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
#     theme(
#       axis.text.x=element_blank(), 
#       axis.ticks.x=element_blank(), 
#       axis.text.y=element_blank(), 
#       axis.ticks.y=element_blank()
#     )
# })

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

library(dplyr)
cellclustersref<-as.data.frame(getCellColData(projMulti_filtered)[c("Clusters_ATAC","Clusters_Combined")])
cellclusterscount<-cellclustersref %>% group_by_all() %>% summarise(COUNT = n())
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

p1<-ggplot(cellclusterscount, aes(fill=Clusters_ATAC, y=COUNT, x=Clusters_Combined)) + 
  geom_bar(stat="identity") + 
  theme_minimal() + scale_fill_manual(values=tableau20)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=10),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.title=element_text(size=0,face="bold"), 
        legend.title = element_blank(),
        strip.text.y =  element_blank(),
        strip.text.x = element_text(size =0, colour = "black"))   + 
  labs(title="Combined-ATAC")

p2<-ggplot(cellclusterscount, aes(fill=Clusters_Combined, y=COUNT, x=Clusters_ATAC)) + 
  geom_bar(stat="identity") + 
  theme_minimal() + scale_fill_manual(values=tableau20)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=10),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.title=element_text(size=0,face="bold"), 
        legend.title = element_blank(),
        strip.text.y =  element_blank(),
        strip.text.x = element_text(size =0, colour = "black")) +
  labs(title="ATAC-Combined")

cellclustersref<-as.data.frame(getCellColData(projMulti_filtered)[c("Clusters_RNA","Clusters_Combined")])
cellclusterscount<-cellclustersref %>% group_by_all() %>% summarise(COUNT = n())
p3<-ggplot(cellclusterscount, aes(fill=Clusters_RNA, y=COUNT, x=Clusters_Combined)) + 
  geom_bar(stat="identity") + 
  theme_minimal() + scale_fill_manual(values=tableau20)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=10),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.title=element_text(size=0,face="bold"), 
        legend.title = element_blank(),
        strip.text.y =  element_blank(),
        strip.text.x = element_text(size =0, colour = "black"))   + 
  labs(title="Combined-RNA")

p4<-ggplot(cellclusterscount, aes(fill=Clusters_Combined, y=COUNT, x=Clusters_RNA)) + 
  geom_bar(stat="identity") + 
  theme_minimal() + scale_fill_manual(values=tableau20)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=10),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.title=element_text(size=0,face="bold"), 
        legend.title = element_blank(),
        strip.text.y =  element_blank(),
        strip.text.x = element_text(size =0, colour = "black")) +
  labs(title="RNA-Combined")

plotPDF(p1, p2, p3,p4,name = "UMAP-scATAC-scRNA-Combined-stackedbar", addDOC = FALSE)

# getAvailableMatrices(projMulti_filtered)
# proj <- addGeneIntegrationMatrix(
#   ArchRProj = projMulti_filtered, 
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "IterativeLSI",
#   seRNA = seRNA,
#   addToArrow = TRUE,
#   groupRNA = "CellType",
#   nameCell = "predictedCell_Un2",
#   nameGroup = "predictedGroup_Un2",
#   nameScore = "predictedScore_Un2",
#   dimsToUse = 1:10,
#   nGenes = 250,
#   force = TRUE
# )
# seurat_object <- ArchRtoSeurat(ArchRProj = projMulti_filtered, useMatrix = "GeneExpressionMatrix")
#mat <- as.matrix(assay(gene_expr_matrix))
#seurat_mat <- CreateSeuratObject(counts = mat, project = "MyArchRProject")


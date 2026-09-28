#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
samplesheet <- args[1]
prefix<-args[2]
outputdir<-args[3]

source("/home/l128405/multiome/src/R/variables/colors.R")
source("/home/l128405/multiome/src/R/functions/violin_plot_qc.R")

#samplesheet<-"/home/l128405/multiome/samplesheet/cellranger-arc_output/samplesheet_S41.csv"
#prefix<-"S41_test"
#outputdir<-"/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen/S41_test"

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
#lapply(x, require, character.only = TRUE)

#library(ArchRtoSignac) # {Link: GitHub https://github.com/swaruplabUCI/ArchRtoSignac}
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.22")
# BiocManager::install("biovizBase")
# devtools::install_github("swaruplabUCI/ArchRtoSignac", dependencies = TRUE)
# BiocManager::install("EnsDb.Mmusculus.v79") 
## Setting default genome to mm10.
#library(EnsDb.Mmusculus.v79)
#library(EnsDb.Hsapiens.v86)
#addArchRGenome("mm10")
#addArchRGenome("hg38")
addArchRThreads(16)
## Setting default number of Parallel threads to 8.
addArchRLocking(locking = TRUE)
## Setting ArchRLocking to TRUE.
set.seed(1)

#BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")
#Full genome sequences for Mus musculus (Mouse) as provided by UCSC (genome mm39, based on assembly GRCm39) and stored in Biostrings objects.
library(BSgenome.Mmusculus.UCSC.mm39)
#BiocManager::install("TxDb.Mmusculus.UCSC.mm39.knownGene")
#TxDb--Bioconductor annotation containers that store genomic locations of transcripts, exons, CDS, and genes
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
#BiocManager::install("org.Mm.eg.db")
#Database contain mappings between Entrez Gene identifiers and GenBank accession number
library(org.Mm.eg.db)

genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Mmusculus.UCSC.mm39)
# Create annotation from TxDb
geneAnnotation <- createGeneAnnotation(
  TxDb = TxDb.Mmusculus.UCSC.mm39.knownGene, 
  OrgDb = org.Mm.eg.db,
)

#mm39gtf<-"/home/l128405/multiome/ref/GRCm39_AAVregen/../GRCm39-2024-A-reference-sources/gencode.vM33.primary_assembly.annotation.gtf"
#geneAnnotation <-createGeneAnnotation(gtf = mm39gtf)
#addArchRGenome(genome = genomeAnnotation)
library(AnnotationHub)
hub <- AnnotationHub()
# Query for Ensembl EnsDb for GRCm39
mm39_edbs <- query(hub, c("EnsDb", "Mus musculus", "GRCm39"))
# Select the most recent (last) entry
edb <- mm39_edbs[[length(mm39_edbs)]]



#inputpath <- c("/lrlhps/scratch/l128405/multiome/cellranger_count/GRCm39_AAVregen/S41A/outs")
#outputdir <- c("/lrlhps/scratch/l128405/multiome/analysis/GRCm39_AAVregen/S41A")

dir.create(outputdir)
setwd(outputdir)
#atacFiles <- inputFiles[grep(pattern = "\\.fragments.tsv.gz$", x = inputFiles)]
#rnaFiles <- inputFiles[grep(pattern = "\\.filtered_feature_bc_matrix.h5$", x = inputFiles)]


# inputpath <- c("/lrlhps/scratch/l128405/multiome/cellranger_count/GRCm39_AAVregen/S38B/outs",
#                 "/lrlhps/scratch/l128405/multiome/cellranger_count/GRCm39_AAVregen/S39B/outs",
#                 "/lrlhps/scratch/l128405/multiome/cellranger_count/GRCm39_AAVregen/S40B/outs",
#                 "/lrlhps/scratch/l128405/multiome/cellranger_count/GRCm39_AAVregen/S41B/outs")

samples <- read.csv(samplesheet,
                    header = FALSE, sep = ",")
inputpath<-samples$V1
atacFiles <- list.files(path=inputpath,pattern = "atac_fragments.tsv.gz$", full.names = TRUE)
rnaFiles <- list.files(path=inputpath,pattern = "filtered_feature_bc_matrix.h5$", full.names = TRUE)
names(atacFiles)<-basename(dirname(inputpath))
names(rnaFiles)<-basename(dirname(inputpath))


#create ArrowFiles from the scATAC-seq fragment files
Multiome_ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  minTSS = 4,
  minFrags = 1000,
  maxFrags = 1e+05,
  minFragSize = 10,
  maxFragSize = 2000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
)


#create an ArchRProject object from those ArrowFiles
projMulti <- ArchRProject(ArrowFiles = Multiome_ArrowFiles,  geneAnnotation = geneAnnotation,
                          genomeAnnotation = genomeAnnotation)
#getCellColData(projMulti)

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
  features = genes(edb)
)
#### seRNA could have duplicate genes!!! 
# sum(duplicated(rownames(seRNA)))
#rownames(seRNA)
#levels(seqnames(seRNA))
#rownames(seRNA)[which(grepl("^mt-", rownames(seRNA)))]
#  add this scRNA-seq data to our ArchRProject via the addGeneExpressionMatrix() function
#length(which(getCellNames(projMulti) %ni% colnames(seRNA)))
# length(colnames(seRNA))
cellsToKeep <- which(getCellNames(projMulti) %in% colnames(seRNA))
# length(cellsToKeep)
#keep cells that pass scRNA-seq quality control and scATAC-seq quality control 
# projMulti_filtered <- subsetArchRProject(ArchRProj = projMulti, cells = getCellNames(projMulti)[cellsToKeep], 
#                                          outputDirectory = outputdir, force = TRUE)
projMulti_filtered<-projMulti[cellsToKeep,]
#getCellColData(projMulti_filtered)
#add the gene expression data to our project
projMulti_filtered <- addGeneExpressionMatrix(input = projMulti_filtered, seRNA = seRNA,chromSizes = getChromSizes(projMulti_filtered),
                                              excludeChr = c("chrM", "chrY"),scaleTo = 10000, strictMatch = TRUE, force = TRUE)
#filter out any doublets
#projMulti_filtered <- addDoubletScores(projMulti_filtered)
projMulti_filtered <- addDoubletScores(
  input = projMulti_filtered,force = FALSE,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

projMulti_filtered <- filterDoublets(projMulti_filtered)
#save project to outputdir
projMulti_filtered <- saveArchRProject(ArchRProj = projMulti_filtered, outputDirectory = outputdir, load = TRUE)

# export scRNA, add QC and save single cell experiment

#sce <- as(seRNA, "SingleCellExperiment")
#subset raw counts from seRNA input by GeneExpressionMatrix geneIDs and cellIDs
#cellsToKeepIDs <-getCellNames(projMulti_filtered)
#sce<-sce[,cellsToKeepIDs]
# Remove duplicate gene names, keeping the first
#sce <- sce[!duplicated(rownames(sce)), ]
# add raw counts from seRNA to sce derived from GeneExpressionMatrix 1) sort by genename and cellIDs 2) add as raw counts
#counts<- assay(sce,"data")
#assay(sce, "counts") <- counts
#get per-cell quality: mito, transcripts, genes based on the counts

# 2. Read all files and create a combined count matrix
library(hdf5r)
rnaFiles <- list.files(path=inputpath,pattern = "filtered_feature_bc_matrix.h5$", full.names = TRUE)
names(rnaFiles)<-basename(dirname(inputpath))
sample_names <-basename(dirname(inputpath))


# 3. Read RNA matrices into a list using an lapply loop
rna_counts_list <- lapply(seq_along(rnaFiles), function(i) {
  
  # Read the full multiome matrix
  full_matrix <- Read10X_h5(rnaFiles[i])
  
  # Extract ONLY the Gene Expression slice
  rna_matrix <- full_matrix$`Gene Expression`
  
  # Prepend the true sample name to the barcodes to prevent duplicates later
  # Example: "AAACCCAAGCGTATGG-1" becomes "SampleName_AAACCCAAGCGTATGG-1"
  colnames(rna_matrix) <- paste(sample_names[i], colnames(rna_matrix), sep = "#")
  
  return(rna_matrix)
})
merged_counts <- Reduce(Matrix::cbind2, rna_counts_list)
# 2. Extract cell barcodes to build the cell metadata table (colData)
cell_names <- colnames(merged_counts)

# 3. Pull the clean sample names back out of the prefixed cell names
# (e.g., converts "Sample1_AAAC..." into "Sample1")
#cell_samples <- sub("_[^_]+$", "", cell_names) 
cell_samples <- sub("#.*", "", cell_names)
# 4. Construct the colData data frame
col_data <- DataFrame(
  Sample = cell_samples,
  row.names = cell_names
)

# 5. Create and save into the SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(counts = merged_counts),
  colData = col_data
)
cellsToKeepIDs <-getCellNames(projMulti_filtered)
sce<-sce[,cellsToKeepIDs]
is.mito <- grep("^mt-", rownames(sce), ignore.case=TRUE)


# counts_mat <- assay(seRNA)
# feature_names <- rowData(seRNA)$name # Extracts gene names
# rownames(counts_mat) <- make.unique(feature_names)
# seurat_obj <- CreateSeuratObject(
#   counts = counts_mat)
# sce <- as.SingleCellExperiment(seurat_obj)
# 
# gene_expr_matrix <- getMatrixFromProject(projMulti_filtered, useMatrix = "GeneExpressionMatrix")
# counts_matrix <- assay(gene_expr_matrix)
# rownames(counts_matrix) <- rowData(gene_expr_matrix)$name
# counts_matrix <- as(counts_matrix, "CsparseMatrix")
# cell_metadata <- as.data.frame(colData(gene_expr_matrix))
# seurat_obj <- CreateSeuratObject(counts = counts_matrix, meta.data = cell_metadata)
# sce <- as.SingleCellExperiment(seurat_obj)
# 
#counts(seRNA)["Tmc1", "S41B#ATTACGTCATGCATAT-1"]
#assay(seRNA["Tmc1", "S41B#ATTACGTCATGCATAT-1"],"data")
#counts(sce)["Tmc1", "S41B#ATTACGTCATGCATAT-1"]


#colData(sce) <- colData(gene_expr_matrix)
sce <- addPerCellQC(sce, subsets=list(Mito=is.mito))
saveRDS(sce, file = file.path(outputdir, paste0(prefix,".unnormalized.sce.rds")))

sampleQC_df_RNA<-as.data.frame(colData(sce))

## add qc prior to filtering###
scRNA_qc_sum<-file.path(outputdir,paste(prefix,"scRNA.sample.qc.sum","csv",sep="."))
sampleQC_df_RNA %>%
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
  write.table(file = scRNA_qc_sum, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)

sampleQC_df_multiome<-as.data.frame(getCellColData(projMulti_filtered))
scmultiome_qc_sum<-file.path(outputdir,paste(prefix,"scmultiome.sample.qc.sum","csv",sep="."))
sampleQC_df_multiome %>%
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
  write.table(file = scmultiome_qc_sum, sep = "\t", quote=FALSE,row.names = FALSE, col.names=TRUE)

df<-sampleQC_df_RNA[,c("sum","detected","subsets_Mito_percent","Sample")]
df_cutoff<-data.frame(variable = c("sum", "detected","subsets_Mito_percent"), cutoff = c(500, 200,20))

#sum_threshold<-20000
sum_threshold<-max(10000,quantile(df$sum, probs = 0.5)*2,quantile(df$sum, probs = 0.9))
#detected_threshold<-10000
detected_threshold<-max(5000,quantile(df$detected, probs = 0.5)*2,quantile(df$detected, probs = 0.9))
subsets_Mito_percent_threshold<-max(20,quantile(df$subsets_Mito_percent, probs = 0.5)*2,quantile(df$subsets_Mito_percent, probs = 0.9))
#df<-df[df$sum<sum_threshold & df$detected<detected_threshold & df$subsets_Mito_percent<subsets_Mito_percent_threshold,]

df<-melt(df, id = c("Sample"))
df <- merge(df,df_cutoff,by="variable")
levels(df$variable) <- c("transcripts", "genes", "mitocondria(%)")
dir.create(file.path(outputdir,"Plots"))
scRNA_qc_sum_violin<-file.path(outputdir,"Plots",paste(prefix,"scRNA.sample.stats.sum","pdf",sep="."))
sampleids<-unique(df$Sample)
pdf(scRNA_qc_sum_violin, height=12,width = 6+1*length(sampleids))
print(violin_plot_qcstats(df))
dev.off()

#sum_threshold<-20000
#TSSEnrichment_threshold<-max(50,quantile(sampleQC_df_RNA$TSSEnrichment, probs = 0.5)*2,quantile(sampleQC_df_RNA$TSSEnrichment, probs = 0.9))
#detected_threshold<-10000
#nFrags_threshold<-max(20000,quantile(sampleQC_df_RNA$nFrags, probs = 0.5)*2,quantile(sampleQC_df_RNA$nFrags, probs = 0.9))
#nGenes_threshold<-max(5000,quantile(sampleQC_df_RNA$Gex_nGenes, probs = 0.5)*2,quantile(sampleQC_df_RNA$Gex_nGenes, probs = 0.9))
#nTranscripts_threshold<-max(10000,quantile(sampleQC_df_RNA$Gex_nUMI, probs = 0.5)*2,quantile(sampleQC_df_RNA$Gex_nUMI, probs = 0.9))

df<-sampleQC_df_multiome[,c("TSSEnrichment","nFrags","Gex_nGenes","Gex_nUMI","Sample")]
df_cutoff<-data.frame(variable = c("TSSEnrichment", "nFrags","Gex_nGenes","Gex_nUMI"), cutoff = c(4,1000, 200,500))
#df<-sampleQC_df_RNA[sampleQC_df_RNA$TSSEnrichment<TSSEnrichment_threshold & sampleQC_df_RNA$nFrags<nFrags_threshold & sampleQC_df_RNA$Gex_nGenes<nGenes_threshold & sampleQC_df_RNA$Gex_nUMI<nTranscripts_threshold,c("Sample","TSSEnrichment","nFrags","Gex_nUMI","Gex_nGenes")]
df_melt<-melt(df, id.vars = c("Sample"))
df <- merge(df_melt,df_cutoff,by="variable")

df$value <- as.numeric(df$value)
sampleids<-unique(df$Sample)
atac_qc_sum_violin<-file.path(outputdir,"Plots",paste(prefix,"scmultiome.sample.stats.sum","pdf",sep="."))
pdf(atac_qc_sum_violin, height=20,width = 6+1*length(sampleids))
print(violin_plot_qcstats(df))
dev.off()

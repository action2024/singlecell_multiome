args <- commandArgs(trailingOnly = TRUE)
inputrds <- args[1]
outputdir<-args[2]
group<-args[3]

library(devtools)
library(ggplot2) # for plotting
library(ggpmisc)
library(reshape2)
library(msigdbr) # for gathering gene sets
library(Seurat)
# library(future) # for parallel computing
# library(future.apply) # for parallel computing
library(dplyr)
#library(patchwork)
#library("Matrix")
#library("readr")
#library("DAPAR")
library(scuttle)
#library(robustbase)
library(scater)
library(scran)
library(scRNAseq)
library(cluster)
library(dendextend)
library(AUCell)
library(clustree)
library(PCAtools)
library(celldex)
library(SC3)
###
library("ggalluvial")
library(future)  # for parallel computing
library(future.apply)  # for parallel computing
#install.packages("anticlust", version = "0.6.1")
library("fastcluster") ## fast clustering function hclust to replace stats::hclust(...)
library("parallelDist") ## Calculates distance matrices in parallel with parDist to replace dist function 
#install.packages("parallelDist")
library(SingleCellExperiment)
#devtools::unload("gsdensity")
#devtools::unload("anticlust")
#remotes::install_version('anticlust', '0.6.1') 
#library(anticlust)
#library(gsdensity)
#library(scDataviz)
library(dittoSeq)
library(gridExtra)
library(bluster)
library(pheatmap)
library(dynamicTreeCut)
customcolor<-c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
               "#D62728", "#FF9896", "#9EDAE5", "#9467BD", "#98DF8A", "#C5B0D5",
               "#8C564B", "#E377C2", "#F7B6D2", "#7F7F7F",
               "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF","#C49C94")

inputrds<-"/cluster/home/liuji/T1D/multiome/analysis/Cyt49_S4D5/scRNAseq/Clusters_Combined.QC_filtered.normalized.rds"
outputdir<-"/cluster/home/liuji/T1D/multiome/analysis/Cyt49_S4D5/scRNAseq"
group<-"Clusters_Combined"
sce<-readRDS(inputrds)
######normalization######
####technical differences in cDNA capture or PCR amplification efficiency across cells due to minimal starting material###
#cells in each cluster are normalized separately and the size factors are 
# rescaled to be comparable across clusters
#quickCluster() will use an approximate algorithm for PCA
clust.sce <- quickCluster(sce) 
#table(clust.sce)
sce <- computeSumFactors(sce, cluster=clust.sce, min.mean=0.1)
sce <- logNormCounts(sce)
#assayNames(sce)

#any(duplicated(rownData(sce)))
normalized.rds<-file.path(outputdir,paste(group,"QC_filtered.normalized","rds",sep="."))
saveRDS(sce, normalized.rds)

#Quantifying per-gene variation
dec.sce <- modelGeneVar(sce)
# Visualizing the fit:
fit.sce <- metadata(dec.sce)

mean_variance_fig<-file.path(outputdir,paste(group,"mean_variance","png",sep="."))
png(mean_variance_fig, width=6, height=6, units="in", res=500)
plot(fit.sce$mean, fit.sce$var, xlab="Mean of log-expression",ylab="Variance of log-expression")
curve(fit.sce$trend(x), col="dodgerblue", add=TRUE, lwd=2)
dev.off()

#remove rowdata to initiate 
rowData(sce) <- NULL

##dec.sce[order(dec.sce$bio,decreasing=TRUE),]
# Taking the top 200 genes here:
hvg.dec.sce.var <- getTopHVGs(dec.sce, n=200)
# remove mitocondria related genes from hvgs
hvg.dec.sce.var <- hvg.dec.sce.var[which(!grepl("^MT-", hvg.dec.sce.var))]
#sce_hvg.dec.sce.var<-sce[hvg.dec.sce.var,]

#select hvgs - top 10% of the whole set
chosen <- getTopHVGs(dec.sce, prop=0.1)
outfile<-file.path(outputdir,paste(group,"hvgs.top10percent.genelist","csv",sep="."))
write.table(chosen, file=outfile, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE, quote=FALSE)
### process the chosen dataset for downstream analysis
#length(chosen)
sce.hvg <- sce[chosen,]
hvg.rds<-file.path(outputdir,paste(group,"hvgs.top10percent","rds",sep="."))
saveRDS(sce.hvg, hvg.rds)

#dim(sce.hvg)
#colData(sce.hvg)

###Dimensionality reduction#####
set.seed(10000)
sce <- runPCA(sce, subset_row=chosen)
#reducedDimNames(sce)

# use the technical component estimates to determine the proportion of variance that should be retained. 
set.seed(001001001)
sce <- denoisePCA(sce, subset.row=chosen, technical=dec.sce)
#ncol(reducedDim(sce, "PCA"))
percent.var <- attr(reducedDim(sce, "PCA"), "percentVar")
chosen.elbow <- findElbowPoint(percent.var)

#plot chosen elbow point for PCA analysis
pca_elbow_fig<-file.path(outputdir,paste(group,"PCA.elbow","png",sep="."))
png(pca_elbow_fig, width=10, height=10, units="in", res=500)
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
abline(v=max(5,chosen.elbow), col="red")
dev.off()

#https://bodenmillergroup.github.io/IMCDataAnalysis/single-cell-visualization.html
#plot the distribution of samples of PCA analysis
pca_fig<-file.path(outputdir,paste(group,"sample.PCA","png",sep="."))
png(pca_fig, width=15, height=10, units="in", res=500)
plotReducedDim(sce, dimred="PCA", colour_by=group,point_size=0.2, text_by=group)+ 
  guides(colour = guide_legend(override.aes = list(size=4)))  + 
  theme(legend.title = element_blank())
dev.off()

# # plot the distribution of samples of PCA analysis, without label for each sample 
# pca_plot<-plotReducedDim(sce, dimred="PCA", colour_by=group,point_size=0.1, text_by=group)+ 
#   guides(colour = guide_legend(override.aes = list(size=4)))  + 
#   theme(legend.title = element_blank())
# pca_fig<-file.path(outputdir,paste(group,"sample.nolabel.PCA","png",sep="."))
# png(pca_fig, width=15, height=10, units="in", res=500)
# pca_plot
# dev.off()

samplenum<-length(unique(colData(sce)[[group]]))

# plot the distribution of samples of PCA analysis separately
# 1) pca plot of all samples 2) loop through each sample 3) arrange plots
sample_plots<-list()
sample_plots<-c(sample_plots,list(plotReducedDim(sce, dimred="PCA", colour_by=group,point_size=0.1, 
                                                 text_by=group,text_size=1)))
for(sample in unique(colData(sce)[[group]])){
  #sample<-'S5D4_R01'
  p<-plotReducedDim(sce[,colData(sce)[[group]]==sample], "PCA", 
                    point_size=0.1)+ 
    guides(colour = guide_legend(override.aes = list(size=4)))  + 
    theme(legend.title = element_blank())
  sample_plots<-c(sample_plots,list(p))
}

pca_fig<-file.path(outputdir,paste(group,"per-sample.PCA","png",sep="."))
png(pca_fig, width=3*round(sqrt(samplenum),0)+9, height=3*round(sqrt(samplenum),0), units="in", res=500)
grid.arrange(grobs = sample_plots, round(sqrt(samplenum),0))## display plot
dev.off()

# plot expression selected genes (top 12 hvgs) - pca plot
# 1) loop through each gene 2) arrange plots
selected_genes<-head(hvg.dec.sce.var,24)
sample_plots<-list()
for(gene in selected_genes){
  #gene<-'DLK1'
  #print(gene)
  p<-plotReducedDim(sce, "PCA", 
                    colour_by=gene,point_size=0.1) 
  sample_plots<-c(sample_plots,list(p))
}

pca_fig<-file.path(outputdir,paste(group,"hvgs.PCA.expression","png",sep="."))
png(pca_fig, width=18, height=12, units="in", res=500)
grid.arrange(grobs = sample_plots, 3,0)## display plot
dev.off()

# plot expression selected genes (top 24 hvgs) -violin plot
hvgs_expression_fig<-file.path(outputdir,paste(group,"hvgs.expression.violin","png",sep="."))
png(hvgs_expression_fig, width=4+samplenum*2, height=10, units="in", res=500)
plotExpression(sce, features=selected_genes,
               x=I(colData(sce)[[group]]),   ncol = 3,point_size=0.01)  + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) 
dev.off()

#save pca sce
PCA.rds<-file.path(outputdir,paste(group,"PCA","rds",sep="."))
saveRDS(sce, PCA.rds)

#sce_count <-  data.frame(table(colnames(sce)))
#dim(sce_count[sce_count$Freq > 1,])
#sce_count <-  data.frame(table(rownames(colData(sce))))

#colnames(sce)<-rownames(colData(sce))
#stopifnot(identical(rownames(colData(sce)), colnames(assay(sce,"counts"))))
#sce_1<-sce



#### remove artifacts caused by overlapping row/cell IDs (not sure of the source) #### 
# rename duplicate rows in the sce 
colnames(sce) = make.names(colnames(sce), unique=TRUE)

# plot expression selected genes (top 200 hvgs) -heatmap plot
hvg_heatmap_fig<-file.path(outputdir,paste(group,"hvgs.heatmap","png",sep="."))
png(hvg_heatmap_fig, width=3+samplenum*1, height=24, units="in", res=500)
plotHeatmap(sce, features=head(hvg.dec.sce.var,200), order_columns_by=group,
            center=TRUE, symmetric=TRUE, zlim=c(-3, 3)) 
dev.off()

# plot expression selected genes (top 12 hvgs) -PCA dot plot
hvgs_pca_fig<-file.path(outputdir,paste(group,"hvgs.expression","png",sep="."))
png(hvgs_pca_fig, width=12, height=6, units="in", res=500)
multi_dittoDimPlot(sce,selected_genes, reduction.use = "PCA",size=0.01)
dev.off()

# plot the gene expression of cell cyle gene
#G1_phase<-c("CCND1","E2F1")
#S_phase<-c("MCM2","PCNA")
#G2_phase<-c("CCND1","CCNB1")
#M_phase<-c("CCND1","E2F1")
#cyclin.genes <- rownames(sce)[which(grepl("^CCN[ABDE][0-9]$", rownames(sce)))]
#plotHeatmap(sce, order_columns_by="ident", cluster_rows=FALSE, features=sort(cyclin.genes))


#inputrds<-'/cluster/home/liuji/T1D/scRNAseq/analysis/R01/S5/S5D4_R01.PCA.rds'
#inputrds<-'/Users/liuji/Downloads/R01/S5D4_R01.PCA.rds'
samplenum<-length(unique(colData(sce)[[group]]))

### https://bioconductor.org/books/3.18/OSCA.advanced/clustering-redux.html
# assess the separation of clusters/ days -- Silhouette width
sil.approx <- approxSilhouette(reducedDim(sce, "PCA"), clusters=sce[[group]])
sil.data <- as.data.frame(sil.approx)
#head(sil.data[sil.data$cluster=="S5D2_R01",],29)
sil.data$closest <- factor(ifelse(sil.data$width > 0, sce[[group]], sil.data$other))
sil.data$cluster <- sce[[group]]
sil.sum<-table(sil.data$cluster, sil.data$closest)


sil.plot<-ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley",size=0.03)+
  theme_classic() +
  theme(axis.text.x=element_text(size=14, angle = 20, hjust=1),
        axis.text.y=element_text(size=14),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size =14, colour = "black"),
        strip.text.y = element_text(size =14, colour = "black"),
        legend.title = element_blank()) + scale_color_manual(values=customcolor) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  stat_summary(fun = median, geom = "point", shape = 23, size =2, color = "red")

sil_plot<-file.path(outputdir,paste(group,"sil","png",sep="."))
png(sil_plot, width=4+samplenum*2, height=6, units="in", res=500)
sil.plot
dev.off()


rmsd <- clusterRMSD(reducedDim(sce, "PCA"), sce[[group]])
rmsd_barplot<-barplot(rmsd, ylab="RMSD", xlab="")
rmsd_plot<-file.path(outputdir,paste(group,"rmsd","png",sep="."))
png(rmsd_plot, width=4+samplenum*2, height=6, units="in", res=500)
rmsd_barplot
dev.off()



### clustering: nearest neighbor graph 
#nn.clusters <- clusterCells(sce, use.dimred="PCA",assay.type="logcounts")
#table(nn.clusters)
#colLabels(sce) <- nn.clusters
#tsne_plot<-file.path(outputdir,paste(group,"TSNE","png",sep="."))
#png(tsne_plot, width=12, height=6, units="in", res=500)
#plotReducedDim(sce, dimred="PCA", colour_by="label", text_by="label")
#dev.off()

#### hierachial clustering 
#sce<-readRDS("/Users/liuji/Downloads/R01/S3.PCA.rds")
#mat <- reducedDim(sce, "PCA")
#dim(mat)
#hp <- HclustParam(clust.fun=fastcluster::hclust,method="ward.D2", cut.dynamic=TRUE,dist.fun=parDist)
#hclust.out <- clusterRows(mat, hp)
#hclust.dyn <- clusterCells(sce, use.dimred="PCA",BLUSPARAM=HclustParam(method="ward.D2", cut.dynamic=TRUE, dist.fun=parDist,clust.fun=fastcluster::hclust,cut.params=list(minClusterSize=10, deepSplit=1)))

#hclust.dyn <- clusterCells(sce, use.dimred="PCA",BLUSPARAM=HclustParam(method="ward.D2", cut.dynamic=TRUE, cut.params=list(minClusterSize=10, deepSplit=1)))
#table(hclust.dyn)
#sce$hclust <- factor(hclust.dyn)
#colLabels(sce) <- hclust.dyn

hclust.dyn <- clusterCells(sce, use.dimred="PCA",BLUSPARAM=HclustParam(method="ward.D2", cut.dynamic=TRUE, 
                                                                           cut.params=list(minClusterSize=10, deepSplit=1)),full=TRUE)
#table(hclust.dyn)
colLabels(sce) <- factor(hclust.dyn$clusters)
sce$hclust <- factor(hclust.dyn$clusters)

#quantify the proportion of clusters
hclust_cluster<-as.data.frame(table(hclust.dyn$clusters))
colnames(hclust_cluster)<-c("cluster","cellnum")
hclust_cluster$cellprop<-round(hclust_cluster$cellnum/sum(hclust_cluster$cellnum)*100,1)
outfile<-file.path(outputdir, paste(group, "hclust.num","csv",sep="."))
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
                 groupLabels=TRUE) %>%
  set("branches_lwd", 1) 
#dend_data <-dendro_data(dend, type = "rectangle")

hclustdend_plot<-file.path(outputdir,paste(group,"hclust.dend","png",sep="."))
png(hclustdend_plot, width=12, height=6, units="in", res=500)
ggd1 <- as.ggdend(dend)
ggplot(ggd1, labels=FALSE) 
dev.off()

#head(colData(sce),5)
hclust_plot<-file.path(outputdir,paste(group,"hclust.PCA","png",sep="."))
png(hclust_plot, width=12, height=6, units="in", res=500)
plotReducedDim(sce, dimred="PCA", colour_by="label",point_size=0.2,text_by="label") +
  guides(colour = guide_legend(override.aes = list(size=4)))  + 
  theme(legend.title = element_blank()) + scale_color_manual(values=customcolor)
dev.off()

#save clustering sce
clust.rds<-file.path(outputdir,paste(group,"clust","rds",sep="."))
saveRDS(sce, clust.rds)

#quantify the distribution of clusters across each day

hclust_ident<-as.data.frame(table(colData(sce)[,c(group,"label")]))
hclust_ident <- hclust_ident %>%
  group_by(!!sym(group)) %>%
  mutate(
    cellnum = sum(Freq),
    cellprop = round((Freq / cellnum) * 100,1)
  )
hclust_ident <- hclust_ident[order(hclust_ident[[group]], decreasing = FALSE),]


outfile<-file.path(outputdir, paste(group, "hclust.distribution","csv",sep="."))
write.csv(hclust_ident, outfile,row.names = FALSE,quote=FALSE)

stacked_barplot<-file.path(outputdir,paste(group,"hclust.distribution","png",sep="."))
p1<-ggplot(hclust_ident, aes(fill=!!sym(group), y=Freq, x=label)) + 
  geom_bar(stat="identity") + 
  theme_minimal() + scale_fill_manual(values=customcolor)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=10),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.title=element_text(size=0,face="bold"), 
        legend.title = element_blank(),
        strip.text.y =  element_blank(),
        strip.text.x = element_text(size =0, colour = "black")) 
p2<-ggplot(hclust_ident, aes(fill=label, y=Freq, x=!!sym(group))) + 
  geom_bar(stat="identity") + 
  geom_text(data = subset(hclust_ident, cellprop > 0), aes(label = paste0(cellprop, "%")),
            position = position_stack(vjust = 0.5),
            color = "white") +
  theme_minimal() + scale_fill_manual(values=customcolor)+
  theme(axis.text.x=element_text(size=15,angle=29, hjust=1),
        axis.text.y=element_text(size=10),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.title=element_text(size=0,face="bold"), 
        legend.title = element_blank(),
        strip.text.y =  element_blank(),
        strip.text.x = element_text(size =0, colour = "black")) 
png(stacked_barplot, width=12, height=6, units="in", res=500)
p1+p2
dev.off()


# assess the separation of clusters -- Silhouette width
sil.approx <- approxSilhouette(reducedDim(sce, "PCA"), clusters=sce$label)
sil.data <- as.data.frame(sil.approx)
#head(sil.data)
sil.data$closest <- factor(ifelse(sil.data$width > 0, sce$label, sil.data$other))
sil.data$cluster <- sce$label
sil.sum<-table(sil.data$cluster, sil.data$closest)


sil.plot<-ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley",size=0.03)+
  theme_classic() +
  theme(axis.text.x=element_text(size=14, angle = 20, hjust=1),
        axis.text.y=element_text(size=14),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size =14, colour = "black"),
        strip.text.y = element_text(size =14, colour = "black"),
        legend.title = element_blank()) + scale_color_manual(values=customcolor) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  stat_summary(fun = median, geom = "point", shape = 23, size =2, color = "red")

sil_plot<-file.path(outputdir,paste(group,"clust.sil","png",sep="."))
png(sil_plot, width=4+samplenum*2, height=6, units="in", res=500)
sil.plot
dev.off()


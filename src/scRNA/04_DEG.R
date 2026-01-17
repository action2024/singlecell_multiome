#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
inputrds <- args[1]
outputdir<-args[2]
group<-args[3]
#inputrds<-"/Users/liuji/Downloads/Cyt49_S4D5/Cyt49_S4D5/sce.rds"
#inputrds<-"/cluster/home/liuji/T1D/multiome/analysis/Cyt49_S4D5/scRNAseq/Clusters_Combined.QC_filtered.normalized.rds"
#group<-"Clusters_Combined"
library(pathfindR)
library(pheatmap)
library(scater)
library(scran)
library(SingleCellExperiment)
library(scuttle) # For logNormCounts
library(dplyr)
library(glmGamPoi)
library(reshape2)

customcolor<-c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
               "#D62728", "#FF9896", "#9EDAE5", "#9467BD", "#98DF8A", "#C5B0D5",
               "#8C564B", "#E377C2", "#F7B6D2", "#7F7F7F",
               "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF","#C49C94")

sce<-readRDS(inputrds)

#is.mito <- any(seqnames(location)=="MT")
#rownames(sce)[is.mito]
#colData(sce)$orig.ident
###### marker detection
# AUC represents the probability that a randomly chosen observation from our cluster of interest is greater than a randomly chosen observation from the other cluster.
# AUC: 1- upregulation, 0.5 - no change, 0 - downregulation
# Cohen’s d: number of standard deviations that separate the means of the two groups
# Cohen’s d: positive - upregulated, negative - downregulation, zero - little difference
# logFC.detected: log-fold change in the proportion of cells with detected expression between clusters

#function to identify pathways in the gene lists
pathfind<-function(selected.vgs,genesets){
  #selected.vgs<-data.frame(up.regulated)
  selected.vgs<-data.frame(selected.vgs)
  output_selected.vgs <- run_pathfindR(selected.vgs,gene_sets=genesets,min_gset_size=3)
  #output_selected.vgs.go.mf <- run_pathfindR(selected.vgs,gene_sets="Reactome",min_gset_size=3)
  #enrichment_chart(result_df = output_selected.vgs,top_terms = 15)
  #head(output_selected.vgs,3)
  input_processed <- input_processing(selected.vgs,p_val_threshold = 0.05)
  if(nrow(output_selected.vgs)>0){
  #gg_list <- visualize_terms(result_df = output_topranked,input_processed = input_processed)  # this function returns a list of ggraph objects (named by Term ID)
  output_topranked_clustered <- cluster_enriched_terms(output_selected.vgs, plot_dend = FALSE, plot_clusters_graph = FALSE)
  #head(output_topranked_clustered)
  #output_topranked_clustered[output_topranked_clustered$ID=="hsa04512",]
  #output_topranked_clustered.mf <- cluster_enriched_terms(output_selected.vgs.go.mf, plot_dend = FALSE, plot_clusters_graph = FALSE)
  #selected_clusters.mf <- output_topranked_clustered.mf
  #selected_clusters <- output_topranked_clustered[output_topranked_clustered$Status=="Representative",]
  return(output_topranked_clustered)
}}

## DE gene analysis by cluster 
genelist_quant<-function(sce,genelist){
  gene_sum<-data.frame()
  for(gene in genelist){
    #gene<-"FOXP3"
    gene_index <- which(rownames(assays(sce)$logcounts) == gene)
    gene_expression <- assays(sce)$logcounts[gene_index, ]
    gene_expression_minquantile <- round(quantile(gene_expression[gene_expression>0], probs = c(0.25))[[1]],1)
    for(cluster in unique(colData(sce)[[group]])){
      # cluster<-"C3"
      cell_index <- rownames(colData(sce)[colData(sce)[[group]]==cluster,])
      # Extract the expression data for 'my_gene'
      gene_expression <- assays(sce)$logcounts[gene_index, cell_index]
      percentage<-round(length(gene_expression[gene_expression>gene_expression_minquantile])/length(cell_index)*100,1)
      gene_expression_sum<-data.frame(t(c(gene,cluster,gene_expression_minquantile,round(mean(gene_expression),1),round(median(gene_expression),1),percentage)))
      gene_sum<-rbind(gene_sum,gene_expression_sum)
    }
  }
  colnames(gene_sum)<-c("gene","cluster","min_cutoff","mean","median","percent")
  gene_sum[, c("min_cutoff","mean","median","percent")] <- sapply(gene_sum[, c("min_cutoff","mean","median","percent")], as.numeric)
  
  gene_order <- unique(gene_sum[order(gene_sum$percent, decreasing=TRUE),]$gene)
  gene_sum_sort<-gene_sum[order(factor(gene_sum$gene, levels = gene_order), gene_sum$cluster,decreasing=FALSE),]
  
  return(gene_sum_sort)}

for(group1 in unique(sce[[group]])){
  group2<-unique(sce[[group]])[!unique(sce[[group]]) %in% group1]
  # Example: Assign cell types manually
  colData(sce)[colData(sce)[[group]] %in% group1,"testgroup"] <- "test"
  colData(sce)[colData(sce)[[group]] %in% group2,"testgroup"] <- "control"
  sce<-sce[,colData(sce)[[group]] %in% group1 | colData(sce)[[group]] %in% group2]
  colData(sce)$rep <- NA
  colData(sce)[colData(sce)[[group]] %in% group1,]$rep<- print(sample(1:2,dim(colData(sce)[colData(sce)[[group]] %in% group1,])[1],replace=TRUE))
  colData(sce)[colData(sce)[[group]] %in% group2,]$rep<- print(sample(1:2,dim(colData(sce)[colData(sce)[[group]] %in% group2,])[1],replace=TRUE))
  reduced_sce <- pseudobulk(sce, group_by = vars(rep,!!sym(group),condition = testgroup),n_cells = n())
  
  # Use DESeq2's size factor calculation procedure
  fit <- glm_gp(reduced_sce, design = ~ condition, size_factor = "ratio", verbose = TRUE)
  res <- test_de(fit, contrast = cond(condition = "test") - cond(condition = "control"),sort_by = "pval")
  #filter by low confidence and low fold change
  res_toprank.up<-res[res$pval<0.05 & res$lfc >2,]
  #sort by most-prevalently expressed genes
  res_toprank.up.sortedsum<-genelist_quant(sce,res_toprank.up$name)
  #res_toprank.up<-res_toprank.up[order(res_toprank.up$lfc, decreasing = TRUE), ]
  res_toprank.up<-res_toprank.up[order(factor(res_toprank.up$name, levels = unique(res_toprank.up.sortedsum$gene))),]
  #reset index 
  rownames(res_toprank.up) <- NULL
  outfile<-file.path(outputdir,paste(group,group1,"upregulated.topranked.DEgene","csv",sep="."))
  write.table(res_toprank.up, file=outfile, row.names =FALSE, col.names = TRUE,sep = ",", append = FALSE, quote=FALSE)
  
  outfile<-file.path(outputdir,paste(group,group1,"upregulated.sorted.DEgene","csv",sep="."))
  write.table(res_toprank.up.sortedsum, file=outfile, row.names =FALSE, col.names = TRUE,sep = ",", append = FALSE, quote=FALSE)
  
  #call pathways
  up.regulated<-res_toprank.up[,c("name","lfc","pval")]
  up.regulated.reactome<-pathfind(up.regulated,"Reactome")
  if(!is.null(up.regulated.reactome)){
    outfile<-file.path(outputdir, paste(group,group1, "upregulated.Reactome.logFC.pathway","csv",sep="."))
    write.csv(up.regulated.reactome, outfile,row.names = FALSE)
  }
   #up.regulated.goall<-pathfind(up.regulated,"GO-All")
  up.regulated.kegg<-pathfind(up.regulated,"KEGG")
 # if(exists(deparse(substitute(up.regulated.kegg)))){
  if(!is.null(up.regulated.kegg)){
  outfile<-file.path(outputdir, paste(group,group1, "upregulated.KEGG.logFC.pathway","csv",sep="."))
  write.csv(up.regulated.kegg, outfile,row.names = FALSE)}
  
  res_toprank.down<-res[res$pval<0.05 & res$lfc < -2,]
  #sort by most-prevalently expressed genes
  res_toprank.down.sortedsum<-genelist_quant(sce,res_toprank.down$name)
  #res_toprank.up<-res_toprank.up[order(res_toprank.up$lfc, decreasing = TRUE), ]
  res_toprank.down<-res_toprank.down[order(factor(res_toprank.down$name, levels = unique(res_toprank.down.sortedsum$gene))),]
  #reset index 
  rownames(res_toprank.down) <- NULL
  #res_toprank.down<-res_toprank.down[order(res_toprank.down$pval, decreasing = FALSE), ]
  outfile<-file.path(outputdir,paste(group,group1,"downregulated.topranked.DEgene","csv",sep="."))
  write.table(res_toprank.down, file=outfile, row.names =FALSE, col.names = TRUE,sep = ",", append = FALSE, quote=FALSE)
  outfile<-file.path(outputdir,paste(group,group1,"downregulated.sorted.DEgene","csv",sep="."))
  write.table(res_toprank.down.sortedsum, file=outfile, row.names =FALSE, col.names = TRUE,sep = ",", append = FALSE, quote=FALSE)
  #pathway
  down.regulated<-res_toprank.down[,c("name","lfc","pval")]
  down.regulated.kegg<-pathfind(down.regulated,"KEGG")
  if(!is.null(down.regulated.kegg)){
    outfile<-file.path(outputdir, paste(group,group1, "downregulated.KEGG.logFC.pathway","csv",sep="."))
  write.csv(down.regulated.kegg, outfile,row.names = FALSE)}
  
  down.regulated.reactome<-pathfind(down.regulated,"Reactome")
  if(!is.null(down.regulated.reactome)){
    outfile<-file.path(outputdir, paste(group,group1, "downregulated.Reactome.logFC.pathway","csv",sep="."))
  write.csv(down.regulated.reactome, outfile,row.names = FALSE)}
  
  #plot heatmap and violin plot
  for(i in 1:round(dim(res_toprank.up)[1]/100)){
    #print(i)
    #i<-1
    #print(c((i-1)*100+1,min(dim(res_toprank.up)[1],i*100)))
    start_subset<-(i-1)*100+1
    end_subset<-min(dim(res_toprank.up)[1],i*100)
    #print(res_toprank.up[start_subset:end_subset,])
    #print(length(na.omit(res_toprank.up[start_subset:end_subset,]$name)))
    genes<-na.omit(res_toprank.up[(i-1)*100+1:min(dim(res_toprank.up)[1],i*100),]$name)
    violin<-plotExpression(sce, features=genes,
                           x=I(colData(sce)[[group]]),   color=I(colData(sce)[[group]]), ncol = 3,point_size=0.01)  + 
      theme(axis.text.x = element_text(angle = 20, hjust = 1)) + facet_wrap(~Feature, scales = "free_y") + 
      scale_color_manual(values=customcolor)
    violinplot<-file.path(outputdir,paste(group,group1,toString(i),"upregulated.gene.violin","png",sep="."))
    png(violinplot, width=30, height=25, units="in", res=200)
    print(violin)
    dev.off()
    
    heatmapplot<-file.path(outputdir,paste(group,group1,toString(i),"upregulated.gene.heatmap","png",sep="."))
    png(heatmapplot, width=18, height=20, units="in", res=200)
    plotHeatmap(sce, features=genes, order_columns_by=group,center=TRUE, symmetric=TRUE, zlim=c(-2, 2)) + scale_color_manual(values=customcolor)
    dev.off()
  }}




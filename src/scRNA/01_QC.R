#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
inputrds <- args[1]
outputdir<-args[2]
group<-args[3]
#inputrds<-"/Users/liuji/Downloads/Cyt49_S4D5/Cyt49_S4D5/sce.rds"
#inputrds<-"/cluster/home/liuji/T1D/multiome/analysis/Cyt49_S4D5/clustering/Cyt49_S4D5.sce.rd"
#group<-"Clusters_Combined"
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

#find peaks for transcripts/genes/mitcondria per cell
ReturnMax<-function(numvec,threshold_lower,threshold_upper){
  # numvec<-df$detected
  numvec<-numvec[numvec>=threshold_lower & numvec<=threshold_upper]
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])}

#inputrds<-'~/Downloads/R01/S5.PCA.rds'
qc_stats<-function(sce){
  #sce<-S5D2_R01_sce
  location <- rowRanges(sce)
  is.mito <- which(grepl("^MT-", rownames(sce)))
  df <- perCellQCMetrics(sce, subsets=list(Mito=is.mito),use.altexps=TRUE)
  df$sample<-colData(sce)[[group]]
  #colData(sce)$Clusters_Combined
  #total amount of transcripts
  transcripts_sum<-summary(df$sum)
  #total amount of genes
  gene_sum<-summary(df$detected)
  #total amount of mitocondria contamination
  mitocondria_sum<-summary(df$subsets_Mito_percent)
  Mito_percent_peak<-ReturnMax(df$subsets_Mito_percent,0,10)
  mitocondria_sum$peak <- Mito_percent_peak
  sum_peak<-ReturnMax(df$sum,1000,10000)
  transcripts_sum$peak<-sum_peak
  detected_peak<-ReturnMax(df$detected,1000,10000)
  gene_sum$peak<-detected_peak
  per_sampleQC<-do.call(rbind, Map(data.frame, gene_sum=gene_sum,transcripts_sum=transcripts_sum,mitocondria_sum=mitocondria_sum))%>% 
    mutate_if(is.numeric, round, digits = 2)
  per_sampleQC$sample<-unique(colData(sce)$ident)
  return(list(per_sampleQC,df))
}

#violin plot for #genes, #transcripts, #percent of 
violin_plot_qcstats<-function(df){
  df<-sampleQC_df[,c("sum","detected","subsets_Mito_percent","sample")]
  sum_threshold<-20000
  detected_threshold<-10000
  subsets_Mito_percent_threshold<-10
  df<-df[df$sum<sum_threshold & df$detected<detected_threshold & df$subsets_Mito_percent<subsets_Mito_percent_threshold,]
  df<-melt(df, id = c("sample"))
  levels(df$variable) <- c("transcripts", "genes", "mitocondria(%)")
  #head(df)
  violin_stat<-ggplot(df, aes(x = sample,y = value,color=variable)) +
    geom_point(position = position_jitterdodge(), shape=16, size=0.1)+
    geom_violin(trim=TRUE,alpha=0.04,color="black",linewidth = 0.2,draw_quantiles = c(0.25, 0.5, 0.75))  + 
    #stat_summary(aes(x=rep_label, y=insert_length), fun.y=median, geom="point", shape=23, size=1,position = position_jitterdodge())+
    theme_minimal() + 
    theme(axis.text.x=element_text(size=14, angle = 20, hjust=1),
          axis.text.y=element_text(size=14),
          plot.title = element_text(size=14, face="bold.italic"),
          axis.title=element_text(size=14,face="bold"), 
          strip.text.x = element_text(size =14, colour = "black"),
          strip.text.y = element_text(size =14, colour = "black"),
          legend.position="none",legend.title = element_blank()) + 
    facet_grid(variable~.,scales = 'free_y') +
    scale_color_manual(values=customcolor) +
    scale_fill_manual(values=customcolor)+
    xlab("samples") + 
    ylab("")
  return(violin_stat)
}

sce<-readRDS(inputrds)
#is.mito <- any(seqnames(location)=="MT")
#rownames(sce)[is.mito]
#colData(sce)$orig.ident
sampleids<-unique(colData(sce)[[group]])
sampleQC<-data.frame()
sampleQC_df<-data.frame()

for(x in sampleids){
  #print(x)
  #x<-"C8"
  sce_subset <- sce[, colData(sce)[[group]] == x]
  location <- rowRanges(sce_subset)
  is.mito <- which(grepl("^MT-", rownames(sce_subset)))
  df <- perCellQCMetrics(sce_subset, subsets=list(Mito=is.mito))
  df$sample<-x
  #total amount of transcripts
  transcripts_sum<-summary(df$sum)
  #total amount of genes
  gene_sum<-summary(df$detected)
  #total amount of mitocondria contamination
  mitocondria_sum<-summary(df$subsets_Mito_percent)
  Mito_percent_peak<-ReturnMax(df$subsets_Mito_percent,0,10)
  mitocondria_sum$peak <- Mito_percent_peak
  sum_peak<-ReturnMax(df$sum,1000,10000)
  transcripts_sum$peak<-sum_peak
  detected_peak<-ReturnMax(df$detected,1000,10000)
  gene_sum$peak<-detected_peak
  per_sampleQC<-do.call(rbind, Map(data.frame, gene_sum=gene_sum,transcripts_sum=transcripts_sum,mitocondria_sum=mitocondria_sum))%>% 
    mutate_if(is.numeric, round, digits = 2)
  per_sampleQC$sample<-x
  sampleQC<-rbind(sampleQC,per_sampleQC)
  sampleQC_df<-rbind(sampleQC_df,data.frame(df))
}
#table(sampleQC_df$sample)
sample_num<-length(sampleids)
#qc_stat_violinplot_file<-''
qcfile<-file.path(outputdir,paste(group,"sample.stats.sum","png",sep="."))
png(qcfile, width=sample_num*1+1.5, height=6, units="in", res=500)
violin_plot_qcstats(sampleQC_df)
dev.off()

#group<-"S1-S5"
qc1<-file.path(outputdir,paste(group,"sample.stats.sum","csv",sep="."))
sampleQC_df %>% 
  group_by(sample) %>% 
  summarise(cells=n(),
            total_transcripts = sum(sum),
            median_transcripts = round(median(sum),0),
            median_genes = round(median(detected),0))  %>% 
  write.table(file = qc1, sep = "\t", quote=FALSE,row.names = FALSE)

qc2<-file.path(outputdir,paste(group,"sample.qc.sum","csv",sep="."))
sampleQC %>% 
  write.table(file = qc2, sep = "\t", quote=FALSE,row.names = TRUE, col.names=NA)


#set filters: transcripts>1000 and <20000, genes>1000, mitocondria percentage<5%
qc.lib <- sampleQC_df$sum < 1000 | sampleQC_df$sum > 20000
qc.nexprs <- sampleQC_df$detected < 1e3 
qc.mito <- df$subsets_Mito_percent >5
discard <- qc.lib | qc.nexprs | qc.mito
summary(discard)

# Summarize the number of cells removed for each reason.
DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
          MitoProp=sum(qc.mito), Total=sum(discard))
sce$discard <- discard
#head(colData(sce),5)

# Keeping the columns we DON'T want to discard.
sce <- sce[,!discard]
#nrow(colData(sce))

outputrds<-file.path(outputdir,paste(group,"QC_filtered","rds",sep="."))
saveRDS(sce, outputrds)

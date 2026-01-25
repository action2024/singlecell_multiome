~/mambaforge/bin/mamba env create --file /cluster/home/liuji/T1D/multiome/src/multiome.yml -n multiome
~/mambaforge/bin/mamba env update --file /cluster/home/liuji/T1D/multiome/src/multiome.yml -n multiome

source ~/mambaforge/etc/profile.d/conda.sh
conda activate multiome

# Clustering
inputpath=/cluster/home/liuji/T1D/multiome/src/batch/$samplename/outs
outputdir=/cluster/home/liuji/T1D/multiome/analysis/$samplename/clustering
mkdir -p $outputdir
Rscript /cluster/home/liuji/T1D/multiome/src/R/01_CombinedClustering.R $inputpath $outputdir

# ATAC peak calling
archRpath=/cluster/home/liuji/T1D/multiome/analysis/$samplename/clustering
outputdir=/cluster/home/liuji/T1D/multiome/analysis/$samplename/peakcalling/Clusters_Combined
mkdir -p $outputdir
Rscript /cluster/home/liuji/T1D/multiome/src/R/03_PeakCalling.R $archRpath $outputdir Clusters_Combined
Rscript /cluster/home/liuji/T1D/multiome/src/R/03_PeakCalling.R $archRpath $outputdir Clusters_ATAC
Rscript /cluster/home/liuji/T1D/multiome/src/R/03_PeakCalling.R $archRpath $outputdir Clusters_RNA

#Marker Identification -RNA
archRpath=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/clustering
outputdir=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/MarkersRNA/Clusters_Combined
mkdir -p $outputdir
Rscript /cluster/home/liuji/T1D/multiome/src/R/02_MarkerDetect-RNA.R $archRpath $outputdir Clusters_Combined

archRpath=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/peakcalling/Clusters_Combined
outputdir=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/MarkersATAC/Clusters_Combined
mkdir -p $outputdir
Rscript /cluster/home/liuji/T1D/multiome/src/R/05_MarkerDetect-ATAC.R $archRpath $outputdir Clusters_Combined

# Differential accessbility/peak calling for individual cluster
archRpath=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/peakcalling/Clusters_Combined
outputdir=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/DAR/Clusters_Combined
mkdir -p $outputdir
Rscript /cluster/home/liuji/T1D/multiome/src/R/04_DAR.R $archRpath $outputdir Clusters_Combined

#positive TF regulator identification
archRpath=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/peakcalling/Clusters_Combined
outputdir=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/positiveTF/Clusters_Combined
mkdir -p $outputdir
Rscript /cluster/home/liuji/T1D/multiome/src/R/07_positiveTFs.R $archRpath $outputdir Clusters_Combined



###scRNA###
source ~/mambaforge/etc/profile.d/conda.sh
conda activate scRNAseq
inputrds=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/clustering/SEM01_S4D5/sce.rds
outputdir=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/scRNAseq/QC
mkdir -p $outputdir
group=Clusters_Combined
Rscript /cluster/home/liuji/T1D/multiome/src/R/multiome_scRNA/01_QC.R $inputrds $outputdir $group

inputrds=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/scRNAseq/Clusters_Combined.QC_filtered.rds
outputdir=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/scRNAseq/Clustering
mkdir -p $outputdir
Rscript /cluster/home/liuji/T1D/multiome/src/R/multiome_scRNA/02_Clustering.R $inputrds $outputdir $group
conda deactivate

source ~/mambaforge/etc/profile.d/conda.sh
conda activate pathfindr
group=Clusters_Combined
inputrds=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/scRNAseq/QC/Clusters_Combined.QC_filtered.rds
outputdir=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/scRNAseq/GeneExpression/$group
mkdir -p $outputdir
Rscript /cluster/home/liuji/T1D/multiome/src/R/multiome_scRNA/03_MarkerDetection.R $inputrds $outputdir $group
Rscript /cluster/home/liuji/T1D/multiome/src/R/multiome_scRNA/04_DEG.R $inputrds $outputdir $group

group=hclust
inputrds=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/scRNAseq/QC/Clusters_Combined.QC_filtered.rds
outputdir=/cluster/home/liuji/T1D/multiome/analysis/SEM01_S4D5/scRNAseq/GeneExpression/$group
mkdir -p $outputdir
Rscript /cluster/home/liuji/T1D/multiome/src/R/multiome_scRNA/03_MarkerDetection.R $inputrds $outputdir $group
Rscript /cluster/home/liuji/T1D/multiome/src/R/multiome_scRNA/04_DEG.R $inputrds $outputdir $group



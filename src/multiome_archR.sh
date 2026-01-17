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

# Marker Identification -RNA
archRpath=/cluster/home/liuji/T1D/multiome/analysis/$samplename/clustering
outputdir=/cluster/home/liuji/T1D/multiome/analysis/$samplename/peakcalling/Clusters_Combined
mkdir -p $outputdir
Rscript /cluster/home/liuji/T1D/multiome/src/R/02_MarkerDetecction-RNA.R $archRpath $outputdir Clusters_Combined

# Differential Accessbility/Peak calling for individual cluster
archRpath=/cluster/home/liuji/T1D/multiome/analysis/$samplename/peakcalling/Clusters_Combined
outputdir=/cluster/home/liuji/T1D/multiome/analysis/$samplename/peakcalling/Clusters_Combined
mkdir -p $outputdir
Rscript /cluster/home/liuji/T1D/multiome/src/R/04_DAR.R $archRpath $outputdir Clusters_Combined

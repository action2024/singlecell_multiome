# install cellranger-arc
wget -O cellranger-arc-2.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-arc/cellranger-arc-2.1.0.tar.gz?Expires=1764835084&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=bLNnGf1c~0J0TFLyB7SAVGrE1BMOf8iu4-r2hJ-33MsUdCffstIFmRHBsDE71ZF6lCr7kUpFt3VT8QRcrfBppddfGXI-CCaIvT2xYQADXi8qr1Yw9HEboMxQk~cQrIlHzNZxLeBx7pRJm3PRc38wE~lt30qh3OwbSAWzF8-fhSZSkcnl1h-TuTr1fOQi5Pz8ZuilrLjgygwspLqI4W~YmUxyS20gvkUW2Qu0n-Sog1TYrB1dv2-nqf8MCkVoa~3MeesBl15gNV9dRQrk-2bMgi9p1d5-V69Z8GGFFlKiTc167oc~eTjkQueG55fyXuF0i-BIzhhGxJJyEs8t1qBpGw__"
# unzip 
tar -xzvf cellranger-arc-2.1.0.tar.gz
# download reference 
wget "https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2024-A.tar.gz"

export PATH=$PATH:~/cellranger-arc-2.1.0/bin

path=~/Downloads
sample=SEM01_S4D5
samplesheet=$path/multiome/src/$sample.fastq.csv
### samplesheet format###
#fastqs,sample,library_type
#/cluster/home/liuji/T1D/multiome/data/fastq/SEM01_S4D5/GEX,SEM01_S4D5_GEX,Gene Expression
#/cluster/home/liuji/T1D/multiome/data/fastq/SEM01_S4D5/ATAC,SEM01_S4D5_ATAC,Chromatin Accessibility

refdir=/cluster/home/liuji/T1D/multiome/ref/refdata-cellranger-arc-GRCh38-2024-A
~/cellranger-arc-2.1.0/bin/cellranger-arc count --id=$sample \
                       --reference=$refdir \
                       --libraries=$samplesheet \
                       --create-bam=false \
                       --localcores=16 \
                       --localmem=64


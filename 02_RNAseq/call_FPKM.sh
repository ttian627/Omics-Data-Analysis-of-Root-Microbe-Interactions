# run on LRZ linux cluster
log in and change to the compute node
srun -p cm4_inter -t 03:00:00 -c 8 --pty bash

#install softwares
conda create -n RNAseq
conda install -c bioconda hisat2
conda install bioconda::stringtie
conda install -c bioconda sra-tools
conda install -c bioconda samtools
conda install fastqc


# download RNA-seq raw data (search from NCBI)
wget -c https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11839278/SRR11839278

# build genome index, download maize genome (Zm-B73-REFERENCE-NAM-5.0.fa) from MaizeGDB
wget -c https://download.maizegdb.org/Genomes/B73/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
gunzip Zm-B73-REFERENCE-NAM-5.0.fa.gz
samtools faidx Zm-B73-REFERENCE-NAM-5.0.fa chr10 > chr10.fa
#hisat2-build Zm-B73-REFERENCE-NAM-5.0.fa B73v5.hisat
hisat2-build chr10.fa B73v5.chr10.hisat

#quality control (skip)
fastqc SRR11839278_1.fastq

#mapping reads to reference genome
hisat2 --summary-file SRR11839278.summary -x B73v5.chr10.hisat -1 SRR11839278_1.fastq -2 SRR11839278_2.fastq | samtools view -bhS -q 30 - |  samtools sort -@ 10 -o SRR11839278.uniq.bam

# call FPKM
stringtie SRR11839278.uniq.bam -G ref/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gtf  -A stringtie_uniq/B73-Z2-2/SRR11839278.abun -o  stringtie_uniq/B73-Z2-2/SRR11839278.gtf -B -e

# merge FPKM of different samples (I have merged the 9 samples and you can download form github with the following link)
https://github.com/ttian627/Omics-Data-Analysis-of-Root-Microbe-Interactions/blob/main/02_RNAseq/B73-exp-matrix-FPKM.txt







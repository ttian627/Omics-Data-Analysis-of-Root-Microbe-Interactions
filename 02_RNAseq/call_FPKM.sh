#install softwares
install conda
conda install -c bioconda hisat2
conda install bioconda::stringtie

# download RNA-seq raw data (search from NCBI)
wget -c https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11839278/SRR11839278
# build genome index, download maize genome (Zm-B73-REFERENCE-NAM-5.0.fa) from MaizeGDB
hisat2-build Zm-B73-REFERENCE-NAM-5.0.fa B73v5.hisat

#mapping reads to reference genome
hisat2 --summary-file SRR11839278.summary -x ref/B73v5.hisat -1 SRR11839278_1.fastq -2 SRR11839278_2.fastq| samtools view -bhS -q 30 - |  samtools sort -@ 10 -o SRR11839278.uniq.bam

# call FPKM
stringtie SRR11839278.uniq.bam -G ref/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gtf  -A stringtie_uniq/SRR11839278/SRR11839278.abun -o  stringtie_uniq/SRR11839278/SRR11839278.gtf -B -e

# merge FPKM of different samples

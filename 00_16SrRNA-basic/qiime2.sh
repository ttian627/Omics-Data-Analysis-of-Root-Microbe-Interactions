
QIIME2 Amplicon Analysis Pipeline (16S / ITS)
1. Install Linux version
# conda env create -n qiime2-amplicon-2024.10 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml

Prepare files:
(1) 01_manifest_file.csv
sample-id   forward-absolute-filepath   reverse-absolute-filepath
wtB73-Z3-3  /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839631_1.fastq   /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839631_2.fastq
wtB73-Z3-2  /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839632_1.fastq   /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839632_2.fastq
wtB73-Z3-1  /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839633_1.fastq   /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839633_2.fastq
wtB73-Z2-3  /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839634_1.fastq   /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839634_2.fastq
wtB73-Z2-2  /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839635_1.fastq   /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839635_2.fastq
wtB73-Z2-1  /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839636_1.fastq   /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839636_2.fastq
wtB73-Z1-3  /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839637_1.fastq   /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839637_2.fastq
wtB73-Z1-2  /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839639_1.fastq   /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839639_2.fastq
wtB73-Z1-1  /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839640_1.fastq   /dss/dsshome1/0D/go82tip2/1microbome/microb-test/16S/SRR11839640_2.fastq
(2) 02_metadata.tsv
#SampleID   Group   Replicate
wtB73-Z3-3  GroupA  1
wtB73-Z3-2  GroupA  2
wtB73-Z3-1  GroupA  3
wtB73-Z2-3  GroupB  1
wtB73-Z2-2  GroupB  2
wtB73-Z2-1  GroupB  3
wtB73-Z1-3  GroupC  1
wtB73-Z1-2  GroupC  2
wtB73-Z1-1  GroupC  3
2. Data import
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path 01_manifest_file.csv \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
3. Inspect imported data
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux-summary.qzv
4. 16S and ITS region extraction #please check the prime sequences used in this study
16S extraction
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-front-f CCTAYGGGRBGCASCAG \
  --p-front-r GGACTACNNGGGTATCTAAT \
  --o-trimmed-sequences 16S-demux.qza
qiime demux summarize \
  --i-data 16S-demux.qza \
  --o-visualization 16S-demux-summary.qzv
ITS extraction
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-front-f CTTGGTCATTTAGAGGAAGTAA \
  --p-front-r GCTGCGTTCTTCATCGATGC \
  --o-trimmed-sequences IST1-demux.qza
qiime demux summarize \
  --i-data ITS1-demux.qza \
  --o-visualization ITS1-demux-summary.qzv
5. Denoising (DADA2)
(Trimming and truncation parameters should be determined based on the quality visualization results in 16S-demux-summary.qzv)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 16S-demux.qza \
  --p-trim-left-f 20 \          # Trim first 20 bp (usually primer region)
  --p-trim-left-r 20 \
  --p-trunc-len-f 250 \         # Truncate forward reads at 250 bp (based on quality profile)
  --p-trunc-len-r 200 \         # Truncate reverse reads
  --p-max-ee 2 \                # Maximum expected errors
  --p-n-threads 4 \             # Number of threads
  --o-representative-sequences rep-seqs.qza \  # Representative ASV sequences (for taxonomy/phylogeny)
  --o-table table.qza \         # ASV abundance table (feature table)
  --o-denoising-stats stats.qza # DADA2 denoising statistics (quality evaluation)
View denoising results
qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv
6. Taxonomic annotation
Download SILVA v138 classifier
(This downloads the SILVA v138 99% similarity OTU reference database for ASV classification)
wget https://data.qiime2.org/2024.2/common/silva-138-99-nb-classifier.qza
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
View taxonomy results
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
7. Remove non-target sequences
Remove mitochondrial and chloroplast sequences
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table filtered-table.qza
ASV filtering
(Keep features present in ≥2 samples and with total reads >10)
qiime feature-table filter-features \
  --i-table filtered-table.qza \
  --p-min-samples 2 \
  --p-min-frequency 10 \
  --o-filtered-table final-table.qza   # Final output table
View filtered table
qiime metadata tabulate \
  --m-input-file final-table.qza \
  --o-visualization final-table.qzv
8. Merge abundance by taxonomic level
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table collapsed-table-genus.qza

9. Normalize microbial abundance using DESeq2 （skipped if do differential analysis）
qiime feature-table normalize \
  --i-table collapsed-table-genus.qza \
  --p-method deseq2 \
  --o-normalized-table collapsed-table-genus-norm.qza

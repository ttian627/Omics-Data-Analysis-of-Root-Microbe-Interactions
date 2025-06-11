library(WGCNA)
library(DESeq2)
library(tidyverse)

#===============================================================================
#
#  Read the gene counts table and plot the sample tree
#
#===============================================================================
#setwd("~/")

# Read the gene counts table
data0=read.table("gene_counts_table_Z1.txt",header=T,row.names=1,sep="\t",check.names = FALSE)
# Normalization with log2(FPKM+1)
sample_metadata = read.csv(file = "sample_info.csv")
#sample_metadata = read.csv(file = "sample_group_info.1.csv")
z1_samples <- sample_metadata[sample_metadata$Zone == "Z1", ]
dataExpr_deseq <- DESeqDataSetFromMatrix(countData = data0[,-61],colData = z1_samples,design = ~ Genotype)
mcols(dataExpr_deseq)$basepairs = data0$geneLength
fpkm_matrix = fpkm(dataExpr_deseq)
datExpr = t(log2(fpkm_matrix+1))


# Calculate sample distance and cluster the samples
sampleTree = hclust(dist(datExpr), method = "average");
# plot sample tree
pdf(file = "0-n-sampleClustering.pdf", width = 10, height = 5);
par(cex = 0.95);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()


#remove outliners

outliners_to_remove <- c("LH205_Z1_4", "LH51_Z1_4")
fpkm_filtered <- fpkm_matrix[, !colnames(fpkm_matrix) %in% outliners_to_remove]


fpkm_mean <- fpkm_filtered %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "fpkm") %>%
  mutate(genotype = sub("_.*", "", sample)) %>%  
  group_by(gene, genotype) %>%
  summarise(mean_fpkm = mean(fpkm), .groups = "drop") %>%  
  pivot_wider(names_from = genotype, values_from = mean_fpkm) %>%
  column_to_rownames("gene") 

datExpr = t(log2(fpkm_mean+1))
# remove genes with the expression of 0 in all samples
zero_genes <- which(colSums(datExpr) == 0)
if (length(zero_genes) > 0) {
  datExpr <- datExpr[, -zero_genes]
}


library(pheatmap)
pdf(file = "1-n-geneCoexpression-heatmap.pdf", width = 10, height = 10);

cor_matrix <- cor(datExpr[,1:50], method = "pearson")
pheatmap(cor_matrix,
         #         clustering_distance_rows = "euclidean",
         #         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Gene co-expression heatmap")

dev.off()

#===============================================================================
#
#  Choose soft threshold parameter
#
#===============================================================================

# Choose a set of soft threshold parameters
powers = c(c(1:20), seq(from = 22, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "2-n-sft.pdf", width = 15, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
	xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=10,col="red") 
dev.off()

#===============================================================================
#
#  Turn data expression into topological overlap matrix
#
#===============================================================================

# Turn data expression into topological overlap matrix
power=sft$powerEstimate
TOM = TOMsimilarityFromExpr(datExpr, power = power)
dissTOM = 1-TOM
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
#pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
#plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
#     labels = FALSE, hang = 0.04);
#dev.off()

#===============================================================================
#
#  Construct modules
#
#===============================================================================

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 30);
table(dynamicMods)
length(table(dynamicMods))
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#pdf(file = "4-module_tree.pdf", width = 10, height = 6);
#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
#                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
#dev.off()

#===============================================================================
#
#  Merge modules
#
#===============================================================================

# Merge close modules
MEDissThres=0.25
#abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  

# display module information
table(mergedColors)
gene_module_df <- data.frame(gene = colnames(datExpr),  module = mergedColors)
blue_genes <- gene_module_df[gene_module_df$module == "blue", "gene"] 


mergedMEs = merge$newMEs  
# Plot merged module tree
pdf(file = "5-merged_Module_Tree.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()
#write.table(merge$oldMEs,file="6-oldMEs.txt");
#write.table(merge$newMEs,file="6-newMEs.txt");
write.table(gene_module_df, file = "6-all_genes_modules.txt", sep = "\t", quote = FALSE)
write.table(blue_genes, file = "6-blue_genes_modules.txt", sep = "\t", quote = FALSE)

#=====================================================================================
#
#  Correlation between gene modules and microbial traits
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)




# Read microbial data as traits
bac_traits_order = read.table("traits_file/f_order_88.txt", header = T, sep = "\t")
bac_traits_filtered <- bac_traits_order %>%
  select(Taxonomy, where(~ is.numeric(.) && max(., na.rm = TRUE) >= 0.1))
bac_traits_Z1 <- bac_traits_filtered[grep("Z1", bac_traits_filtered$Taxonomy), ]

# mean fromn3 repeats
bac_traits_Z1$genotype <- sub("_.*", "", bac_traits_Z1$Taxonomy)
bac_traits_Z1_means <- bac_traits_Z1 %>%
  as.data.frame() %>%
  select(-Taxonomy) %>%  
  group_by(genotype) %>%
  summarise(across(everything(), ~mean(., na.rm = TRUE))) %>%
  as.data.frame()
bac_traits<-bac_traits_Z1_means


rownames(bac_traits) = bac_traits[, 1]
bac_traits = bac_traits[, -1]
#rownames(MEs) = paste(substr(rownames(MEs), 1, nchar(rownames(MEs))-1), rep(c("1", "2", "3"), 60), sep = "")
# sample names should be consistent in eigen genes and traits !!!!
bac_traits = bac_traits[match(rownames(MEs), rownames(bac_traits)), ]
table(rownames(MEs) == rownames(bac_traits))
# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, bac_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# rm NA
write.table(moduleTraitCor,file="7-moduleTrait_correlation.txt");
write.table(moduleTraitPvalue,file="7-moduleTrait_pValue.txt");


#=====================================================================================
#
#  Plot heatmap of module-traits relationship
#
#=====================================================================================

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("7-module-traits-bacteria-order.pdf", width = 20, height = 10)
par(mar = c(10, 10, 2, 2));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(bac_traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = NULL,
               setStdMargins = FALSE,
               cex.text = 1.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()



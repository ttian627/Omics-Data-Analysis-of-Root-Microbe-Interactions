#run in Rstudio
#PCA analysis(R script)
#setwd("change to your working directory")
setwd("/Users/ttian627163.com/Documents/TUM/03_Course/Omics-Data-Analysis-of-Root-Microbe-Interactions-main/02_RNAseq")
library(ggplot2)
library(ggrepel)
library(plotly)
#read expression matrix
expr_matrix_B73 <- read.table("B73-exp-matrix-FPKM.txt", row.names=1,header = T,check.names = FALSE)
#expr_matrix_B73 <- expr_matrix[,91:99]
#expr_matrix_other <- expr_matrix[,1:9]
expr_matrix_B73_t <- t(expr_matrix_B73)
expr_matrix_B73_filtered <- expr_matrix_B73_t[, apply(expr_matrix_B73_t, 2, function(x) sd(x) > 10)]
pca_result <- prcomp(expr_matrix_B73_filtered, scale. = TRUE)

#extact PCA result       
pca_data <- data.frame(Sample = rownames(pca_result$x),
                       PC1 = pca_result$x[,1],
                       PC2 = pca_result$x[,2],
                       PC3 = pca_result$x[,3])

#add group information
group_info <- read.csv("sample_group_info.csv")
pca_data <- merge(pca_data, group_info, by="Sample")

#2D PCA
ggplot(pca_data, aes(x=PC1, y=PC2, color=Group)) +
     geom_point(size=3) +
     geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 100) +
     scale_color_manual(values = c("Z1" = "#f5c242", "Z2" = "#4fad5b", "Z3" = "#4fadea")) + 
     theme_minimal() +
     labs(title="PCA of Transcriptome Data",
          x=paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
          y=paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"))

# 3D PCA
plot_ly(pca_data, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Group,
        colors = c("Z1" = "#f5c242", "Z2" = "#4fad5b", "Z3" = "#4fadea"),
        type = "scatter3d", mode = "markers+text", text = ~Sample,textposition = "top center", marker = list(size = 5)) %>%
  plotly::layout(title = "3D PCA of Transcriptome Data",
         scene = list(xaxis = list(title = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100,1), "%)")),
                      yaxis = list(title = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100,1), "%)")),
                      zaxis = list(title = paste0("PC3 (", round(summary(pca_result)$importance[2,3]*100,1), "%)"))))
                      
#heatmap analysis(R script)
library(pheatmap)
#cor_matrix <- cor(expr_matrix_B73, method = "pearson")
cor_matrix <- cor(expr_matrix_B73, method = "spearman")
pheatmap(cor_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Sample Correlation Heatmap")


                            

# DEGs
library(DESeq2)
#all( colnames(expr_matrix)[-181] == group_info$Sample)   # 应为 TRUE
#colData$Group <- factor(colData$Group)

expr_matrix=read.table("../03_WGCNA/gene_counts_table_WGCNA.txt",header=T,row.names=1,sep="\t")
sample_metadata = read.csv(file = "sample_group_info.csv")

sub_metadata <- sample_metadata[sample_metadata$Genotype %in% c("B73"), ]
sub_expr_matrix <- expr_matrix[, sub_metadata$Sample]
dataExpr_deseq <- DESeqDataSetFromMatrix(countData = sub_expr_matrix,colData = sub_metadata,design = ~ Group)


dataExpr_deseq <- DESeq(dataExpr_deseq)
res <- results(dataExpr_deseq, contrast = c("Group", "Z2", "Z1"))
res <- na.omit(res)
sig_res <- res[(res$padj < 0.05) & (abs(res$log2FoldChange) > 1), ]


# Volcano Plot
df <- as.data.frame(res)

# mark up and down regulated genes
df$significance <- "Not Sig"
df$significance[df$padj < 0.05 & df$log2FoldChange > 1] <- "Up"
df$significance[df$padj < 0.05 & df$log2FoldChange < -1] <- "Down"

ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Volcano plot",
       x = "log2 Fold Change",
       y = "-log10 Adjusted p-value") +
  theme_minimal() +
  theme(legend.title = element_blank())


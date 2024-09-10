setwd("C:/Users/tobij/Desktop/HB Essays, Papers and Tasks/R Tasks/R Course Tasks") # Setting my current working directory
dres <- read.csv('dres.csv') # Loading the DESeq results
str(dres) # Inspecting the sturcture of the dataframe
sum(is.na(dres)) # Checking for null values

library(ggplot2) # Loading ggplot2

pdf('dres_volcano.pdf', width = 10, height = 10) # Opening a pdf device for a volcano plot

# Setting thresholds and creating a column to show whether a gene is significant or not
abs_lfc_thresh <- 1
significance_threshold_value <- 0.01
dres$thresh <- with(dres, abs(log2FoldChange) > abs_lfc_thresh & pvalue < significance_threshold_value)

# Plotting the volcano plot
ggplot(dres, aes(x = log2FoldChange, y = -log10(pvalue))) + geom_point(aes(colour = thresh), size=1.5) + theme_minimal() + xlab('Log2 Fold Change') + ylab('-Log10 P-Value') + ggtitle('Volcano Plot of -log10pvalue vs log2FoldChange') + scale_colour_manual(values = c('orange', 'navy'), name = 'Significant?')

dev.off() # Closing the device

# Selecting upregulated and downregulated genes based on the pvalue and fold change thresholds
upregulated <- subset(dres, log2FoldChange > abs_lfc_thresh & pvalue < significance_threshold_value) 
downregulated <- subset(dres, log2FoldChange < -abs_lfc_thresh & pvalue < significance_threshold_value)

# Ranking the genes based on the pvalue and foldchange metrics. It is the best to show the magnitude of expression to select the most significantly over or underexpressed genes
upregulated$rank_score <- upregulated$log2FoldChange * -log10(upregulated$pvalue) 
downregulated$rank_score <- downregulated$log2FoldChange * -log10(downregulated$pvalue)

# Sorting the genes based on the rank score
upregulated <- upregulated[order(upregulated$rank_score, decreasing = T), ] # Sorted in decreasing order
downregulated <- downregulated[order(downregulated$rank_score), ] # Sorted in increasing order

# Picking the top genes
top_5_upregulated_genes <- head(upregulated, n=5)$gene
top_5_downregulated_genes <- head(downregulated, n=5)$gene


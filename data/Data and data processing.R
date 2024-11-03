# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install missing Bioconductor packages
BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("tximport", force = TRUE)
BiocManager::install("pheatmap", force = TRUE)
BiocManager::install("tximportData", force = TRUE)

# Load necessary libraries
library(DESeq2)
library(tximport)
library(pheatmap)
library(tidyverse)


# Check if tximport has been loaded correctly
if ("tximport" %in% loadedNamespaces()) {
  print("tximport loaded successfully.")
} else {
  stop("Failed to load tximport. Please reinstall the package.")
}


# Load metadata from a CSV or tab-delimited file
metadata <- read.csv("SRP518774_metadata.txt", sep = ",", header = TRUE)

# Print column names to verify structure
print(colnames(metadata))

# Ensure the 'Run' column matches the folder names
print(metadata$Run)

quant_files <- list.files(path = "/Users/putriramadani/Documents/Drug repurposing for AML/data/quants", pattern = "quant.sf", full.names = TRUE, recursive = TRUE)

# Check if quant files are detected
print(quant_files)


sample_names <- basename(dirname(quant_files))
names(quant_files) <- sample_names

# Print sample names to check if they match with 'Run'
print(sample_names)

# Filter metadata to include only samples present in the quant files
metadata_filtered <- metadata[metadata$Run %in% sample_names, ]
print(nrow(metadata_filtered))  # Check if this matches the number of quant files

# Add a 'Condition' column based on 'Sample Name' or other logic
metadata_filtered$Condition <- ifelse(grepl("HD", metadata_filtered$Sample.Name, ignore.case = TRUE), "Healthy", "AML")

# Print the filtered metadata to verify
print(metadata_filtered)

# Prepare coldata for DESeq2
coldata <- data.frame(
  row.names = metadata_filtered$Run,
  condition = metadata_filtered$Condition
)

# Print coldata to verify
print(coldata)

# Check if coldata matches quant_files
if (!all(rownames(coldata) %in% names(quant_files))) {
  stop("Mismatch between coldata row names and quant file names")
}

# Import quant.sf files using tximport
txi <- tximport(quant_files, type = "salmon", txOut = TRUE)




# Print filtered metadata to ensure it has rows
print("Filtered metadata:")
print(metadata_filtered)

# Ensure rownames are set properly
coldata <- data.frame(
  condition = metadata_filtered$Condition,
  row.names = metadata_filtered$Run
)

# Print coldata to verify it has row names
print("Coldata with row names:")
print(rownames(coldata))

# Check if coldata matches the quant_files
if (!all(rownames(coldata) %in% names(quant_files))) {
  stop("Mismatch between coldata row names and quant file names")
} else {
  print("All row names in coldata match the quant file names.")
}

















# Load necessary libraries
library(DESeq2)
library(pheatmap)
library(ggplot2)

# Proceed with DESeq2 pipeline
# Create DESeq2 dataset from imported data and coldata
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~condition)

# Run DESeq2 normalization and differential expression analysis
dds <- DESeq(dds)

# Variance stabilizing transformation for visualization
vsd <- vst(dds, blind = FALSE)

# PCA plot to visualize data distribution
pca_plot <- plotPCA(vsd, intgroup = "condition")
print(pca_plot)

# Get differential expression results
res <- results(dds)

# Order results by adjusted p-value
resOrdered <- res[order(res$padj),]

# Print summary of results
summary(res)

# Visualize significant DEGs with a volcano plot
volcano_data <- as.data.frame(res)
volcano_data$significant <- volcano_data$padj < 0.05 & abs(volcano_data$log2FoldChange) > 1
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")

# Generate heatmap of top differentially expressed genes
top_genes <- head(order(res$padj, na.last = NA), 50)
pheatmap(assay(vsd)[top_genes, ], cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, annotation_col = coldata)


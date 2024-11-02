Exploring transcript data profile of healthy and AML patient samples
================

## Loading the required libraries

``` r
library(tximport)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
```

## Loading the data

Prepare the sample table of the experiment and load quant data generated
by Salmon into R

``` r
sample_table <- read_csv("../data/SRP518774_metadata.txt") %>%
  filter(`Assay Type` == "RNA-Seq") %>%
  select(`Run`, `isolate`) %>%
  mutate(condition = if_else(str_detect(isolate, "health"), "Control/Healthy", "AML Patient")) %>%
  select(`Run`, condition) %>%
  mutate(sample_name = Run) %>%
  select(sample_name, condition)
```

    ## Rows: 17 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (28): Run, AGE, Assay Type, BIOMATERIAL_PROVIDER, BioProject, BioSample...
    ## dbl   (4): AvgSpotLen, Bases, Bytes, version
    ## dttm  (2): ReleaseDate, create_date
    ## date  (1): Collection_Date
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
sample_files <- paste0(pull(sample_table, 
            sample_name), '_quant/quant.sf')

sample_files <- paste0("../data/quants/", sample_files)

paste0(sample_files, '_quant/quant.sf')
```

    ##  [1] "../data/quants/SRR29752624_quant/quant.sf_quant/quant.sf"
    ##  [2] "../data/quants/SRR29752625_quant/quant.sf_quant/quant.sf"
    ##  [3] "../data/quants/SRR29752626_quant/quant.sf_quant/quant.sf"
    ##  [4] "../data/quants/SRR29752627_quant/quant.sf_quant/quant.sf"
    ##  [5] "../data/quants/SRR29752628_quant/quant.sf_quant/quant.sf"
    ##  [6] "../data/quants/SRR29752629_quant/quant.sf_quant/quant.sf"
    ##  [7] "../data/quants/SRR29752630_quant/quant.sf_quant/quant.sf"
    ##  [8] "../data/quants/SRR29752631_quant/quant.sf_quant/quant.sf"
    ##  [9] "../data/quants/SRR29752632_quant/quant.sf_quant/quant.sf"
    ## [10] "../data/quants/SRR29752633_quant/quant.sf_quant/quant.sf"
    ## [11] "../data/quants/SRR29752634_quant/quant.sf_quant/quant.sf"
    ## [12] "../data/quants/SRR29752635_quant/quant.sf_quant/quant.sf"
    ## [13] "../data/quants/SRR29752636_quant/quant.sf_quant/quant.sf"
    ## [14] "../data/quants/SRR29752639_quant/quant.sf_quant/quant.sf"
    ## [15] "../data/quants/SRR29752640_quant/quant.sf_quant/quant.sf"

``` r
gene_map <- read_csv('../data/gene_map.csv', 
                     col_names = c('enstid', 'ensgid'))
```

    ## Rows: 495632 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): enstid, ensgid
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
names(sample_files) <- pull(sample_table, sample_name)


count_data <- tximport(files = sample_files,
                       type = 'salmon', 
                       tx2gene = gene_map,
                       ignoreTxVersion = TRUE)
```

    ## reading in files with read_tsv
    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 
    ## removing duplicated transcript rows from tx2gene
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

## Load the quant data and sample information into DESeq2 datatype

``` r
sample_table$condition <- factor((c('AML', 'AML', 'AML', 'AML', 'Healthy',
                                         'Healthy', 'AML', 'AML', 'Healthy', 'Healthy',
                                         'Healthy', 'Healthy', 'Healthy', 'AML', 'AML')),
                                      levels = c('AML', 'Healthy'))

dds_aml <- DESeqDataSetFromTximport(txi = count_data,
                                     colData = sample_table,
                                     design = ~condition)
```

## Principal Component Analysis

``` r
# Transform the data so it is suitable for PCA
vst_aml <- varianceStabilizingTransformation(dds_aml)
# create matrix
vst_mat <- assay(vst_aml)

pca <- prcomp(t(vst_mat))

# plot out the data frame to produce PCA
df <- as.data.frame(pca$x)
df$condition <- sample_table$condition

pve <- round(pca$sdev^2/sum(pca$sdev^2) * 100, 2)

rownames_to_column(df, var = "sample_name") %>% as.tibble() %>%
  ggplot(., aes(x=PC1, y=PC2, color = condition)) +
  geom_point(size = 4) +
  #geom_text(aes(label = sample_name, color = condition), vjust = -1, size = 4) +
  xlab(label = paste0("PC1 (", pve[1], "%)")) +
  ylab(label = paste0("PC2 (", pve[2], "%)")) +
  theme_classic() +
  ggtitle("Principal Component Analysis")
```

![](eda_files/figure-gfm/pca-1.png)<!-- -->

``` r
ggsave("../figures/PCA.png")
```

Looks like we have two dots that consist of two overlap healthy samples,
so that we only see five dots for healthy

## Heatmap

``` r
distance <- dist(t(assay(vst_aml)))
distance_matrix <- as.matrix(distance)
rownames(distance_matrix) <- vst_aml$condition
colnames(distance_matrix) <- vst_aml$condition
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)

pheatmap(distance_matrix,
         clustering_distance_rows=distance,
         clustering_distance_cols=distance,
         col=colors)
```

![](eda_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggsave("../figures/heatmap.png")
```

comnfirmed with the heatmap, two dots consist of two samples that very
very similiar. Is that come from same person but different processed
samples (technical replicates)?

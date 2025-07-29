# Lib ----
library(readr)
library(readxl)
library(tidyverse)
library(tidyr)
library(dplyr)

library(DESeq2)
library(limma)
library(edgeR)

# data analysis ----
merged_gene_counts <- read_delim("data/salmon.merged.gene_counts.tsv", delim = "\t")
d_se_filter_corrected <- readRDS("data/d_se_filter_corrected.rds")
test <- as.data.frame(colData(d_se_filter_corrected))

counts <-  d_se_filter_corrected@assays@data@listData[["counts"]]

gene_name <- merged_gene_counts %>%
  select(gene_id, gene_name)



counts <- as.data.frame(counts)
counts <- tibble::rownames_to_column(counts, var = "gene_id")

counts <- counts %>%
  left_join(gene_name, by = "gene_id")
counts <- counts %>%
  relocate(gene_name, .after = gene_id)

t1_cols <- counts %>% select(ends_with("T1"))
t4_cols <- counts %>% select(ends_with("T4"))

ncol(t1_cols)
ncol(t4_cols)

# Stats ----
# basic stata
counts2 <- counts %>%
  mutate(
    mean_T1 = rowMeans(select(., ends_with("T1"))),
    mean_T4 = rowMeans(select(., ends_with("T4"))),
    log2FC = log2(mean_T4 + 1) - log2(mean_T1 + 1) 
  )

# DESeq2
count_matrix <- counts %>%
  select(ends_with("T1") | ends_with("T4")) %>%
  round() %>%
  as.matrix()
rownames(count_matrix) <- counts$gene_name


sample_ids <- colnames(count_matrix)
condition <- ifelse(grepl("T1$", sample_ids), "T1", "T4")

col_data <- data.frame(
  row.names = sample_ids,
  condition = factor(condition, levels = c("T1", "T4"))
)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "T4", "T1"))

res_df <- as.data.frame(res)
res_df$gene_name <- rownames(res_df)

counts_deseq2 <- counts %>%
  left_join(res_df, by = "gene_name") %>%
  arrange(padj)


# limma
group <- factor(condition)

dge <- DGEList(counts = count_matrix, group = group)
dge <- calcNormFactors(dge)  

design <- model.matrix(~ group)
v <- voom(dge, design, plot = TRUE) 

fit <- lmFit(v, design)
fit <- eBayes(fit)

res_limma <- topTable(fit, coef = 2, number = Inf)
res_limma$gene_name <- rownames(res_limma)

counts_limma <- counts %>%
  left_join(res_limma, by = "gene_name") %>%
  arrange(adj.P.Val)

#filtering for significance ----
# For DESeq2
sig_deseq2 <- counts_deseq2 %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

# For limma
sig_limma <- counts_limma %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1)


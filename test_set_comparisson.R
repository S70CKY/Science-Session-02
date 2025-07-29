library(readr)
library(dplyr)
library(tibble)
library(DESeq2)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)

counts_raw <- read_delim("data/salmon.merged.gene_counts.tsv", delim = "\t")

gene_names <- counts_raw %>% select(gene_id, gene_name)
count_matrix <- counts_raw %>%
  select(-gene_name) %>%
  column_to_rownames("gene_id") %>%
  round() %>%
  as.matrix()

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
res_deseq2 <- results(dds, contrast = c("condition", "T4", "T1"))

res_deseq2_df <- as.data.frame(res_deseq2) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_names, by = "gene_id") %>%
  arrange(padj)

volcano_data <- res_deseq2_df %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up (T4)",
      padj < 0.05 & log2FoldChange < -1 ~ "Down (T4)",
      TRUE ~ "Not significant"
    )
  )

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up (T4)" = "firebrick", "Down (T4)" = "dodgerblue", "Not significant" = "grey70")) +
  labs(
    title = "Volcano Plot: T4 vs T1",
    x = "log2 Fold Change (T4 / T1)",
    y = "-log10 Adjusted p-value",
    color = "Direction"
  ) +
  theme_minimal()

res_df_clean <- as.data.frame(res_deseq2)
res_df_clean$gene_id <- rownames(res_df_clean)

res_df_clean <- res_df_clean[!is.na(res_df_clean$padj), ]

res_df_clean <- res_df_clean[order(res_df_clean$padj), ]

top_genes <- res_df_clean$gene_id[1:30]


norm_counts <- counts(dds, normalized = TRUE)
norm_top <- norm_counts[top_genes, ]

rownames(norm_top) <- res_deseq2_df %>%
  filter(gene_id %in% top_genes) %>%
  arrange(match(gene_id, top_genes)) %>%
  pull(gene_name)

pheatmap(norm_top,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Top 30 DE Genes: T4 vs T1")



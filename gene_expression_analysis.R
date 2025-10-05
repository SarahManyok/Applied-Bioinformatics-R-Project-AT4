# gene_expression_analysis.R
# Analysis of RNA-seq count data from gene_expression.tsv

library(dplyr)
library(ggplot2)

# ---- Load data ----
gene_expr <- read.delim("gene_expression.tsv", row.names = 1, stringsAsFactors = FALSE)

# 1. Show first six genes
head(gene_expr, 6)

# 2. Add column: mean of other columns
gene_expr$mean_expr <- rowMeans(gene_expr, na.rm = TRUE)
head(gene_expr, 6)

# 3. List 10 genes with highest mean expression
top10_genes <- head(gene_expr[order(-gene_expr$mean_expr), ], 10)
print(top10_genes)

# 4. Number of genes with mean < 10
low_expr_count <- sum(gene_expr$mean_expr < 10)
cat("Number of genes with mean expression < 10:", low_expr_count, "\n")

# 5. Histogram of mean values
p <- ggplot(gene_expr, aes(x = mean_expr)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Mean Gene Expression",
       x = "Mean Expression",
       y = "Count")
print(p)

ggsave("gene_expression_histogram.png", p, width = 8, height = 5)

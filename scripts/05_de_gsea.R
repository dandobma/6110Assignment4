# =============================================================================
# Assignment 4 - scRNA-seq Analysis
# Script 05: Differential Expression & GSEA
# =============================================================================

# --- Install packages (run once, then comment out) ---------------------------
# BiocManager::install("DESeq2")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db")
# install.packages("ggrepel")

library(Seurat)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggrepel)
library(enrichplot)

# Load checkpoint from script 04
seu <- readRDS("data/seurat_checkpoint_04.rds")

# =============================================================================
# PART 1: PSEUDOBULK DIFFERENTIAL EXPRESSION
# =============================================================================

# --- 1. Subset to cell type and conditions of interest -----------------------
# Focus: Macrophages, RM tissue, Naive vs D05
# Exclude unassigned mouse_id cells — DESeq2 needs replicate labels

mac_rm <- subset(seu, subset = cell_type    == "Macrophages" &
                               organ_custom == "RM"          &
                               time         %in% c("Naive", "D05") &
                               replicate    != "unassigned")

cat("Cells in analysis:\n")
table(mac_rm$time, mac_rm$mouse_id)

# --- 2. Create pseudobulk counts ----------------------------------------------
# Sum raw counts across all cells per sample (mouse_id × timepoint)
# This converts single-cell data into bulk-like counts for DESeq2

# Extract raw count matrix
counts_mat <- GetAssayData(mac_rm, assay = "RNA", layer = "counts")

# Create a grouping variable: replicate + timepoint
mac_rm$pseudobulk_id <- paste(mac_rm$mouse_id, mac_rm$time, sep = "_")

# Aggregate counts per pseudobulk sample
pseudobulk_counts <- sapply(unique(mac_rm$pseudobulk_id), function(id) {
  cells <- WhichCells(mac_rm, expression = pseudobulk_id == id)
  rowSums(counts_mat[, cells, drop = FALSE])
})

cat("\nPseudobulk sample sizes (cells per sample):\n")
print(table(mac_rm$pseudobulk_id))

# --- 3. Create sample metadata for DESeq2 ------------------------------------
sample_meta <- data.frame(
  pseudobulk_id = colnames(pseudobulk_counts)
) %>%
  mutate(
    time      = str_extract(pseudobulk_id, "Naive|D05"),
    mouse_id  = str_remove(pseudobulk_id, "_Naive|_D05"),
    condition = factor(time, levels = c("Naive", "D05"))  # Naive = reference
  )
rownames(sample_meta) <- sample_meta$pseudobulk_id

cat("\nSample metadata:\n")
print(sample_meta)

# --- 4. Run DESeq2 ------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts,
  colData   = sample_meta,
  design    = ~ condition
)

# Filter lowly expressed genes: keep genes with ≥10 counts in ≥2 samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
cat("\nGenes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds, test = "Wald")

# --- 5. Extract results -------------------------------------------------------
# D05 vs Naive: positive log2FC = upregulated at D05
res <- results(dds,
               contrast  = c("condition", "D05", "Naive"),
               alpha     = 0.05)

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

cat("\nDifferential expression summary:\n")
summary(res)

cat("\nTop 20 DEGs:\n")
print(head(res_df, 20))

write.csv(res_df, "results/de_macrophages_RM_NaiveVsD05.csv", row.names = FALSE)

# --- 6. Volcano plot ----------------------------------------------------------
res_df <- res_df %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up in D05",
      padj < 0.05 & log2FoldChange < -1 ~ "Down in D05",
      TRUE ~ "Not significant"
    ),
    label = ifelse(padj < 0.05 & abs(log2FoldChange) > 2, gene, "")
  )

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj),
                                 color = significance, label = label)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Up in D05"       = "#E34234",
                                "Down in D05"     = "#4169E1",
                                "Not significant" = "grey60")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title    = "Macrophages (RM): D05 vs Naive",
       x        = "log2 fold change",
       y        = "-log10(adjusted p-value)",
       color    = "") +
  theme_classic() +
  theme(legend.position = "top")

ggsave("figures/05_volcano_macrophages_RM.png",
       p_volcano, width = 8, height = 7, dpi = 150, bg = "white")

# =============================================================================
# PART 2: GSEA
# =============================================================================

# --- 7. Prepare ranked gene list for GSEA ------------------------------------
# Rank all genes by their log2FC × -log10(pvalue) score
# This preserves direction and significance simultaneously

gsea_input <- res_df %>%
  filter(!is.na(log2FoldChange), !is.na(pvalue)) %>%
  mutate(rank_score = log2FoldChange * -log10(pvalue + 1e-300)) %>%
  arrange(desc(rank_score))

# Convert gene symbols to Entrez IDs (required by clusterProfiler)
gene_ids <- bitr(gsea_input$gene,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Mm.eg.db)

gsea_input <- gsea_input %>%
  left_join(gene_ids, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

# Named numeric vector of rank scores
gene_list        <- gsea_input$rank_score
names(gene_list) <- gsea_input$ENTREZID

# --- 8. Run GSEA on GO Biological Process terms ------------------------------
gsea_go <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",           # Biological Process
  keyType      = "ENTREZID",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

cat("\nTop GSEA GO terms:\n")
print(head(as.data.frame(gsea_go)[, c("Description", "NES", "pvalue", "p.adjust")], 20))

write.csv(as.data.frame(gsea_go),
          "results/gsea_macrophages_RM_NaiveVsD05.csv", row.names = FALSE)

# --- 9. GSEA dot plot ---------------------------------------------------------
p_gsea <- dotplot(gsea_go, showCategory = 20, split = ".sign") +
          facet_grid(. ~ .sign) +
          ggtitle("GSEA GO:BP — Macrophages (RM): D05 vs Naive") +
          theme(axis.text.y = element_text(size = 8))

ggsave("figures/05_gsea_dotplot_macrophages_RM.png",
       p_gsea, width = 14, height = 8, dpi = 150, bg = "white")

# --- 10. GSEA enrichment plot for top terms -----------------------------------
# Plot enrichment curve for the single most significant term
top_term <- gsea_go@result$ID[1]

p_enrich <- gseaplot2(gsea_go,
                       geneSetID = top_term,
                       title     = gsea_go@result$Description[1])

ggsave("figures/05_gsea_enrichment_plot.png",
       p_enrich, width = 8, height = 6, dpi = 150, bg = "white")

# --- 11. Save checkpoint ------------------------------------------------------
saveRDS(seu, "data/seurat_checkpoint_05.rds")

cat("\nScript 05 complete.\n")

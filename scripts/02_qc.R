# =============================================================================
# Assignment 4 - scRNA-seq Analysis
# Script 02: Quality Control & Filtering
# =============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)

# Load checkpoint from script 01
seu <- readRDS("data/seurat_checkpoint_01.rds")

# --- 1. Calculate mitochondrial percentage ------------------------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

summary(seu@meta.data$percent.mt)

# --- 2. Visualise QC metrics before filtering ---------------------------------
# Group by orig.ident (one bar per sample) to spot any outlier samples

p_vln_before <- VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by  = "orig.ident",
  pt.size   = 0,
  ncol      = 3
)
ggsave("figures/02_qc_violin_before.png",
       p_vln_before, width = 16, height = 5, dpi = 150)

# Scatter: nCount vs nFeature
# Scatter: nCount vs percent.mt (high mt = dying cells)
p_scatter <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                             group.by = "organ_custom") +
             FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt",
                             group.by = "organ_custom")
ggsave("figures/02_qc_scatter.png",
       p_scatter, width = 12, height = 5, dpi = 150)

# --- 3. Inspect the distribution to choose thresholds ------------------------
# Look at these quantiles to inform cutoffs
quantile(seu@meta.data$nFeature_RNA, probs = c(0.01, 0.05, 0.95, 0.99))
quantile(seu@meta.data$nCount_RNA,   probs = c(0.01, 0.05, 0.95, 0.99))
quantile(seu@meta.data$percent.mt,   probs = c(0.90, 0.95, 0.99))

# --- 4. Filter cells ----------------------------------------------------------
# Thresholds to justify in introduction:
#   nFeature_RNA > 200  : removes empty droplets (too few genes = no real cell)
#   nFeature_RNA < 6000 : removes likely doublets (two cells in one droplet)
#   nCount_RNA   < 40000: removes high-count outliers / doublets
#   percent.mt   < 15   : removes dying/damaged cells

seu_filtered <- subset(
  seu,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 6000 &
           nCount_RNA   < 40000 &
           percent.mt   < 15
)

# How many cells were removed?
cat("Cells before filtering:", ncol(seu),          "\n")
cat("Cells after filtering: ", ncol(seu_filtered), "\n")
cat("Cells removed:         ", ncol(seu) - ncol(seu_filtered), "\n")

# --- 5. Visualise QC metrics after filtering ----------------------------------
p_vln_after <- VlnPlot(
  seu_filtered,
  features  = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by  = "orig.ident",
  pt.size   = 0,
  ncol      = 3
)
ggsave("figures/02_qc_violin_after.png",
       p_vln_after, width = 16, height = 5, dpi = 150)

# --- 6. Handle cells with missing mouse_id ------------------------------------
table(seu_filtered@meta.data$organ_custom[seu_filtered@meta.data$mouse_id == ""])
table(seu_filtered@meta.data$time[seu_filtered@meta.data$mouse_id == ""])

# For pseudobulk DE analysis, need replicate labels.
# Add a clean replicate column that marks blanks explicitly:
seu_filtered@meta.data$replicate <- ifelse(
  seu_filtered@meta.data$mouse_id == "",
  "unassigned",
  seu_filtered@meta.data$mouse_id
)

# --- 7. Save checkpoint -------------------------------------------------------
saveRDS(seu_filtered, "data/seurat_checkpoint_02.rds")

cat("\nScript 02 complete. Check figures in results/figures/\n")
cat("Review threshold choices and update if needed before proceeding.\n")

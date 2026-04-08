# =============================================================================
# Assignment 4 - scRNA-seq Analysis
# Script 03: Normalization, Integration, Clustering & UMAP
# =============================================================================

# --- Install packages (run once, then comment out) ---------------------------
# install.packages("harmony")
# BiocManager::install("harmony")  # if above doesn't work

library(Seurat)
library(tidyverse)
library(ggplot2)
library(harmony)

# Load checkpoint from script 02
seu <- readRDS("data/seurat_checkpoint_02.rds")

# --- 1. Normalize with SCTransform --------------------------------------------
# SCTransform is preferred over basic LogNormalize (lecture 17):
# - Uses negative binomial regression to normalize for sequencing depth
# - Automatically identifies highly variable features
# - Regress out percent.mt to remove dying-cell signal from clustering
#
# Note: SCTransform replaces the need to separately run:
#   NormalizeData() -> FindVariableFeatures() -> ScaleData()
# This will take a few minutes on 156k cells.

seu <- SCTransform(
  seu,
  vars.to.regress = "percent.mt",
  verbose         = TRUE
)

# --- 2. PCA -------------------------------------------------------------------
seu <- RunPCA(seu, verbose = FALSE)

# Elbow plot: helps choose how many PCs to use downstream
# Look for where the curve flattens — that's your cutoff
p_elbow <- ElbowPlot(seu, ndims = 50)
ggsave("results/figures/03_elbow_plot.png", p_elbow, width = 7, height = 5, dpi = 150)

# Visually inspect the elbow plot before proceeding.
# A common choice is 30 PCs for large heterogeneous datasets.
# We use 30 below — adjust if your elbow is clearly earlier (e.g. 20).
N_PCS <- 30

# --- 3. Check for batch effects before integration ----------------------------
# Run a preliminary UMAP on PCA (no correction) to see if cells cluster
# by sample/tissue rather than by cell type — if so, Harmony is needed.

seu <- RunUMAP(seu, dims = 1:N_PCS, reduction = "pca",
               reduction.name = "umap.uncorrected", verbose = FALSE)

p_batch_check <- DimPlot(seu, reduction = "umap.uncorrected",
                         group.by = "orig.ident", pt.size = 0.1) +
  DimPlot(seu, reduction = "umap.uncorrected",
          group.by = "organ_custom", pt.size = 0.1)
ggsave("results/figures/03_umap_precorrection.png",
       p_batch_check, width = 14, height = 6, dpi = 150)

# If cells cluster strongly by orig.ident (sample) rather than biology,
# Harmony correction is warranted. With 3 tissues and 5 timepoints across
# 3 mice, batch effects are very likely — we apply Harmony regardless.

# --- 4. Harmony batch correction ----------------------------------------------
# Integrate across mouse_id (biological replicate) to remove per-mouse
# technical variation while preserving biological signal.
# We use orig.ident as the batch variable since it captures both
# tissue and timepoint per sample — this is the finest-grain correction.
#
# Note: cells with mouse_id == "unassigned" are still included here;
# Harmony operates on the PCA embedding, not on metadata labels.

seu <- RunHarmony(
  seu,
  group.by.vars   = "orig.ident",
  reduction       = "pca",
  reduction.save  = "harmony",
  verbose         = TRUE
)

# --- 5. UMAP on Harmony embedding ---------------------------------------------
seu <- RunUMAP(seu, dims = 1:N_PCS, reduction = "harmony",
               reduction.name = "umap", verbose = FALSE)

# Visual check: do cells now mix by sample but separate by biology?
p_post_harmony <- DimPlot(seu, reduction = "umap",
                          group.by = "orig.ident", pt.size = 0.1) +
  DimPlot(seu, reduction = "umap",
          group.by = "organ_custom", pt.size = 0.1)
ggsave("results/figures/03_umap_postcorrection.png",
       p_post_harmony, width = 14, height = 6, dpi = 150)

# --- 6. Clustering ------------------------------------------------------------
# FindNeighbors builds the shared nearest-neighbour graph on the
# Harmony embedding (not raw PCA) so batch correction carries through.
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)

# Resolution controls cluster granularity (lecture 17):
#   Lower (0.3-0.5) = fewer, broader clusters
#   Higher (0.8-1.2) = more, finer clusters
# For a large heterogeneous dataset (nasal mucosa + immune cells),
# 0.5 is a reasonable starting point. We also try 0.3 and 0.8 for comparison.
# You may wish to revisit this after annotation.

seu <- FindClusters(seu, resolution = 0.3, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.8, verbose = FALSE)

# Set 0.5 as the active clustering for now
Idents(seu) <- "SCT_snn_res.0.5"

# --- 7. UMAP plots at each resolution -----------------------------------------
p_res <- DimPlot(seu, reduction = "umap", group.by = "SCT_snn_res.0.3",
                 label = TRUE, pt.size = 0.1) + ggtitle("Resolution 0.3") +
  DimPlot(seu, reduction = "umap", group.by = "SCT_snn_res.0.5",
          label = TRUE, pt.size = 0.1) + ggtitle("Resolution 0.5") +
  DimPlot(seu, reduction = "umap", group.by = "SCT_snn_res.0.8",
          label = TRUE, pt.size = 0.1) + ggtitle("Resolution 0.8")
ggsave("results/figures/03_umap_resolutions.png",
       p_res, width = 18, height = 6, dpi = 150)

# Final UMAP coloured by tissue and timepoint for biological context
p_meta <- DimPlot(seu, reduction = "umap", group.by = "organ_custom",
                  pt.size = 0.1, label = FALSE) + ggtitle("Tissue") +
  DimPlot(seu, reduction = "umap", group.by = "time",
          pt.size = 0.1, label = FALSE) + ggtitle("Timepoint")
ggsave("results/figures/03_umap_metadata.png",
       p_meta, width = 14, height = 6, dpi = 150)

# --- 8. Save checkpoint -------------------------------------------------------
saveRDS(seu, "data/seurat_checkpoint_03.rds")

cat("\nScript 03 complete.\n")
cat("Check results/figures/03_elbow_plot.png — if the elbow is clearly\n")
cat("before PC 30, update N_PCS and re-run from step 5 onward.\n")
cat("Check results/figures/03_umap_resolutions.png and choose a resolution\n")
cat("that gives biologically meaningful cluster separation.\n")
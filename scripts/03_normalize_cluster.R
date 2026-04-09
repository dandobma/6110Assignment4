# =============================================================================
# Assignment 4 - scRNA-seq Analysis
# Script 03: Normalization, Integration, Clustering & UMAP
# =============================================================================

# --- Install packages (run once, then comment out) ---------------------------
# BiocManager::install("harmony")  
# BiocManager::install("glmGamPoi")

library(Seurat)
library(tidyverse)
library(ggplot2)
library(harmony)
library(glmGamPoi)

# Load checkpoint from script 02
seu <- readRDS("data/seurat_checkpoint_02.rds")

# --- 1. Normalize with SCTransform --------------------------------------------
# Did not have enough memory to run this locally, switched to log-normalization
#gc()
#options(future.globals.maxSize = 8000 * 1024^2)  # sets limit to 8 GB
#seu <- SCTransform(
#  seu,
#  vars.to.regress = "percent.mt",
#  method = "glmGamPoi",
#  variable.features.n = 2000,
#  verbose         = TRUE
#)
seu <- NormalizeData(seu, verbose = FALSE)

seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)

seu <- ScaleData(seu, vars.to.regress = "percent.mt", verbose = FALSE)


# --- 2. PCA -------------------------------------------------------------------
seu <- RunPCA(seu, verbose = FALSE)

# Elbow plot: helps choose how many PCs to use downstream
p_elbow <- ElbowPlot(seu, ndims = 50)
ggsave("figures/03_elbow_plot.png", p_elbow, 
       width = 7, height = 5, dpi = 150, bg = "white")

# Visually inspect the elbow plot before proceeding.
N_PCS <- 30

# --- 3. Check for batch effects before integration ----------------------------
# Run a preliminary UMAP on PCA (no correction) to see if cells cluster by sample/tissue rather than by cell type

seu <- RunUMAP(seu, dims = 1:N_PCS, reduction = "pca",
               reduction.name = "umap.uncorrected", verbose = FALSE)

p_batch_check <- DimPlot(seu, reduction = "umap.uncorrected",
                         group.by = "orig.ident", pt.size = 0.1) +
  DimPlot(seu, reduction = "umap.uncorrected",
          group.by = "organ_custom", pt.size = 0.1)
ggsave("figures/03_umap_precorrection.png",
       p_batch_check, width = 14, height = 6, dpi = 150)

# --- 4. Harmony batch correction ----------------------------------------------
# Integrate across mouse_id (biological replicate) to remove per-mouse
# technical variation while preserving biological signal.
# Use orig.ident as the batch variable since it captures both
# tissue and timepoint per sample

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
ggsave("figures/03_umap_postcorrection.png",
       p_post_harmony, width = 14, height = 6, dpi = 150)

# --- 6. Clustering ------------------------------------------------------------
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)

# Resolution controls cluster granularity:
#   Lower (0.3-0.5) = fewer, broader clusters
#   Higher (0.8-1.2) = more, finer clusters
seu <- FindClusters(seu, resolution = 0.3, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.8, verbose = FALSE)

# Set 0.5 as the active clustering for now

# --- 7. UMAP plots at each resolution -----------------------------------------
Idents(seu) <- "RNA_snn_res.0.5"

p_res <- DimPlot(seu, reduction = "umap", group.by = "RNA_snn_res.0.3",
                 label = TRUE, pt.size = 0.1) + ggtitle("Resolution 0.3") +
  DimPlot(seu, reduction = "umap", group.by = "RNA_snn_res.0.5",
          label = TRUE, pt.size = 0.1) + ggtitle("Resolution 0.5") +
  DimPlot(seu, reduction = "umap", group.by = "RNA_snn_res.0.8",
          label = TRUE, pt.size = 0.1) + ggtitle("Resolution 0.8")

ggsave("figures/03_umap_resolutions.png",
       p_res, width = 18, height = 6, dpi = 150)

# Final UMAP coloured by tissue and timepoint for biological context
p_meta <- DimPlot(seu, reduction = "umap", group.by = "organ_custom",
                  pt.size = 0.1, label = FALSE) + ggtitle("Tissue") +
  DimPlot(seu, reduction = "umap", group.by = "time",
          pt.size = 0.1, label = FALSE) + ggtitle("Timepoint")
ggsave("03_umap_metadata.png",
       p_meta, width = 14, height = 6, dpi = 150)

# --- 8. Save checkpoint -------------------------------------------------------
saveRDS(seu, "data/seurat_checkpoint_03.rds")

cat("\nScript 03 complete.\n")

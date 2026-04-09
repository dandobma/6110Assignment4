# =============================================================================
# Assignment 4 - scRNA-seq Analysis
# Script 07: Trajectory Analysis (Slingshot)
# =============================================================================

# --- Install packages (run once, then comment out) ---------------------------
# BiocManager::install("slingshot")
# BiocManager::install("SingleCellExperiment")
# install.packages("viridis")

library(Seurat)
library(tidyverse)
library(ggplot2)
library(slingshot)
library(SingleCellExperiment)
library(viridis)

# Load checkpoint from script 06
seu <- readRDS("data/seurat_checkpoint_06.rds")

# =============================================================================
# PART 1: SUBSET TO MYELOID CELLS IN RM
# =============================================================================

# --- 1. Subset myeloid cells --------------------------------------------------
# Include: Macrophages, Monocytes, Dendritic cells, Neutrophils, Immature neutrophils
# These are the cells most likely to show differentiation dynamics during influenza infection

myeloid_types <- c("Macrophages", "Monocytes", "Dendritic cells",
                   "Neutrophils", "Immature neutrophils")

myeloid <- subset(seu,
                  subset = cell_type    %in% myeloid_types &
                           organ_custom == "RM")

cat("Myeloid cells in RM:\n")
print(table(myeloid$cell_type))
cat("\nBy timepoint:\n")
print(table(myeloid$cell_type, myeloid$time))

# --- 2. Re-cluster myeloid cells ----------------------------------------------
# Re-run PCA and UMAP on just the myeloid subset for better resolution of trajectories within this compartment

myeloid <- NormalizeData(myeloid, verbose = FALSE)
myeloid <- FindVariableFeatures(myeloid, nfeatures = 2000, verbose = FALSE)
myeloid <- ScaleData(myeloid, verbose = FALSE)
myeloid <- RunPCA(myeloid, verbose = FALSE)

# Elbow plot for myeloid subset
p_elbow_m <- ElbowPlot(myeloid, ndims = 30)
ggsave("figures/07_elbow_myeloid.png",
       p_elbow_m, width = 6, height = 4, dpi = 150, bg = "white")

# Use 15 PCs for the myeloid subset (fewer cells, less complexity)
myeloid <- RunUMAP(myeloid, dims = 1:15, verbose = FALSE)

# UMAP coloured by cell type
p_myeloid_umap <- DimPlot(myeloid, reduction = "umap",
                           group.by = "cell_type",
                           label = TRUE, repel = TRUE, pt.size = 0.5) +
                  ggtitle("Myeloid cells in RM tissue") +
                  theme(legend.position = "right")

ggsave("figures/07_umap_myeloid.png",
       p_myeloid_umap, width = 9, height = 6, dpi = 150, bg = "white")

# UMAP coloured by timepoint
p_myeloid_time <- DimPlot(myeloid, reduction = "umap",
                           group.by = "time",
                           pt.size = 0.5) +
                  ggtitle("Myeloid cells by timepoint") +
                  scale_color_viridis_d(option = "plasma") +
                  theme(legend.position = "right")

ggsave("figures/07_umap_myeloid_time.png",
       p_myeloid_time, width = 9, height = 6, dpi = 150, bg = "white")

# =============================================================================
# PART 2: SLINGSHOT TRAJECTORY ANALYSIS
# =============================================================================

# --- 3. Convert to SingleCellExperiment for Slingshot ------------------------
sce <- as.SingleCellExperiment(myeloid)

# --- 4. Run Slingshot ---------------------------------------------------------
# Start lineage from Monocytes and allow Slingshot to infer connections to other cell types
# clusterLabels: use cell type annotations
# reducedDim: use UMAP embedding

sce <- slingshot(
  sce,
  clusterLabels = sce$cell_type,
  reducedDim    = "UMAP",
  start.clus    = "Monocytes"   # monocytes as root
)

cat("\nSlingshot lineages:\n")
print(slingLineages(sce))

# --- 5. Extract pseudotime values ---------------------------------------------
# slingPseudotime returns a matrix with one column per lineage
pt <- slingPseudotime(sce)
cat("\nNumber of lineages identified:", ncol(pt), "\n")
print(head(pt))

# Add pseudotime to metadata
# If multiple lineages, use the mean pseudotime across lineages
myeloid$pseudotime <- rowMeans(pt, na.rm = TRUE)

# Also add per-lineage pseudotime if more than one lineage
if (ncol(pt) > 1) {
  for (i in seq_len(ncol(pt))) {
    myeloid[[paste0("pseudotime_L", i)]] <- pt[, i]
  }
}

# --- 6. UMAP coloured by pseudotime ------------------------------------------
p_pseudotime <- FeaturePlot(myeloid, reduction = "umap",
                             features = "pseudotime",
                             pt.size = 0.8) +
               scale_color_viridis_c(option = "magma", name = "Pseudotime") +
               ggtitle("Myeloid pseudotime (RM tissue)")

ggsave("figures/07_umap_pseudotime.png",
       p_pseudotime, width = 8, height = 6, dpi = 150, bg = "white")

# --- 7. Overlay Slingshot curves on UMAP -------------------------------------
# Extract curve coordinates for plotting
curves     <- slingCurves(sce, as.df = TRUE)
umap_coords <- reducedDim(sce, "UMAP")

umap_df <- data.frame(
  UMAP1     = umap_coords[, 1],
  UMAP2     = umap_coords[, 2],
  cell_type = sce$cell_type,
  pseudotime = rowMeans(pt, na.rm = TRUE)
)

p_curves <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point(size = 0.6, alpha = 0.7) +
  geom_path(data = curves,
            aes(x = umap_1, y = umap_2, group = Lineage),
            color = "black", linewidth = 1.2, inherit.aes = FALSE) +
  labs(title = "Slingshot trajectories — Myeloid cells (RM)",
       color = "Cell type") +
  theme_classic() +
  theme(legend.position = "right")

ggsave("figures/07_trajectories.png",
       p_curves, width = 9, height = 6, dpi = 150, bg = "white")

# --- 8. Pseudotime by timepoint violin plot -----------------------------------

myeloid$time <- factor(myeloid$time,
                       levels = c("Naive", "D02", "D05", "D08", "D14"))

p_pt_violin <- ggplot(myeloid@meta.data,
                      aes(x = time, y = pseudotime,
                          fill = time)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  scale_fill_viridis_d(option = "plasma") +
  labs(title = "Pseudotime distribution by timepoint (Myeloid, RM)",
       x     = "Timepoint",
       y     = "Pseudotime") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("figures/07_pseudotime_by_timepoint.png",
       p_pt_violin, width = 8, height = 5, dpi = 150, bg = "white")

# --- 9. Save checkpoint -------------------------------------------------------
saveRDS(seu,     "data/seurat_checkpoint_07.rds")
saveRDS(myeloid, "data/seurat_myeloid_trajectory.rds")

cat("\nScript 07 complete.\n")

# =============================================================================
# Assignment 4 - scRNA-seq Analysis
# Script 04: Cluster Annotation
# =============================================================================

# --- Install packages (run once, then comment out) ---------------------------
# BiocManager::install("SingleR")
# BiocManager::install("celldex")
# BiocManager::install("scrapper")
# install.packages("devtools")
# devtools::install_github("immunogenomics/presto")

library(Seurat)
library(tidyverse)
library(ggplot2)
library(SingleR)
library(celldex)
library(scrapper)

# Load checkpoint from script 03
seu <- readRDS("data/seurat_checkpoint_03.rds")

# Set resolution 0.5 as active clustering
Idents(seu) <- "RNA_snn_res.0.5"

# =============================================================================
# PART 1: AUTOMATED ANNOTATION WITH SingleR
# =============================================================================

# --- 1. Load reference dataset ------------------------------------------------
ref_immgen <- ImmGenData()        # mouse immune cell atlas
ref_mouse  <- MouseRNAseqData()   # broader mouse cell type reference

# --- 2. Run SingleR -----------------------------------------------------------
# Extract normalized expression matrix
expr_matrix <- GetAssayData(seu, assay = "RNA", layer = "data")

# Run against immune reference
singler_immgen <- SingleR(
  test     = expr_matrix,
  ref      = ref_immgen,
  labels   = ref_immgen$label.main,
  clusters = Idents(seu)
)

# Run against broader mouse reference
singler_mouse <- SingleR(
  test     = expr_matrix,
  ref      = ref_mouse,
  labels   = ref_mouse$label.main,
  clusters = Idents(seu)
)

# --- 3. Inspect SingleR results -----------------------------------------------
# Compare what each reference calls each cluster side by side
singler_comparison <- data.frame(
  cluster     = rownames(singler_immgen),
  ImmGen      = singler_immgen$labels,
  MouseRNAseq = singler_mouse$labels
)
print(singler_comparison)

# Save to file for reference during manual annotation
write.csv(singler_comparison, "results/singler_labels.csv", row.names = FALSE)

# Diagnostic heatmap: shows confidence of each label assignment
# Low-confidence assignments (pale colours) need extra manual scrutiny
png("figures/04_singler_heatmap_immgen.png",
    width = 1400, height = 900, res = 120)
plotScoreHeatmap(singler_immgen)
dev.off()

# =============================================================================
# PART 2: MANUAL ANNOTATION WITH MARKER GENES
# =============================================================================

markers <- FindAllMarkers(
  seu, 
  only.pos            = TRUE,   # only upregulated markers
  min.pct             = 0.25,   # gene must be expressed in 25%+ of cluster cells
  logfc.threshold     = 0.5,
  test.use            = "wilcox",
  max.cells.per.ident = 500,   # downsample to 500 cells per cluster
  verbose             = TRUE
)

# Save full marker table
write.csv(markers, "results/cluster_markers.csv", row.names = FALSE)

# Top 5 markers per cluster for quick inspection
top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

print(top5, n = Inf)
write.csv(top5, "results/top5_markers.csv", row.names = FALSE)


# =============================================================================
# PART 3: APPLY ANNOTATIONS
# =============================================================================

cluster_annotations <- c(
  "0"  = "Olfactory sensory neurons",
  "1"  = "Macrophages",
  "2"  = "Basal epithelial cells",
  "3"  = "B cells",
  "4"  = "Monocytes",
  "5"  = "Neutrophils",
  "6"  = "Endothelial cells",
  "7"  = "Sustentacular cells",
  "8"  = "NK cells",
  "9"  = "Neuronal progenitors",
  "10" = "Olfactory epithelium",
  "11" = "Fibroblasts",
  "12" = "Dendritic cells",
  "13" = "Olfactory neurons",
  "14" = "Ductal epithelial cells",
  "15" = "Secretory epithelial cells",
  "16" = "Olfactory epithelium",
  "17" = "LNG secretory cells",
  "18" = "Immature neutrophils",
  "19" = "Smooth muscle cells",
  "20" = "Goblet cells",
  "21" = "Mesothelial cells",
  "22" = "Serous acinar cells",
  "23" = "LNG secretory cells",
  "24" = "Osteoblasts",
  "25" = "Tuft cells",
  "26" = "Schwann cells",
  "27" = "Ciliated epithelial cells",
  "28" = "Proliferating cells",
  "29" = "Chondrocytes",
  "30" = "Mature olfactory neurons",
  "31" = "Squamous epithelial cells",
  "32" = "Neuronal progenitors",
  "33" = "Proliferating cells",
  "34" = "T cells",
  "35" = "Vomeronasal sensory neurons"
)

# Check what the assignment actually produces
test <- cluster_annotations[as.character(Idents(seu))]
head(test)
length(test)
length(Cells(seu))
names(test) <- Cells(seu)
seu$cell_type <- test
head(seu$cell_type)
table(seu$cell_type)

# --- Final annotated UMAP ----------------------------------------------------
 p_annotated <- DimPlot(seu, reduction = "umap", group.by = "cell_type",
                         label = TRUE, repel = TRUE, pt.size = 0.1) +
                ggtitle("Cell type annotation") +
                theme(legend.position = "right")
 ggsave("figures/04_umap_annotated.png",
        p_annotated, width = 12, height = 8, dpi = 150, bg = "white")

# --- Feature plots of canonical markers ---------------------------------------
# Visualise key marker genes on the UMAP to validate annotation

 p_feat <- FeaturePlot(seu, reduction = "umap",
                       features = c("Cd3d", "Cd19", "Cd68",
                                    "Epcam", "Col1a1", "S100a8"),
                       ncol = 3, pt.size = 0.3)
 ggsave("figures/04_feature_plots_markers.png",
        p_feat, width = 15, height = 10, dpi = 150, bg = "white")

# --- Save checkpoint ----------------------------------------------------------
 saveRDS(seu, "data/seurat_checkpoint_04.rds")
 
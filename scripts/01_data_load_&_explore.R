# =============================================================================
# Assignment 4 - scRNA-seq Analysis
# Script 01: Load Seurat Object & Explore Metadata
# =============================================================================

# Install packages 
# install.packages("Seurat")
# install.packages("tidyverse")
# install.packages("ggplot2")

library(Seurat)
library(tidyverse)
library(ggplot2)

# --- 1. Load the Seurat object ------------------------------------------------
# wget -O data/seurat_ass4.rds "https://aacgenomicspublic.blob.core.windows.net/public/seurat_ass4.rds"

seu <- readRDS("data/seurat_ass4.rds")

# --- 2. Data exploration -----------------------------------------------------
print(seu)

# Check which assays are available
Assays(seu)

# Check which dimensionality reductions are already present
Reductions(seu)

# --- 3. Explore the metadata --------------------------------------------------
head(seu@meta.data)

# Get all column names in the metadata
colnames(seu@meta.data)

# How many cells are in each tissue type?
table(seu@meta.data$organ_custom)

# How many cells at each timepoint?
table(seu@meta.data$time)

# Cross-tabulate tissue × timepoint to confirm experimental design
table(seu@meta.data$organ_custom, seu@meta.data$time)

# How many cells per sample (biological replicate)?
table(seu@meta.data$mouse_id)        

table(seu@meta.data$organ_custom[seu@meta.data$mouse_id == ""])
table(seu@meta.data$time[seu@meta.data$mouse_id == ""])

# --- 4. Check gene counts and basic QC metrics --------------------------------
summary(seu@meta.data$nCount_RNA)   # total UMI counts per cell
summary(seu@meta.data$nFeature_RNA) # number of genes detected per cell

# Check if mitochondrial % is already calculated
"percent.mt" %in% colnames(seu@meta.data)

# If NOT already present, calculate it
if (!"percent.mt" %in% colnames(seu@meta.data)) {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
}

# --- 5. Quick QC violin plots -------------------------------------------------
# Visualise the three main QC metrics across samples before any filtering
# Save to results/figures/ for the README

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

p_qc <- VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "orig.ident",   
  pt.size = 0,           # hide individual points for legibility at scale
  ncol = 3
)

  ggsave("results/figures/01_qc_violin_raw.png", p_qc, width = 14, height = 5, dpi = 150)

# --- 6. Save a checkpoint -----------------------------------------------------
saveRDS(seu, "data/seurat_checkpoint_01.rds")

cat("\nStep 1 complete. Check the metadata column names above and update\n")
cat("'tissue', 'timepoint', and 'sample' references in subsequent scripts.\n")

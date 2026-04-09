# Assignment 4: scRNA-seq Analysis of Influenza A Virus Infection in the Nasal Mucosa
## Matei Dan-Dobre
### 1407506

# Introduction

Single-cell RNA sequencing (scRNA-seq) has emerged as a transformative approach in transcriptomics, enabling genome-wide gene expression profiling at the resolution of individual cells. Unlike bulk RNA sequencing, which measures population-averaged expression, scRNA-seq captures the transcriptional heterogeneity present within complex tissues, making it especially powerful for dissecting immune responses where diverse cell populations act in concert (Jovic et al., 2022). The technique involves isolating individual cells, barcoding their RNA, and sequencing the resulting library, with downstream analysis enabling clustering of cells by transcriptional similarity, annotation of cell types, and comparison of gene expression between conditions (Hao et al., 2021).

Influenza A Virus (IAV) infection of the upper respiratory tract triggers a well-characterized innate and adaptive immune response involving coordinated recruitment of myeloid and lymphoid cells to the site of infection (Iwasaki & Pillai, 2014). The nasal mucosa serves as the primary site of IAV entry and replication, and comprises three anatomically and functionally distinct regions: the respiratory mucosa (RM), the olfactory mucosa (OM), and the lateral nasal gland (LNG). These tissues differ substantially in their cellular composition, including differences in epithelial subtypes, neuronal populations, and resident immune cells, making them an ideal system for studying tissue-specific immune responses to viral infection (Chen et al., 2024). Understanding how infection alters cell type composition and transcriptional programs within each tissue compartment provides insight into both local antiviral defense mechanisms and potential pathological consequences such as olfactory dysfunction.

Tissue-resident macrophages play a central role in the innate immune response to respiratory viral infection. In the nasal mucosa, macrophages serve as sentinel cells that detect and respond to viral pathogens, producing cytokines and recruiting additional immune cells to the site of infection (Hashimoto et al., 2013). During influenza infection, blood-derived monocytes are recruited to infected tissues and can differentiate into monocyte-derived macrophages, replenishing tissue-resident populations that may be depleted during active viral replication (Guilliams et al., 2013). This monocyte-to-macrophage differentiation trajectory has been observed in lung tissue during IAV infection but has been less thoroughly characterized in the nasal mucosa (Misharin et al., 2017).

The dataset analyzed in this study was derived from a published single-cell atlas of the murine nasal mucosa across five timepoints of IAV infection: naive (uninfected), 2, 5, 8, and 14 days post-infection (dpi), spanning key phases of viral replication and immune clearance (Chen et al., 2024). Samples were collected from three tissue regions (RM, OM, LNG) with three biological replicates per condition. The primary analytical focus of this study is the respiratory mucosa, with a comparison of naive and D05 timepoints representing the transition from homeostasis to peak myeloid recruitment, using macrophages as the focal cell type for differential expression analysis. This timepoint was selected because D05 corresponds to peak myeloid cell recruitment as described in the original experimental design, making it a biologically relevant point at which to interrogate macrophage transcriptional reprogramming.

Computational analysis was performed using Seurat for clustering and dimensionality reduction (Hao et al., 2021), DESeq2 for pseudobulk differential expression (Love et al., 2014), clusterProfiler for gene set enrichment analysis (Yu et al., 2012), SingleR for automated cell type annotation (Aran et al., 2019), and Slingshot for trajectory inference (Street et al., 2018). Cell type composition was assessed across tissues and timepoints to characterize infection-driven changes in cellular abundance. Together, these analyses provide a comprehensive characterization of the transcriptional and compositional landscape of the murine nasal mucosa during IAV infection.

---

# Methods

### Computational Environment

All analyses were performed in R (version 4.4.x; R Core Team, 2024) within RStudio on Windows 11. All scripts are provided in the `scripts/` directory of this repository.

### Data Acquisition

A pre-processed Seurat object containing scRNA-seq data and metadata from Chen et al. (2024) was downloaded directly from a public repository hosted by the authors. The dataset encompasses 156,572 cells from three nasal tissue compartments (RM, OM, LNG) across five timepoints (naive, D02, D05, D08, D14) with three biological replicates per condition (n = 3 mice). The original data were generated using droplet-based single-cell sequencing and aligned to the mouse genome.

### Quality Control

Quality control metrics were assessed using the Seurat package (Hao et al., 2021). Mitochondrial gene content was calculated using `PercentageFeatureSet()` with the pattern `^mt-`, reflecting the lowercase mitochondrial gene naming convention of the mouse (*Mus musculus*) genome. Cells were retained based on the following thresholds: number of detected features (nFeature_RNA) between 200 and 6,000; total UMI counts (nCount_RNA) below 15,000; and mitochondrial fraction (percent.mt) below 15%. These thresholds were selected based on inspection of empirical quantile distributions: the 99th percentile of nFeature_RNA was 3,438, of nCount_RNA was 6,203, and of percent.mt was 13.5%, indicating the data had been substantially pre-filtered by the original authors. The applied thresholds therefore served as conservative safeguards to remove residual low-quality cells. A total of 435 cells (0.28%) were removed, yielding 156,137 cells for downstream analysis.

### Normalization

Data normalization was performed using `NormalizeData()` with default log-normalization parameters, followed by identification of the top 3,000 highly variable features using `FindVariableFeatures()`. Data were scaled and mitochondrial fraction was regressed out using `ScaleData()` with `vars.to.regress = "percent.mt"`. SCTransform (Hafemeister & Satija, 2019) was attempted but was not feasible due to memory constraints on a dataset of this size; log-normalization was used as the alternative.

### Dimensionality Reduction and Batch Correction

Principal component analysis (PCA) was performed using `RunPCA()`. Thirty principal components were retained based on inspection of the elbow plot, which showed a gradual decline in variance explained beyond PC 15–20 with no clear inflection point, making 30 PCs a conservative and well-supported choice. Batch correction was applied using Harmony (Korsunsky et al., 2019), integrating across sample identity (`orig.ident`) to remove per-sample technical variation while preserving biological signal. Uniform Manifold Approximation and Projection (UMAP) was computed on the Harmony-corrected embedding using `RunUMAP()` with 30 dimensions (McInnes et al., 2018).

### Clustering

Shared nearest-neighbour graphs were constructed using `FindNeighbors()` on the Harmony embedding with 30 dimensions. Clustering was performed at resolutions 0.3, 0.5, and 0.8 using `FindClusters()`. Resolution 0.5, yielding 36 clusters, was selected as the optimal balance between cluster granularity and biological interpretability for a dataset comprising three anatomically diverse tissue types.

### Cell Type Annotation

Automated annotation was performed using SingleR (Aran et al., 2019) with two reference datasets: ImmGenData for immune cell types and MouseRNAseqData for broader cell type coverage, both accessed via the celldex package. Annotation was performed at the cluster level by passing `clusters = Idents(seu)` to `SingleR()`. Manual annotation was performed in parallel using marker genes identified by `FindAllMarkers()` with a Wilcoxon rank-sum test (`min.pct = 0.25`, `logfc.threshold = 0.5`, `max.cells.per.ident = 500`). Canonical marker genes used for manual annotation included: *Cd3d* and *Cd3e* (T cells), *Cd19* and *Iglc2* (B cells), *Cd68* and *Trem2* (macrophages), *Ly6g* and *Cxcr2* (neutrophils), *Klrb1c* and *Ncr1* (NK cells), *Flt3* and *Cd209a* (dendritic cells), *Epcam* and *Krt15* (epithelial cells), *Col1a1* and *Apod* (fibroblasts), *Emcn* and *Ptprb* (endothelial cells), *Omp* and *Nrn1l* (olfactory neurons), *Neurog1* and *Neurod1* (neuronal progenitors), *Mpz* and *Gpr37l1* (Schwann cells), *Tmem212* and *Sntn* (ciliated epithelial cells), *Trpm5* and *Gnat3* (tuft cells), *Bglap* and *Ibsp* (osteoblasts), and *Col2a1* and *Cytl1* (chondrocytes). Final cluster labels were assigned by reconciling SingleR predictions with manual marker gene evidence.

### Differential Expression Analysis

Pseudobulk differential expression analysis was performed on macrophages from the RM tissue comparing Naive and D05 timepoints using DESeq2 (Love et al., 2014). Cells with unassigned biological replicates (mouse_id = "") were excluded, as replicate identity is required for pseudobulk aggregation. Raw counts were summed across all cells within each sample (mouse_id × timepoint combination) using `rowSums()`, yielding six pseudobulk samples (three Naive, three D05). Genes with fewer than 10 counts in fewer than two samples were excluded prior to model fitting. A Wald test was performed using the design formula `~ condition`, with Naive as the reference level. Genes with an adjusted p-value < 0.05 and |log₂ fold change| ≥ 1 were considered significantly differentially expressed. Multiple testing correction was applied using the Benjamini-Hochberg procedure (Benjamini & Hochberg, 1995).

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

### Gene Set Enrichment Analysis

Gene set enrichment analysis (GSEA) was performed using the `gseGO()` function from clusterProfiler (Yu et al., 2012) against Gene Ontology Biological Process (GO:BP) terms. Genes were ranked by a composite score of log₂ fold change × −log₁₀(p-value), preserving both effect size and statistical significance. Gene symbols were converted to Entrez IDs using `bitr()` with the `org.Mm.eg.db` annotation database. GSEA was run with `minGSSize = 15`, `maxGSSize = 500`, and `pvalueCutoff = 0.05`. Enrichment plots were generated using `gseaplot2()` from the enrichplot package.

### Cell Composition Analysis

Cell type proportions were calculated for the RM tissue across all five timepoints by aggregating cell counts per cell type per biological replicate and computing proportions relative to total RM cells per replicate. Immune cell proportions over time were visualized as a line plot with mean ± standard error. Macrophage proportions were modeled as a function of timepoint using a linear model (`lm(proportion ~ time_num)`), where timepoints were encoded numerically (0, 2, 5, 8, 14 dpi). Cell type composition was also compared across all three tissue types.

### Trajectory Analysis

Trajectory analysis was performed on myeloid cells (macrophages, monocytes, dendritic cells, neutrophils, and immature neutrophils) from the RM tissue using Slingshot (Street et al., 2018). Myeloid cells were re-normalized and re-clustered independently of the full dataset using log-normalization, 2,000 variable features, and UMAP with 15 principal components. Slingshot was applied to the UMAP embedding using cell type annotations as cluster labels and monocytes as the designated root cluster, based on the biological hypothesis that infiltrating monocytes differentiate into macrophages during infection resolution. Pseudotime values were extracted using `slingPseudotime()` and visualized on the UMAP and as violin plots stratified by timepoint.

---

## Results

A total of 156,137 cells passed quality control filtering and were retained for analysis. The UMAP of the full dataset, coloured by annotated cell type, reveals 36 distinct clusters spanning a broad range of cell lineages expected in the nasal mucosa (Figure 1). Immune cell populations including macrophages, monocytes, neutrophils, B cells, T cells, NK cells, and dendritic cells cluster together in the upper-left region of the UMAP. A large central region is occupied by olfactory sensory neurons, which represent the most abundant cell type in the dataset (37,258 cells), consistent with the inclusion of olfactory mucosa tissue. Stromal and structural cell types including fibroblasts, endothelial cells, chondrocytes, osteoblasts, and smooth muscle cells occupy the lower portion of the UMAP. Glandular epithelial populations including LNG secretory cells, serous acinar cells, and goblet cells cluster in the lower-left, reflecting the LNG tissue contribution. The spatial organization of cell types within the UMAP is biologically coherent and consistent with the known anatomy of the three tissue compartments included in this dataset.

![Figure 1](results/figures/04_umap_annotated.png)

**Figure 1:** UMAP of all 156,137 cells coloured by annotated cell type. Annotation was performed by combining SingleR automated labelling with two mouse reference datasets (ImmGenData and MouseRNAseqData) and manual marker gene validation. Thirty-six clusters were identified at resolution 0.5. Cell types are indicated by colour and labeled directly on the plot.

Feature plots of canonical lineage marker genes confirm the accuracy of cell type annotation (Figure 2). *Cd3d* expression is concentrated in a small discrete cluster corresponding to T cells within the immune compartment. *Cd19* marks the B cell cluster in the upper-centre of the UMAP. *Cd68*, a pan-macrophage marker, is expressed in the macrophage cluster in the left immune region. *Epcam* shows broad expression across multiple epithelial clusters spanning the central and lower UMAP. *Col1a1* is restricted to the fibroblast cluster in the lower stromal region. *S100a8* marks a tight cluster of neutrophils in the upper-left immune compartment. The spatial concordance between marker gene expression and cluster annotation supports the reliability of the annotated cell identities.

![Figure 2](results/figures/04_feature_plots_markers.png)

**Figure 2:** Feature plots showing normalized expression of six canonical cell type marker genes on the UMAP embedding. Expression level is indicated by colour intensity (purple = high). Marker genes shown are *Cd3d* (T cells), *Cd19* (B cells), *Cd68* (macrophages), *Epcam* (epithelial cells), *Col1a1* (fibroblasts), and *S100a8* (neutrophils).

Pseudobulk differential expression analysis comparing macrophages in the RM tissue between Naive and D05 timepoints identified 779 significantly differentially expressed genes (adjusted p-value < 0.05, |log₂FC| ≥ 1), of which 100 were upregulated and 679 were downregulated at D05 relative to Naive (Figure 3). The most significantly downregulated genes included multiple ribosomal protein genes (*Rps21*, *Rps27*, *Rpl39*, *Rpl37a*, *Rps28*, *Rpl34*) and immediate early response genes (*Fos*, *Ier2*, *Dusp1*, *Egr1*, *Jund*). Among the most significantly upregulated genes were *Ly6c2*, a marker of inflammatory monocyte-derived cells, *Isg15*, an interferon-stimulated gene with established antiviral function, and *Gm42418*, a long non-coding RNA. The volcano plot illustrates the asymmetry of the response, with the downregulated signal being quantitatively dominant but the upregulated genes including highly biologically relevant antiviral effectors.

![Figure 3](results/figures/05_volcano_macrophages_RM.png)

**Figure 3:** Volcano plot of pseudobulk differential expression in RM macrophages comparing D05 to Naive (n = 3 biological replicates per condition). The x-axis shows log₂ fold change and the y-axis shows −log₁₀ adjusted p-value. Red points indicate genes significantly upregulated at D05 (adjusted p-value < 0.05, log₂FC > 1); blue points indicate significantly downregulated genes (log₂FC < −1). Selected genes with |log₂FC| > 2 are labeled. Dashed lines indicate significance thresholds. Differential expression was performed using DESeq2 with the Benjamini–Hochberg multiple testing correction.

Gene set enrichment analysis of RM macrophages at D05 versus Naive revealed that all significantly enriched GO Biological Process terms were in the suppressed direction, with no pathways significantly enriched among upregulated genes at the GO:BP level (Figure 4). The top enriched terms included cytoplasmic translation (NES = −1.51, adjusted p = 2.51 × 10⁻⁷), translation (NES = −1.48, adjusted p = 2.51 × 10⁻⁷), ribonucleoprotein complex biogenesis (NES = −1.43, adjusted p = 1.42 × 10⁻⁴), and ribosome biogenesis (NES = −1.44, adjusted p = 1.56 × 10⁻³). All enriched terms were related to translational machinery and ribosomal biology. The GSEA enrichment curve for the top-ranked term, cytoplasmic translation, shows a steadily declining running enrichment score that reaches its minimum near the end of the ranked gene list, confirming that ribosomal and translation genes are concentrated among the most strongly downregulated genes at D05 (Figure 5).

![Figure 4](results/figures/05_gsea_dotplot_macrophages_RM.png)

**Figure 4:** GSEA dot plot of GO Biological Process terms enriched in RM macrophages at D05 versus Naive. All significantly enriched terms (adjusted p-value < 0.05) fall in the suppressed category, indicating downregulation of these pathways at D05. The x-axis represents the gene ratio (proportion of genes in the leading edge), dot size indicates the number of genes, and colour indicates the adjusted p-value. Analysis was performed using clusterProfiler with the org.Mm.eg.db mouse annotation database.

![Figure 5](results/figures/05_gsea_enrichment_plot.png)

**Figure 5:** GSEA enrichment plot for the top-ranked GO Biological Process term, cytoplasmic translation, in RM macrophages at D05 versus Naive. The running enrichment score (green line) declines continuously across the ranked gene list, with the minimum enrichment score occurring near the right end of the distribution. Tick marks indicate positions of cytoplasmic translation gene set members in the ranked list. The colour bar at the bottom indicates the ranking metric (red = upregulated at D05, blue = downregulated).

Cell type composition analysis of the RM tissue across all five timepoints reveals dynamic changes in immune cell proportions over the course of infection (Figure 6). The stacked bar plot shows that the relative abundance of immune cell types increases substantially from Naive to D05 and D08, with a corresponding relative decrease in structural and epithelial populations. The overall compositional shift is most pronounced at D05 and D08, consistent with the experimental design's designation of these timepoints as periods of peak myeloid recruitment and T cell infiltration, respectively.

![Figure 6](results/figures/06_composition_stacked_RM.png)

**Figure 6:** Stacked bar plot showing mean cell type proportions in the RM tissue at each of the five timepoints (Naive, D02, D05, D08, D14). Each bar represents the mean proportion across three biological replicates. Cell types are indicated by colour according to the legend. Proportions were calculated relative to total RM cells per replicate.

The dynamics of individual immune cell populations over the course of infection are shown in Figure 7. Monocytes exhibit a pronounced increase from Naive through D05 and peak at D08, before declining toward D14. NK cells show a similar trajectory, peaking at D08. In contrast, macrophages are most abundant in the Naive state and decline through D05 and D08, with partial recovery at D14. Neutrophils increase modestly from D02 onward. B cells and T cells, which are adaptive immune cell types, remain at low proportions across all timepoints, consistent with the relatively early timepoints captured by this dataset.

![Figure 7](results/figures/06_immune_proportions_RM.png)

**Figure 7:** Line plot showing mean proportions (± standard error) of eight immune cell types in the RM tissue across five timepoints. Proportions were calculated per biological replicate (n = 3) and summarized as mean ± SE. Each line represents one cell type as indicated by colour.

The proportion of macrophages per biological replicate across all five timepoints is shown in Figure 8. While the mean trend (dashed line) suggests a decline in macrophage proportion from Naive through D05–D08 followed by partial recovery at D14, there is substantial inter-replicate variability. A linear model of macrophage proportion as a function of timepoint was not statistically significant (p = 0.945, R² < 0.001), reflecting the non-linear nature of the response and the limited statistical power of three biological replicates per condition.

![Figure 8](results/figures/06_macrophage_proportion_RM.png)

**Figure 8:** Macrophage proportion per biological replicate in RM tissue across five timepoints. Each point represents one replicate, coloured by mouse identity. The dashed black line indicates the mean proportion across replicates at each timepoint. Proportions are expressed as a fraction of total RM cells per replicate.

To investigate whether monocytes differentiate into macrophages during the course of influenza infection, Slingshot trajectory analysis was performed on the myeloid compartment of the RM tissue. The UMAP of myeloid cells coloured by pseudotime shows a continuous gradient from low pseudotime (dark purple) in the monocyte region to high pseudotime (orange/yellow) in the macrophage and neutrophil regions (Figure 9). The pseudotime gradient transitions smoothly across the monocyte-macrophage boundary, indicating a transcriptional continuum rather than a discrete boundary between the two cell types. The neutrophil cluster on the right side of the UMAP shows a separate pseudotime progression reflecting the independent granulocyte maturation trajectory. Slingshot identified three lineages, all originating from the monocyte cluster: Lineage 1 (Monocytes → Dendritic cells → Immature neutrophils), Lineage 2 (Monocytes → Dendritic cells → Macrophages), and Lineage 3 (Monocytes → Dendritic cells → Neutrophils). The inferred trajectories, overlaid as black curves on the cell type UMAP, confirm the branched myeloid differentiation structure with dendritic cells occupying an intermediate transcriptional state (Figure 10).

![Figure 9](results/figures/07_umap_pseudotime.png)

**Figure 9:** UMAP of myeloid cells from RM tissue coloured by Slingshot pseudotime. Pseudotime was computed using monocytes as the root cluster. Colour indicates pseudotime value (dark purple = low, yellow = high) using the magma colour scale. Myeloid cells comprise macrophages, monocytes, dendritic cells, neutrophils, and immature neutrophils (n = 14,037 cells).

![Figure 10](results/figures/07_trajectories.png)

**Figure 10:** Slingshot trajectory curves overlaid on the myeloid cell UMAP coloured by cell type. Black curves indicate the three inferred differentiation lineages: Lineage 1 (Monocytes → Dendritic cells → Immature neutrophils), Lineage 2 (Monocytes → Dendritic cells → Macrophages), and Lineage 3 (Monocytes → Dendritic cells → Neutrophils). All lineages originate from the monocyte cluster.

---

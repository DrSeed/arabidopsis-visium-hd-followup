Spatial transcriptomics of Arabidopsis seedlings with Visium HD
================
Dr. Seed
November 2026

- [About this project](#about-this-project)
- [Step 1: Install the packages
  (once)](#step-1-install-the-packages-once)
- [Step 2: Download the data (once, 768
  MB)](#step-2-download-the-data-once-768-mb)
- [Step 3: Load the data into Seurat](#step-3-load-the-data-into-seurat)
- [Step 4: Quality control](#step-4-quality-control)
- [Step 5: Normalise, then check known marker
  genes](#step-5-normalise-then-check-known-marker-genes)
- [Step 6: Cluster the bins](#step-6-cluster-the-bins)
- [Step 7: Find marker genes and name the
  clusters](#step-7-find-marker-genes-and-name-the-clusters)
  - [Marker genes](#marker-genes)
  - [Heat map: a bug fix that matters
    later](#heat-map-a-bug-fix-that-matters-later)
  - [Naming the clusters](#naming-the-clusters)
- [Step 8: Find spatially variable genes (Moran’s
  I)](#step-8-find-spatially-variable-genes-morans-i)
  - [Spatially patterned does not mean specific to one
    cluster](#spatially-patterned-does-not-mean-specific-to-one-cluster)
- [Step 9: Differential expression between two
  tissues](#step-9-differential-expression-between-two-tissues)
- [Step 10: Pathway analysis (GO
  enrichment)](#step-10-pathway-analysis-go-enrichment)
  - [ORA: which processes characterise each
    tissue?](#ora-which-processes-characterise-each-tissue)
  - [GSEA: Cotyledon vs Germinating
    seed](#gsea-cotyledon-vs-germinating-seed)
- [Step 11: Gradients from the shoot apical
  meristem](#step-11-gradients-from-the-shoot-apical-meristem)
  - [Separating the seedling pieces](#separating-the-seedling-pieces)
  - [Fitting a smooth curve (GAM) for each
    gene](#fitting-a-smooth-curve-gam-for-each-gene)
- [Summary](#summary)
  - [This run compared with the
    original](#this-run-compared-with-the-original)
  - [Bugs fixed in the original
    scripts](#bugs-fixed-in-the-original-scripts)
  - [Problems hit along the way, and how they were
    solved](#problems-hit-along-the-way-and-how-they-were-solved)
  - [Limitations](#limitations)
  - [Ideas for next steps](#ideas-for-next-steps)
- [Session information](#session-information)

# About this project

This notebook reproduces, one step at a time, the Visium HD analysis of
*Arabidopsis thaliana* seedlings published by Anita Shrestha in
[ashrest58-commits/arabidopsis-visium-hd-spatial-transcriptomics](https://github.com/ashrest58-commits/arabidopsis-visium-hd-spatial-transcriptomics)
(MIT licence). The data are the public 10x Genomics dataset *Visium HD
3′ Arabidopsis* (Space Ranger 4.0.1 outputs).

I re-ran the whole analysis on a Mac in RStudio with Seurat 5.6. Along
the way I fixed several bugs in the original scripts and checked every
number against the original. **This notebook reports the numbers my run
actually produced.** Where they differ from the original, the text says
so and explains why.

**The questions the analysis asks:**

1.  Can we recover the tissues of a seedling (cotyledon, root, shoot
    meristem, …) from gene expression alone, and do they fall in the
    right places on the tissue section?
2.  Which genes have a spatial pattern?
3.  Which genes and biological processes differ between tissues?
4.  Does gene expression change gradually with distance from the shoot
    apical meristem (SAM)?

**Results at a glance (this run):**

- After quality control, 1,516 bins of 8 × 8 µm and 17,852 genes remain.
- Clustering gives 7 groups. Their marker genes and their positions on
  the slide match 7 seedling tissues: cotyledon, germinating seed, root
  elongation zone, root, meristematic tissue, hypocotyl and SAM.
- The most strongly patterned gene is **EXT3** (Moran’s I = 0.79), an
  extensin expressed along the whole root.
- Cotyledon vs germinating seed: 331 genes tested, but only **50** pass
  both significance and a two-fold change. One gene (CAB1,
  photosynthesis) is up in cotyledon and 49 are up in germinating seed,
  many of them ribosome and chloroplast-building genes.
- GO enrichment recovers the expected biology: photosynthesis in
  cotyledon, ribosome biogenesis in meristematic tissue, cell wall in
  root, phyllome (leaf) development in the SAM.
- **AFO** and **ANT** are expressed only in a small area around the SAM,
  and both fall off significantly with distance from it (ANT p =
  0.00026, the same as the original). **CAB1**, a photosynthesis gene,
  shows no significant trend (p = 0.14).

**How to run this notebook**

1.  Open the project’s `.Rproj` file in RStudio. This sets the working
    directory to the project folder, so the relative paths (`data/`,
    `results/`) work.
2.  Run the chunks in **Step 1** and **Step 2** once, by clicking the
    green ▶ on each chunk. They are marked `eval = FALSE`, so pressing
    **Knit** does not reinstall the packages or download the 768 MB
    dataset again every time.
3.  Press **Knit** (or *Run All*). The whole analysis takes about one
    minute on a laptop.

------------------------------------------------------------------------

# Step 1: Install the packages (once)

| Package | What it is used for |
|----|----|
| **Seurat** | The main toolkit: stores the data, normalises, clusters, finds markers |
| **hdf5r** | Reads the 10x count matrix (`.h5` file) |
| **arrow** | Reads the bin positions (`tissue_positions.parquet`) |
| **Rfast2** | Fast Moran’s I with permutation p-values (Step 8) |
| **ggplot2, dplyr, patchwork, ggrepel** | Plotting and data handling |
| **mgcv** | Smooth curves (GAMs) for the gradient analysis (Step 11) |
| **clusterProfiler, org.At.tair.db, enrichplot** | GO enrichment for Arabidopsis (Step 10). These come from Bioconductor, not CRAN, so they are installed with `BiocManager` |

``` r
options(timeout = 3600)   # R's default 60-second download limit is too short for big packages
install.packages(c("Seurat", "hdf5r", "arrow", "Rfast2", "ggplot2", "dplyr",
                   "patchwork", "ggrepel", "mgcv", "BiocManager"))
BiocManager::install(c("clusterProfiler", "org.At.tair.db", "enrichplot"),
                     update = FALSE, ask = FALSE)
```

> **Problem I hit:** the first installation stopped with a timeout. R
> gives each download 60 seconds by default, and large packages such as
> `arrow` need longer, especially while the dataset is downloading at
> the same time. Setting `options(timeout = 3600)` and installing
> *before* starting the data download fixed it.

> **Why Rfast2 matters:** when Rfast2 is installed, Seurat computes
> Moran’s I p-values by shuffling the data about 1,000 times. Without
> it, Seurat uses a formula from the `ape` package instead, and the
> p-values come out different (many become exactly 0). The Moran’s I
> values themselves are the same.

# Step 2: Download the data (once, 768 MB)

**What Visium HD measures:** the slide is covered in a continuous grid
of 2 × 2 µm squares. Each square has its own DNA barcode, so every mRNA
molecule captured on it can be traced back to its position. Space Ranger
(10x Genomics’ processing software) then adds neighbouring squares
together into **bins** of 8 × 8 µm and 16 × 16 µm.

**Why 8 µm bins?** A single 2 µm square holds very few molecules, too
few to work with. 16 µm bins have more molecules but blur the borders
between tissues. 8 µm is the usual compromise, and it is what the
original analysis used.

**A bin is not a cell.** Plant cells are large and irregular, and
segmenting cells from the tissue image did not work for this sample. So
the whole analysis treats each bin as if it were a cell. A bin on the
border between two tissues can contain a mixture of both.

``` r
options(timeout = 3600)
dir.create("data", showWarnings = FALSE)
url <- paste0("https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Arabidopsis/",
              "Visium_HD_3prime_Arabidopsis_binned_outputs.tar.gz")
download.file(url, "data/binned_outputs.tar.gz", mode = "wb")
untar("data/binned_outputs.tar.gz", exdir = "data")
file.exists("data/binned_outputs/square_008um/filtered_feature_bc_matrix.h5")   # should be TRUE
```

After extraction, the folder Seurat reads looks like this:

    data/binned_outputs/
    ├── square_002um/   2 µm squares
    ├── square_008um/   8 µm bins   <- the ones used here
    │   ├── filtered_feature_bc_matrix.h5   gene × bin count matrix (bins under tissue only)
    │   └── spatial/                        tissue image, bin positions, scale factors
    └── square_016um/   16 µm bins

# Step 3: Load the data into Seurat

`Load10X_Spatial()` reads the count matrix, the tissue image and the
position of every bin into one Seurat object. `bin.size = 8` tells it to
take the 8 µm folder. The data end up in an *assay* called
`Spatial.008um`, which is why later column names end in
`_Spatial.008um`.

``` r
library(Seurat); library(ggplot2); library(dplyr); library(patchwork)

obj <- Load10X_Spatial("data", bin.size = 8)
obj
```

    ## An object of class Seurat 
    ## 31109 features across 7356 samples within 1 assay 
    ## Active assay: Spatial.008um (31109 features, 0 variable features)
    ##  1 layer present: counts
    ##  1 spatial field of view present: slice1.008um

``` r
n_genes_raw <- nrow(obj)
n_bins_raw  <- ncol(obj)
head(obj@meta.data)
```

    ##                       orig.ident nCount_Spatial.008um nFeature_Spatial.008um
    ## s_008um_00262_00403-1          s                   26                     25
    ## s_008um_00437_00287-1          s                    5                      5
    ## s_008um_00233_00391-1          s                   76                     68
    ## s_008um_00456_00333-1          s                  192                    169
    ## s_008um_00531_00393-1          s                   50                     46
    ## s_008um_00463_00324-1          s                  121                    112

Each row of `meta.data` is one bin:

- **nCount_Spatial.008um** is the number of mRNA molecules (UMIs)
  captured in the bin.
- **nFeature_Spatial.008um** is the number of different genes detected
  in it.

**Check:** 31,109 genes × 7,356 bins. This is identical to the original
(31,109 × 7,356).

> When you run this chunk, Seurat warns *“Feature names cannot have
> underscores, replacing with dashes”*. This is harmless: a few gene
> names contain `_`, which Seurat does not allow, so it renames them.

``` r
SpatialDimPlot(obj) + NoLegend() + ggtitle("All 8 µm bins under tissue")
```

![](arabidopsis_analysis_files/figure-gfm/load-plot-1.png)<!-- -->

The capture area holds **several separate seedling pieces**. Keep this
in mind for Step 11.

# Step 4: Quality control

**Why:** many bins sit at the edges of the tissue or on background and
captured very little RNA. Such bins are mostly noise and would form
meaningless clusters, so they are removed before anything else.

``` r
summary(obj$nCount_Spatial.008um)   # UMIs per bin before filtering
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       0      14      39     138      84    5500

``` r
median_umi_raw <- median(obj$nCount_Spatial.008um)

VlnPlot(obj, "nCount_Spatial.008um", pt.size = 0) + NoLegend() + ggtitle("UMIs per bin") |
  SpatialFeaturePlot(obj, "nCount_Spatial.008um", pt.size.factor = 6)
```

![](arabidopsis_analysis_files/figure-gfm/qc-before-1.png)<!-- -->

The filter keeps bins with **more than 100 UMIs and more than 50
detected genes**. Then it drops genes that are not detected in any of
the remaining bins, because they carry no information.

``` r
obj <- subset(obj, nCount_Spatial.008um > 100 & nFeature_Spatial.008um > 50)
obj <- obj[rowSums(LayerData(obj, layer = "counts")) > 0, ]
dim(obj)   # genes, bins
```

    ## [1] 17852  1516

**Check:** 17,852 genes × 1,516 bins, identical to the original (17,852
× 1,516). That means only 21% of the bins are kept. This is a big cut,
but the median bin had only 39 UMIs before filtering, so most bins
really are nearly empty.

**Why no mitochondrial filter?** In single-cell studies, a high share of
mitochondrial reads usually marks dying cells. In plant tissue,
metabolically active cells naturally have many mitochondrial and
chloroplast transcripts, so this filter could remove real biology. The
original analysis skipped it, and so do I.

``` r
SpatialFeaturePlot(obj, "nCount_Spatial.008um", pt.size.factor = 6) +
  ggtitle("UMIs per bin after filtering")
```

![](arabidopsis_analysis_files/figure-gfm/qc-after-1.png)<!-- -->

> **Try it yourself:** change the cut-off from 100 to 50 UMIs and see
> how many more bins you keep.

# Step 5: Normalise, then check known marker genes

**Why normalise?** Some bins capture more RNA than others for purely
technical reasons. To compare bins fairly, `NormalizeData()` divides
each gene’s count by the bin’s total, multiplies by 10,000 and takes the
log (`log1p`). The log stops a few very highly expressed genes from
dominating.

``` r
obj <- NormalizeData(obj, verbose = FALSE)
```

**Why check markers now?** Before trusting any clustering, it is worth
confirming that known biology shows up where it should. **PHB**
(PHABULOSA) and **ERL2** (ERECTA-LIKE 2) are both expressed in and
around the shoot apical meristem. If the data are good, they should
light up in small, specific spots, not everywhere.

``` r
SpatialFeaturePlot(obj, c("PHB", "ERL2"), pt.size.factor = 6)
```

![](arabidopsis_analysis_files/figure-gfm/known-markers-1.png)<!-- -->

# Step 6: Cluster the bins

Clustering groups together bins with similar expression. **It uses
expression only, never position.** So if the clusters later turn out to
form coherent regions on the tissue, that is independent evidence that
they are real.

The steps, and why each one is needed:

| Step | What it does | Why |
|----|----|----|
| `FindVariableFeatures` | Picks the 3,000 genes that vary most between bins | Genes that are the same everywhere (housekeeping genes) add noise and no information |
| `ScaleData` | Centres each gene on 0 and scales it to the same spread | Otherwise highly expressed genes would dominate the next step |
| `RunPCA` | Summarises the 3,000 genes as 50 axes of co-varying genes | Removes noise and makes the next steps faster |
| `ElbowPlot` | Shows how much each axis explains | Helps choose how many axes to keep |
| `FindNeighbors` | Links each bin to the bins most similar to it (first 15 axes) | Builds the network the clustering works on |
| `FindClusters` | Finds tightly connected groups in that network (Louvain algorithm) | `resolution = 0.5` controls how many groups: higher gives more, smaller clusters |
| `RunUMAP` | Squashes the 15 axes into a 2D picture | For looking only: distances on a UMAP are not reliable |

``` r
obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
```

``` r
ElbowPlot(obj, ndims = 50)
```

![](arabidopsis_analysis_files/figure-gfm/elbow-1.png)<!-- -->

The curve flattens at around 10 axes. The original used 15. Keeping a
few more axes than the elbow is safer than keeping too few, because a
small tissue may only show up in a later axis.

``` r
obj <- FindNeighbors(obj, dims = 1:15, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:15, verbose = FALSE)

cluster_sizes <- table(cluster = Idents(obj))
knitr::kable(t(as.matrix(cluster_sizes)), caption = "Bins per cluster")
```

|   0 |   1 |   2 |   3 |   4 |   5 |   6 |
|----:|----:|----:|----:|----:|----:|----:|
| 567 | 440 | 183 | 145 |  83 |  57 |  41 |

Bins per cluster

**Check:** 7 clusters, the same number as the original. `FindClusters`
and `RunUMAP` use fixed random seeds by default, so you get the same
result every time you run this.

``` r
DimPlot(obj, label = TRUE) + ggtitle("UMAP") |
  SpatialDimPlot(obj, pt.size.factor = 10) + ggtitle("Clusters on the tissue")
```

![](arabidopsis_analysis_files/figure-gfm/clusters-plot-1.png)<!-- -->

# Step 7: Find marker genes and name the clusters

## Marker genes

`FindAllMarkers()` compares each cluster with all the other bins using a
Wilcoxon rank-sum test. `only.pos = TRUE` keeps only genes that are
*higher* in the cluster, because those are what define it. For each
cluster I list the top five markers that are at least two-fold higher
(`avg_log2FC > 1`).

``` r
markers <- FindAllMarkers(obj, only.pos = TRUE, verbose = FALSE)

top5 <- markers %>% group_by(cluster) %>% filter(avg_log2FC > 1) %>% slice_head(n = 5)
top5 %>% summarise(top_markers = paste(gene, collapse = ", ")) %>% knitr::kable()
```

| cluster | top_markers                              |
|:--------|:-----------------------------------------|
| 0       | CAB1, RBCS1A, LHB1B2, LHCB3, AT2G10940   |
| 1       | ICL, FSD1, MLS, PORB, CP29               |
| 2       | AGP31, BGLU21, AT2G43610, GLL22, JAL34   |
| 3       | EXT3, UCC2, AGP4, XYP2, AGP14            |
| 4       | GRP4, NOP10, HDA3, AT5G27120, DI21       |
| 5       | ACO2, AGP12, SIAR1, AT5G61660, AT2G04850 |
| 6       | AT4G29030, LTP6, CYP78A5, ERL2, AFO      |

## Heat map: a bug fix that matters later

The heat map needs *scaled* values for every gene it shows, and some top
markers are not among the 3,000 variable genes scaled in Step 6. So they
must be scaled now.

**The original script scaled only the ~35 heat-map genes.** That
replaces Seurat’s whole table of scaled values, so the 3,000 variable
genes lose theirs. Moran’s I (Step 8) runs on the scaled values, so in
the original it could only see those 35 genes. The fix is to scale the
**union** of the variable genes and the heat-map genes.

``` r
obj <- ScaleData(obj, features = union(VariableFeatures(obj), top5$gene), verbose = FALSE)
DoHeatmap(obj, features = top5$gene, size = 4)
```

![](arabidopsis_analysis_files/figure-gfm/heatmap-1.png)<!-- -->

## Naming the clusters

Cluster numbers can change between Seurat versions, so **check the
markers before naming**. The dot plot below shows one key marker per
tissue. The dot size is the share of bins in which the gene is detected,
and the colour is the average expression.

``` r
DotPlot(obj, features = c("CAB1", "ICL", "EXPA9", "EXT3", "GRP4", "XYP2", "CYP78A5")) + RotatedAxis()
```

![](arabidopsis_analysis_files/figure-gfm/dotplot-1.png)<!-- -->

Each marker has its highest average expression in the cluster expected
from the original mapping, so I use the same names. Two markers are less
clean than the others:

- **EXPA9** is highest on average in cluster 2, but it is *detected* in
  more bins of cluster 4 (64% vs 51%).
- **XYP2** is highest on average in cluster 5, but it is detected in
  only 11% of that cluster’s bins, and it is also expressed in root
  (cluster 3).

| Cluster | Name | Key evidence |
|----|----|----|
| 0 | Cotyledon | **CAB1**, RBCS1A, LHCB3: light harvesting and carbon fixation in photosynthesis |
| 1 | Germinating seed | **ICL**, MLS: the glyoxylate cycle, which turns stored seed oil into sugar in the first days after germination |
| 2 | Root elongation zone | **AGP31**, **EXPA9** (an expansin that loosens cell walls so cells can elongate) |
| 3 | Root | **EXT3**, XYP2, AGP4: root cell-wall and vascular proteins |
| 4 | Meristematic tissue | **GRP4**, NOP10, HDA3: RNA, ribosome and chromatin genes typical of dividing cells |
| 5 | hypocotyl | ACO2, AGP12 and **XYP2** (see above). This is the **weakest label** in the analysis |
| 6 | SAM | **CYP78A5** (KLUH), ERL2 and AFO, which mark the shoot apical meristem and young leaf primordia around it |

``` r
ids <- c("0" = "Cotyledon", "1" = "Germinating seed", "2" = "Root elongation zone",
         "3" = "Root", "4" = "Meristematic tissue", "5" = "hypocotyl", "6" = "SAM")
obj <- RenameIdents(obj, ids)
obj$tissue <- Idents(obj)

# The marker table was made before renaming, so it still uses cluster numbers.
# Add the tissue names to it as well, because Step 10 selects genes by tissue name.
markers$tissue <- unname(ids[as.character(markers$cluster)])
table(markers$tissue)   # number of marker genes per tissue
```

    ## 
    ##            Cotyledon     Germinating seed            hypocotyl 
    ##                 1564                 1008                  265 
    ##  Meristematic tissue                 Root Root elongation zone 
    ##                 7663                 2693                 3687 
    ##                  SAM 
    ##                 1112

``` r
saveRDS(obj, "results/annotated.rds")
```

> **Problems I hit at this step:**
>
> - `WhichCells(obj, idents = "Root")` failed with *“Cannot find the
>   following identities in the object: Root”*. `levels(Idents(obj))`
>   showed `"0" ... "6"`: the renaming chunk had never been run. Run
>   `levels(Idents(obj))` whenever you are unsure which names an object
>   is using.
> - Later, Step 10 failed with *“object ‘tissue’ not found”*. The
>   `markers` table was created before the renaming, so it had no
>   `tissue` column. That is why `markers$tissue` is added right here.

**The real test is the tissue section:** do the named tissues sit where
they should?

``` r
SpatialDimPlot(obj, cells.highlight = CellsByIdentities(obj), facet.highlight = TRUE,
               pt.size.factor = 4, ncol = 4)
```

![](arabidopsis_analysis_files/figure-gfm/tissues-spatial-1.png)<!-- -->

# Step 8: Find spatially variable genes (Moran’s I)

**What Moran’s I measures:** for each gene, it asks whether neighbouring
bins have similar expression.

- **+1** means the gene forms patches (neighbours are alike).
- **0** means random scatter.
- **−1** means a checkerboard (neighbours are opposite).

It uses **only positions and expression, never the cluster labels**, so
it is an independent way of finding spatial structure.

``` r
obj <- FindSpatiallyVariableFeatures(obj, features = VariableFeatures(obj),
                                     selection.method = "moransi", verbose = FALSE)
svg <- SVFInfo(obj, method = "moransi")
svg <- svg[!is.na(svg$MoransI_p.value), ]                  # genes that could not be tested
svg$p.adj <- p.adjust(svg$MoransI_p.value, method = "BH")  # correct for testing 3,000 genes
svg <- svg[order(svg$p.adj, -svg$MoransI_observed), ]
n_svg <- sum(svg$p.adj < 0.05)

knitr::kable(signif(head(svg, 10), 3), caption = "Top 10 spatially variable genes")
```

|           | MoransI_observed | MoransI_p.value |   p.adj |
|:----------|-----------------:|----------------:|--------:|
| EXT3      |            0.791 |        0.000976 | 0.00205 |
| AT1G26380 |            0.740 |        0.000976 | 0.00205 |
| AT4G18205 |            0.686 |        0.000976 | 0.00205 |
| IGMT2     |            0.587 |        0.000976 | 0.00205 |
| AGP31     |            0.545 |        0.000976 | 0.00205 |
| DI21      |            0.501 |        0.000976 | 0.00205 |
| RGF3      |            0.481 |        0.000976 | 0.00205 |
| RGF2      |            0.465 |        0.000976 | 0.00205 |
| AT1G02360 |            0.463 |        0.000976 | 0.00205 |
| AT1G56680 |            0.454 |        0.000976 | 0.00205 |

Top 10 spatially variable genes

**Why sort by p-value and then by Moran’s I?** With about 1,000
shuffles, the smallest p-value the test can return is about 0.001, and
many genes hit that floor. So the p-value alone cannot rank the top
genes; Moran’s I breaks the ties. In total, 1,966 of 3,000 tested genes
are significant (adjusted p \< 0.05).

**Check:** EXT3 is on top with Moran’s I = 0.7908, and AGP31 has 0.5449.
Both are identical to the original.

> **Honest correction:** in the original top-20 list, AGP31 is second.
> That list was drawn from only the ~35 heat-map genes, because of the
> scaling bug explained in Step 7. Ranked over all 3,000 variable genes,
> AGP31 is in **5th** place. AT1G26380, AT4G18205 and IGMT2 come before
> it. RGF2 and RGF3, two root-meristem growth factors, also appear in
> the top 10, which is a good biological sign.

``` r
SpatialFeaturePlot(obj, rownames(svg)[1:6], pt.size.factor = 6, ncol = 3)
```

![](arabidopsis_analysis_files/figure-gfm/moran-plot-1.png)<!-- -->

## Spatially patterned does not mean specific to one cluster

The table below shows where the bins expressing each gene fall: for each
gene, the percentage of its expressing bins that lie in each tissue.

``` r
overlap <- sapply(c("EXT3", "XYP2", "AGP31"), function(g) {
  expressing <- FetchData(obj, vars = g)[, 1] > 0
  round(100 * prop.table(table(Idents(obj)[expressing])), 1)
})
knitr::kable(overlap, caption = "% of each gene's expressing bins found in each tissue")
```

|                      | EXT3 | XYP2 | AGP31 |
|:---------------------|-----:|-----:|------:|
| Cotyledon            |  0.8 |  0.0 |   2.1 |
| Germinating seed     |  1.2 |  0.0 |   5.4 |
| Root elongation zone | 43.1 |  0.0 |  43.4 |
| Root                 | 46.7 | 82.9 |  15.7 |
| Meristematic tissue  |  8.2 |  2.4 |  23.8 |
| hypocotyl            |  0.0 | 14.6 |   2.4 |
| SAM                  |  0.0 |  0.0 |   7.2 |

% of each gene’s expressing bins found in each tissue

EXT3 has the strongest spatial pattern of all genes, yet only 46.7% of
its bins are in the Root cluster. Most of the rest are in the root
elongation zone: EXT3 marks the **whole root**, which spans two
clusters. XYP2 is more cluster-specific (82.9% in Root). A high Moran’s
I means “spatially clumped”, not “specific to one cluster”.

# Step 9: Differential expression between two tissues

`FindMarkers()` compares two groups of bins directly, here **Cotyledon
vs Germinating seed**, using a Wilcoxon rank-sum test.

- `min.pct = 0.1` only tests genes detected in at least 10% of the bins
  of one of the two groups.
- `logfc.threshold = 0.25` only tests genes with at least a small
  difference between the groups.

Both settings are pre-filters that save time. They do **not** decide
significance.

``` r
de <- FindMarkers(obj, ident.1 = "Cotyledon", ident.2 = "Germinating seed",
                  min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
de$gene <- rownames(de)
de$sig <- ifelse(de$p_val_adj < 0.05 & de$avg_log2FC >  1, "Up in Cotyledon",
          ifelse(de$p_val_adj < 0.05 & de$avg_log2FC < -1, "Up in Germinating seed",
                 "Not significant"))
write.csv(de, "results/DE_cot_vs_germ.csv", row.names = FALSE)

c(genes_tested = nrow(de), adjusted_p_below_0.05 = sum(de$p_val_adj < 0.05))
```

    ##          genes_tested adjusted_p_below_0.05 
    ##                   331                   208

``` r
table(de$sig)
```

    ## 
    ##        Not significant        Up in Cotyledon Up in Germinating seed 
    ##                    281                      1                     49

> **Honest correction:** the original reports this comparison as “330
> DEGs”. That is the number of genes **tested** (rows returned by
> `FindMarkers`), not the number of differentially expressed genes. In
> this run, 331 genes were tested and 208 have an adjusted p \< 0.05.
> Only 50 also differ at least two-fold (\|log2FC\| \> 1). The
> difference between 330 and 331 tested genes comes from the newer
> Seurat version.

**The largest fold changes are not the most reliable.** The genes with
the biggest fold changes are detected in only a few bins, so they are
not significant:

``` r
de %>% arrange(-avg_log2FC) %>% head(5) %>%
  dplyr::select(gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)),
         p_val_adj = formatC(p_val_adj, format = "g", digits = 2)) %>% knitr::kable(row.names = FALSE)
```

| gene      | avg_log2FC | pct.1 | pct.2 | p_val_adj |
|:----------|-----------:|------:|------:|:----------|
| CA1       |       2.50 | 0.109 | 0.041 | 0.49      |
| GRP3S     |       2.42 | 0.152 | 0.080 | 1         |
| AT2G45180 |       1.46 | 0.226 | 0.205 | 1         |
| PSBQ-2    |       1.44 | 0.150 | 0.132 | 1         |
| LHB1B1    |       1.44 | 0.150 | 0.134 | 1         |

For example, CA1 is 5.7 times higher in cotyledon, but it is detected in
only 11% of cotyledon bins against 4% of germinating-seed bins, and its
adjusted p is 0.5. The genes that pass both cut-offs:

``` r
de %>% filter(sig != "Not significant") %>% arrange(p_val_adj) %>% head(15) %>%
  dplyr::select(gene, sig, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)),
         p_val_adj = formatC(p_val_adj, format = "g", digits = 2)) %>% knitr::kable(row.names = FALSE)
```

| gene      | sig                    | avg_log2FC | pct.1 | pct.2 | p_val_adj |
|:----------|:-----------------------|-----------:|------:|------:|:----------|
| AT3G49910 | Up in Germinating seed |      -1.69 | 0.044 | 0.307 | 1e-22     |
| RPL12-C   | Up in Germinating seed |      -1.27 | 0.074 | 0.316 | 5.6e-16   |
| AT5G21020 | Up in Germinating seed |      -1.09 | 0.078 | 0.323 | 1.3e-15   |
| AT3G62870 | Up in Germinating seed |      -1.60 | 0.048 | 0.257 | 2.1e-15   |
| emb2171   | Up in Germinating seed |      -1.01 | 0.093 | 0.352 | 2.7e-15   |
| AT3G60245 | Up in Germinating seed |      -1.14 | 0.060 | 0.270 | 7.8e-14   |
| PORA      | Up in Germinating seed |      -1.23 | 0.046 | 0.241 | 8.9e-14   |
| AT4G27090 | Up in Germinating seed |      -1.04 | 0.065 | 0.282 | 9.8e-14   |
| AT1G52300 | Up in Germinating seed |      -1.04 | 0.065 | 0.270 | 3.1e-12   |
| NUC-L1    | Up in Germinating seed |      -1.00 | 0.074 | 0.277 | 2.3e-11   |
| BBC1      | Up in Germinating seed |      -1.02 | 0.060 | 0.250 | 2.4e-11   |
| AT3G08520 | Up in Germinating seed |      -1.53 | 0.026 | 0.177 | 3.4e-11   |
| AT2G27720 | Up in Germinating seed |      -1.05 | 0.053 | 0.232 | 8.2e-11   |
| AT5G56710 | Up in Germinating seed |      -1.47 | 0.034 | 0.184 | 2.9e-10   |
| AT2G47780 | Up in Germinating seed |      -1.14 | 0.048 | 0.214 | 4.3e-10   |

``` r
to_label <- de %>% filter(sig != "Not significant") %>% arrange(p_val_adj) %>% head(12)
ggplot(de, aes(avg_log2FC, -log10(p_val_adj + 1e-300), colour = sig)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60") +
  ggrepel::geom_text_repel(data = to_label, aes(label = gene), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c("Not significant" = "grey70", "Up in Cotyledon" = "#1b9e77",
                                 "Up in Germinating seed" = "#d95f02")) +
  labs(x = "log2 fold change (positive = higher in Cotyledon)", y = "-log10 adjusted p",
       colour = NULL, title = "Cotyledon vs Germinating seed") +
  theme_classic()
```

![](arabidopsis_analysis_files/figure-gfm/volcano-1.png)<!-- -->

`+ 1e-300` avoids `log10(0)`. Some p-values are so small that R stores
them as exactly 0, and without this those points would disappear from
the plot.

**What the result means:** only **CAB1** (light harvesting) is
significantly higher in cotyledon, and it is detected at a similar rate
in both tissues (62% vs 63%). The cotyledon differs mainly in *how much*
photosynthesis mRNA it has. The 49 genes higher in the germinating seed
include many ribosomal proteins (RPL12-C, RPL16B, RPS11-BETA), nucleolar
genes (NUC-L1) and chloroplast-building genes (PORA, ALB3, CHL-CPN10,
CLPP3). This fits young tissue that is still building its
protein-synthesis and chloroplast machinery. The glyoxylate-cycle genes
ICL and MLS are also strongly significant (adjusted p below 1e-11), but
just under the two-fold cut-off.

``` r
top_germ <- de %>% filter(sig == "Up in Germinating seed") %>% arrange(p_val_adj) %>% head(4)
SpatialFeaturePlot(obj, c("CAB1", "ICL", top_germ$gene), pt.size.factor = 6, ncol = 3)
```

![](arabidopsis_analysis_files/figure-gfm/de-spatial-1.png)<!-- -->

> **Caveat:** bins are not independent biological replicates. They all
> come from one section of one sample, so these p-values describe *this
> slide*, not Arabidopsis seedlings in general. The p-values are also
> over-confident, because neighbouring bins are not independent.

# Step 10: Pathway analysis (GO enrichment)

There are two flavours of GO (Gene Ontology) enrichment, and they answer
different questions.

- **ORA (over-representation analysis)** takes a *list* of genes, for
  example the strong markers of one tissue. It asks which GO terms
  appear in that list more often than expected by chance.
- **GSEA (gene set enrichment analysis)** takes *all* genes of a
  comparison, *ranked* by fold change. It asks whether the genes of a GO
  term pile up at the top or the bottom of the ranking.

## ORA: which processes characterise each tissue?

**Choosing each tissue’s gene list:** markers with adjusted p \< 0.05,
at least two-fold higher (`avg_log2FC > 1`) and detected in more than
25% of the tissue’s bins. Two tissues get special cut-offs:

- **Meristematic tissue** has hundreds of markers, so it gets a stricter
  cut-off (log2FC \> 1.5, detected in more than 40%) to keep only the
  strongest ones.
- **Hypocotyl** has very few strong markers, so it gets a looser one
  (log2FC \> 0.5, detected in more than 10%).

**The background (`universe`)** is every gene detected on this slide.
Only those genes could ever have become markers, so using the whole
genome as background would exaggerate enrichment.

``` r
library(clusterProfiler); library(org.At.tair.db); library(enrichplot)

sel <- markers %>% filter(p_val_adj < 0.05) %>%
  filter(case_when(tissue == "Meristematic tissue" ~ avg_log2FC > 1.5 & pct.1 > 0.40,
                   tissue == "hypocotyl"           ~ avg_log2FC > 0.5 & pct.1 > 0.10,
                   TRUE                            ~ avg_log2FC > 1.0 & pct.1 > 0.25))
table(sel$tissue)   # genes per tissue going into ORA
```

    ## 
    ##            Cotyledon     Germinating seed            hypocotyl 
    ##                   13                   25                   15 
    ##  Meristematic tissue                 Root Root elongation zone 
    ##                  313                  143                  192 
    ##                  SAM 
    ##                   21

``` r
ora <- compareCluster(split(sel$gene, sel$tissue), fun = "enrichGO",
                      OrgDb = org.At.tair.db, keyType = "SYMBOL", ont = "BP",
                      universe = rownames(obj), pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ora <- simplify(ora, cutoff = 0.7)   # drop near-duplicate GO terms (similarity > 0.7)
```

**Why `simplify()`?** GO terms are nested and overlap heavily. For
example, “rRNA processing” sits inside “ribosome biogenesis”.
`simplify()` keeps the most significant term from each group of very
similar terms.

``` r
dotplot(ora, showCategory = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("GO biological processes enriched in each tissue's markers")
```

![](arabidopsis_analysis_files/figure-gfm/ora-plot-1.png)<!-- -->

The numbers under each tissue are the genes that could be matched to GO
annotations.

``` r
as.data.frame(ora) %>% group_by(Cluster) %>% slice_min(p.adjust, n = 3, with_ties = FALSE) %>%
  dplyr::select(Cluster, Description, GeneRatio, p.adjust) %>%
  mutate(p.adjust = formatC(p.adjust, format = "g", digits = 2)) %>% knitr::kable()
```

| Cluster | Description | GeneRatio | p.adjust |
|:---|:---|:---|:---|
| Cotyledon | photosynthesis, light harvesting in photosystem I | 7/11 | 6.7e-15 |
| Cotyledon | response to high light intensity | 2/11 | 0.017 |
| Germinating seed | translation | 6/19 | 0.0012 |
| Meristematic tissue | translation | 24/121 | 2e-12 |
| Meristematic tissue | ribosome biogenesis | 18/121 | 5.4e-12 |
| Meristematic tissue | rRNA processing | 14/121 | 1.3e-09 |
| Root | nucleoside phosphate metabolic process | 16/96 | 2.5e-06 |
| Root | cell wall organization or biogenesis | 17/96 | 5e-05 |
| Root | carbohydrate derivative metabolic process | 16/96 | 0.00024 |
| Root elongation zone | translation | 17/111 | 2.3e-06 |
| Root elongation zone | protein folding | 12/111 | 1.6e-05 |
| Root elongation zone | protein maturation | 12/111 | 0.0005 |
| SAM | phyllome development | 6/15 | 0.015 |
| SAM | cytoplasmic translation | 2/15 | 0.042 |

> **Why `dplyr::select` and not just `select`?** Loading clusterProfiler
> also loads AnnotationDbi and IRanges. These packages have their own
> `select()`, `slice()` and `rename()` functions, which *mask* dplyr’s
> versions. Writing `dplyr::select()` makes R use the dplyr one.

**Interpretation:** each tissue gets the processes you would expect.

- Cotyledon: photosynthesis.
- Meristematic tissue: translation, ribosome biogenesis and rRNA
  processing, the protein-making machinery of dividing cells.
- Root: cell wall organisation and cellulose.
- SAM: phyllome development (*phyllome* means leaves and leaf-like
  organs).

> **Honest correction: hypocotyl has no enriched terms here, but 5 in
> the original.** In the original script, the special hypocotyl and
> meristematic cut-offs never took effect. The markers file stored
> cluster *numbers* (0–6), while the filter looked for tissue *names*.
> So every tissue got the default cut-off, which left hypocotyl with
> only 2 annotated genes. With so few genes, one or two hits can make a
> term look significant. Here the looser cut-off does apply: hypocotyl
> gets 15 genes, and none of its GO terms reaches adjusted p \< 0.05
> (the best is 0.15). Its top terms, such as “response to starvation”,
> are the same as in the original, but they are not significant. This is
> consistent with hypocotyl being the weakest cluster.

``` r
mapped <- suppressWarnings(bitr(rownames(obj), fromType = "SYMBOL", toType = "ENTREZID",
                                OrgDb = org.At.tair.db))
c(detected_genes = nrow(obj), with_a_GO_usable_symbol = length(unique(mapped$SYMBOL)))
```

    ##          detected_genes with_a_GO_usable_symbol 
    ##                   17852                    7519

> **A limitation:** only 42% of the detected genes can be matched by
> gene symbol. The others only have a locus code, such as AT2G45180, and
> are silently dropped. This is the reason for the *“…% of input gene
> IDs are fail to map”* warnings. Matching by TAIR locus ID
> (`keyType = "TAIR"`) instead would keep many more genes. The original
> analysis has the same limitation.

## GSEA: Cotyledon vs Germinating seed

The ranking is the fold change from Step 9, so **positive means higher
in Cotyledon**. In the plot:

- **activated** terms sit at the top of the ranking (higher in
  Cotyledon);
- **suppressed** terms sit at the bottom (higher in Germinating seed).

GSEA estimates its p-values by shuffling, so a fixed seed makes the
result repeatable.

``` r
ranked <- sort(setNames(de$avg_log2FC, de$gene), decreasing = TRUE)
set.seed(123)
gsea <- gseGO(ranked, OrgDb = org.At.tair.db, keyType = "SYMBOL", ont = "BP",
              minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, seed = TRUE, verbose = FALSE)
n_gsea_raw <- nrow(gsea)
gsea <- simplify(gsea, cutoff = 0.7)

as.data.frame(gsea)[, c("Description", "setSize", "NES", "p.adjust")] %>%
  mutate(NES = round(NES, 2), p.adjust = formatC(p.adjust, format = "g", digits = 2)) %>%
  knitr::kable(row.names = FALSE)
```

| Description                    | setSize |  NES | p.adjust |
|:-------------------------------|--------:|-----:|:---------|
| photosynthesis, light reaction |      32 | 3.54 | 6.2e-09  |
| response to light stimulus     |      37 | 2.45 | 0.00018  |
| cellular response to stress    |      18 | 2.10 | 0.0068   |
| biological_process             |     184 | 1.62 | 0.011    |
| cellular process               |     176 | 1.59 | 0.043    |
| metabolic process              |     158 | 1.55 | 0.015    |

``` r
dotplot(gsea, split = ".sign") + facet_grid(. ~ .sign)
```

![](arabidopsis_analysis_files/figure-gfm/gsea-plot-1.png)<!-- -->

``` r
gseaplot2(gsea, geneSetID = 1, title = gsea$Description[1])
```

![](arabidopsis_analysis_files/figure-gfm/gsea-plot-2.png)<!-- -->

**Interpretation:** the clearest result is **“photosynthesis, light
reaction”** at the top of the ranking (NES ≈ 3.5): photosynthesis genes
are higher in cotyledon, as expected. GSEA finds 19 terms before
simplifying and 6 after (the original: 22 and 8). All of them are
“activated”. Terms such as “biological_process” and “cellular process”
sit at the root of the GO tree and say little on their own.

> **Why are there no “suppressed” terms, even though 49 genes are up in
> germinating seed?** GSEA here only sees the 331 genes that passed the
> pre-filters in Step 9, and many of them have no GO symbol. That gives
> it little power. Running `FindMarkers` with `logfc.threshold = 0` and
> `min.pct = 0` and ranking all genes would be a better input. That is a
> natural next step.

# Step 11: Gradients from the shoot apical meristem

**The question:** does gene expression change gradually with distance
from the SAM (a gradient), or switch sharply (a boundary), or not depend
on distance at all?

## Separating the seedling pieces

The capture area holds several separate pieces. A distance between bins
on different pieces means nothing, so I first group bins by **position
only** (k-means on x and y, 8 groups). Then I keep the piece that
contains the most SAM bins.

``` r
library(mgcv)
coords <- GetTissueCoordinates(obj)
if ("cell" %in% colnames(coords)) rownames(coords) <- coords$cell

genes_grad <- c("CYP78A5", "AFO", "ANT", "CAB1", "RBCS1A")
df <- FetchData(obj, vars = genes_grad)
df$x <- coords[rownames(df), "x"]
df$y <- coords[rownames(df), "y"]
df$tissue <- Idents(obj)[rownames(df)]

set.seed(123)   # k-means starts from random centres
df$piece <- factor(kmeans(df[, c("x", "y")], centers = 8)$cluster)
table(df$piece, df$tissue)
```

    ##    
    ##     Cotyledon Germinating seed Root elongation zone Root Meristematic tissue
    ##   1        80               45                    0    0                   0
    ##   2        83                1                    0    0                   0
    ##   3         0                5                   90    1                  83
    ##   4       138              325                   93  135                   0
    ##   5       148               39                    0    0                   0
    ##   6        64               24                    0    0                   0
    ##   7         0                0                    0    9                   0
    ##   8        54                1                    0    0                   0
    ##    
    ##     hypocotyl SAM
    ##   1        51  37
    ##   2         0   0
    ##   3         0   0
    ##   4         3   4
    ##   5         1   0
    ##   6         2   0
    ##   7         0   0
    ##   8         0   0

``` r
ggplot(df, aes(x, y, colour = piece)) + geom_point(size = 0.8) +
  coord_fixed() + scale_y_reverse() + theme_minimal() +
  ggtitle("Bins grouped into pieces by position (k-means)")
```

![](arabidopsis_analysis_files/figure-gfm/pieces-plot-1.png)<!-- -->

``` r
p <- names(which.max(table(df$piece[df$tissue == "SAM"])))   # the piece containing the SAM
shoot <- droplevels(subset(df, piece == p))

# Coordinates are in image pixels. The scale-factor file says how many µm one pixel is.
# Check: neighbouring 8 µm bins are 29.2 px apart, and 29.2 px × 0.274 µm/px = 8 µm.
um_per_px <- jsonlite::fromJSON(
  "data/binned_outputs/square_008um/spatial/scalefactors_json.json")$microns_per_pixel

centre <- colMeans(shoot[shoot$tissue == "SAM", c("x", "y")])
shoot$dist_um <- sqrt((shoot$x - centre[1])^2 + (shoot$y - centre[2])^2) * um_per_px

shoot %>% group_by(tissue) %>%
  summarise(bins = n(), mean_distance_um = round(mean(dist_um))) %>%
  arrange(mean_distance_um) %>% knitr::kable()
```

| tissue           | bins | mean_distance_um |
|:-----------------|-----:|-----------------:|
| SAM              |   37 |               19 |
| Germinating seed |   45 |               49 |
| Cotyledon        |   80 |              116 |
| hypocotyl        |   51 |              148 |

**Check:** the SAM piece holds 80 cotyledon, 45 germinating-seed, 51
hypocotyl and 37 SAM bins (the original: 80, 44, 51, 38). The tissues
line up by distance from the SAM in the same order as in the original:
SAM, then germinating seed, cotyledon and hypocotyl. I report distances
in µm instead of pixels so they mean something physically.

``` r
ggplot(shoot, aes(x, y, colour = tissue)) + geom_point(size = 2) +
  coord_fixed() + scale_y_reverse() + theme_minimal() +
  ggtitle(paste("Piece", p, ": the shoot"))
```

![](arabidopsis_analysis_files/figure-gfm/shoot-plot-1.png)<!-- -->

## Fitting a smooth curve (GAM) for each gene

A GAM (generalised additive model) fits a smooth curve of expression
against distance, without assuming that the relationship is a straight
line.

- `k = 4` limits how wiggly the curve can be, so it cannot chase noise.
- The **p-value** tests whether the curve is flat.
- **Deviance explained** is the share of the bin-to-bin variation that
  distance accounts for.

``` r
gam_results <- do.call(rbind, lapply(genes_grad, function(g) {
  fit <- summary(gam(as.formula(paste(g, "~ s(dist_um, k = 4)")), data = shoot))
  data.frame(gene = g,
             p_value = format.pval(fit$s.table[1, "p-value"], digits = 2),
             deviance_explained = paste0(round(100 * fit$dev.expl, 1), "%"))
}))
knitr::kable(gam_results)
```

| gene    | p_value | deviance_explained |
|:--------|:--------|:-------------------|
| CYP78A5 | 0.0011  | 8.2%               |
| AFO     | \<2e-16 | 18.8%              |
| ANT     | 0.00026 | 9.6%               |
| CAB1    | 0.14    | 3.2%               |
| RBCS1A  | 0.52    | 1%                 |

``` r
plot_gene <- function(g) {
  ggplot(shoot, aes(dist_um, .data[[g]])) +
    geom_point(aes(colour = tissue), alpha = 0.6) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 4), colour = "black") +
    labs(title = g, x = "Distance from SAM centre (µm)", y = "log-normalised expression") +
    theme_minimal()
}
(plot_gene("AFO") | plot_gene("ANT") | plot_gene("CAB1")) + plot_layout(guides = "collect")
```

![](arabidopsis_analysis_files/figure-gfm/gam-plot-1.png)<!-- -->

**How to read the curves:** a cliff means a sharp boundary, a slope
means a gradient, and a flat line means distance from the SAM does not
matter.

- **AFO** (FILAMENTOUS FLOWER / YABBY1) is expressed only in bins within
  about 100 µm of the SAM centre, mostly in SAM bins, and is absent
  further away. This is a sharp, local pattern, as expected for a gene
  that marks the young leaf primordia beside the meristem.
- **ANT** (AINTEGUMENTA) shows the same kind of pattern, expressed only
  within about 60 µm of the SAM centre, in SAM and germinating-seed
  bins. The curve is significant (p = 0.00026, the same value as the
  original) but weaker than for AFO, because ANT is detected in fewer
  bins. The original describes ANT as a “modest gradient”. In this plot
  it looks more like a small local domain than a gentle slope.
- **CAB1** has no significant relationship with distance (p = 0.14; the
  original: 0.135). It is lower in SAM bins than in cotyledon and
  hypocotyl bins, but it is zero in many bins at every distance, so the
  curve stays flat overall. Photosynthesis genes follow light and
  chloroplasts rather than the developmental axis. RBCS1A, another
  photosynthesis gene, behaves the same way.

> **Honest note:** distance explains less than 20% of the variation for
> every gene. Look at the row of points along zero at the bottom of each
> plot: in 8 µm bins most genes are not detected in most bins. The GAM
> treats those zeros as real measurements, which pulls the curves down
> and weakens the fit. The patterns are real, but they are trends, not
> tight curves. Modelling whether a gene is detected at all (a binomial
> GAM), or using 16 µm bins, would handle this better.

------------------------------------------------------------------------

# Summary

## This run compared with the original

| What | Original | This run |
|----|----|----|
| Genes × bins loaded | 31,109 × 7,356 | 31,109 × 7,356 |
| Genes × bins after QC | 17,852 × 1,516 | 17,852 × 1,516 |
| Clusters | 7 | 7 |
| Top spatially variable gene (Moran’s I) | EXT3 (0.791) | EXT3 (0.791) |
| AGP31 rank among spatially variable genes | 2nd (from 35 genes only) | 5th (from all 3,000) |
| Cotyledon vs germinating seed: genes tested | 330 (reported as “DEGs”) | 331 |
| … significant with \|log2FC\| \> 1 | not reported | 50 |
| GO terms for hypocotyl (ORA) | 5 (from 2 genes) | 0 |
| GSEA terms after simplifying | 8 | 6 |
| ANT gradient p-value | 0.00026 | 0.00026 |
| CAB1 gradient p-value | 0.135 | 0.14 |

Everything up to clustering is identical. The later differences come
from two sources: bugs in the original scripts (explained below), and
small changes between software versions.

## Bugs fixed in the original scripts

| Problem in the original | Effect | Fix used here |
|----|----|----|
| `06_SVG_detection.R` contains pasted console output (line 41) | The script does not run (syntax error) | Removed |
| `06` uses `top_svg` without ever creating it | Error | Defined from the sorted Moran’s I table |
| `06` writes to a hard-coded path on the author’s computing cluster | Error on any other computer | Relative path `results/` |
| `05` re-runs `ScaleData()` on only the ~35 heat-map genes | Moran’s I can then only test those 35 genes | Scale variable genes + heat-map genes |
| The markers file stores cluster numbers, but `08` filters by tissue names | Tissue-specific ORA cut-offs silently do nothing | Add a `tissue` column and filter on it |
| `09` assumes k-means “piece 1” is the shoot | Breaks when k-means numbers the pieces differently | Pick the piece with the most SAM bins |
| Volcano plots use `-log10(0)` | Points with p = 0 vanish | Add `1e-300` before taking the log |
| `FindMarkers` row counts reported as “DEGs” | Overstates the number of differentially expressed genes | Report tested, significant, and significant + two-fold separately |

## Problems hit along the way, and how they were solved

| Symptom | Cause | Fix |
|----|----|----|
| `cannot open file` / `No such file or directory` | Commands were run from the wrong folder | Open the `.Rproj` file so the working directory is the project folder |
| Package installation stopped with a timeout | R’s 60-second download limit | `options(timeout = 3600)`, and install before downloading the data |
| `Cannot find the following identities in the object: Root` | The clusters had not been renamed yet | Run the renaming chunk; check with `levels(Idents(obj))` |
| `object 'tissue' not found` in Step 10 | `markers` was made before renaming | Add `markers$tissue` right after `RenameIdents()` |
| `unable to find an inherited method for function 'select'` | clusterProfiler’s dependencies mask dplyr’s `select()` | Write `dplyr::select()` |
| Knitting reinstalled packages and re-downloaded 768 MB | The setup chunks ran on every knit | `eval = FALSE` on those chunks; run them once by hand |

## Limitations

- **One section of one sample.** Bins are not biological replicates, so
  the statistics describe this slide, not seedlings in general.
- **Bins are not cells.** A bin on a tissue border can mix cell types.
- **The tissue names rest on a handful of marker genes.** “hypocotyl” is
  the weakest label: it has few specific markers and no enriched GO
  terms.
- **GO analysis loses genes.** Only about 42% of detected genes can be
  matched by gene symbol.

## Ideas for next steps

All four are tried in the follow-up notebook,
[follow_up_analyses.md](follow_up_analyses.md).

- Repeat the analysis with 16 µm bins (`bin.size = 16`). Are the same
  tissues found, with less noise?
- Try `SCTransform()` for normalisation, or BANKSY for clustering that
  also uses position.
- Run GSEA on **all** genes, and match genes by TAIR locus ID to keep
  more of them.
- Fit the gradient along the **root** (meristematic tissue → elongation
  zone → root) instead of the shoot.

# Session information

The exact software versions used for this run:

``` r
sessionInfo()
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: aarch64-apple-darwin23
    ## Running under: macOS Golden Gate 27.2
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_GB/en_GB/en_GB/C/en_GB/en_GB
    ## 
    ## time zone: Africa/Lagos
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] mgcv_1.9-4             nlme_3.1-169           enrichplot_1.32.1     
    ##  [4] org.At.tair.db_3.22.0  AnnotationDbi_1.74.0   IRanges_2.46.0        
    ##  [7] S4Vectors_0.50.3       Biobase_2.72.0         BiocGenerics_0.58.1   
    ## [10] generics_0.1.4         clusterProfiler_4.20.3 patchwork_1.3.2       
    ## [13] dplyr_1.2.1            ggplot2_4.0.3          Seurat_5.6.0          
    ## [16] SeuratObject_5.4.0     sp_2.2-4              
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RcppAnnoy_0.0.23        splines_4.6.1           later_1.4.8            
    ##   [4] ggplotify_0.1.4         tibble_3.3.1            polyclip_1.10-7        
    ##   [7] enrichit_0.2.7          fastDummies_1.7.6       lifecycle_1.0.5        
    ##  [10] httr2_1.3.0             globals_0.20.0          processx_3.9.0         
    ##  [13] Rnanoflann_0.0.3        lattice_0.22-9          hdf5r_1.3.16           
    ##  [16] MASS_7.3-65             magrittr_2.0.5          plotly_4.12.1          
    ##  [19] rmarkdown_2.32          yaml_2.3.12             ggtangle_0.1.5         
    ##  [22] httpuv_1.6.17           otel_0.2.0              sctransform_0.4.3      
    ##  [25] spam_2.11-4             spatstat.sparse_3.2-0   reticulate_1.47.0      
    ##  [28] cowplot_1.2.0           pbapply_1.7-5           DBI_1.3.0              
    ##  [31] RColorBrewer_1.1-3      abind_1.4-8             Rtsne_0.17             
    ##  [34] purrr_1.2.2             presto_1.1.0            yulab.utils_0.2.5      
    ##  [37] tweenr_2.0.3            rappdirs_0.3.4          aisdk_1.4.12           
    ##  [40] gdtools_0.5.2           ggrepel_0.9.8           irlba_2.4.1            
    ##  [43] listenv_1.1.0           spatstat.utils_3.2-5    tidytree_0.4.8         
    ##  [46] goftest_1.2-3           RSpectra_0.16-2         spatstat.random_3.5-2  
    ##  [49] fitdistrplus_1.2-6      parallelly_1.48.0       codetools_0.2-20       
    ##  [52] ggforce_0.5.0           DOSE_4.6.0              tidyselect_1.2.1       
    ##  [55] aplot_0.3.2             farver_2.1.2            matrixStats_1.5.0      
    ##  [58] spatstat.explore_3.8-3  Seqinfo_1.2.0           jsonlite_2.0.0         
    ##  [61] progressr_1.0.0         ggridges_0.5.7          survival_3.8-6         
    ##  [64] systemfonts_1.3.2       ggnewscale_0.5.2        tools_4.6.1            
    ##  [67] treeio_1.36.1           ica_1.0-3               Rcpp_1.1.2             
    ##  [70] glue_1.8.1              gridExtra_2.3.1         xfun_0.62              
    ##  [73] qvalue_2.44.0           withr_3.0.3             fastmap_1.2.0          
    ##  [76] callr_3.8.0             digest_0.6.39           gridGraphics_0.5-1     
    ##  [79] R6_2.6.1                mime_0.13               scattermore_1.2        
    ##  [82] GO.db_3.23.1            tensor_1.5.1            dichromat_2.0-1        
    ##  [85] spatstat.data_3.1-9     RSQLite_3.53.3          Rfast2_0.1.5.6         
    ##  [88] tidyr_1.3.2             fontLiberation_0.1.0    data.table_1.18.6.1    
    ##  [91] httr_1.4.9              htmlwidgets_1.6.4       scatterpie_0.2.6       
    ##  [94] uwot_0.2.5              pkgconfig_2.0.3         gtable_0.3.6           
    ##  [97] blob_1.3.0              lmtest_0.9-40           S7_0.2.2               
    ## [100] XVector_0.52.0          htmltools_0.5.9         fontBitstreamVera_0.1.1
    ## [103] dotCall64_1.2           zigg_0.0.2              scales_1.4.0           
    ## [106] png_0.1-9               spatstat.univar_3.2-0   ggfun_0.2.2            
    ## [109] knitr_1.52              reshape2_1.4.5          zoo_1.9-1              
    ## [112] cachem_1.1.0            stringr_1.6.0           KernSmooth_2.23-26     
    ## [115] parallel_4.6.1          miniUI_0.1.2            vipor_0.4.7            
    ## [118] arrow_26.0.0            ggrastr_1.0.2           pillar_1.11.1          
    ## [121] grid_4.6.1              vctrs_0.7.3             RANN_2.6.3             
    ## [124] promises_1.5.0          tidydr_0.0.7            xtable_1.8-8           
    ## [127] cluster_2.1.8.2         beeswarm_0.4.0          evaluate_1.0.5         
    ## [130] cli_3.6.6               compiler_4.6.1          rlang_1.3.0            
    ## [133] crayon_1.5.3            future.apply_1.21.0     labeling_0.4.3         
    ## [136] plyr_1.8.9              fs_2.1.0                ggbeeswarm_0.7.3       
    ## [139] ggiraph_0.9.6           stringi_1.8.9           viridisLite_0.4.3      
    ## [142] deldir_2.0-4            assertthat_0.2.1        Biostrings_2.80.2      
    ## [145] lazyeval_0.2.3          spatstat.geom_3.8-3     fontquiver_0.2.1       
    ## [148] GOSemSim_2.38.3         Matrix_1.7-5            RcppHNSW_0.7.0         
    ## [151] bit64_4.8.6             future_1.80.0           KEGGREST_1.52.2        
    ## [154] shiny_1.14.0            ROCR_1.0-12             Rfast_2.1.5.2          
    ## [157] igraph_2.3.4            memoise_2.0.1           RcppParallel_6.2.1     
    ## [160] ggtree_4.2.1            bit_4.6.0               gson_0.2.2             
    ## [163] ape_5.8-1

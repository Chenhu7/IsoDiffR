# IsoDiffR

IsoDiffR is a computational tool specifically designed to identify RNA isoforms whose expression patterns differ from those of their corresponding genes or major isoforms across distinct cell types, using single-cell RNA sequencing (scRNA-seq) data. It supports both pairwise and multi-cell-type comparisons, offering flexibility for diverse experimental designs.

The analysis of RNA isoforms in long-read single-cell RNA-seq data represents a cutting-edge direction in transcriptomics, providing insights that transcend conventional gene-level analyses. Despite the promise of this approach, dedicated computational methods remain scarce, underscoring an urgent need for methodological innovation in this rapidly evolving field.

To address this gap, we developed IsoDiffR as a versatile and accessible tool that facilitates isoform-level differential expression analysis. By enabling more nuanced exploration of transcriptional diversity and cell-type-specific splicing regulation, IsoDiffR contributes to a deeper understanding of gene expression complexity at single-cell resolution.

## Installation

IsoDiffR can be installed via this command

```
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}
devtools::install_github("Eveqian98/IsoDiffR", build_vignettes = TRUE)
```

or 

```
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
BiocManager::install(version='devel')
BiocManager::install("DiffIsoR")
```

## Data Preparation

IsoDiffR is designed to work with the following required inputs:

- A Seurat object that contains:
  - Isoform-level count matrix, stored in a dedicated assay (e.g., `"isoform"`).
  - Cell type annotations, stored as a column in the `meta.data` slot.
- A GTF file that includes:
  - Transcript IDs (`transcript_id`)
  - Corresponding gene IDs (`gene_id`)

Both the annotated Seurat object and the properly formatted GTF file are essential for running IsoDiffR and cannot be omitted.

The `meta.data` matrix of the Seurat object must include a column with cell type annotations.

```
> unique(isoform.0h@meta.data$cluster)
[1] "Corneal epithelial cells" "Basal cells" "Conjunctival cells" "Limbal stem cells"  
```

IsoDiffR expects the GTF file to include the mapping between isoforms and their corresponding genes, with the attributes named "gene_id" and "transcript_id" respectively.

```
> gtf[sample(nrow(gtf), 5), ]
        seqnames PacBio type     start       end strand           transcript_id width gene_id
1343338        7 PacBio exon 137591288 137591387      +  PB.58172.31(DCAF4~ISM)    99   DCAF4
1441262        9 PacBio exon  56276808  56276896      -  PB.62348.202(PPIF~NIC)    88    PPIF
640775        19 PacBio exon  49142503  49142605      +  PB.25747.221(EMP3~NNC)   102    EMP3
538191        17 PacBio exon  28938710  28938879      - PB.20816.480(KPNA3~ISM)   169   KPNA3
1318949        7 PacBio exon  94837377  94837458      - PB.55305.179(STRN3~NNC)    81   STRN3
```

## Exmaple Usage

### getDEIso

The function `getDEIso()` returns the identified differentially expressed isoforms (DEIs) based on the input Seurat object and GTF file.

```
DEIs = getDEIso(seurat_obj = isoform.0h,gtf = gtf,subset_ident = unique(isoform.0h@meta.data$cluster),cluster_column = "cluster")
```

The returned result includes four data frames.

- `multi.gene` DEIs compared with their corresponding genes across multiple cell types
- `multi.major` DEIs compared with their corresponding major isoforms across multiple cell types
- `pair.gene` DEIs compared with their corresponding genes between two cell types
- `pair.major` DEIs compared with their corresponding major isoforms between two cell types

The returned data frames contain various pieces of information In the `multi.gene` and `multi.major` data frames

- `isoform`  Transcript identifiers to each DEI
- `gene_id` The gene identifier corresponding to each DEI
- `major_isoform`  The major isoform corresponding to each DEI.
- `ratio_range` The range of the isoform expression relative to the corresponding gene expression across clusters.
- `adj.r.squared`  Adjusted R-squared value of the DEI.

```
> DEIs$multi.gene[sample(nrow(DEIs$multi.gene), 5), ]
                       isoform  gene_id              major_isoform ratio_range adj.r.squared
963      PB.17036.10(KLF4~NNC)     KLF4       PB.17036.3(KLF4~NNC)       0.093        -0.053
680 PB.23661.186(C18orf32~NNC) C18orf32 PB.23661.186(C18orf32~NNC)       0.188        -0.349
996     PB.11125.216(RTN4~FSM)     RTN4     PB.11125.158(RTN4~NNC)       0.213        -0.259
545     PB.23991.6(MYL12A~FSM)   MYL12A   PB.23991.140(MYL12A~FSM)       0.075        -0.414
730      PB.10371.19(RHOB~FSM)     RHOB      PB.10371.19(RHOB~FSM)       0.136        -0.292
```

The `pair.gene` and `pair.major` data frames return several distinct items.

- `condition`  The two cell types between which the DEI is identified
- `ratio.1`/`ratio.2` The proportion of the DEI's expression relative to the gene's total expression in each of the two cell types
- `ratio_difference`  The difference in the proportion of the DEI's expression relative to the gene's total expression between the two cell types
- `r×adj.R²` The product of the DEI's Pearson correlation coefficient and cosine similarity

```
> DEIs$pair.gene[sample(nrow(DEIs$pair.gene), 5), 4:8]
# A tibble: 5 × 5
  condition                            ratio.1 ratio.2 ratio_difference `r×adj.R²`
  <chr>                                  <dbl>   <dbl>            <dbl>      <dbl>
1 Basal cells_Conjunctival cells         0.401   0.232            0.169     -0.964
2 Conjunctival cells_Limbal stem cells   0.604   0.71             0.106     -0.997
3 Conjunctival cells_Limbal stem cells   0.519   0.639            0.12      -0.995
4 Basal cells_Conjunctival cells         0.373   0.505            0.132     -0.989
5 Basal cells_Limbal stem cells          0.42    0.213            0.207     -0.945
```

### plotDEIso

The function `plotDEIso()` generates a line plot showing the expression pattern of a specific isoform, along with its corresponding gene and major isoform, across specified cell types.

```
p <- plotDEIso(seurat_obj = isoform.0h,gtf = gtf,subset_ident = unique(isoform.0h@meta.data$cluster),cluster_column = "cluster",transcript_id = "PB.17036.10 (KLF4~NNC)")
```
 ![Image text](https://github.com/Eveqian98/IsoDiffR/blob/master/image/p1.png)

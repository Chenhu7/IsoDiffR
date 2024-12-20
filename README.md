### IsoDiffR: a tool designed to identify RNA isoforms with expression patterns distinct from their corresponding genes or major isoforms across different cell types in single-cell RNA-seq data, enabling both pairwise and multi-cell-type comparisons.

##### The analysis of long-read single-cell RNA sequencing(scRNA-seq) data at the RNA isoform level represents a frontier in current research, offering a transformative perspective beyond traditional gene expression analyses. Despite its potential, analytical methods specifically designed for this technology remain limited, highlighting the critical need for further methodological advancements to keep pace with its rapid development.

##### In response, we developed IsoDiffR,provides a significant resource for advancing our understanding of the complexity of single-cell gene expression and the diverse mechanisms underlying its regulation.
# Installation

##### IsoDiffR can be installed via this command

```
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}
devtools::install_github("Eveqian98/IsoDiffR", build_vignettes = TRUE)
```

##### or

```
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
BiocManager::install(version='devel')
BiocManager::install("DiffIsoR")
```


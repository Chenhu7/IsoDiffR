#' Differential Isoform Expression Patterns Across Cell Types in Single-Cell RNA-seq Data
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' This object should include isoform expression data, and ideally, have cell type annotations or clustering information.
#' The function will use the expression patterns across cells to identify differentially expressed isoforms
#' based on their variation across different cell types or clusters.
#' @param cluster_column A character string specifying the name of the column in `seurat_obj@meta.data`
#' that contains cell type annotations or clustering information. This column is used to group cells
#' based on their cell type or cluster for differential isoform expression analysis.
#' @param subset_ident A character vector specifying the cell types or clusters of interest to be extracted
#' from the `cluster_column` in `seurat_obj@meta.data`. Only the cells corresponding to the specified cell types
#' or clusters will be included in the differential isoform expression analysis.
#' @param gtf A GTF file containing the isoform IDs and corresponding gene IDs.
#' This file is used to map the isoforms in the single-cell RNA-seq data to their respective genes
#' for differential expression analysis of isoforms across cell types or clusters.
#' @param min.pct A numeric value specifying the minimum percentage of cells in a cluster
#' that must express an isoform for it to be considered in the differential expression analysis.
#' This parameter is passed to the `min.pct` argument in Seurat's `FindAllMarkers` function.
#' Isoforms with expression below this threshold across cell types will be filtered out.
#' The default value is 0.25 (i.e., isoforms expressed in at least 25% of cells in a cluster).
#' @param rate_cutoff A numeric value specifying the minimum ratio of isoform expression to gene expression
#' required for an isoform to be considered in the differential expression analysis.
#' Isoforms with a ratio (rate) below this threshold will be filtered out.
#' The default value is 0.1, meaning isoforms must account for at least 10% of the total gene expression to be included.
#' @param multi_cutoff A numeric value specifying the maximum adjusted R-squared value allowed for an isoform
#' to be considered as showing differential expression patterns across more than two cell types.
#' The linear regression function used is adjusted R-squared (adj.R.squared).
#' Isoforms with an adjusted R-squared value greater than this threshold will be filtered out.
#' The default value is 0.2, meaning only isoforms with an adjusted R-squared value less than 0.2,
#' indicating weak or less significant expression pattern differences across multiple cell types, will be retained.
#' @param ratio_range A numeric value specifying the minimum change in the isoform's expression as a proportion of the total gene expression
#' across different cell types. This parameter filters isoforms based on the variability in their expression ratio across multiple cell types.
#' The default value is 0.05, meaning isoforms must show at least a 5% difference in their expression ratio across cell types to be considered.
#' @param pair_cutoff A numeric vector of length 2 specifying the range for the product of Pearson correlation coefficient
#' and cosine similarity used to filter isoforms showing differential expression patterns between two cell types.
#' Isoforms with a product of these two values outside of this range will be excluded.
#' The default value is c(-1, 0.9), meaning only isoforms with a product of Pearson correlation coefficient and cosine similarity
#' greater than -1 and less than 0.9 will be retained.
#' @param ratio_difference A numeric value specifying the minimum difference in the isoform's expression ratio as a proportion of the total gene expression
#' between two cell types. This parameter filters isoforms based on the absolute difference in their expression ratios across the two cell types.
#' Isoforms with a difference smaller than this threshold will be excluded. The default value is 0.1, meaning isoforms must show at least a 10% difference
#' in their expression ratios between the two cell types to be considered.
#'
#' @return A list containing four data frames:
#' 1. **multi.gene**: A data frame with isoforms showing differential expression patterns across more than two cell types, relative to gene expression.
#' 2. **multi.major**: A data frame with isoforms showing differential expression patterns across more than two cell types, relative to major isoform expression.
#' 3. **pair.gene**: A data frame with isoforms showing differential expression patterns between two cell types, relative to gene expression.
#' 4. **pair.major**: A data frame with isoforms showing differential expression patterns between two cell types, relative to major isoform expression.
#'
#' Each data frame contains isoform IDs, gene IDs, and other statistical metrics used for filtering (e.g., correlation values, expression ratios).
#' Only isoforms that meet the specified thresholds for differential expression patterns across the relevant cell types will be included in the respective data frames.
#'
#' @export
getDEIso <- function(seurat_obj,cluster_column,subset_ident,gtf,
                     min.pct = 0.25,rate_cutoff = 0.1,
                     multi_cutoff = 0.2,ratio_range = 0.05,
                     pair_cutoff = c(-1,0.9),ratio_difference = 0.1){

  if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
  library(Seurat)
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  library(dplyr)
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
  library(tidyr)

  # ident using metadata column
  if (!cluster_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: The column", cluster_column, "does not exist in seurat_obj@meta.data."))
  }
  Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_column]]

  # reset ident metadata_column using subset_ident
  if (!cluster_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: The column", cluster_column, "does not exist in the meta.data of the Seurat object."))
  }
  seurat_obj <- subset(seurat_obj , idents = subset_ident)

  # Set the default assay to RNA
  DefaultAssay(seurat_obj) = "RNA"

  # Check which isoforms are present in Seurat object
  existing_genes = row.names(seurat_obj@assays$RNA)
  gtf_exist = gtf[gtf$transcript_id %in% existing_genes,]
  gtf_trans = gtf_exist[gtf_exist$type == "transcript",]

  #Find marker isoforms as candidates for diff isoforms
  # Filter the percentage of isoform expression in clusters
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = min.pct,logfc.threshold = 0)
  markers <- data.frame(gene = unique(markers$gene))
  colnames(markers)[1] <- "isoform"

  # Match isoforms with their corresponding genes
  match_rows <- match(markers$isoform,gtf_trans$transcript_id)
  markers$gene_id <- gtf_trans$gene_id[match_rows]
  markers <- na.omit(markers)

  # Remove genes with only one isoform
  markers$relevant <- NA
  for (i in 1:nrow(markers)) {
    gene_id <- markers$gene_id[i]
    relevant <- length(gtf_trans[gtf_trans$gene_id == gene_id, "transcript_id"])
    markers$relevant[i] <- relevant
  }
  markers <- subset(markers,markers[,3] != 1)
  markers <- markers[,-3]

  # Match isoforms with their corresponding major isoforms
  n_iter <- nrow(markers)
  cat("\033[31mMatch isoforms with their corresponding major isoforms\033[0m\n")
  pb <- txtProgressBar(min = 0,max = n_iter, style = 3,width = 60,char = "+")
  for(i in 1:n_iter) {
    raw_counts_isoform <- GetAssayData(seurat_obj,assay = "RNA",slot <-"counts")
    gene_id <- markers$gene_id[i]
    relevant_transcripts <- unlist(na.omit(as.data.frame(gtf_trans[gtf_trans$gene_id == gene_id, "transcript_id"])))
    total_expressions <- apply(raw_counts_isoform[relevant_transcripts, ], 1, sum)
    major_isoform <- names(which.max(total_expressions))
    markers$major_isoform[i] <- major_isoform
    setTxtProgressBar(pb, i)
  }
  ## Compute gene major isoform and marker isoform counts (in CPM) for each gene across idents
  idents <- levels(seurat_obj@active.ident)
  total_iterations <- length(idents) * nrow(markers)
  pb <- txtProgressBar(min = 0, max = total_iterations, style = 3, width = 60, char = "+")
  iteration_counter <- 0

  for (ident in idents) {
    subset_data <- subset(seurat_obj, idents = ident)
    raw_counts_subset <- subset_data@assays$RNA$counts
    total_counts <- sum(raw_counts_subset)

    isoform_cpm_column <- numeric(nrow(markers))
    major_cpm_column <- numeric(nrow(markers))
    gene_cpm_column <- numeric(nrow(markers))
    ratio_column <- numeric(nrow(markers))

    for (i in 1:nrow(markers)) {
      isoform_name <- markers[i, "isoform"]
      isoform_cpm_column[i] <- round((sum(raw_counts_subset[isoform_name,]) / total_counts) * 1e6,3)

      gene_id <- markers[i, "gene_id"]
      gene_of_isoform <- unlist(na.omit(as.data.frame(gtf_trans$transcript_id[gtf_trans$gene_id == gene_id])))
      gene_cpm_column[i] <- round((sum(raw_counts_subset[gene_of_isoform,]) / total_counts) * 1e6,3)

      major_name <- markers[i,"major_isoform"]
      major_cpm_column[i] <- round((sum(raw_counts_subset[major_name,]) / total_counts) * 1e6,3)

      ratio_column[i] <- round(isoform_cpm_column[i] / gene_cpm_column[i],3)

      iteration_counter <- iteration_counter + 1
      setTxtProgressBar(pb, iteration_counter)

      if (iteration_counter %% nrow(markers) == 0) {
        cat(sprintf("\033[31mProcessing ident: %s (%d/%d)\033[0m\n", ident, iteration_counter / nrow(markers), length(idents)))
        flush.console()
      }
    }

    markers[, paste(ident, "CPM", sep = ".")] <- isoform_cpm_column
    markers[, paste(ident, "gene_CPM", sep = ".")] <- gene_cpm_column
    markers[, paste(ident, "major_CPM", sep = ".")] <- major_cpm_column
    markers[,paste(ident,"ratio",sep = ".")] <- ratio_column
  }


  # Remove isoforms with expression ratios below the rate_cutoff
  ratio_vector <- nrow(markers)
  for (i in 1:nrow(markers)) {
    ratio_max <-  max(markers[i,grep("ratio",names(markers))], na.rm = TRUE)
    ratio_vector[i] <- ratio_max
  }
  markers$ratio_max <- ratio_vector

  markers <- markers[markers$ratio_max > rate_cutoff,]

  markers <- markers[,-ncol(markers)]

  #The expression patterns of isoforms that show significant differences in multiple clusters compared to genes or major isoforms
  multi <- markers
  # Remove isoforms with expression ratio ranges below the ratio_range
  range_vector <- numeric(nrow(multi))
  for (i in 1:nrow(multi)){
    ratio_vector <- as.numeric(multi[i,grepl(".ratio$",names(multi))])
    range_value <- round(max(ratio_vector, na.rm = TRUE) - min(ratio_vector, na.rm = TRUE), 3)
    range_vector[i] <- range_value
  }
  multi$ratio_range <- range_vector

  multi <- multi[multi$ratio_range > ratio_range,]

  # Calculate the adj.R.squared between isoform and major isoform / gene expression vectors
  adj.r.squared.vector <- numeric(nrow(multi))
  adj.r.squared_major.vector <- numeric(nrow(multi))
  for (i in 1:nrow(multi)) {
    isoform_cpm_values <- multi[i, grepl(".CPM$", names(multi)) &
                                  !grepl(".gene_CPM$", names(multi))&
                                  !grepl(".major_CPM$", names(multi))]
    gene_cpm_values <- multi[i, grepl(".gene_CPM$", names(multi))]
    major_cpm_values <- multi[i, grepl(".major_CPM$", names(multi))]

    isoform_cpm_values = as.numeric(isoform_cpm_values)
    gene_cpm_values = as.numeric(gene_cpm_values)
    major_cpm_values = as.numeric(major_cpm_values)

    model <- lm(isoform_cpm_values ~ gene_cpm_values)
    summary_model <- summary(model)
    adj_r_squared <- summary_model$adj.r.squared

    model_major <- lm(isoform_cpm_values ~ major_cpm_values)
    summary_model_major <- summary(model_major)
    adj_r_squared_major <- summary_model_major$adj.r.squared

    adj.r.squared.vector[i] <- round(adj_r_squared,3)
    adj.r.squared_major.vector[i] <- round(adj_r_squared_major,3)
  }
  multi$adj.r.squared <- adj.r.squared.vector
  multi$adj.r.squared_major <- adj.r.squared_major.vector

  # Filter isoforms with an adjusted R-squared value less than the multi_cutoff in relation to gene and major isoform
  multi_res <- subset(multi,multi$adj.r.squared < multi_cutoff |
                        multi$adj.r.squared_major < multi_cutoff)


  # Remove major isoforms from multi.major where the major isoform itself is a diff isoform

  multi.gene = multi_res[multi_res$adj.r.squared < multi_cutoff,]
  multi.major = multi_res[multi_res$adj.r.squared_major < multi_cutoff,]

  multi.major$dima = ifelse(multi.major$major_isoform %in% multi.gene$isoform,"yes","no")
  multi.major = multi.major[multi.major$dima == "no",]
  multi.major = multi.major[,-ncol(multi.major)]

  #The expression patterns of isoforms that show significant differences in two clusters compared to genes or major isoforms
  # Determine the expression differences of isoforms between pairs of clusters using cosine similarity
  pair.gene <- markers
  pair.major <- markers
  # Pairwise combine the clusters
  levels_x <- levels(Idents(seurat_obj))
  combinations <- combn(levels_x, 2)

  combinations_names <- apply(combinations, 2, function(combo) {
    paste(combo, collapse = "_")
  })

  for (name in combinations_names) {
    pair.gene[[name]] <- NA
  }

  # Calculate the cosine similarity between each pair of groups
  start_col <- which(names(pair.gene) == combinations_names[1])
  end_col <- which(names(pair.gene) == combinations_names[length(combinations_names)])

  for (i in start_col:end_col) {
    col = colnames(pair.gene)[i]
    a =  strsplit(col, "_")[[1]][1]
    b =  strsplit(col, "_")[[1]][2]

    cosine_similarity_values <- numeric(nrow(pair.gene))
    for (j in 1:nrow(pair.gene)) {
      isoform_1 <- pair.gene[j,paste(a, ".CPM", sep = "")]
      isoform_2 <- pair.gene[j,paste(b, ".CPM", sep = "")]
      x = c(isoform_1,isoform_2)

      gene_1 <- pair.gene[j,paste(a, ".gene_CPM", sep = "")]
      gene_2 <- pair.gene[j,paste(b, ".gene_CPM", sep = "")]
      y = c(gene_1,gene_2)

      dot_product <- sum(x * y)
      norm_x <- sqrt(sum(x^2))
      norm_y <- sqrt(sum(y^2))
      cosine_similarity <- dot_product / (norm_x * norm_y)
      correlation <- cor(x,y)
      cosine_similarity_values[j] <- round(sign(correlation) * cosine_similarity,3)

    }
    pair.gene[,col] <- cosine_similarity_values
  }

  negative <- pair_cutoff[1]
  positive <- pair_cutoff[2]

  pair.gene <- pair.gene%>%
    rowwise() %>%
    mutate(condition = list(names(pair.gene)[start_col:end_col][
      (c_across(colnames(pair.gene)[start_col]:colnames(pair.gene)[end_col]) > negative & c_across(colnames(pair.gene)[start_col]:colnames(pair.gene)[end_col]) < 0) |
        (c_across(colnames(pair.gene)[start_col]:colnames(pair.gene)[end_col]) > 0 & c_across(colnames(pair.gene)[start_col]:colnames(pair.gene)[end_col]) < positive)
    ])) %>%
    ungroup() %>%
    filter(lengths(condition) > 0) %>%
    unnest(condition)

  # ratio difference
  pair.gene <- na.omit(pair.gene)

  ratio_vector.1 <- numeric(nrow(pair.gene))
  ratio_vector.2 <- numeric(nrow(pair.gene))
  for (i in 1:nrow(pair.gene)) {
    condition <- as.character(pair.gene[i,"condition"])
    cluster.1 <- paste0(sub("_.*", "", condition),".ratio")
    cluster.2 <- paste0(sub(".*_", "", condition),".ratio")
    ratio_vector.1[i] <- as.numeric(pair.gene[i,cluster.1])
    ratio_vector.2[i] <- as.numeric(pair.gene[i,cluster.2])
  }
  pair.gene$ratio.1 <- ratio_vector.1
  pair.gene$ratio.2 <- ratio_vector.2

  pair.gene$ratio_difference <- abs(pair.gene$ratio.1 - pair.gene$ratio.2)

  pair.gene <- pair.gene[pair.gene$ratio_difference > ratio_difference,]


  cs_vector <- nrow(pair.gene)
  for (i in 1:nrow(pair.gene)) {
    condition <- as.character(pair.gene[i,"condition"])
    cs <- as.numeric(pair.gene[i,condition])
    cs_vector[i] <- cs
  }
  pair.gene$cosine_similarity <- cs_vector

  # Calculate the cosine similarity with the major isoform ####
  combinations_names <- apply(combinations, 2, function(combo) {
    paste(combo, collapse = "_major_")
  })

  for (name in combinations_names) {
    pair.major[[name]] <- NA
  }

  start_col <- which(names(pair.major) == combinations_names[1])
  end_col <- which(names(pair.major) == combinations_names[length(combinations_names)])

  for (i in start_col:end_col) {
    col <- colnames(pair.major)[i]
    a <- strsplit(col, "_")[[1]][1]
    b <- strsplit(col, "_")[[1]][3]

    cosine_similarity_values <- numeric(nrow(pair.major))
    for (j in 1:nrow(pair.major)) {
      isoform_1 <- pair.major[j,paste(a, ".CPM", sep = "")]
      isoform_2 <- pair.major[j,paste(b, ".CPM", sep = "")]
      x = c(isoform_1,isoform_2)

      gene_1 <- pair.major[j,paste(a, ".major_CPM", sep = "")]
      gene_2 <- pair.major[j,paste(b, ".major_CPM", sep = "")]
      y = c(gene_1,gene_2)

      dot_product <- sum(x * y)
      norm_x <- sqrt(sum(x^2))
      norm_y <- sqrt(sum(y^2))
      cosine_similarity <- dot_product / (norm_x * norm_y)

      correlation <- cor(x,y)
      cosine_similarity_values[j] <- round(sign(correlation) * cosine_similarity,3)
    }
    pair.major[,col] <- cosine_similarity_values
  }

  pair.major <- pair.major%>%
    rowwise() %>%
    mutate(condition = list(names(pair.major)[start_col:end_col][
      (c_across(colnames(pair.major)[start_col]:colnames(pair.major)[end_col]) > negative & c_across(colnames(pair.major)[start_col]:colnames(pair.major)[end_col]) < 0) |
        (c_across(colnames(pair.major)[start_col]:colnames(pair.major)[end_col]) > 0 & c_across(colnames(pair.major)[start_col]:colnames(pair.major)[end_col]) < positive)
    ])) %>%
    ungroup() %>%
    filter(lengths(condition) > 0) %>%
    unnest(condition)

  # ratio difference
  pair.major <- na.omit(pair.major)

  ratio_vector.1 <- numeric(nrow(pair.major))
  ratio_vector.2 <- numeric(nrow(pair.major))
  for (i in 1:nrow(pair.major)) {
    condition = as.character(pair.major[i,"condition"])
    cluster.1 = paste0(sub("_.*", "", condition),".ratio")
    cluster.2 <- paste0(sub(".*_", "", condition),".ratio")
    ratio_vector.1[i] = as.numeric(pair.major[i,cluster.1])
    ratio_vector.2[i] = as.numeric(pair.major[i,cluster.2])
  }
  pair.major$ratio.1 <- ratio_vector.1
  pair.major$ratio.2 <- ratio_vector.2

  pair.major$ratio_difference <- abs(pair.major$ratio.1 - pair.major$ratio.2)

  pair.major <- pair.major[pair.major$ratio_difference > ratio_difference,]


  cs_vector <- nrow(pair.major)
  for (i in 1:nrow(pair.major)) {
    condition <- as.character(pair.major[i,"condition"])
    cs <- as.numeric(pair.major[i,condition])
    cs_vector[i] <- cs
  }
  pair.major$cosine_similarity <- cs_vector

  # Remove major isoforms from pair.major where the major isoform itself is a diff isoform
  pair.gene$ isoform.condition <- paste(pair.gene$isoform,pair.gene$condition,sep = ".")
  pair.major $ major.condition <- paste(pair.major$major_isoform,pair.major$condition,sep = ".")
  pair.major$madi <- ifelse(pair.major$major.condition %in% pair.gene$isoform.condition,"yes","no")
  pair.major <- pair.major[pair.major$madi == "no",]
  pair.major <- pair.major[,-ncol(pair.major)]

  # Organization of the output results of diff isoforms
  multi.gene <- multi.gene[,c(1:3,ncol(multi.gene)-2,ncol(multi.gene)-1)]
  multi.major <- multi.major[,c(1:3,ncol(multi.major)-1,ncol(multi.major))]

  pair.gene <- pair.gene[,c(1:3,ncol(pair.gene)-5,ncol(pair.gene)-4,ncol(pair.gene)-3,ncol(pair.gene)-2,ncol(pair.gene)-1)]
  pair.major <- pair.major[,c(1:3,ncol(pair.major)-5,ncol(pair.major)-4,ncol(pair.major)-3,ncol(pair.major)-2,ncol(pair.major)-1)]

  colnames(pair.gene)[ncol(pair.gene)] <- "r×adj.R²"
  colnames(pair.major)[ncol(pair.major)] <- "r×adj.R²"

  DEIso <- setNames(list(multi.gene, multi.major, pair.gene, pair.major),
                    c("multi.gene", "multi.major", "pair.gene", "pair.major"))

  return(DEIso)
}

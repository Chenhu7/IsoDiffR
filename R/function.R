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
process_data <- function(seurat_obj, cluster_column, subset_ident, gtf,
                         min.pct = 0.25, rate_cutoff = 0.1){
  
  if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
  library(Seurat)
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  library(dplyr)
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
  library(tidyr)
  
  if (!cluster_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: The column", cluster_column, "does not exist in seurat_obj@meta.data."))
  }
  Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_column]]
  
  if (!cluster_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: The column", cluster_column, "does not exist in the meta.data of the Seurat object."))
  }
  seurat_obj <- subset(seurat_obj , idents = subset_ident)
  
  DefaultAssay(seurat_obj) = "RNA"
  
  existing_genes = row.names(seurat_obj@assays$RNA)
  gtf_exist = gtf[gtf$transcript_id %in% existing_genes,]
  gtf_trans = gtf_exist[gtf_exist$type == "transcript",]
  
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = min.pct,logfc.threshold = 0)
  markers <- data.frame(isoform = unique(markers$gene))
  cat("After FindAllMarkers ->", nrow(markers), "\n")
  
  match_rows <- match(markers$isoform,gtf_trans$transcript_id)
  markers$gene_id <- gtf_trans$gene_id[match_rows]
  markers <- na.omit(markers)
  cat("After match with GTF ->", sum(!is.na(match_rows)), "\n")
  
  markers$relevant <- sapply(markers$gene_id, function(gid) {
    length(unique(gtf_trans$transcript_id[gtf_trans$gene_id == gid]))
  })
  markers <- subset(markers, relevant != 1)
  cat("After remove single isoform genes ->", nrow(markers), "\n")
  markers <- markers[, !names(markers) %in% "relevant"]
  
  if (nrow(markers) == 0) {
    stop("No candidate isoforms found after filtering. Check FindAllMarkers and GTF matching.")
  }
  
  # Major isoform matching
  cat("\033[31mOptimized: Match isoforms with their corresponding major isoforms\033[0m\n")
  
  raw_counts_all <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  genes_to_process <- unique(markers$gene_id)
  relevant_gtf <- gtf_trans[gtf_trans$gene_id %in% genes_to_process, ]
  
  transcripts_to_sum <- relevant_gtf$transcript_id
  transcripts_to_sum <- intersect(transcripts_to_sum, rownames(raw_counts_all))
  
  total_expressions <- apply(raw_counts_all[transcripts_to_sum, , drop = FALSE], 1, sum)
  
  expr_df <- data.frame(
    transcript_id = names(total_expressions),
    total_expr = total_expressions
  )
  expr_df <- left_join(expr_df, relevant_gtf[, c("transcript_id", "gene_id")], by = "transcript_id")
  expr_df <- na.omit(expr_df)
  
  major_isoforms_map <- expr_df %>%
    group_by(gene_id) %>%
    slice_max(order_by = total_expr, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(gene_id, major_isoform = transcript_id)
  
  markers <- markers %>%
    left_join(major_isoforms_map, by = "gene_id") %>%
    na.omit() 
  
  cat("Major isoform matching complete. Markers count:", nrow(markers), "\n")
  
  
  # CPM normalization
  cat("\033[31mOptimized: Compute gene major isoform and marker isoform counts (in CPM)\033[0m\n")
  
  raw_counts_all <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  cell_idents <- Idents(seurat_obj)
  idents <- levels(cell_idents)
  results_list <- list()
  
  markers_for_join <- markers %>% 
    select(isoform, gene_id, major_isoform) %>% 
    mutate(row_id = row_number())
  
  all_relevant_transcripts <- relevant_gtf$transcript_id
  all_relevant_transcripts <- intersect(all_relevant_transcripts, rownames(raw_counts_all))
  
  if(length(all_relevant_transcripts) == 0) {
    warning("No relevant transcripts found in raw counts.")
  }
  
  n_idents = length(idents)
  pb <- txtProgressBar(min = 0, max = n_idents, style = 3, width = 60, char = "+")
  
  for (k in 1:n_idents) {
    ident = idents[k]
    cells_in_ident <- names(cell_idents[cell_idents == ident])
    raw_counts_subset <- raw_counts_all[, cells_in_ident, drop = FALSE]
    total_counts <- sum(raw_counts_subset)
    
    # Isoform and Major CPM
    transcripts_to_sum_in_cluster <- c(markers_for_join$isoform, markers_for_join$major_isoform)
    transcripts_to_sum_in_cluster <- unique(intersect(transcripts_to_sum_in_cluster, rownames(raw_counts_subset)))
    
    iso_major_sums <- rowSums(raw_counts_subset[transcripts_to_sum_in_cluster, , drop = FALSE])
    
    # Map Isoform CPM
    isoform_sums <- iso_major_sums[markers_for_join$isoform]
    isoform_sums[is.na(isoform_sums)] <- 0
    isoform_cpm_column <- round((isoform_sums / total_counts) * 1e6, 3)
    
    # Map Major CPM
    major_sums <- iso_major_sums[markers_for_join$major_isoform]
    major_sums[is.na(major_sums)] <- 0
    major_cpm_column <- round((major_sums / total_counts) * 1e6, 3)
    
    # Gene CPM
    gene_cpm_column <- numeric(nrow(markers_for_join))
    
    if (length(all_relevant_transcripts) > 0) {
      all_relevant_subset_counts <- raw_counts_subset[all_relevant_transcripts, , drop = FALSE]
      counts_df <- data.frame(
        transcript_id = rownames(all_relevant_subset_counts),
        counts_sum = rowSums(all_relevant_subset_counts)
      )
      
      counts_df <- left_join(counts_df, relevant_gtf[, c("transcript_id", "gene_id")], by = "transcript_id")
      counts_df <- na.omit(counts_df) 
      
      gene_sums <- counts_df %>%
        group_by(gene_id) %>%
        summarise(gene_sum = sum(counts_sum, na.rm = TRUE), .groups = 'drop')
      
      gene_cpm_map <- setNames(gene_sums$gene_sum, gene_sums$gene_id)
      gene_sums_ordered <- gene_cpm_map[markers_for_join$gene_id]
      gene_sums_ordered[is.na(gene_sums_ordered)] <- 0 
      
      gene_cpm_column <- round((gene_sums_ordered / total_counts) * 1e6, 3)
    }
    
    ratio_column <- round(isoform_cpm_column / gene_cpm_column, 3)
    
    results_list[[ident]] <- data.frame(
      row_id = markers_for_join$row_id,
      ident = ident,
      isoform_cpm = isoform_cpm_column,
      gene_cpm = gene_cpm_column,
      major_cpm = major_cpm_column,
      ratio = ratio_column
    )
    
    setTxtProgressBar(pb, k)
    
    cat(sprintf("\n\033[31mProcessing ident: %s (%d/%d)\033[0m\n", ident, k, n_idents))
    flush.console()
  }
  
  # 5. Combine and transform results (Corrected names_glue)
  combined_results <- bind_rows(results_list)
  
  wide_results <- combined_results %>% 
    pivot_wider(
      names_from = ident, 
      values_from = c(isoform_cpm, gene_cpm, major_cpm, ratio),
      names_glue = "{ident}.{.value}"
    ) 
  
  new_names <- names(wide_results)
  new_names <- sub("\\.isoform_cpm$", ".CPM", new_names)
  new_names <- sub("\\.gene_cpm$", ".gene_CPM", new_names)
  new_names <- sub("\\.major_cpm$", ".major_CPM", new_names)
  
  names(wide_results) <- new_names
  
  # Final join back to markers
  markers <- markers %>%
    mutate(row_id = row_number()) %>% 
    left_join(wide_results, by = "row_id") %>%
    select(-row_id)  
  
  # Remove isoforms with expression ratios below the rate_cutoff
  ratio_vector <- numeric(nrow(markers))
  for (i in 1:nrow(markers)) {
    ratio_max <-  max(markers[i,grep("\\.ratio$", names(markers))], na.rm = TRUE) 
    ratio_vector[i] <- ratio_max
  }
  markers$ratio_max <- ratio_vector
  
  markers <- markers[markers$ratio_max > rate_cutoff,]
  
  markers <- markers[,-ncol(markers)]
  
  return(markers)
}


## identify multil.gene multil.major
identify_multi_DiffIso = function(markers,
                                  multi_cutoff = 0.2,ratio_range = 0.05){
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
  multi_res <- subset(multi,multi$adj.r.squared <= multi_cutoff |
                        multi$adj.r.squared_major <= multi_cutoff)
  
  
  # Remove major isoforms from multi.major where the major isoform itself is a diff isoform
  multi.gene = multi_res[multi_res$adj.r.squared <= multi_cutoff,]
  multi.major = multi_res[multi_res$adj.r.squared_major <= multi_cutoff,]
  
  # print(length(row.names(multi.major)))
  multi.major$dima = ifelse(multi.major$major_isoform %in% multi.gene$isoform,"yes","no")
  multi.major = multi.major[multi.major$dima == "no",]
  multi.major = multi.major[,-ncol(multi.major)]
  
  keep_cols = c("isoform", "gene_id", "major_isoform", "ratio_range", "adj.r.squared")
  multi.gene = multi.gene[, keep_cols]
  keep_cols = c("isoform", "gene_id", "major_isoform", "ratio_range", "adj.r.squared_major")
  multi.major = multi.major[, keep_cols]
  
  cat(paste('Identified multi DiffIso: ', length(row.names(multi.gene)),'(gene);',length(row.names(multi.major)),'(major)\n'))
  return(list(multi.gene = multi.gene, multi.major = multi.major))
}


## identify pair.gene pair.major
identify_pair_DiffIso = function(markers, compare_ident,
                                 pair_cutoff = c(-1,0.9),ratio_difference = 0.1){
  
  pair.gene <- markers
  pair.major <- markers
  combinations <- combn(compare_ident, 2)
  print(combinations)
  
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
      (c_across(colnames(pair.gene)[start_col]:colnames(pair.gene)[end_col]) > negative &
         c_across(colnames(pair.gene)[start_col]:colnames(pair.gene)[end_col]) < 0) |
        (c_across(colnames(pair.gene)[start_col]:colnames(pair.gene)[end_col]) > 0 &
           c_across(colnames(pair.gene)[start_col]:colnames(pair.gene)[end_col]) <= positive)
    ])) %>%
    ungroup() %>%
    # filter(lengths(condition) > 0) %>%
    filter(!is.na(condition), lengths(condition) > 0) %>%
    unnest(condition)
  
  # ratio difference
  pair.gene <- na.omit(pair.gene)
  
  if (nrow(pair.gene) > 0) {
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
  }else{
    message("[Warning] pair.gene is empty after ratio_difference filtering.")
    pair.gene <- markers[0,]
  }
  
  if (nrow(pair.gene) > 0) {
    cs_vector <- nrow(pair.gene)
    for (i in 1:nrow(pair.gene)) {
      condition <- as.character(pair.gene[i,"condition"])
      cs <- as.numeric(pair.gene[i,condition])
      cs_vector[i] <- cs
    }
    pair.gene$cosine_similarity <- cs_vector
    
  }else{
    message("[Warning] pair.gene is empty after ratio_difference filtering.")
    pair.gene <- markers[0,]
  }
  
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
        (c_across(colnames(pair.major)[start_col]:colnames(pair.major)[end_col]) > 0 & c_across(colnames(pair.major)[start_col]:colnames(pair.major)[end_col]) <= positive)
    ])) %>%
    ungroup() %>%
    # filter(lengths(condition) > 0) %>%
    filter(!is.na(condition), lengths(condition) > 0) %>%
    unnest(condition)
  
  # ratio difference
  pair.major <- na.omit(pair.major)
  
  if (nrow(pair.major) == 0) {
    message("[Warning] pair.major is empty after ratio_difference filtering.")
    pair.major <- markers[0,]
    return(list(pair.gene = pair.gene, pair.major = pair.major))
  }
  
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
  
  if (nrow(pair.major) == 0) {
    message("[Warning] pair.major is empty after ratio_difference filtering.")
    pair.major <- markers[0,]
    return(list(pair.gene = pair.gene, pair.major = pair.major))
  }
  
  cs_vector <- numeric(nrow(pair.major))
  for (i in seq_len(nrow(pair.major))) {
    condition <- as.character(pair.major[i, "condition"])
    
    # skipping NA
    if (is.na(condition) || condition == "" || !(condition %in% colnames(pair.major))) {
      cs_vector[i] <- NA
      next
    }
    
    cs <- suppressWarnings(as.numeric(pair.major[i, condition]))
    
    if (is.na(cs)) {
      cs_vector[i] <- NA
    } else {
      cs_vector[i] <- cs
    }
  }
  pair.major$cosine_similarity <- cs_vector
  
  #
  pair.major <- pair.major %>% 
    mutate(cosine_similarity = as.numeric(cosine_similarity)) %>%
    filter(cosine_similarity <= as.numeric(pair_cutoff[2]))
  
  pair.gene <- pair.gene %>% 
    mutate(cosine_similarity = as.numeric(cosine_similarity)) %>%
    filter(cosine_similarity <= as.numeric(pair_cutoff[2]))
  
  # Remove major isoforms from pair.major where the major isoform itself is a diff isoform
  pair.gene$isoform.condition <- paste(pair.gene$isoform,pair.gene$condition,sep = ".")
  pair.major$major.condition <- paste(pair.major$major_isoform,pair.major$condition,sep = ".")
  pair.major$madi <- ifelse(pair.major$major.condition %in% pair.gene$isoform.condition,"yes","no")
  pair.major <- pair.major[pair.major$madi == "no",]
  pair.major <- pair.major[,-ncol(pair.major)]
  
  
  keep_cols = c("isoform", "gene_id","major_isoform", "condition", "ratio.1", "ratio.2","cosine_similarity")
  pair.gene = pair.gene[, keep_cols]
  pair.major = pair.major[, keep_cols]
  
  names(pair.gene)[names(pair.gene) == "cosine_similarity"] = "r×adj.R²"
  names(pair.major)[names(pair.major) == "cosine_similarity"] = "r×adj.R²"
  
  cat(paste('Identified pair DiffIso: ', length(row.names(pair.gene)),'(gene);',length(row.names(pair.major)),'(major)\n'))
  return(list(pair.gene = pair.gene, pair.major = pair.major))
  
}

## output
output_DEIs <- function(multi_result,pair_result){
  multi.gene = multi_result$multi.gene
  multi.major = multi_result$multi.major
  pair.gene = pair_result$pair.gene
  pair.major = pair_result$pair.major
  
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


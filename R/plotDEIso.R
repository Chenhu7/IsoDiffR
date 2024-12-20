#' Plot Isoform, Major Isoform, and Gene Expression Across Cell Types
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' This object should include isoform expression data, and ideally, have cell type annotations or clustering information.
#' The function will use the expression patterns across cells to identify differentially expressed isoforms
#' based on their variation across different cell types or clusters.
#' @param gtf A GTF file containing the isoform IDs and corresponding gene IDs.
#' This file is used to map the isoforms in the single-cell RNA-seq data to their respective genes
#' for differential expression analysis of isoforms across cell types or clusters.
#' @param cluster_column A character string specifying the name of the column in `seurat_obj@meta.data`
#' that contains cell type annotations or clustering information. This column is used to group cells
#' based on their cell type or cluster for differential isoform expression analysis.
#' @param subset_ident  A character vector specifying the cell types or clusters of interest to be extracted
#' from the `cluster_column` in `seurat_obj@meta.data`. Only the cells corresponding to the specified cell types
#' or clusters will be included in the differential isoform expression analysis.
#' @param transcript_id A character string representing the isoform of interest.
#'   It must be a valid transcript ID present in the GTF file provided.
#'
#' @return A line plot displaying the expression levels of the input `transcript_id` isoform,
#'   its corresponding major isoform, and the gene across different cell types.
#' @export
plotDEIso <- function(seurat_obj,gtf,cluster_column,subset_ident,transcript_id){
  if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
  library(Seurat)
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  library(dplyr)
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("dplyr")
  library(ggplot2)

  # ident using metadata column
  Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_column]]

  # reset ident metadata_column using subset_ident
  seurat_obj <- subset(seurat_obj , idents = subset_ident)

  # Set the default assay to RNA
  DefaultAssay(seurat_obj) = "RNA"

  # Check which isoforms are present in Seurat object
  existing_genes = row.names(seurat_obj@assays$RNA)
  gtf_exist = gtf[gtf$transcript_id %in% existing_genes,]
  gtf_trans = gtf_exist[gtf_exist$type == "transcript",]

  # check if the isoform is major isoform
  isoform = as.data.frame(transcript_id)
  colnames(isoform)[1] = "isoform"

  # match gene_id
  match_rows = match(isoform$isoform,gtf_trans$transcript_id)

  if (any(is.na(match_rows))){
    stop("Error:Unable to match gene ID, please check your input")
  }

  isoform$gene_id = gtf_trans$gene_id[match_rows]

  # Obtaining expression levels through Seurat objects
  raw_counts_isoform = GetAssayData(seurat_obj,assay = "RNA",layer = "counts")

  # Check if this isoform is the only isoform of the corresponding gene
  gene_id = isoform$gene_id
  isoform$revelant = length(gtf_trans[gtf_trans$gene_id == gene_id , "transcript_id"])

  if (isoform[1,3] == 1) {
    stop("The isoform input is the only isoform of the corresponding gene")
  }

  # The corresponding major isoform for the input isoform
  relevant_transcripts = gtf_trans[gtf_trans$gene_id == gene_id , "transcript_id"]
  relevant_transcripts = unlist(na.omit(as.data.frame(relevant_transcripts)))
  total_expressions = apply(raw_counts_isoform[relevant_transcripts,],1,sum)
  isoform$major = names(which.max(total_expressions))

  # Check the expression levels of the isoform, the major isoform, and the gene separately
  df = as.data.frame(raw_counts_isoform[transcript_id,])
  colnames(df)[1] = isoform$isoform
  df$major = raw_counts_isoform[isoform$major,]
  colnames(df)[2] = isoform$major
  Gene <- as.data.frame(raw_counts_isoform[relevant_transcripts, ])
  Gene = colSums(Gene)
  Gene_df = as.data.frame(Gene)
  matched_rows <- match(rownames(df), rownames(Gene_df))
  df$Gene <- Gene_df$Gene[matched_rows]

  #Differentiate between major isoform and diff isoform
  if (isoform$isoform != isoform$major) {
    colnames(df)[1] <- paste("Specified", colnames(df)[1], sep = "_")
    colnames(df)[2] <- paste("Major", colnames(df)[2], sep = "_")
  } else {
    colnames(df)[1] <- paste("Major/Specified", colnames(df)[1], sep = "_")
    df <- df[,-2]
  }

  # Assign cells to identities
  match_rows <- match(rownames(df), rownames(seurat_obj@meta.data))
  df$idents <- seurat_obj@meta.data[["cluster"]][match_rows]

  # Convert counts to CPM
  col_num = ncol(df)

  idents <- levels(seurat_obj@active.ident)
  for (ident in idents) {
    subset_data <- subset(seurat_obj, idents = ident)
    raw_count_subset <- subset_data@assays$RNA$counts
    total_counts <- sum(raw_count_subset)
    rows_to_change <- which(df[, col_num] == ident)
    df[rows_to_change, 1:(col_num-1)] <- (df[rows_to_change, 1:(col_num-1)] / total_counts) * 1e6
  }
  # Calculate and add error bars.
  dff = aggregate(. ~ idents, data = df, sum)

  dff_long <- tidyr::pivot_longer(dff,
                                  cols = -idents,
                                  names_to = "Gene",
                                  values_to = "Expression")
  # Calculate the standard error.
  df_long <- tidyr::pivot_longer(df,
                                 cols = -idents,
                                 names_to = "Gene",
                                 values_to = "Expression")
  df_se <- df_long %>%
    group_by(idents, Gene) %>%
    summarise(SE = sd(Expression, na.rm = TRUE) / sqrt(n()), .groups = 'drop')  # 计算标准误差

  #### se + CPM
  dff_long = dff_long %>% left_join(df_se,by = c("idents" = "idents" ,"Gene" = "Gene"))
  dff_long$Expression = log(dff_long$Expression + 1)
  dff_long$SE = log(dff_long$SE + 1)


  # Check the number of values and set the number of legend rows
  num_values <- length(unique(dff_long$Gene))
  if (num_values == 2) {
    mycol1 <- c('#1A0841','#FF5959')
  } else if (num_values == 3) {
    mycol1 <- c('#1A0841','#4F9DA6','#FF5959')
  }

  legend_rows <- num_values

  # Plot a line chart.
  p <- ggplot(dff_long, aes(x = idents, y = Expression, group = Gene, color = Gene)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Expression - SE, ymax = Expression + SE), width = 0.25, size = 0.75) +
    scale_color_manual(values = mycol1) +
    scale_fill_manual(values = mycol1) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
      axis.text.y = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      panel.grid.major = element_line(colour = "gray80"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(color = 'black', size = 16),
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(linewidth = 1, color = 'black'),
      legend.key.size = unit(20, "pt")
    ) +
    guides(color = guide_legend(reverse = F, nrow = legend_rows, byrow = TRUE)) +
    labs(x = "Cluster", y = "log(CPM+1)")
}


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
plotDEIso <- function(seurat_obj,gtf,cluster_column,subset_ident,transcript_id,plot_nrow){
  if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
  library(Seurat)
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  library(dplyr)
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("dplyr")
  library(ggplot2)
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("dplyr")
  library(patchwork)

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
  # gene_id = isoform$gene_id
  # isoform$relevant = length(gtf_trans[gtf_trans$gene_id == gene_id , "transcript_id"])
  gene_id <- unique(isoform$gene_id)
  transcript_count <- gtf_trans %>%
    group_by(gene_id) %>%
    summarise(transcript_count = n(), .groups = 'drop')
  isoform <- isoform %>%
    left_join(transcript_count, by = "gene_id") %>%
    rename(relevant = transcript_count)

  ## plot
  results_list <- list()

  for (i in 1:nrow(isoform)) {

    current_gene_id <- isoform$gene_id[i]
    current_isoform_id <- isoform$isoform[i]
    current_relevant <- isoform$relevant[i]

    # check
    if (current_relevant == 1) {
      message_text <- paste0(
        "The isoform '",
        current_isoform_id,
        "' is the only isoform of the corresponding gene (Gene ID: ",
        current_gene_id,
        ")"
      )
      print(message_text)
      next
    }

    # relevant_transcripts
    relevant_transcripts <- gtf_trans$transcript_id[gtf_trans$gene_id == current_gene_id]
    relevant_transcripts <- na.omit(relevant_transcripts)
    relevant_transcripts <- as.character(relevant_transcripts)
    available_transcripts <- relevant_transcripts[relevant_transcripts %in% rownames(raw_counts_isoform)]

    if (length(available_transcripts) == 0) {
      warning(paste("No expression data for gene", current_gene_id))
      next
    }

    # isoform
    total_expressions <- rowSums(raw_counts_isoform[available_transcripts, , drop = FALSE])
    major_isoform <- names(which.max(total_expressions))

    # df
    df <- data.frame(
      Specified_Isoform = as.numeric(raw_counts_isoform[current_isoform_id, ]),
      Major_Isoform = as.numeric(raw_counts_isoform[major_isoform, ])
    )
    rownames(df) <- colnames(raw_counts_isoform)
    gene_expression <- colSums(raw_counts_isoform[available_transcripts, , drop = FALSE])
    df$Gene_Total <- gene_expression

    # colnames
    if (current_isoform_id != major_isoform) {
      colnames(df)[1] <- paste0("Specified_", current_isoform_id)
      colnames(df)[2] <- paste0("Major_", major_isoform)
    } else {
      df <- df[, c(1, 3), drop = FALSE]
      colnames(df)[1] <- paste0("Major_Specified_", current_isoform_id)
    }

    # annotation
    match_rows <- match(rownames(df), rownames(seurat_obj@meta.data))
    df$idents <- seurat_obj@meta.data[[cluster_column]][match_rows]

    # CPM
    col_num <- ncol(df)
    idents <- levels(seurat_obj@active.ident)

    for (ident in idents) {
      subset_data <- subset(seurat_obj, idents = ident)
      raw_count_subset <- subset_data@assays$RNA$counts
      total_counts <- sum(raw_count_subset)
      rows_to_change <- which(df[, col_num] == ident)
      df[rows_to_change, 1:(col_num-1)] <- (df[rows_to_change, 1:(col_num-1)] / total_counts) * 1e6
    }

    results_list[[current_gene_id]] <- list(
      data = df,
      gene_id = current_gene_id,
      isoform_id = current_isoform_id,
      major_isoform = major_isoform
    )
  }

  # polt function
  plot_expression_profiles <- function(results_list) {
    plots <- list()

    for (i in seq_along(results_list)) {
      result <- results_list[[i]]
      df <- result$data
      dff <- aggregate(. ~ idents, data = df, sum)
      dff_long <- tidyr::pivot_longer(dff,
                                      cols = -idents,
                                      names_to = "Gene",
                                      values_to = "Expression")
      df_long <- tidyr::pivot_longer(df,
                                     cols = -idents,
                                     names_to = "Gene",
                                     values_to = "Expression")

      df_se <- df_long %>%
        group_by(idents, Gene) %>%
        summarise(SE = sd(Expression, na.rm = TRUE) / sqrt(n()), .groups = 'drop')

      dff_long <- dff_long %>%
        left_join(df_se, by = c("idents", "Gene"))
      dff_long$Expression <- log(dff_long$Expression + 1)
      dff_long$SE <- log(dff_long$SE + 1)

      # color
      num_values <- length(unique(dff_long$Gene))
      if (num_values == 2) {
        mycol1 <- c('#1A0841', '#FF5959')
      } else if (num_values == 3) {
        mycol1 <- c('#1A0841', '#4F9DA6', '#FF5959')
      } else {
        mycol1 <- scales::hue_pal()(num_values)
      }

      # plot
      p <- ggplot(dff_long, aes(x = idents, y = Expression, group = Gene, color = Gene)) +
        geom_line(linewidth = 1.5) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = Expression - SE, ymax = Expression + SE), width = 0.25, size = 0.75) +
        scale_color_manual(values = mycol1) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
          axis.text.y = element_text(size = 12),
          legend.position = 'top',
          legend.justification = c(0, 1),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          panel.grid.major = element_line(colour = "gray80"),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = 'black', size = 16),
          axis.line = element_line(linewidth = 1),
          axis.ticks = element_line(linewidth = 1, color = 'black'),
          legend.key.size = unit(20, "pt")
        ) +
        guides(color = guide_legend(nrow = num_values, byrow = TRUE)) +
        labs(
          x = "Cluster",
          y = "log(CPM+1)",
          title = paste("Gene:", result$gene_id)
        )

      plots[[result$gene_id]] <- p
    }

    return(plots)
  }

  all_plots <- plot_expression_profiles(results_list)

  combined_plot <- wrap_plots(all_plots, nrow = plot_nrow)
  return(combined_plot)

}

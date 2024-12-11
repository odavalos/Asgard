#' @title Prepare Drug Reference using a single updated dataset.
#' @description Prepare tissue-specific drug reference profiles from a single L1000 drug response dataset.
#' @details This function converts L1000 data to the tissue specific drug rank matrix for a single dataset.
#' @param cell.info The local path and name of the cell.info text file (from LINCS L1000).
#' @param gene.info The local path and name of the gene.info text file.
#' @param sig.info The local path and name of the signature info (sig.info) text file for the updated dataset.
#' @param gctx The local path and name of the GCTX file for the updated dataset.
#' @param Output.Dir The output directory for the generated files.
#' @import cmapR
#' @export
PrepRefv2 <- function(cell.info = NULL,
                             gene.info = NULL,
                             sig.info = NULL,
                             gctx = NULL,
                             Output.Dir = "./") {
  # Load cell information
  cell_data <- read.table(file = cell.info, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
  celldata$cell_id <- celldata$cell_iname
  
  # Identify tissues
  tissues <- unique(as.character(cell_data$cell_lineage))
  tissues <- tissues[which(tissues != "-666")]
  
  # Load gene information
  gene_data <- read.table(file = gene.info, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
  my_gene_info <- gene_data[, 1:2]
  colnames(my_gene_info) <- c("ID", "Gene.Symbol")
  
  # Load signature information
  col_meta <- read.delim(sig.info, sep = "\t", stringsAsFactors = FALSE)
  col_meta$cell_id <- col_meta$cell_iname
  
  for (tissue in tissues) {
    message("Processing tissue: ", tissue)
    
    # Subset cell lines for this tissue
    cell_ids <- which(cell_data$cell_lineage == tissue)
    cell_names <- cell_data$cell_id[cell_ids]
    
    # Filter signatures for this tissue
    idx <- which(col_meta$cell_id %in% cell_names & col_meta$pert_type == "trt_cp")
    sig_ids <- col_meta$sig_id[idx]
    # Remove replicate suffixes if any (e.g., remove REP.)
    rm.ids <- grep('REP\\.', sig_ids)
    if (length(rm.ids) > 0) {
      sig_ids <- sig_ids[-rm.ids]
    }
    
    # If no signatures are found for this tissue, skip
    if (length(sig_ids) == 0) {
      message("No drug treatments found for tissue: ", tissue)
      next
    }
    
    # Parse GCTX file for selected signatures
    my_ds <- parse_gctx(gctx, cid = sig_ids)
    
    # Rank transformation function
    myrank <- function(x) {
      temp <- rank(-x, ties.method = "min")
      return(temp)
    }
    
    # Create rank matrix
    rank_matrix <- apply(my_ds@mat, 2, myrank)
    rank_matrix <- as.data.frame(rank_matrix)
    
    # Clean column names
    colnames(rank_matrix) <- gsub(":", "_", colnames(rank_matrix))
    cnames <- colnames(rank_matrix)
    # Assign integer column names for the rank matrix
    colnames(rank_matrix) <- 1:length(cnames)
    dcnames <- colnames(rank_matrix)
    
    # Add probe_id column
    rank_matrix$probe_id <- row.names(rank_matrix)
    rank_matrix <- rank_matrix[, c('probe_id', dcnames)]
    
    # Write out rank matrix
    filename <- file.path(Output.Dir, paste0(gsub(" ", "-", tissue), "_rankMatrix.txt"))
    write.table(rank_matrix, file = filename, quote = FALSE, row.names = FALSE, sep = "\t")
    
    # Write out gene info
    gene_filename <- file.path(Output.Dir, paste0(gsub(" ", "-", tissue), "_gene_info.txt"))
    write.table(my_gene_info, file = gene_filename, quote = FALSE, row.names = FALSE, sep = "\t")
    
    # Create drug info file
    # Keep only selected sig_ids
    sub_meta <- col_meta[col_meta$sig_id %in% cnames, ]
    # If no matching signatures are found in sub_meta, skip
    if (nrow(sub_meta) == 0) {
      message("No matching sig_ids found in metadata for tissue: ", tissue)
      next
    }
    
    # Generate drug info table
    # We combine pert_iname, pert_id, concentration (pert_idose), duration (pert_itime), and cell_id
    my_drug_info <- data.frame(
      instance_id = 1:length(sub_meta$sig_id),
      cmap_name = paste(sub_meta$pert_iname, sub_meta$pert_id, sep = "_"),
      concentration..M = sub_meta$pert_idose,
      duration..h = sub_meta$pert_itime,
      cell2 = sub_meta$cell_id,
      catalog_name = sub_meta$pert_id,
      treatment = paste(sub_meta$pert_iname, "_", sub_meta$sig_id, sep = "")
    )
    
    # Write out drug info
    drug_filename <- file.path(Output.Dir, paste0(gsub(" ", "-", tissue), "_drug_info.txt"))
    write.table(my_drug_info, file = drug_filename, quote = FALSE, row.names = FALSE, sep = "\t")
  }
}

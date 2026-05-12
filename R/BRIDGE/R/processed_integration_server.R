#' @export
processed_integration <- function(input, output, session, rv) {
    observeEvent(input$process_integrate_data, {
        req(input$processed_integration, input$heatmap_k)

        withProgress(message = "Running processed integration", value = 0, {

        selected_tables <- input$processed_integration
        filtered_ids <- list()
        dim_info <- list()
        species_vec <- unique(tolower(unlist(rv$species[selected_tables])))
        cross_species_mode <- length(species_vec[!is.na(species_vec) & nzchar(species_vec)]) > 1

        make_match_key <- function(df, cross_species = FALSE) {
            gene_id <- trimws(as.character(df$Gene_ID))
            gene_name <- tolower(trimws(as.character(df$Gene_Name)))
            xid <- tolower(trimws(as.character(df$XID)))

            gene_id[is.na(gene_id)] <- ""
            gene_name[is.na(gene_name)] <- ""
            xid[is.na(xid)] <- ""

            if (cross_species) {
                # Cross-species: prefer orthology/xref, then case-insensitive gene symbol.
                dplyr::case_when(
                    nzchar(xid) ~ paste0("xid:", xid),
                    nzchar(gene_name) ~ paste0("gn:", gene_name),
                    TRUE ~ paste0("gid:", tolower(gene_id))
                )
            } else {
                # Same-species: keep original behavior preference (Gene_ID first).
                dplyr::case_when(
                    nzchar(gene_id) ~ paste0("gid:", tolower(gene_id)),
                    nzchar(xid) ~ paste0("xid:", xid),
                    TRUE ~ paste0("gn:", gene_name)
                )
            }
        }

        extract_features <- function(dep, datatype, contrast) {
            rd <- as.data.frame(SummarizedExperiment::rowData(dep))
            # Add Gene_Name/Gene_ID consistently
            if (datatype == "rnaseq") {
                rd <- as.data.frame(dep@test_result)
                rd$Gene_ID <- gsub("^.*_", "", rownames(rd))
                rd$Gene_Name <- gsub("_.*$", "", rownames(rd))
            } else {
                rd <- as.data.frame(SummarizedExperiment::rowData(dep))
            }

            if (!"XID" %in% colnames(rd)) {
                rd$XID <- NA_character_
            }

            rownames(rd) <- NULL
            rd <- rd %>%
                dplyr::select(where(~ !is.numeric(.)), where(is.numeric)) %>%
                dplyr::mutate(Gene_Name = stringr::str_to_title(Gene_Name))
            # Grab contrast‑specific columns
            lfc_col <- paste0(contrast, "_diff")
            padj_col <- paste0(contrast, "_p.adj")

            if (!all(c(lfc_col, padj_col) %in% colnames(rd))) {
                stop("contrast columns missing")
            }

            rd <- rd |>
                dplyr::transmute(Gene_ID, Gene_Name, XID, diff = .data[[lfc_col]], p.adj = .data[[padj_col]])

            rd$Match_Key <- make_match_key(rd, cross_species = cross_species_mode)
            # message("Extracted features for ", datatype, ": ", nrow(rd), " rows\n COLS: ", paste(colnames(rd), collapse = ", "), "\nHEAD: ", head(rd, n = 3))
            rd
        }

        normalize_ids <- function(dep, datatype) {
            if (methods::is(dep, "DEGdata") || datatype == "rnaseq") {
                rd <- as.data.frame(dep@test_result)
                # If gene IDs/names not separately present, split rownames once
                if (!("Gene_ID" %in% colnames(rd))) {
                    rd$Gene_ID <- sub("^.*_", "", rownames(rd))
                    rd$Gene_Name <- sub("_.*$", "", rownames(rd))
                }
                # message("Normalizing IDs for datatype ", datatype, ": ", nrow(rd), " rows\n COLS: ", paste(colnames(rd), collapse = ", "), "\nHEAD: ", head(rd, n = 1))
            } else {
                rd <- as.data.frame(SummarizedExperiment::rowData(dep))
                # message("Normalizing IDs for datatype ", datatype, ": ", nrow(rd), " rows\n COLS: ", paste(colnames(rd), collapse = ", "), "\nHEAD: ", head(rd, n = 1))
                # Often rowData already has Gene_ID & Gene_Name
                if (!("Gene_ID" %in% colnames(rd))) {
                    stop("No Gene_ID found in rowData for datatype: ", datatype)
                }
            }

            # Ensure correctly typed/clean
            rd <- rd %>%
                dplyr::mutate(
                    Gene_ID   = as.character(Gene_ID),
                    Gene_Name = stringr::str_to_title(as.character(Gene_Name))
                ) %>%
                dplyr::select(where(~ !is.numeric(.)), where(is.numeric))

            rd
        }

        strip_sig <- function(dep_obj) {
            # message("Stripping old sig cols")
            if (methods::is(dep_obj, "DEGdata")) {
                # message("DEGdata: ", str(dep_obj))
                tr <- dep_obj@test_result
                tr <- tr[, !grepl("_significant$|^significant$", colnames(tr))]
                dep_obj@test_result <- tr
            } else {
                # message("SummarizedExperiment: ", class(dep_obj))
                rd <- SummarizedExperiment::rowData(dep_obj)
                rd <- rd[, !grepl("_significant$|^significant$", colnames(rd))]
                # message("RD: ", class(rd))
                SummarizedExperiment::rowData(dep_obj) <- rd
            }
            # message("Stripped old sig cols")
            dep_obj
        }

        get_depflt_matrix <- function(dep_obj, lfc_cut, p_cut, cids) {
            dep_obj <- strip_sig(dep_obj)
            sig <- DEP2::add_rejections(dep_obj, alpha = p_cut, lfc = lfc_cut)
            # If there are no significant rows, or the object is empty, exit early
            if (is.null(sig) || nrow(sig) == 0L) {
                return(list(mat = NULL, mat_scaled = NULL))
            }
            mat <- as.data.frame(SummarizedExperiment::assay(sig))

            if (methods::is(sig, "DEGdata")) {
                df <- as.data.frame(sig@test_result)
                df$id_name <- rownames(df)
                # Rownames are expected as <Gene_Name>_<Gene_ID>; split on LAST underscore
                df$Gene_Name <- sub("_[^_]+$", "", df$id_name)
                df$Gene_ID   <- sub("^.*_", "", df$id_name)
                df$XID <- "" # There are no XID in rnaseq data
                df <- cbind(df, mat)
            } else {
                df <- cbind(as.data.frame(SummarizedExperiment::rowData(sig)), mat)  
            }

            if (!"XID" %in% colnames(df)) {
                df$XID <- NA_character_
            }
            df$Match_Key <- make_match_key(df, cross_species = cross_species_mode)

            keep_idx <- which(df$Match_Key %in% cids)
            mat <- mat[keep_idx, , drop = FALSE]
            # If assay is NULL or empty, exit early
            mat <- as.matrix(mat)
            df_filt <- df[keep_idx, , drop = FALSE]
            rownames(mat) <- paste0(df_filt$Gene_ID, "_", df_filt$Gene_Name)
            if (is.null(mat) || length(mat) == 0L) {
                return(list(mat = NULL, mat_scaled = NULL))
            }
            # sanitize matrix (rows/cols need finite variance and ≥2 finite values)
            row_ok <- rowSums(is.finite(mat)) >= 2 & apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
            col_ok <- colSums(is.finite(mat)) >= 2 & apply(mat, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
            mat <- mat[row_ok, col_ok, drop = FALSE]            
            mat_scaled <- safe_row_scale(mat)
            #message("Scaled matrix: ", nrow(mat_scaled), " rows, ", ncol(mat_scaled), " cols", "\n IDS: ", paste(colnames(mat_scaled), collapse = ", "), "\nGENES: ", head(rownames(mat_scaled), n = 3), "\n")
            mat_scaled
        }

        get_df_from_dep <- function(dep) {            
            if (methods::is(dep, "DEGdata")) {
                res <- as.data.frame(dep@test_result) 
                res$id_name <- rownames(res)
                # Rownames are expected as <Gene_Name>_<Gene_ID>; split on LAST underscore
                res$Gene_Name <- sub("_[^_]+$", "", res$id_name)
                res$Gene_ID   <- sub("^.*_", "", res$id_name)
                res$XID <- "" # There are no XID in rnaseq data
            } else {
                res <- as.data.frame(SummarizedExperiment::rowData(dep))
            }
            if (!"XID" %in% colnames(res)) {
                res$XID <- NA_character_
            }
            return(res)
        }

        get_matrix_from_dep <- function(dep) {    
            mat <- as.data.frame(SummarizedExperiment::assay(dep))
            if (methods::is(dep, "DEGdata")) {
                df <- as.data.frame(dep@test_result)
                df$id_name <- rownames(df)
                # Rownames are expected as <Gene_Name>_<Gene_ID>; split on LAST underscore
                df$Gene_Name <- sub("_[^_]+$", "", df$id_name)
                df$Gene_ID   <- sub("^.*_", "", df$id_name)
                df$XID <- "" # There are no XID in rnaseq data
            } else {
                df <- as.data.frame(SummarizedExperiment::rowData(dep))
            }       
            #mat$names <- rownames(mat)
            mat$Gene_Name <- stringr::str_to_title(df$Gene_Name)  # gsub("_.*", "", mat$names, perl = TRUE)
            mat$Gene_ID <- df$Gene_ID  # gsub(".*_", "", mat$names, perl = TRUE)
            if("pepG" %in% colnames(df)) {
                mat$pepG <- df$pepG
            }
            if("Protein_ID" %in% colnames(df)) {
                mat$Protein_ID <- df$Protein_ID
            }
            if ("XID" %in% colnames(df)) {
                mat$XID <- df$XID
            } else {
                mat$XID <- NA_character_
            }
            return(mat)
        }

        incProgress(0.1, detail = "Filtering significant features")
        # Filter each table by thresholds
        for (tbl in selected_tables) {
            contrast <- input[[paste0("pi_comparison_selected_", tbl)]]
            dep <- rv$dep_output[[tbl]]
            datatype <- rv$datatype[[tbl]]

            feats <- extract_features(dep, datatype, contrast)            
            keep <- with(
                feats,
                !is.na(Match_Key) &
                    nzchar(Match_Key) &
                    !is.na(diff) &
                    !is.na(p.adj) &
                    abs(diff) >= input$lfc_thresh_pi &
                    p.adj <= input$pval_thresh_pi
            )

            filtered_ids[[tbl]] <- unique(feats$Match_Key[keep])

            dim_info[[tbl]] <- list(
                original   = dim(SummarizedExperiment::assay(dep)),
                filtered   = c(length(filtered_ids[[tbl]]), ncol(SummarizedExperiment::assay(dep)))
            )
            incProgress(0.2 / max(1, length(selected_tables)), detail = paste("Filtered", tbl))
        }

        incProgress(0.1, detail = "Intersecting feature IDs")
        # Intersect IDs
        common_ids <- Reduce(intersect, filtered_ids)
        if (length(common_ids) < input$heatmap_k) {
            if (length(common_ids) == 0) {
                showNotification("Not enough intersecting significant IDs found, consider changing the parameters.", type = "error")
                rv$intersected_tables_processed <- NULL
                rv$integration_preview_dims <- NULL
                return()
            }
        }

        incProgress(0.15, detail = "Building intersected tables")
        # Subset SummarizedExperiment::assays by intersected gene names
        intersected_list <- lapply(selected_tables, function(tbl) {
            dep <- rv$dep_output[[tbl]]

            res <- get_df_from_dep(dep)
            mat <- get_matrix_from_dep(dep)
            res$Match_Key <- make_match_key(res, cross_species = cross_species_mode)
            mat$Match_Key <- make_match_key(mat, cross_species = cross_species_mode)

            dep_flt <- res[res$Match_Key %in% common_ids, ]
            mat_flt <- mat[mat$Match_Key %in% common_ids, ]

            join_keys <- "Match_Key"
            
            if ("pepG" %in% colnames(mat_flt) && "pepG" %in% colnames(dep_flt)) {
                join_keys <- c(join_keys, "pepG")
            } else if ("Protein_ID" %in% colnames(mat_flt) && "Protein_ID" %in% colnames(dep_flt)) {
                join_keys <- c(join_keys, "Protein_ID")
            }

            data <- dplyr::inner_join(mat_flt, dep_flt, by = join_keys, keep=FALSE, suffix=c("",".y")) %>%
                dplyr::select(-ends_with(".y"))

            # Generate unique IDs
            if ("pepG" %in% colnames(data)) {
                data$unique_id <- paste(data$Gene_Name, data$pepG, sep = "_")
            } else if ("Protein_ID" %in% colnames(data)) {
                data$unique_id <- paste(data$Gene_Name, data$Protein_ID, sep = "_")
            } else {
                data$unique_id <- data$Gene_Name
            }            
            data
        })

        names(intersected_list) <- selected_tables

        incProgress(0.15, detail = "Extracting intersected matrices")
        # Extract intersected expression matrices per table
        intersected_matrix <- lapply(selected_tables, function(tbl) {
            dep <- rv$dep_output[[tbl]]
            datatype <- rv$datatype[[tbl]]
            mat <- get_depflt_matrix(dep, input$lfc_thresh_pi, input$pval_thresh_pi, common_ids)
            mat
        })
        names(intersected_matrix) <- selected_tables
        
        #message("Intersected matrices extracted: ", paste(sapply(selected_tables, function(tbl) {
        #    if (is.null(intersected_matrix[[tbl]])) {
        #        paste0(tbl, "NULL")
        #    } else {
        #        paste(tbl, dim(intersected_matrix[[tbl]]), collapse = " x ")
        #    }
        #}), collapse = ", "), "\n")

        # Record dimensions
        for (tbl in selected_tables) {
            dim_info[[tbl]]$intersected <- dim(intersected_matrix[[tbl]])
        }

        incProgress(0.1, detail = "Preparing scatter data")
        # Prepare scatter data: wide table by Gene_ID
        scatter_data <- NULL
        for (tbl in selected_tables) {
            contrast <- input[[paste0("pi_comparison_selected_", tbl)]]
            dep <- rv$dep_output[[tbl]]
            datatype <- rv$datatype[[tbl]]
            feats <- extract_features(dep, datatype, contrast)
            #message("Features for scatter from ", tbl, ": ", nrow(feats), " rows\n COLS: ", paste(colnames(feats), collapse = ", "), "\nHEAD: ", head(feats$Gene_ID, n = 3))
            feats <- feats[feats$Match_Key %in% common_ids, ]
            feats <- feats |>
                dplyr::select(Match_Key, Gene_Name, diff) |>
                dplyr::rename(!!paste0("LFC_", tbl) := diff)

            if (is.null(scatter_data)) {
                scatter_data <- feats
            } else {
                feats_small <- feats |> dplyr::select(Match_Key, !!paste0("LFC_", tbl))
                scatter_data <- dplyr::full_join(scatter_data, feats_small, by = "Match_Key")
            }
        }

        incProgress(0.1, detail = "Building scatter plots")
        # Build pairwise scatterplots
        lfc_cols <- grep("^LFC_", names(scatter_data), value = TRUE)
        plot_list <- list()
        for (i in 1:(length(lfc_cols) - 1)) {
            for (j in (i + 1):length(lfc_cols)) {
                df_plot <- scatter_data |>
                    dplyr::select(Gene_Name, x = !!lfc_cols[i], y = !!lfc_cols[j]) |>
                    na.omit()

                if (nrow(df_plot) > 0) {
                    axis_rng <- range(c(df_plot$x, df_plot$y), na.rm = TRUE)
                    if (!all(is.finite(axis_rng))) next
                    if (diff(axis_rng) == 0) {
                        axis_rng <- axis_rng + c(-0.5, 0.5)
                    }

                    p <- ggplot(df_plot, aes(x, y, text = Gene_Name)) +
                        geom_point(alpha = .7, color = BRIDGE_COLORS$accent) +
                        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#AAAAAA") +
                        coord_equal(xlim = axis_rng, ylim = axis_rng, expand = TRUE) +
                        labs(title = paste(lfc_cols[i], "vs", lfc_cols[j]), x = paste(lfc_cols[i]), y = paste(lfc_cols[j])) +
                        theme_minimal()
                    plot_list[[paste(lfc_cols[i], lfc_cols[j], sep = "_vs_")]] <- p
                }
            }
        }

        message("Processed integration with ", length(common_ids), " common IDs across ", length(selected_tables), " tables.\nIntersect: ", paste(selected_tables, collapse = ", "), "\nDimensions: ", paste(sapply(selected_tables, function(tbl) {
            dim_info[[tbl]]$original[1]
        }), collapse = " x "), " original rows; ", paste(sapply(selected_tables, function(tbl) {
            dim_info[[tbl]]$filtered[1]
        }), collapse = " x "), " filtered rows; ", paste(sapply(selected_tables, function(tbl) {
            dim_info[[tbl]]$intersected[1]
        }), collapse = " x "), " intersected rows.\n")

        incProgress(0.1, detail = "Saving results")
        # Save results
        rv$scatter_plots <- plot_list
        rv$intersected_tables_processed <- intersected_list
        rv$intersected_matrix_processed <- intersected_matrix
        rv$integration_preview_dims <- lapply(selected_tables, function(tbl) dim_info[[tbl]])
        names(rv$integration_preview_dims) <- selected_tables
        # Estimate cluster number on stacked data
        # Check if all matrices have the same number of columns
        col_counts <- sapply(intersected_matrix, ncol)
        all_same_cols <- length(unique(col_counts)) == 1

        # Conditionally create stacked
        if (all_same_cols) {
            stacked <- do.call(rbind, intersected_matrix)  # genes x samples 
            # If same number of columns stacking is possible
        } else {
            stacked <- intersected_matrix[[1]]  # use only the first table
            # If differing number of columns use only first matrix for optimal k computation
        }
        if (nrow(stacked) < 2) {
            showNotification("Not enough genes for clustering.", type = "error")
            rv$optimal_k <- NULL
            return()
        }
        if (length(common_ids) < input$heatmap_k) {
            optimal_k = 1
        } else {
            optimal_k <- safe_nbclust(stacked, k_min = 2, k_max = 10)
        }
        rv$optimal_k <- ifelse(is.na(optimal_k), NULL, optimal_k)

        incProgress(0.1, detail = "Done")

        })

    })
}

#' @export
heatmap_server <- function(id, rv) {
    moduleServer(id, function(input, output, session) {        
        shiny::observe({
            lapply(rv$table_names, function(tbl_name) {
                output[[paste0("raw_ht_", tbl_name)]] <- shiny::renderPlot({
                    ComplexHeatmap::Heatmap(rv$ht_matrix[[tbl_name]], show_row_dend = FALSE, show_row_names = FALSE)
                })
            })
        })
    })
}

#' @export
RawHeatmapServer <- function(id, rv, tbl_name) {
    moduleServer(id, function(input, output, session) {
        ready <- reactiveVal(FALSE)

        task <- ExtendedTask$new(function(raw_input) {
            promises::future_promise({
                raw_data <- raw_input$raw_data
                datapoint_cols <- raw_input$datapoint_cols
                datatype <- raw_input$datatype

                clean <- as.matrix(raw_data[datapoint_cols])
                rn <- switch(datatype,
                    proteomics        = raw_data$Gene_Name,
                    rnaseq            = raw_data$Gene_Name,
                    phosphoproteomics = raw_data$pepG
                )
                rownames(clean) <- rn
                list(clean_matrix = clean)
            })
        })


        observeEvent(input$compute, {
            ready(TRUE)
            task$invoke(list(
                raw_data       = isolate(rv$tables[[tbl_name]]),
                datapoint_cols = isolate(rv$data_cols[[tbl_name]]),
                datatype       = isolate(rv$datatype[[tbl_name]])
            ))
        }, ignoreInit = TRUE)

        # Mount spinner only after "Compute" is pressed
        output$plot_slot <- renderUI({
            if (!ready()) {
                return(div(
                    style = "padding:10px; color:#777;",
                    "Click “Compute Raw Heatmap” to generate the plot."
                ))
            }
            shinycssloaders::withSpinner(
                plotOutput(session$ns("plot"), height = "520px"),
                type = 8, color = "#2b8cbe", caption = "Loading..."
            )
        })

        output$plot <- renderPlot({
            res <- task$result()
            req(res)
            ComplexHeatmap::Heatmap(res$clean_matrix,
                show_row_dend = FALSE,
                show_row_names = FALSE
            )
        })
    })
}


#### THIS FUNCTION CONTAINS THE SERVER FOR THE DEP HEATMAP
#' @export
DepHeatmapServer <- function(id, rv, cache, tbl_name) {
    moduleServer(id, function(input, output, session) {
        # Initialize per-dataset reactive storage if it doesn't exist
        if (is.null(rv$heatmap_state)) {
            rv$heatmap_state <- reactiveValues()
        }
        
        # Create isolated state for each dataset
        if (is.null(rv$heatmap_state[[tbl_name]])) {
            rv$heatmap_state[[tbl_name]] <- reactiveValues(
                heatmap_ready      = FALSE,
                cluster_ready      = FALSE,
                depflt_cache       = list(),
                last_params        = NULL,
                last_mat           = NULL,
                roword             = NULL,
                heatmap_rows       = NULL,
                cluster_table_data = NULL
            )
        }
        
        # ref to this dataset's state
        state <- rv$heatmap_state[[tbl_name]]
        ns <- session$ns

        strip_sig <- function(dep_obj) {
            if (methods::is(dep_obj, "DEGdata")) {
                tr <- dep_obj@test_result
                tr <- tr[, !grepl("_significant$|^significant$", colnames(tr))]
                dep_obj@test_result <- tr
            } else {
                rd <- SummarizedExperiment::rowData(dep_obj)
                rd <- rd[, !grepl("_significant$|^significant$", colnames(rd))]
                SummarizedExperiment::rowData(dep_obj) <- rd
            }
            dep_obj
        }

        get_dep_result <- function() {
            key <- rv$current_dep_heatmap_key[[tbl_name]]
            if (!is.null(key) && cache$exists(key)) {
                return(cache$get(key))
            }
            heatmap_task$result()
        }

        get_depflt <- function(params) {
            key <- paste(
                params$tbl_name, params$columns_key,
                sprintf("pcut=%.4f", params$p_cut),
                sprintf("lfc=%.3f", params$lfc_cut),
                "depflt",
                sep = "_"
            )

            if (!is.null(state$depflt_cache[[key]])) {
                return(state$depflt_cache[[key]])
            }

            dep_output <- params$dep_output
            dep_output <- strip_sig(dep_output)
            dep_flt    <- DEP2::add_rejections(dep_output, alpha = params$p_cut, lfc = params$lfc_cut)

            mats <- get_depflt_matrix(dep_flt, datatype = params$datatype)
            mat  <- mats$mat
            mat_scaled <- mats$mat_scaled
            meta <- mats$meta
            unique_id <- mats$unique_id

            state$last_mat <- mat_scaled

            ret_list <- list(
                mat       = mat,
                mat_scaled = mat_scaled,
                dep_flt    = dep_flt,
                meta       = meta,
                unique_id  = unique_id
            )
            state$depflt_cache[[key]] <- ret_list
            ret_list
        }


        get_depflt_matrix <- function(dep_obj, datatype) {
            sig <- DEP2::get_signicant(dep_obj)
            if (is.null(sig) || nrow(sig) == 0L) {
                #message("get_depflt_matrix: sig is NULL or empty")
                state$last_mat <- NULL
                return(list(mat = NULL, mat_scaled = NULL, unique_id = NULL, meta = NULL))
            }

            sig_assay <- SummarizedExperiment::assay(sig)
            mat <- as.matrix(sig_assay)
            if (is.null(mat) || length(mat) == 0L) {
                #message("get_depflt_matrix: mat is NULL or empty")
                state$last_mat <- NULL
                return(list(mat = NULL, mat_scaled = NULL, unique_id = NULL, meta = NULL))
            }

            row_ok <- rowSums(is.finite(mat)) >= 2 &
                apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
            col_ok <- colSums(is.finite(mat)) >= 2 &
                apply(mat, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)

            if (!any(row_ok) || !any(col_ok)) {
                #message("get_depflt_matrix: all rows or cols filtered out")
                state$last_mat <- NULL
                return(list(mat = NULL, mat_scaled = NULL, unique_id = NULL, meta = NULL))
            }

            mat <- mat[row_ok, col_ok, drop = FALSE]
            #message("get_depflt_matrix: after filtering, nrow(mat)=", nrow(mat),
            #        " ncol(mat)=", ncol(mat), " datatype=", datatype)

            mat_scaled <- safe_row_scale(mat)
            state$last_mat <- mat_scaled

            rd <- SummarizedExperiment::rowData(sig)[row_ok, , drop = FALSE]

            # Generate unique_id based on datatype
            unique_id <- switch(datatype,
                proteomics        = paste0(stringr::str_to_title(rd$Gene_Name), "_", rd$Protein_ID),
                phosphoproteomics = paste0(stringr::str_to_title(rd$Gene_Name), "_", rd$Protein_ID, "_", rd$pepG),
                rnaseq            = {
                    rn <- rownames(mat)
                    if (is.null(rn) || all(rn == "")) {
                        paste0("gene_", seq_len(nrow(mat)))
                    } else {
                        rn
                    }
                },
                {
                    # Fallback if datatype is somehow missing
                    rn <- rownames(mat)
                    if (is.null(rn) || all(rn == "")) {
                        paste0("feat_", seq_len(nrow(mat)))
                    } else {
                        rn
                    }
                }
            )

            #message("get_depflt_matrix: length(unique_id)=", length(unique_id))

            meta <- rd
            list(mat = mat, mat_scaled = mat_scaled, unique_id = unique_id, meta = meta)
        }


        heatmap_task <- ExtendedTask$new(function(args) {
            promises::future_promise(
                {
                    dep_output <- args$dep_output
                    p_cut <- args$p_cut
                    lfc_cut <- args$lfc_cut
                    datatype <- args$datatype

                    if (methods::is(dep_output, "DEGdata")) {
                        tr <- dep_output@test_result
                        tr <- tr[, !grepl("_significant$|^significant$", colnames(tr))]
                        dep_output@test_result <- tr
                    } else {
                        rd <- SummarizedExperiment::rowData(dep_output)
                        rd <- rd[, !grepl("_significant$|^significant$", colnames(rd))]
                        SummarizedExperiment::rowData(dep_output) <- rd
                    }

                    dep_flt <- DEP2::add_rejections(dep_output, alpha = p_cut, lfc = lfc_cut)
                    
                    # Use the SAME filtering logic as get_depflt_matrix
                    sig <- DEP2::get_signicant(dep_flt)
                    mat <- as.matrix(SummarizedExperiment::assay(sig))
                    row_ok <- rowSums(is.finite(mat)) >= 2 & apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
                    col_ok <- colSums(is.finite(mat)) >= 2 & apply(mat, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
                    mat <- mat[row_ok, col_ok, drop = FALSE]
                    mat_scaled <- safe_row_scale(mat)
                    
                    tmp <- tempfile(fileext = ".pdf")
                    pdf(tmp)
                    on.exit({
                        dev.off()
                        unlink(tmp)
                    }, add = TRUE)

                    optimal_k <- safe_nbclust(mat_scaled, k_min = 2, k_max = 10)

                    # Build df from the FILTERED sig object (we may still want it even if k is NA)
                    rd_filtered <- SummarizedExperiment::rowData(sig)[row_ok, , drop = FALSE]
                    df <- cbind(as.data.frame(mat), as.data.frame(rd_filtered))

                    # Generate unique_id matching the filtered rows
                    unique_id <- switch(datatype,
                        proteomics        = paste0(stringr::str_to_title(df$Gene_Name), "_", df$Protein_ID),
                        phosphoproteomics = paste0(stringr::str_to_title(df$Gene_Name), "_", df$Protein_ID, "_", df$pepG),
                        rnaseq            = rownames(df),
                        rownames(df)
                    )
                    df$unique_id <- unique_id

                    if (datatype == "rnaseq") {
                        df$Gene_Name <- gsub("_.*", "", df$unique_id, perl = TRUE)
                        df$Gene_ID   <- gsub(".*_", "", df$unique_id, perl = TRUE)
                    }

                    df <- df %>% dplyr::select(Gene_ID, Gene_Name, unique_id, dplyr::everything())

                    # Always return a list; if optimal_k is NA this just means “no suggestion”
                    list(optimal_k = optimal_k, df = df, datatype = datatype)
                },
                seed = TRUE
            )
        })

        observeEvent(input$recompute_heatmap, {
            state$cluster_ready <- FALSE
            state$roword <- NULL
            state$heatmap_rows <- NULL
            state$cluster_table_data <- NULL
            
            req(rv$dep_output[[tbl_name]])
            params <- list(
                tbl_name    = tbl_name,
                p_cut       = input$heatmap_pcutoff,
                lfc_cut     = input$heatmap_fccutoff,
                clustering  = isTRUE(input$clustering),
                k           = input$num_clusters,
                dep_output  = isolate(rv$dep_output[[tbl_name]]),
                datatype    = isolate(rv$datatype[[tbl_name]]),
                columns_key = paste(isolate(rv$data_cols[[tbl_name]]), collapse = "_"),
                dummy       = NULL
            )

            state$last_params <- params

            key <- paste(
                tbl_name, params$columns_key,
                sprintf("pcut=%.4f", params$p_cut),
                sprintf("lfc=%.3f", params$lfc_cut),
                "dep_heatmap_task",
                sep = "_"
            )
            rv$current_dep_heatmap_key[[tbl_name]] <- key

            if (!cache$exists(key)) {
                heatmap_task$invoke(list(
                    dep_output = params$dep_output,
                    p_cut      = params$p_cut,
                    lfc_cut    = params$lfc_cut,
                    datatype   = params$datatype
                ))
            }

            # No prewarming of depflt cache anymore
            #depflt_key <- paste(
            #    tbl_name, params$columns_key,
            #    sprintf("pcut=%.4f", params$p_cut),
            #    sprintf("lfc=%.3f", params$lfc_cut),
            #    "depflt",
            #    sep = "_"
            #)

            #if (is.null(state$depflt_cache[[depflt_key]])) {
            #    pr <- promises::future_promise({
            #        dep_output <- strip_sig(params$dep_output)
            #        dep_flt <- DEP2::add_rejections(dep_output, alpha = params$p_cut, lfc = params$lfc_cut)

            #        sig <- DEP2::get_signicant(dep_flt)
            #        if (is.null(sig) || nrow(sig) == 0L) {
            #            return(list(mat = NULL, mat_scaled = NULL, dep_flt = dep_flt))
            #        }

            #        mat <- as.matrix(SummarizedExperiment::assay(sig))
            #        row_ok <- rowSums(is.finite(mat)) >= 2 &
            #            apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
            #        col_ok <- colSums(is.finite(mat)) >= 2 &
            #            apply(mat, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
            #        if (!any(row_ok) || !any(col_ok)) {
            #            return(list(mat = NULL, mat_scaled = NULL, dep_flt = dep_flt))
            #        }

            #        mat <- mat[row_ok, col_ok, drop = FALSE]
            #        mat_scaled <- safe_row_scale(mat)

            #        list(mat = mat, mat_scaled = mat_scaled, dep_flt = dep_flt)
            #    }, seed = TRUE)

            #    pr <- promises::then(pr, function(ret) {
            #        state$depflt_cache[[depflt_key]] <- ret
            #        state$last_mat <- ret$mat_scaled
            #        NULL
            #    })

            #    pr <- promises::catch(pr, function(e) {
            #        NULL
            #    })
            #}
            state$heatmap_ready <- TRUE
        }, ignoreInit = TRUE)


        output$ht_slot <- renderUI({
            if (!state$heatmap_ready) {
                return(div(
                    style = "padding:10px; color:#777;",
                    "Choose settings and click \"Compute heatmap\"."
                ))
            }
            shinycssloaders::withSpinner(
                plotOutput(session$ns("ht"), height = "800px", width = "100%"),
                type = 8, color = "#2b8cbe", caption = "Loading..."
            )
        })


        output$ht <- renderPlot({
            req(isTRUE(state$heatmap_ready))
            key    <- rv$current_dep_heatmap_key[[tbl_name]]
            res    <- get_dep_result()
            req(res)
            params <- req(state$last_params)

            if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

            dep_flt_list <- get_depflt(params)
            dep_flt      <- dep_flt_list$dep_flt
            meta         <- dep_flt_list$meta
            unique_id    <- dep_flt_list$unique_id

            # store ids used in the heatmap
            state$heatmap_rows <- unique_id

            k <- NULL
            if (params$clustering) {
                k <- as.integer(params$k)
                k_max <- max(2L, nrow(dep_flt) - 1L)
                if (!is.finite(k) || k < 2L) k <- 2L
                if (k > k_max) {
                    showNotification(sprintf(
                        "Reducing k from %d to %d (only %d rows).",
                        input$num_clusters, k_max, nrow(dep_flt)
                    ), type = "warning")
                    k <- k_max
                }
            }

            #message("DepHeatmapServer[", tbl_name, "]: datatype=", params$datatype,
            #" clustering=", params$clustering,
            #" k=", params$k)

            local_mat <- dep_flt_list$mat_scaled
                        
            if (is.null(local_mat) || nrow(local_mat) == 0L || ncol(local_mat) == 0L) {
                # Still try to build a minimal cluster table if possible
                isolate({
                    if (is.null(state$cluster_table_data)) {
                        # No valid matrix, so no cluster table either
                        state$cluster_table_data <- data.frame(
                            Message = "No significant features with the current FDR/log2FC cutoffs.",
                            stringsAsFactors = FALSE
                        )
                        state$cluster_ready <- TRUE
                    }
                })
                
                # Show warning plot
                plot.new()
                par(mar = c(1, 1, 1, 1))
                text(0.5, 0.5,
                    paste0(
                        "Not enough differentially expressed genes found\n",
                        "with current settings:\n\n",
                        sprintf("FDR cutoff: %.4f\n", params$p_cut),
                        sprintf("log2FC cutoff: %.3f\n\n", params$lfc_cut),
                        "Try relaxing the cutoffs to see more genes."
                    ),
                    cex = 1.2,
                    col = "#d95f02",
                    font = 2
                )
                return(invisible())
            }

            # proceed with heatmap
            y_lim <- range(local_mat, na.rm = TRUE)
            if (!all(is.finite(y_lim)) || diff(y_lim) == 0) y_lim <- c(-1, 1)

            if (!is.null(k)) {
                ht <- ComplexHeatmap::Heatmap(
                    local_mat,
                    row_km         = k,
                    show_row_names = FALSE,
                    cluster_columns = FALSE
                ) + ComplexHeatmap::rowAnnotation(
                    profile = ComplexHeatmap::anno_empty(
                        which  = "row",
                        width  = grid::unit(55, "mm"),
                        height = grid::unit(30, "mm")
                    ),
                    annotation_name_gp = grid::gpar(col = NA)
                )
            } else {
                ht <- ComplexHeatmap::Heatmap(
                    local_mat,
                    cluster_rows   = FALSE,
                    show_row_names = FALSE,
                    cluster_columns = FALSE
                )
            }

            ht <- ComplexHeatmap::draw(ht, merge_legend = TRUE, newpage = TRUE)

            rord_val <- ComplexHeatmap::row_order(ht)
            state$roword <- rord_val

            # Build and freezw the cluster table data here
            isolate({
                ids_in_heatmap <- unique_id
                cluster.all    <- NULL

                #message("DepHeatmapServer[", tbl_name, "]: building cluster table, n_ids=",
                #        length(ids_in_heatmap), " is.list(rord)=", is.list(rord_val))

                if (!length(ids_in_heatmap) || !length(rord_val)) {
                    state$cluster_table_data <- NULL
                    state$cluster_ready      <- FALSE
                    return(NULL)
                }

                if (is.list(rord_val)) {
                    for (i in seq_along(rord_val)) {
                        idx <- rord_val[[i]]
                        if (!length(idx)) next

                        clu_ids <- as.character(ids_in_heatmap[idx])
                        if (!length(clu_ids)) next

                        cluster.all <- rbind(
                            cluster.all,
                            data.frame(
                                unique_id = clu_ids,
                                Cluster   = rep(paste0("cluster ", i), length(clu_ids)),
                                stringsAsFactors = FALSE
                            )
                        )
                    }
                } else {
                    if (length(rord_val)) {
                        clu_ids <- as.character(ids_in_heatmap[rord_val])
                        if (length(clu_ids)) {
                            cluster.all <- data.frame(
                                unique_id = clu_ids,
                                Cluster   = rep("cluster 1", length(clu_ids)),
                                stringsAsFactors = FALSE
                            )
                        }
                    }
                }

                if (is.null(cluster.all) || !nrow(cluster.all)) {
                    #message("DepHeatmapServer[", tbl_name, "]: cluster.all is empty, no table.")
                    state$cluster_table_data <- NULL
                    state$cluster_ready      <- FALSE
                    return(NULL)
                }

                # Build df from meta
                df <- as.data.frame(meta)
                df$unique_id <- as.character(ids_in_heatmap)

                if (!"Gene_ID" %in% colnames(df)) {
                    df$Gene_ID <- stringr::str_to_upper(gsub("[a-zA-Z0-9]+_", "", rownames(df)))
                }
                if (!"Gene_Name" %in% colnames(df)) {
                    alt_name_cols <- intersect(
                        c("Gene_Name", "gene_name", "GENE_NAME", "Symbol", "SYMBOL"),
                        colnames(df)
                    )
                    if (length(alt_name_cols) > 0) {
                        df$Gene_Name <- df[[alt_name_cols[1]]]
                    } else {
                        df$Gene_Name <- stringr::str_to_title(gsub("_[a-zA-Z0-9]+", "", rownames(df)))
                    }
                }
                df$Gene_Name <- stringr::str_to_title(df$Gene_Name)

                # Join cluster mapping
                df_clustered <- dplyr::left_join(cluster.all, df, by = "unique_id")
                if ("ID" %in% colnames(df_clustered)) {
                    df_clustered <- df_clustered %>% dplyr::select(-ID)
                }

                leading_cols <- c("Gene_ID", "Gene_Name", "unique_id", "Cluster")
                leading_cols <- intersect(leading_cols, colnames(df_clustered))
                other_cols   <- setdiff(colnames(df_clustered), leading_cols)
                df_clustered <- df_clustered[, c(leading_cols, other_cols), drop = FALSE]
                rownames(df_clustered) <- NULL

                #message("DepHeatmapServer[", tbl_name, "]: cluster table rows=", nrow(df_clustered))

                state$cluster_table_data <- df_clustered
                state$cluster_ready      <- TRUE
            })

            # --- Decorate profile annotation (only if clustering is enabled) ---
            if (!is.null(k)) {
                if (is.list(rord_val)) {
                    n_slices <- length(rord_val)
                } else {
                    rord_val <- list(rord_val)
                    n_slices <- 1L
                }

                cord <- ComplexHeatmap::column_order(ht)
                if (is.list(cord)) cord <- cord[[1L]]

                labs <- colnames(local_mat)[cord]
                xs   <- seq_along(cord)

                for (s in seq_len(n_slices)) {
                    idx <- rord_val[[s]]
                    sub <- local_mat[idx, cord, drop = FALSE]
                    if (!is.matrix(sub) || nrow(sub) == 0L) next

                    med <- apply(sub, 2, stats::median,   na.rm = TRUE)
                    q1  <- apply(sub, 2, stats::quantile, probs = 0.25, na.rm = TRUE)
                    q3  <- apply(sub, 2, stats::quantile, probs = 0.75, na.rm = TRUE)

                    ComplexHeatmap::decorate_annotation("profile", slice = s, {
                        grid::pushViewport(grid::viewport(
                            xscale = c(0.5, length(xs) + 0.5),
                            yscale = y_lim,
                            clip   = "on"
                        ))
                        on.exit(grid::popViewport(), add = TRUE)

                        # midline at 0
                        grid::grid.lines(
                            x = grid::unit(c(1, length(xs)), "native"),
                            y = grid::unit(c(0, 0),         "native"),
                            gp = grid::gpar(col = "blue")
                        )

                        # IQR ribbon
                        ok <- is.finite(q1) & is.finite(q3)
                        if (any(ok)) {
                            x_ok <- xs[ok]
                            grid::grid.polygon(
                                x = grid::unit(c(x_ok, rev(x_ok)), "native"),
                                y = grid::unit(c(q1[ok], rev(q3[ok])), "native"),
                                gp = grid::gpar(fill = "#00000022", col = NA)
                            )
                        }

                        # median trend
                        okm <- is.finite(med)
                        if (any(okm)) {
                            runs <- split(xs[okm], cumsum(c(1, diff(xs[okm]) != 1)))
                            for (r in runs) {
                                grid::grid.lines(
                                    x = grid::unit(r,        "native"),
                                    y = grid::unit(med[r],   "native"),
                                    gp = grid::gpar(col = "black", lwd = 2)
                                )
                            }
                            grid::grid.points(
                                x = grid::unit(xs[okm], "native"),
                                y = grid::unit(med[okm],"native"),
                                pch = 16, size = grid::unit(0.7, "mm")
                            )
                        } else {
                            grid::grid.text("no data", gp = grid::gpar(cex = 0.6, col = "grey50"))
                        }

                        # x‑axis labels only on bottom slice
                        if (s == n_slices) {
                            grid::grid.xaxis(at = xs, label = FALSE)

                            for (i in xs) {
                                tg    <- grid::textGrob(labs[i], gp = grid::gpar(cex = 0.55), rot = 90)
                                w_nat <- grid::convertX(grid::grobWidth(tg),  "native")
                                h_mm  <- grid::convertY(grid::grobHeight(tg), "mm", valueOnly = TRUE)

                                x_pos <- grid::unit(i, "native")
                                if (i == xs[1]) {
                                    x_pos <- x_pos + 0.2 * w_nat
                                } else if (i == xs[length(xs)]) {
                                    x_pos <- x_pos - 0.2 * w_nat
                                }

                                y_pos <- grid::unit(h_mm / 2 + 1.2, "mm")

                                grid::grid.text(
                                    labs[i],
                                    x   = x_pos,
                                    y   = y_pos,
                                    rot = 90,
                                    just = c(0.5, 0.5),
                                    gp  = grid::gpar(cex = 0.55)
                                )
                            }
                        }
                    })
                }
            }
        })


        output$ht_panel <- renderUI({
            DT::DTOutput(session$ns("ht_sig"), height = "300px")
        })

        output$ht_sig <- DT::renderDT({
            req(isTRUE(state$cluster_ready))
            df_clustered <- req(state$cluster_table_data)

            DT::datatable(
                df_clustered %>%
                    dplyr::select(where(~ !is.numeric(.)), where(is.numeric)),
                extensions = "Buttons",
                filter     = "top",
                options    = list(
                    scrollX     = TRUE,
                    processing  = TRUE,
                    pageLength  = 10,
                    lengthMenu  = c(5, 10, 25, 50, 100),
                    dom         = "Blfrtip",
                    buttons     = c("copy", "csv", "excel", "pdf", "print")
                )
            )
        })


        output$optimal_k <- renderUI({
            res <- get_dep_result()
            req(res)
            key <- rv$current_dep_heatmap_key[[tbl_name]]
            if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

            if (is.null(res$optimal_k) || is.na(res$optimal_k)) {
                return(tags$ul(
                    tags$li("Suggested k: not available (not enough clean data after filtering).")
                ))
            }

            tags$ul(tags$li(paste("Suggested k (from filtered matrix):", res$optimal_k)))
        })
    })
}
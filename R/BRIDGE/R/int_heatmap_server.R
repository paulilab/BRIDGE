#' @export
int_heatmap_server <- function(input, output, session, rv) {
  fallback_notified <- shiny::reactiveVal(character(0))

  shiny::observeEvent(input$process_integrate_data, {
    fallback_notified(character(0))
  }, ignoreInit = TRUE)

  output$integrated_heatmaps <- shiny::renderUI({
    shiny::req(rv$intersected_matrix_processed, input$heatmap_k)
    htmltools::tagList(
      shiny::div(
        style = "display:flex; gap:10px; align-items:end; flex-wrap:wrap; margin-bottom:10px;",
        shiny::selectInput(
          "integration_heatmap_download_table",
          "Heatmap dataset",
          choices = names(rv$intersected_matrix_processed),
          selected = names(rv$intersected_matrix_processed)[1],
          width = "260px"
        ),
        shiny::selectInput(
          "integration_heatmap_download_download_format",
          "Download format",
          choices = c("pdf", "svg"),
          selected = "pdf",
          width = "150px"
        ),
        shiny::downloadButton("integration_heatmap_download_download", "Download plot")
      ),
      lapply(names(rv$intersected_matrix_processed), function(tbl) {
        plotOutput(outputId = paste0("heatmap_", tbl), width = "80%")
      }),
      shiny::br(),
      shinydashboard::box(
        title = "Cluster Table",
        width = 12,
        solidHeader = TRUE,
        status = "info",
        collapsible = TRUE,
        collapsed = FALSE,
        DT::DTOutput(outputId = "Cluster_table")
      )
    )
  })

  # Helper: build the exact ID used as heatmap rownames for a given table
  make_unique_id <- function(meta_df, mat, datatype) {
    if (!is.null(datatype) &&
        datatype == "phosphoproteomics" &&
        all(c("Gene_Name", "pepG", "Protein_ID") %in% colnames(meta_df))) {
      as.character(paste(meta_df$Gene_Name, meta_df$pepG, meta_df$Protein_ID, sep = "_"))
    } else if (!is.null(datatype) &&
               datatype == "proteomics" &&
               all(c("Gene_Name", "Protein_ID") %in% colnames(meta_df))) {
      as.character(paste(meta_df$Gene_Name, meta_df$Protein_ID, sep = "_"))
    } else if (!is.null(datatype) &&
               datatype == "rnaseq" &&
               "Gene_Name" %in% colnames(meta_df)) {
      as.character(meta_df$Gene_Name)
    } else {
      # fallback: keep matrix rownames if present
      rn <- rownames(mat)
      if (is.null(rn)) as.character(seq_len(nrow(mat))) else as.character(rn)
    }
  }

  # Ensure ID vector length always matches matrix row count
  coerce_unique_id <- function(meta_df, mat, datatype, tbl_name = NULL) {
    uid <- make_unique_id(meta_df, mat, datatype)
    n_mat <- nrow(mat)

    if (length(uid) != n_mat) {
      n_keep <- min(length(uid), n_mat)
      if (is.null(n_keep) || !is.finite(n_keep) || n_keep < 1) {
        return(rep_len(seq_len(max(1L, n_mat)), n_mat))
      }

      # keep only aligned prefix; this avoids rownames<- length mismatch
      uid <- uid[seq_len(n_keep)]
      if (n_keep < n_mat) {
        uid <- c(uid, paste0("row_", seq.int(n_keep + 1L, n_mat)))
      }

      message(sprintf(
        "[integration heatmap] ID length mismatch%s: meta ids=%d, matrix rows=%d; coercing to %d.",
        if (!is.null(tbl_name)) paste0(" for ", tbl_name) else "",
        length(make_unique_id(meta_df, mat, datatype)),
        n_mat,
        length(uid)
      ))
    }

    uid <- as.character(uid)
    uid[is.na(uid) | !nzchar(uid)] <- paste0("row_", which(is.na(uid) | !nzchar(uid)))
    uid
  }

  # One reactive that computes everything needed for heatmaps + cluster table
  heatmap_state <- shiny::reactive({
    shiny::req(rv$intersected_matrix_processed, rv$intersected_tables_processed, rv$datatype, input$heatmap_k)

    all_mats  <- rv$intersected_matrix_processed
    all_meta  <- rv$intersected_tables_processed

    k <- suppressWarnings(as.integer(input$heatmap_k))
    if (length(k) < 1L || !is.finite(k[[1]]) || k[[1]] < 2L) {
      k <- 2L
    } else {
      k <- as.integer(k[[1]])
    }

    safe_kmeans_cluster <- function(mat, centers) {
      n_rows <- nrow(mat)
      if (!is.finite(n_rows) || n_rows < 1L) {
        return(list(cluster = integer(0), centers_used = 0L, requested = as.integer(centers)))
      }

      requested <- min(max(1L, as.integer(centers)), as.integer(n_rows))
      c_try <- requested
      while (c_try >= 1L) {
        km <- tryCatch(
          stats::kmeans(mat, centers = c_try, nstart = 10, iter.max = 100),
          error = function(e) NULL
        )
        if (!is.null(km) && !is.null(km$cluster) && length(km$cluster) == n_rows) {
          return(list(
            cluster = as.integer(km$cluster),
            centers_used = as.integer(c_try),
            requested = as.integer(requested)
          ))
        }
        c_try <- c_try - 1L
      }

      # Ultimate fallback: one cluster for all rows
      list(
        cluster = rep.int(1L, as.integer(n_rows)),
        centers_used = 1L,
        requested = as.integer(requested)
      )
    }

    res <- list()   # per-table results
    tbls <- intersect(names(all_mats), names(all_meta))

    shiny::req(length(tbls) > 0)

    for (tbl_name in tbls) {
      mat <- all_mats[[tbl_name]]
      meta_df <- all_meta[[tbl_name]]
      datatype <- rv$datatype[[tbl_name]]

      if (is.null(mat) || is.null(meta_df)) next

      # Make sure mat is a numeric matrix (kmeans requires numeric)
      mat <- as.matrix(mat)
      if (length(dim(mat)) != 2L || nrow(mat) < 1L || ncol(mat) < 1L) next
      storage.mode(mat) <- "double"

      # Drop rows with any non-finite values (avoids kmeans hanging/failing silently)
      keep <- apply(mat, 1, function(x) all(is.finite(x)))
      if (!length(keep)) next
      mat <- mat[keep, , drop = FALSE]
      if (nrow(mat) < 1L) next

      if (!(is.data.frame(meta_df) || is.matrix(meta_df))) {
        meta_df <- as.data.frame(meta_df, stringsAsFactors = FALSE)
      }
      if (nrow(meta_df) < 1L) next

      # meta_df can drift from mat row count when integration keys collapse/expand
      if (nrow(meta_df) == length(keep)) {
        meta_df <- meta_df[keep, , drop = FALSE]
      } else {
        n_align <- min(nrow(meta_df), nrow(mat))
        if (n_align < 1) next
        mat <- mat[seq_len(n_align), , drop = FALSE]
        meta_df <- meta_df[seq_len(n_align), , drop = FALSE]
      }

      # Need at least 2 rows to cluster
      if (nrow(mat) < 2) next

      unique_id <- coerce_unique_id(meta_df, mat, datatype, tbl_name = tbl_name)
      rownames(mat) <- unique_id

      # Bound requested k to valid range [1, nrow(mat)] and compute clusters safely.
      k_requested <- as.integer(k)
      k_bounded <- min(max(1L, k_requested), as.integer(nrow(mat)))
      km_res <- safe_kmeans_cluster(mat, centers = k_bounded)
      cluster_assign <- km_res$cluster
      if (!length(cluster_assign)) next

      k_used <- as.integer(km_res$centers_used)
      if (!is.finite(k_used) || k_used < 1L) k_used <- 1L

      if (k_used < k_requested) {
        note_key <- paste(tbl_name, nrow(mat), k_requested, k_bounded, k_used, sep = "|")
        already <- fallback_notified()
        if (!(note_key %in% already)) {
          shiny::showNotification(
            sprintf(
              "Integrated heatmap for '%s': requested k=%d is not feasible for current data (n=%d). Using k=%d.",
              tbl_name, k_requested, nrow(mat), k_used
            ),
            type = "warning",
            duration = 8
          )
          fallback_notified(c(already, note_key))
        }
      }

      present_levels <- sort(unique(as.integer(cluster_assign)))
      cluster_vec <- factor(as.integer(cluster_assign), levels = present_levels)
      ord <- order(cluster_vec)

      mat_ordered <- mat[ord, , drop = FALSE]
      cluster_vec <- cluster_vec[ord]

      cluster_df <- data.frame(
        dataset  = tbl_name,
        unique_id = rownames(mat_ordered),
        Cluster  = as.integer(cluster_vec),
        stringsAsFactors = FALSE
      )

      res[[tbl_name]] <- list(
        mat_ordered = mat_ordered,
        cluster_vec = cluster_vec,
        cluster_df  = cluster_df,
        meta_df     = meta_df,     # filtered meta (same rows as mat before ordering)
        datatype    = datatype,
        k_use       = k_used
      )
    }

    shiny::req(length(res) > 0)
    res
  })

  # Render each heatmap from computed state
  shiny::observe({
    state <- heatmap_state()

    lapply(names(state), function(tbl_name) {
      local({
        nm <- tbl_name
        output[[paste0("heatmap_", nm)]] <- shiny::renderPlot({
          x <- state[[nm]]
          mat_ordered <- x$mat_ordered
          cluster_vec <- x$cluster_vec
          k_use <- x$k_use

          # Trend lines per cluster (safe)
          cluster_ids <- split(seq_len(nrow(mat_ordered)), cluster_vec)
          line_profiles <- t(vapply(cluster_ids, function(idxs) {
            colMeans(mat_ordered[idxxs <- idxs, , drop = FALSE], na.rm = TRUE)
          }, FUN.VALUE = numeric(ncol(mat_ordered))))

          if (nrow(line_profiles) == 0 || ncol(line_profiles) == 0) return(NULL)

          line_profiles_norm <- t(apply(line_profiles, 1, function(v) {
            rng <- range(v, na.rm = TRUE)
            if (!all(is.finite(rng)) || diff(rng) == 0) rep(0.5, length(v)) else (v - rng[1]) / diff(rng)
          }))

          trend_anno <- ComplexHeatmap::rowAnnotation(
            trend = ComplexHeatmap::anno_link(
              align_to = cluster_vec,
              which = "row",
              panel_fun = function(index, nm2) {
                grid::grid.rect()
                grid::grid.lines(
                  x = seq_len(ncol(line_profiles_norm)) / ncol(line_profiles_norm),
                  y = line_profiles_norm[as.integer(nm2), ],
                  gp = ggfun::gpar(col = BRIDGE_COLORS$accent, lwd = 1)
                )
              },
              side = "right",
              size = unit(3, "cm"),
              width = unit(5, "cm")
            )
          )

          ComplexHeatmap::Heatmap(
            mat_ordered,
            name = nm,
            col = bridge_heatmap_col(max(abs(mat_ordered), na.rm = TRUE)),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = FALSE,
            show_column_dend = TRUE,
            show_row_names = FALSE,
            row_split = cluster_vec,
            row_names_gp = ggfun::gpar(fontsize = 6),
            column_names_gp = ggfun::gpar(fontsize = 8),
            heatmap_legend_param = list(title = "Expression"),
            right_annotation = trend_anno
          )
        })
      })
    })
  })

  # Render cluster table (directly from unique_id + cluster, joined to metadata if desired)
  output$Cluster_table <- DT::renderDT({
    state <- heatmap_state()

    # Build a combined table; join meta by unique_id (no gene mapping)
    df_list <- lapply(names(state), function(tbl_name) {
      x <- state[[tbl_name]]
      meta_df <- x$meta_df
      datatype <- x$datatype

      # Recompute unique_id for meta_df to match clustering ids (same function, same filtering)
      unique_id <- coerce_unique_id(meta_df, x$mat_ordered, datatype, tbl_name = tbl_name)
      meta_df$unique_id <- as.character(unique_id)

      joined <- dplyr::inner_join(
        meta_df,
        x$cluster_df,
        by = c("unique_id" = "unique_id")
      )

      # prefix cluster with dataset (your original UI convention)
      joined$Cluster <- paste0(tbl_name, "_", joined$Cluster)

      # Columns to show
      if (datatype == "rnaseq") {
        cols_to_show <- c("Gene_Name", "Gene_ID", "Cluster")
      } else if (datatype == "phosphoproteomics") {
        cols_to_show <- c("Gene_Name", "Gene_ID", "pepG", "Protein_ID", "Cluster")
      } else if (datatype == "proteomics") {
        cols_to_show <- c("Gene_Name", "Gene_ID", "Protein_ID", "Cluster")
      } else {
        cols_to_show <- c("unique_id", "Cluster")
      }

      joined %>%
        dplyr::select(dplyr::any_of(cols_to_show)) %>%
        dplyr::distinct()
    })

    final_df <- dplyr::bind_rows(df_list) %>%
      dplyr::arrange(Cluster) %>%
      dplyr::distinct()

    DT::datatable(
      final_df,
      extensions = "Buttons",
      filter = "top",
      options = list(
        scrollX = TRUE,
        pageLength = 5,
        lengthMenu = c(5, 10, 25, 50, 100),
        dom = "Blfrtip",
        buttons = list(
          list(extend = "copy", exportOptions = list(modifier = list(page = "all"))),
          list(extend = "csv",  exportOptions = list(modifier = list(page = "all"))),
          list(extend = "excel",exportOptions = list(modifier = list(page = "all"))),
          list(extend = "pdf",  exportOptions = list(modifier = list(page = "all"))),
          "print"
        )
      )
    )
  })

  output$lfc_scatter_selector <- shiny::renderUI({
    shiny::req(rv$scatter_plots)
    shiny::selectInput(
      "scatter_comparisons",
      "Select comparison:",
      choices = names(rv$scatter_plots)[!grepl("phosphoproteomics", names(rv$scatter_plots))]
    )
  })

  output$lfc_scatter_ui <- shiny::renderUI({
    shiny::req(rv$intersected_tables_processed, rv$scatter_plots)
    if (length(names(rv$scatter_plots)[!grepl("phosphoproteomics", names(rv$scatter_plots))]) == 0) {
      shiny::div(
        style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
        "Scatter plot is not available for this integration."
      )
    } else {
      plotly::plotlyOutput("lfc_scatter_plot", height = "400px", width = "100%")
    }
  })

  build_native_lfc_scatter <- function(df_plot, comparison_name = NULL) {
    shiny::req(is.data.frame(df_plot), all(c("x", "y") %in% colnames(df_plot)))

    plot_df <- as.data.frame(df_plot)
    plot_df <- plot_df[is.finite(plot_df$x) & is.finite(plot_df$y), , drop = FALSE]
    shiny::req(nrow(plot_df) > 0)

    hover_text <- if ("Gene_Name" %in% colnames(plot_df)) {
      as.character(plot_df$Gene_Name)
    } else {
      sprintf("x: %.3f<br>y: %.3f", plot_df$x, plot_df$y)
    }

    axis_rng <- range(c(plot_df$x, plot_df$y), na.rm = TRUE)
    if (!all(is.finite(axis_rng))) axis_rng <- c(-1, 1)
    if (diff(axis_rng) == 0) axis_rng <- axis_rng + c(-0.5, 0.5)

    axis_labels <- c("LFC X", "LFC Y")
    if (!is.null(comparison_name) && grepl("_vs_", comparison_name, fixed = TRUE)) {
      parts <- strsplit(comparison_name, "_vs_", fixed = TRUE)[[1]]
      if (length(parts) == 2) axis_labels <- parts
    }

    p_native <- plotly::plot_ly(
      data = plot_df,
      x = ~x,
      y = ~y,
      type = "scattergl",
      mode = "markers",
      text = hover_text,
      hoverinfo = "text",
      marker = list(
        color = BRIDGE_COLORS$accent,
        size = 7,
        opacity = 0.7
      )
    )

    p_native <- plotly::add_segments(
      p_native,
      x = axis_rng[1], xend = axis_rng[2],
      y = axis_rng[1], yend = axis_rng[2],
      inherit = FALSE,
      line = list(color = "#AAAAAA", dash = "dash"),
      hoverinfo = "skip",
      showlegend = FALSE
    )

    p_native |>
      plotly::layout(
        title = comparison_name,
        xaxis = list(title = axis_labels[1], range = axis_rng),
        yaxis = list(title = axis_labels[2], range = axis_rng, scaleanchor = "x", scaleratio = 1)
      )
  }

  output$lfc_scatter_plot <- plotly::renderPlotly({
    shiny::req(rv$scatter_plots, input[["scatter_comparisons"]])

    comparison <- input[["scatter_comparisons"]]
    p_obj <- rv$scatter_plots[[comparison]]
    shiny::req(!is.null(p_obj))

    tryCatch({
      plotly::ggplotly(p_obj, tooltip = "text", dynamicTicks = FALSE)
    }, error = function(e) {
      df_plot <- p_obj$data
      if (!is.data.frame(df_plot) || !all(c("x", "y") %in% colnames(df_plot))) {
        shiny::showNotification(
          paste("Failed to render scatter:", conditionMessage(e)),
          type = "error",
          duration = NULL
        )
        return(plotly::plot_ly())
      }

      shiny::showNotification(
        "Using native scatter renderer due to ggplotly conversion error.",
        type = "message",
        duration = 4
      )
      build_native_lfc_scatter(df_plot, comparison_name = comparison)
    })
  })

  output$lfc_scatter_download_ui <- shiny::renderUI({
    if (is.null(rv$scatter_plots) || !length(rv$scatter_plots)) return(NULL)
    available <- names(rv$scatter_plots)[!grepl("phosphoproteomics", names(rv$scatter_plots))]
    if (!length(available)) return(NULL)
    plot_download_controls(session$ns, "lfc_scatter")
  })

  register_plot_download(
    input = input,
    output = output,
    session = session,
    id_prefix = "lfc_scatter",
    filename_prefix = "integration_lfc_scatter",
    plot_fun = function() {
      shiny::req(rv$scatter_plots, input[["scatter_comparisons"]])
      rv$scatter_plots[[input[["scatter_comparisons"]]]]
    },
    width = 8,
    height = 6
  )

  output$integration_heatmap_download_download <- shiny::downloadHandler(
    filename = function() {
      fmt <- input$integration_heatmap_download_download_format
      if (is.null(fmt) || !(fmt %in% c("pdf", "svg"))) fmt <- "pdf"
      tbl <- input$integration_heatmap_download_table
      if (is.null(tbl) || !nzchar(tbl)) tbl <- "integration"
      paste0("integration_heatmap_", tbl, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", fmt)
    },
    content = function(file) {
      fmt <- input$integration_heatmap_download_download_format
      if (is.null(fmt) || !(fmt %in% c("pdf", "svg"))) fmt <- "pdf"
      if (identical(fmt, "svg")) {
        grDevices::svg(file, width = 10, height = 8)
      } else {
        grDevices::pdf(file, width = 10, height = 8, onefile = FALSE)
      }
      on.exit(grDevices::dev.off(), add = TRUE)

      state <- heatmap_state()
      tbl <- input$integration_heatmap_download_table
      shiny::req(tbl, tbl %in% names(state))

      x <- state[[tbl]]

      cluster_ids <- split(seq_len(nrow(x$mat_ordered)), x$cluster_vec)
      line_profiles <- t(vapply(cluster_ids, function(idxs) {
        colMeans(x$mat_ordered[idxs, , drop = FALSE], na.rm = TRUE)
      }, FUN.VALUE = numeric(ncol(x$mat_ordered))))

      trend_anno <- NULL
      if (nrow(line_profiles) > 0 && ncol(line_profiles) > 0) {
        line_profiles_norm <- t(apply(line_profiles, 1, function(v) {
          rng <- range(v, na.rm = TRUE)
          if (!all(is.finite(rng)) || diff(rng) == 0) rep(0.5, length(v)) else (v - rng[1]) / diff(rng)
        }))
        trend_anno <- ComplexHeatmap::rowAnnotation(
          trend = ComplexHeatmap::anno_link(
            align_to = x$cluster_vec,
            which = "row",
            panel_fun = function(index, nm2) {
              grid::grid.rect()
              grid::grid.lines(
                x = seq_len(ncol(line_profiles_norm)) / ncol(line_profiles_norm),
                y = line_profiles_norm[as.integer(nm2), ],
                gp = ggfun::gpar(col = BRIDGE_COLORS$accent, lwd = 1)
              )
            },
            side = "right",
            size = unit(3, "cm"),
            width = unit(5, "cm")
          )
        )
      }

      ComplexHeatmap::draw(
        ComplexHeatmap::Heatmap(
          x$mat_ordered,
          name = tbl,
          col = bridge_heatmap_col(max(abs(x$mat_ordered), na.rm = TRUE)),
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_dend = FALSE,
          show_column_dend = TRUE,
          show_row_names = FALSE,
          row_split = x$cluster_vec,
          row_names_gp = ggfun::gpar(fontsize = 6),
          column_names_gp = ggfun::gpar(fontsize = 8),
          heatmap_legend_param = list(title = "Expression"),
          right_annotation = trend_anno
        ),
        merge_legend = TRUE,
        newpage = FALSE
      )
    }
  )
}
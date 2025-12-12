#' @export
int_heatmap_server <- function(input, output, session, rv) {

  output$integrated_heatmaps <- shiny::renderUI({
    shiny::req(rv$intersected_matrix_processed, input$heatmap_k)
    htmltools::tagList(
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

  # One reactive that computes everything needed for heatmaps + cluster table
  heatmap_state <- shiny::reactive({
    shiny::req(rv$intersected_matrix_processed, rv$intersected_tables_processed, rv$datatype, input$heatmap_k)

    all_mats  <- rv$intersected_matrix_processed
    all_meta  <- rv$intersected_tables_processed

    k <- suppressWarnings(as.integer(input$heatmap_k))
    if (!is.finite(k) || k < 2) k <- 2L

    res <- list()   # per-table results
    tbls <- intersect(names(all_mats), names(all_meta))

    shiny::req(length(tbls) > 0)

    for (tbl_name in tbls) {
      mat <- all_mats[[tbl_name]]
      meta_df <- all_meta[[tbl_name]]
      datatype <- rv$datatype[[tbl_name]]

      # Make sure mat is a numeric matrix (kmeans requires numeric)
      mat <- as.matrix(mat)
      storage.mode(mat) <- "double"

      # Drop rows with any non-finite values (avoids kmeans hanging/failing silently)
      keep <- apply(mat, 1, function(x) all(is.finite(x)))
      mat <- mat[keep, , drop = FALSE]
      meta_df <- meta_df[keep, , drop = FALSE]

      # Need at least 2 rows to cluster
      if (nrow(mat) < 2) next

      unique_id <- make_unique_id(meta_df, mat, datatype)
      rownames(mat) <- unique_id

      k_max <- max(2L, nrow(mat) - 1L)
      k_use <- min(k, k_max)
      if (k_use < 2L) next

      km <- stats::kmeans(mat, centers = k_use, nstart = 10, iter.max = 100)

      cluster_vec <- factor(km$cluster, levels = 1:k_use)
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
        k_use       = k_use
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
                  gp = ggfun::gpar(col = "#2b8cbe", lwd = 1)
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
      unique_id <- make_unique_id(meta_df, x$mat_ordered, datatype)
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
        buttons = c("copy", "csv", "excel", "pdf", "print")
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

  output$lfc_scatter_plot <- plotly::renderPlotly({
    shiny::req(rv$scatter_plots, input[["scatter_comparisons"]])
    plotly::ggplotly(rv$scatter_plots[[input[["scatter_comparisons"]]]], tooltip = "text")
  })
}
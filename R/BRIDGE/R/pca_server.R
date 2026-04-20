#' @export
pcaServer <- function(id, rv, tbl_name) {
    moduleServer(id, function(input, output, session) {
        pca_task <- ExtendedTask$new(function(dep) {
            promises::future_promise({
                mat <- SummarizedExperiment::assay(dep)
                vars <- matrixStats::rowVars(mat)
                select <- head(order(vars, decreasing = TRUE), min(500, length(vars)))

                pca <- stats::prcomp(t(mat[select, , drop = FALSE]), center = TRUE, scale. = FALSE)
                var_expl <- pca$sdev^2 / sum(pca$sdev^2)

                # PCA score data for plotting
                score_df <- as.data.frame(pca$x)
                score_df$sample <- rownames(score_df)

                cd <- as.data.frame(SummarizedExperiment::colData(dep))
                if ("condition" %in% colnames(cd)) {
                    score_df$group <- as.factor(cd$condition)
                } else if (ncol(cd) > 0) {
                    score_df$group <- as.factor(cd[[1]])
                } else {
                    score_df$group <- as.factor(score_df$sample)
                }

                # Top-contribution table for PC loadings
                loadings_df <- as.data.frame(pca$rotation) |>
                    rownames_to_column("gene") |>
                    pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loading") |>
                    mutate(
                        abs_loading = abs(loading),
                        pc_idx = as.integer(sub("PC", "", PC)),
                        var_explained = var_expl[pc_idx]
                    )

                contrib <- loadings_df %>%
                    group_by(PC) %>%
                    mutate(contribution = (loading^2) / sum(loading^2)) %>%
                    ungroup()

                top_contrib <- contrib %>%
                    filter(PC == "PC1" | PC == "PC2") %>%
                    arrange(desc(contribution))

                list(
                    score_df = score_df,
                    var_expl = var_expl,
                    top_contrib = as.data.frame(top_contrib)
                )
            })
        })

        observeEvent(input$compute,
            {
                req(rv$dep_output[[tbl_name]])
                pca_task$invoke(isolate(rv$dep_output[[tbl_name]]))
            },
            ignoreInit = TRUE
        )

        pca_plot <- reactive({
            res <- pca_task$result()
            req(res)
            score_df <- res$score_df
            req(is.data.frame(score_df), "PC1" %in% colnames(score_df), "PC2" %in% colnames(score_df))

            var_expl <- res$var_expl
            x_lab <- sprintf("PC1 (%.1f%%)", 100 * var_expl[1])
            y_lab <- sprintf("PC2 (%.1f%%)", 100 * var_expl[2])

            ggplot2::ggplot(score_df, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data$group, text = .data$sample)) +
                ggplot2::geom_point(size = 3, alpha = 0.9) +
                ggplot2::labs(x = x_lab, y = y_lab, color = "Group") +
                ggplot2::theme_minimal()
        })

        output$plot_slot <- renderUI({
            if (isTRUE(input$interactive)) {
                plotly::plotlyOutput(session$ns("plotly"), height = "480px")
            } else {
                plotOutput(session$ns("plot"), height = "480px")
            }
        })

        output$plot <- renderPlot({
            pca_plot() +
                ggplot2::ggtitle(paste("PCA for", tbl_name)) +
                ggplot2::theme_minimal()
        })

        output$plotly <- plotly::renderPlotly({
            plotly::ggplotly(
                pca_plot() + ggplot2::ggtitle(paste("PCA for", tbl_name)),
                tooltip = c("text", "x", "y", "colour")
            )
        })

        output$plot_download_ui <- renderUI({
            plot_download_controls(session$ns, "pca")
        })

        register_plot_download(
            input = input,
            output = output,
            session = session,
            id_prefix = "pca",
            filename_prefix = paste0("pca_", tbl_name),
            plot_fun = function() {
                pca_plot() + ggplot2::ggtitle(paste("PCA for", tbl_name))
            },
            width = 9,
            height = 6
        )

        output$pcs_panel <- renderUI({
            DT::DTOutput(session$ns("pcs"), height = "300px")
        })

        output$pcs <- DT::renderDT({
            res <- pca_task$result()
            top_contrib <- if (!is.null(res)) res$top_contrib else NULL
            if (is.null(top_contrib)) {
                validate(need(FALSE, htmltools::tagList(
                    tags$span(class = "lds-ring", tags$div(), tags$div(), tags$div(), tags$div()),
                    " Loading PCA loadings…"
                )))
            }
            req(!is.null(top_contrib))
            DT::datatable(top_contrib %>% dplyr::select(where(~!is.numeric(.)), where(is.numeric)),
                extensions = "Buttons",
                filter = "top",
                options = list(
                    scrollX = TRUE, processing = TRUE, pageLength = 10, lengthMenu = c(5, 10, 25, 50, 100), dom = "Blfrtip", buttons = c("copy", "csv", "excel", "pdf", "print")
                )
            )
        })
    })
}

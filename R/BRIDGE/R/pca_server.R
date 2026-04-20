#' @export
pcaServer <- function(id, rv, cache, tbl_name) {
    moduleServer(id, function(input, output, session) {
        pca_task <- ExtendedTask$new(function(mat, group, ntop) {
            promises::future_promise({
                vars <- matrixStats::rowVars(mat)
                ntop <- as.integer(ntop)
                if (!is.finite(ntop) || ntop < 2L) ntop <- 500L
                select <- head(order(vars, decreasing = TRUE), min(ntop, length(vars)))
                xmat <- t(mat[select, , drop = FALSE])

                # We never need more than 10 PCs in the UI; cap for speed.
                max_rank <- min(10L, max(1L, nrow(xmat) - 1L), ncol(xmat))
                pca <- stats::prcomp(xmat, center = TRUE, scale. = FALSE, rank. = max_rank)
                var_expl <- pca$sdev^2 / sum(pca$sdev^2)

                # PCA score data for plotting
                score_df <- as.data.frame(pca$x)
                score_df$sample <- rownames(score_df)

                if (is.null(group) || length(group) != nrow(score_df)) {
                    score_df$group <- as.factor(score_df$sample)
                } else {
                    score_df$group <- as.factor(group)
                }

                list(
                    score_df = score_df,
                    var_expl = var_expl,
                    rotation = as.data.frame(pca$rotation)
                )
            })
        })

        pca_result_cache <- reactiveVal(NULL)
        pca_cache_key <- reactiveVal(NULL)
        loadings_cache <- reactiveVal(NULL)

        get_top_n <- function() {
            n <- suppressWarnings(as.integer(input$top_n_features))
            if (!is.finite(n) || n < 2L) n <- 500L
            n
        }

        get_pca_key <- function(n_top) {
            cols <- isolate(rv$data_cols[[tbl_name]])
            cols_key <- if (is.null(cols) || !length(cols)) "no_cols" else paste(sort(trimws(cols)), collapse = "_")
            paste(
                tbl_name,
                cols_key,
                paste0("ntop=", as.integer(n_top)),
                "rank=10",
                "pca_v1",
                sep = "_"
            )
        }

        get_pca_result <- function() {
            res_cached <- pca_result_cache()
            if (!is.null(res_cached)) return(res_cached)

            res_task <- pca_task$result()
            if (!is.null(res_task)) {
                key <- pca_cache_key()
                if (!is.null(cache) && !is.null(key) && !cache$exists(key)) {
                    cache$set(key, res_task)
                }
                pca_result_cache(res_task)
            }
            res_task
        }

        observeEvent(input$compute,
            {
                req(rv$dep_output[[tbl_name]])
                n_top <- get_top_n()
                key <- get_pca_key(n_top)
                pca_cache_key(key)
                dep <- isolate(rv$dep_output[[tbl_name]])
                mat <- SummarizedExperiment::assay(dep)
                cd <- as.data.frame(SummarizedExperiment::colData(dep))
                grp <- if ("condition" %in% colnames(cd)) {
                    cd$condition
                } else if (ncol(cd) > 0) {
                    cd[[1]]
                } else {
                    colnames(mat)
                }

                if (!is.null(cache) && cache$exists(key)) {
                    pca_result_cache(cache$get(key))
                } else {
                    pca_result_cache(NULL)
                    pca_task$invoke(mat, grp, n_top)
                }
                loadings_cache(NULL)
            },
            ignoreInit = TRUE
        )

        output$pc_axes_ui <- renderUI({
            res <- get_pca_result()
            if (is.null(res) || is.null(res$score_df)) {
                return(div(style = "padding-top: 6px; color: #777;", "Compute PCA to select axes."))
            }
            pc_choices <- grep("^PC[0-9]+$", colnames(res$score_df), value = TRUE)
            if (length(pc_choices) < 2) {
                return(div(style = "padding-top: 6px; color: #777;", "Need at least 2 PCs to select axes."))
            }

            x_sel <- if (!is.null(input$x_pc) && input$x_pc %in% pc_choices) input$x_pc else pc_choices[1]
            y_default <- if (length(pc_choices) >= 2) pc_choices[2] else pc_choices[1]
            y_sel <- if (!is.null(input$y_pc) && input$y_pc %in% pc_choices) input$y_pc else y_default

            fluidRow(
                column(6, selectInput(session$ns("x_pc"), "X axis", choices = pc_choices, selected = x_sel)),
                column(6, selectInput(session$ns("y_pc"), "Y axis", choices = pc_choices, selected = y_sel))
            )
        })

        pca_plot <- reactive({
            res <- get_pca_result()
            req(res)
            score_df <- res$score_df
            req(is.data.frame(score_df))

            x_pc <- req(input$x_pc)
            y_pc <- req(input$y_pc)
            validate(need(x_pc != y_pc, "Please choose two different PCs for X and Y axes."))
            req(x_pc %in% colnames(score_df), y_pc %in% colnames(score_df))

            var_expl <- res$var_expl
            x_idx <- as.integer(sub("PC", "", x_pc))
            y_idx <- as.integer(sub("PC", "", y_pc))
            req(is.finite(x_idx), is.finite(y_idx), x_idx <= length(var_expl), y_idx <= length(var_expl))
            x_lab <- sprintf("%s (%.1f%%)", x_pc, 100 * var_expl[x_idx])
            y_lab <- sprintf("%s (%.1f%%)", y_pc, 100 * var_expl[y_idx])

            ggplot2::ggplot(score_df, ggplot2::aes(x = .data[[x_pc]], y = .data[[y_pc]], color = .data$group, text = .data$sample)) +
                ggplot2::geom_point(size = 3, alpha = 0.9) +
                ggplot2::labs(x = x_lab, y = y_lab, color = "Group") +
                ggplot2::theme_minimal()
        })

        output$plot <- renderPlot({
            pca_plot() +
                ggplot2::ggtitle(paste("PCA for", tbl_name)) +
                ggplot2::theme_minimal()
        })

        output$plot_download_ui <- renderUI({
            res <- get_pca_result()
            if (is.null(res) || is.null(res$score_df)) return(NULL)
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
            top_contrib <- loadings_cache()
            if (is.null(top_contrib)) {
                res <- get_pca_result()
                req(res, res$rotation, res$var_expl)

                loadings_df <- as.data.frame(res$rotation) |>
                    rownames_to_column("gene") |>
                    pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loading") |>
                    mutate(
                        abs_loading = abs(loading),
                        pc_idx = as.integer(sub("PC", "", PC)),
                        var_explained = res$var_expl[pc_idx]
                    )

                top_contrib <- loadings_df %>%
                    group_by(PC) %>%
                    mutate(contribution = (loading^2) / sum(loading^2)) %>%
                    arrange(desc(contribution), .by_group = TRUE) %>%
                    slice_head(n = 200) %>%
                    ungroup() %>%
                    as.data.frame()

                loadings_cache(top_contrib)
            }

            if (is.null(top_contrib)) {
                validate(need(FALSE, htmltools::tagList(
                    tags$span(class = "lds-ring", tags$div(), tags$div(), tags$div(), tags$div()),
                    " Loading PCA loadings…"
                )))
            }
            req(!is.null(top_contrib))
            selected_pcs <- unique(stats::na.omit(c(input$x_pc, input$y_pc)))
            if (length(selected_pcs)) {
                top_contrib <- top_contrib %>% dplyr::filter(PC %in% selected_pcs)
            }
            top_contrib <- top_contrib %>% dplyr::arrange(PC, dplyr::desc(contribution))

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

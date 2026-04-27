#' @export
EnrichmentServer <- function(id, rv, tbl_name) {
    moduleServer(id, function(input, output, session) {
        enrichment_plot_obj <- reactiveVal(NULL)
        enrichment_dep_ready <- reactiveVal(FALSE)
        enrichment_dep_error <- reactiveVal(NULL)

        output$plot_download_ui <- renderUI({
            if (is.null(enrichment_plot_obj())) return(NULL)
            plot_download_controls(session$ns, "enrichment")
        })

        register_plot_download(
            input = input,
            output = output,
            session = session,
            id_prefix = "enrichment",
            filename_prefix = paste0("enrichment_", tbl_name),
            plot_fun = function() {
                req(!is.null(enrichment_plot_obj()))
                enrichment_plot_obj()
            },
            width = 9,
            height = 6
        )

        # Populate contrast choices when available
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

        check_sig <- function(dep_obj) {
            if (methods::is(dep_obj, "DEGdata")) {
                tr <- as.data.frame(dep_obj@test_result)
                check <- nrow(tr %>% dplyr::filter(significant == TRUE))
            } else {
                rd <- as.data.frame(SummarizedExperiment::rowData(dep_obj))
                check <- nrow(rd %>% dplyr::filter(significant == TRUE))
            }
            # message("Significant hits: ", check)
            check
        }

        observe({
            req(rv$contrasts[[tbl_name]])
            updateSelectInput(session, "contrasts_enrichment",
                choices = rv$contrasts[[tbl_name]]
            )
        })

        # Render contrast select UI (so we can safely update it)
        output$contrast_ui <- renderUI({
            selectInput(session$ns("contrasts_enrichment"),
                "Contrast",
                choices = rv$contrasts[[tbl_name]]
            )
        })

        # Preflight dependencies when species context is available.
        # This keeps package/org-db installation out of the compute click path.
        observeEvent(rv$species[[tbl_name]], {
            req(rv$species[[tbl_name]])

            species <- stringr::str_to_title(rv$species[[tbl_name]])
            if (!(species %in% c("Zebrafish", "Human", "Mouse"))) {
                species <- "Human"
            }

            enrichment_dep_ready(FALSE)
            enrichment_dep_error(NULL)

            tryCatch({
                DEP2::check_enrichment_depends()
                DEP2::check_organismDB_depends(species, install = TRUE)
                enrichment_dep_ready(TRUE)
            }, error = function(e) {
                enrichment_dep_error(conditionMessage(e))
                enrichment_dep_ready(FALSE)
                showNotification(
                    paste0(
                        "Enrichment dependencies are not ready for species '", species,
                        "'. Please check package/org-db installation.\nDetails: ", conditionMessage(e)
                    ),
                    type = "error",
                    duration = NULL
                )
            })
        }, ignoreInit = FALSE)

        observeEvent(input$compute_enrichment,
            {
                req(rv$dep_output[[tbl_name]], rv$datatype[[tbl_name]], rv$species[[tbl_name]])

                # Phospho: show message and stop
                if (identical(rv$datatype[[tbl_name]], "phosphoproteomics")) {
                    enrichment_plot_obj(NULL)
                    output$enrichment <- renderUI({
                        div(
                            style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                            "No enrichment plot available for phosphoproteomics datasets."
                        )
                    })
                    return(invisible())
                }

                if (!isTRUE(enrichment_dep_ready())) {
                    detail <- enrichment_dep_error()
                    output$enrichment <- renderUI({
                        div(
                            style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                            if (!is.null(detail) && nzchar(detail)) {
                                paste0("Enrichment dependencies are not ready yet: ", detail)
                            } else {
                                "Enrichment dependencies are still being prepared. Please retry in a moment."
                            }
                        )
                    })
                    return(invisible())
                }

                species <- stringr::str_to_title(rv$species[[tbl_name]])
                enrichment_type <- input$comparison_db
                contrast <- input$contrasts_enrichment
                dep_output <- rv$dep_output[[tbl_name]]

                # Fallback for unsupported organisms
                Gene_Names <- stringr::str_to_lower(gsub("_.*$", "", rownames(dep_output)))
                #rownames(dep_output) <- stringr::str_to_upper(rownames(dep_output))
                rownames(dep_output) <- Gene_Names
                #message("NAMES: ", paste(rownames(dep_output)[1:5], collapse = ", "))

                if (!(species %in% c("Zebrafish", "Human", "Mouse"))) {                 
                    species <- "Human"
                }

                dep_output <- strip_sig(dep_output) # Remove old sig cols if present
                # Sig filtering used in your original code before ORA
                dep_pg <- DEP2::add_rejections(dep_output, alpha = input$enrichment_pcutoff, lfc = input$enrichment_fccutoff)                
                sig_count <- check_sig(dep_pg)
                #message("Significant hits for enrichment: ", sig_count)
                if (sig_count < 10) {                
                    enrichment_plot_obj(NULL)
                    output$enrichment <- renderUI({
                        div(
                            style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                            paste0("Not enough significant hits (n < 10) for enrichment analysis with contrast '", contrast, "'.")
                        )
                    })
                    return(invisible())
                } else{
                    message("Performing ORA enrichment analysis with contrast '", contrast, "' species '", species, "' and type '", enrichment_type, "' on '", sig_count, "' hits like '",paste(rownames(dep_pg)[1:5], collapse = ',') ,"' ...")
                    
                    # Try/catch for DEP2::test_ORA
                    res_ora <- tryCatch({
                        withCallingHandlers(
                            DEP2::test_ORA(
                                dep_pg,
                                type = enrichment_type,
                                species = species,
                                contrasts = contrast
                            ),
                            warning = function(w) {
                                if (grepl("no applicable method for ", conditionMessage(w))) {
                                    output$enrichment <- renderUI({
                                        div(
                                            style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                                            "No genes matched with annotation in the selected organism database! Please check your cutoffs or try another contrast."
                                        )
                                    })
                                    invokeRestart("muffleWarning")
                                }
                            }
                        )
                    }, error = function(e) {
                        if (grepl("no applicable method for ", conditionMessage(e))) {
                            output$enrichment <- renderUI({
                                div(
                                    style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                                    "No genes matched with annotation in the selected organism database! Please check your cutoffs or try another contrast."
                                )
                            })
                            return(NULL)
                        } else {
                            stop(e)
                        }
                    })

                    if (!is.null(res_ora)) {
                        enrichment_plot_obj(enrichplot::dotplot(res_ora))
                        output$enrichment <- renderUI({
                            plotOutput(session$ns("enrichment_plot"), height = "520px")
                        })
                        output$enrichment_plot <- renderPlot({
                            req(!is.null(enrichment_plot_obj()))
                            enrichment_plot_obj()
                        })
                    }
                }
            },
            ignoreInit = TRUE
        )
    })
}

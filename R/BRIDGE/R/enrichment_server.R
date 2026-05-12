#' @export
EnrichmentServer <- function(id, rv, cache, tbl_name) {
    moduleServer(id, function(input, output, session) {
        enrichment_plot_obj <- reactiveVal(NULL)
        enrichment_dep_ready <- reactiveVal(FALSE)
        enrichment_dep_error <- reactiveVal(NULL)

        make_cache_key <- function(...) {
            paste("enrich", tbl_name, ..., sep = "|")
        }

        # Helper function to auto-detect best gene name casing for any species
        detect_best_gene_case <- function(gene_names, species) {
            # Map species to organism database name
            org_db_map <- list(
                "Human" = "org.Hs.eg.db",
                "Mouse" = "org.Mm.eg.db",
                "Zebrafish" = "org.Dr.eg.db"
            )
            org_db_name <- org_db_map[[species]]
            
            # Fallback for unknown species
            if (is.null(org_db_name)) {
                message("Species '", species, "' not in predefined list; attempting to detect best case...")
                org_db_name <- "org.Hs.eg.db"  # Use human as last resort
            }
            
            # cache case-detection result per table/species/data-cols fingerprint
            data_cols_fp <- paste(sort(rv$data_cols[[tbl_name]]), collapse = ",")
            case_cache_key <- make_cache_key("gene_case", species, data_cols_fp)
            if (!is.null(cache) && cache$exists(case_cache_key)) {
                best_case <- cache$get(case_cache_key)
                message("Using cached gene name case: '", best_case, "'.")
                return(switch(best_case,
                    uppercase = function(x) stringr::str_to_upper(x),
                    lowercase = function(x) stringr::str_to_lower(x),
                    titlecase = function(x) stringr::str_to_title(x),
                    function(x) x
                ))
            }

            # Sample genes to test (up to 100)
            sample_n <- min(100, length(gene_names))
            if (sample_n == 0) {
                return(function(x) x)
            }
            sample_genes <- gene_names[seq_len(sample_n)]
            
            # Try different cases
            cases_to_test <- list(
                uppercase = stringr::str_to_upper(sample_genes),
                lowercase = stringr::str_to_lower(sample_genes),
                titlecase = stringr::str_to_title(sample_genes)
            )
            
            # Test each case by attempting ENTREZID mapping
            match_rates <- sapply(cases_to_test, function(genes) {
                tryCatch({
                    org_db <- get(org_db_name)
                    mapped <- AnnotationDbi::mapIds(org_db, keys = genes, column = "ENTREZID",
                                                    keytype = "SYMBOL", multiVals = "first")
                    sum(!is.na(mapped))
                }, error = function(e) 0)
            })
            
            # Return the case that gives best match rate
            best_case <- names(which.max(match_rates))
            best_rate <- max(match_rates)
            message("Auto-detected best gene name case: '", best_case, "' (", best_rate, " / ", length(sample_genes), " genes matched)")

            if (!is.null(cache) && !cache$exists(case_cache_key)) {
                cache$set(case_cache_key, best_case)
            }
            
            # Return the transformation function
            switch(best_case,
                uppercase = function(x) stringr::str_to_upper(x),
                lowercase = function(x) stringr::str_to_lower(x),
                titlecase = function(x) stringr::str_to_title(x),
                function(x) x  # fallback: no transformation
            )
        }

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

                if (is.null(contrast) || !nzchar(contrast)) {
                    enrichment_plot_obj(NULL)
                    output$enrichment <- renderUI({
                        div(
                            style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                            "Please select a contrast before running enrichment."
                        )
                    })
                    return(invisible())
                }

                output$enrichment <- renderUI({
                    div(
                        style = "padding: 20px; color: #337ab7; font-weight: bold; text-align: center;",
                        "Running enrichment analysis..."
                    )
                })

                # Extract gene names and auto-detect best case for this organism
                Gene_Names <- gsub("_.*$", "", rownames(dep_output))
                
                # Auto-detect which case (uppercase, lowercase, titlecase) gives best matches
                case_transform <- detect_best_gene_case(Gene_Names, species)
                Gene_Names <- case_transform(Gene_Names)
                
                rownames(dep_output) <- Gene_Names
                #message("NAMES: ", paste(rownames(dep_output)[1:5], collapse = ", "))

                if (!(species %in% c("Zebrafish", "Human", "Mouse"))) {                 
                    species <- "Human"
                }

                dep_output <- strip_sig(dep_output) # Remove old sig cols if present

                withProgress(message = "Computing enrichment", value = 0, {
                    incProgress(0.2, detail = "Filtering significant hits")
                    data_cols_fp <- paste(sort(rv$data_cols[[tbl_name]]), collapse = ",")
                    dep_pg_cache_key <- make_cache_key(
                        "dep_pg",
                        species,
                        contrast,
                        paste0("pcut=", format(input$enrichment_pcutoff, scientific = FALSE)),
                        paste0("fccut=", format(input$enrichment_fccutoff, scientific = FALSE)),
                        paste0("cols=", data_cols_fp)
                    )

                    if (!is.null(cache) && cache$exists(dep_pg_cache_key)) {
                        dep_pg <- cache$get(dep_pg_cache_key)
                    } else {
                        dep_pg <- DEP2::add_rejections(dep_output, alpha = input$enrichment_pcutoff, lfc = input$enrichment_fccutoff)
                        if (!is.null(cache) && !cache$exists(dep_pg_cache_key)) {
                            cache$set(dep_pg_cache_key, dep_pg)
                        }
                    }
                    sig_count <- check_sig(dep_pg)

                    if (sig_count < 10) {
                        enrichment_plot_obj(NULL)
                        output$enrichment <- renderUI({
                            div(
                                style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                                paste0("Not enough significant hits (n < 10) for enrichment analysis with contrast '", contrast, "'.")
                            )
                        })
                        return(invisible())
                    }

                    message("Performing ORA enrichment analysis with contrast '", contrast, "' species '", species, "' and type '", enrichment_type, "' on '", sig_count, "' hits like '", paste(rownames(dep_pg)[1:5], collapse = ','), "' ...")

                    incProgress(0.6, detail = "Running ORA")
                    ora_cache_key <- make_cache_key(
                        "ora",
                        species,
                        enrichment_type,
                        contrast,
                        paste0("pcut=", format(input$enrichment_pcutoff, scientific = FALSE)),
                        paste0("fccut=", format(input$enrichment_fccutoff, scientific = FALSE)),
                        paste0("cols=", data_cols_fp)
                    )

                    if (!is.null(cache) && cache$exists(ora_cache_key)) {
                        message("Using cached ORA result for key: ", ora_cache_key)
                        res_ora <- cache$get(ora_cache_key)
                    } else {
                        res_ora <- tryCatch({
                            DEP2::test_ORA(
                                dep_pg,
                                type = enrichment_type,
                                species = species,
                                contrasts = contrast
                            )
                        }, error = function(e) {
                            output$enrichment <- renderUI({
                                div(
                                    style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                                    paste0("Enrichment failed: ", conditionMessage(e))
                                )
                            })
                            return(NULL)
                        })

                        if (!is.null(res_ora) && !is.null(cache) && !cache$exists(ora_cache_key)) {
                            cache$set(ora_cache_key, res_ora)
                        }
                    }

                    if (is.null(res_ora)) return(invisible())

                    ora_df <- tryCatch(as.data.frame(res_ora), error = function(e) NULL)
                    if (is.null(ora_df) || nrow(ora_df) == 0) {
                        enrichment_plot_obj(NULL)
                        output$enrichment <- renderUI({
                            div(
                                style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                                "No enriched terms found for current cutoffs/contrast."
                            )
                        })
                        return(invisible())
                    }

                    incProgress(0.2, detail = "Preparing enrichment plot")
                    plot_cache_key <- make_cache_key("ora_plot", ora_cache_key)
                    if (!is.null(cache) && cache$exists(plot_cache_key)) {
                        plot_obj <- cache$get(plot_cache_key)
                    } else {
                        plot_obj <- tryCatch(enrichplot::dotplot(res_ora), error = function(e) {
                            output$enrichment <- renderUI({
                                div(
                                    style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                                    paste0("Enrichment completed but plotting failed: ", conditionMessage(e))
                                )
                            })
                            NULL
                        })
                    }

                    if (is.null(plot_obj)) return(invisible())

                    if (!is.null(cache) && !cache$exists(plot_cache_key)) {
                        cache$set(plot_cache_key, plot_obj)
                    }

                    enrichment_plot_obj(plot_obj)
                    output$enrichment <- renderUI({
                        plotOutput(session$ns("enrichment_plot"), height = "520px")
                    })
                    output$enrichment_plot <- renderPlot({
                        req(!is.null(enrichment_plot_obj()))
                        enrichment_plot_obj()
                    })
                })
            },
            ignoreInit = TRUE
        )
    })
}

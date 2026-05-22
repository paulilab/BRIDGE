#' @export
server_function <- function(input, output, session, db_path) {
    #### DATA RETRIEVING SERVER ####

    con <- connect_db(db_path)
    session$onSessionEnded(function() {
        try(DBI::dbDisconnect(con), silent = TRUE)
        gc()
    })

    all_tables <- DBI::dbListTables(con)

    # Memory management: configurable cap in MB (default 8GB)
    mem_cap_mb <- as.numeric(Sys.getenv("BRIDGE_MEMORY_CAP_MB", "8192"))

    get_session_mem_mb <- function() {
        # Use R's gc() report for current session memory (in MB)
        gc_info <- gc(verbose = FALSE, reset = FALSE)
        sum(gc_info[, 2])  # "Mb" column = current usage
    }

    trim_module_caches <- function() {
        # Trim depflt_cache in heatmap_state to most recent 3 entries per table
        for (tbl in rv$table_names) {
            state <- rv$heatmap_state[[tbl]]
            if (!is.null(state) && length(state$depflt_cache) > 3) {
                keys <- names(state$depflt_cache)
                to_drop <- keys[seq_len(length(keys) - 3)]
                for (k in to_drop) state$depflt_cache[[k]] <- NULL
            }
        }
        gc()
    }

    # Periodic memory check (every 30 seconds)
    mem_check_timer <- shiny::reactiveTimer(30000)
    shiny::observe({
        mem_check_timer()
        mem_mb <- get_session_mem_mb()
        if (mem_mb > mem_cap_mb) {
            trim_module_caches()
            message(sprintf("BRIDGE: memory pressure (%.0f MB > %.0f MB cap), caches trimmed.", mem_mb, mem_cap_mb))
        }
    })

    # creation of cache starr object
    cache <- storr::storr_dbi(
        tbl_data = "storr_data",
        tbl_keys = "storr_keys",
        con = con
    )

    # Database cache maintenance: prune stale entries and cap per-table entries
    local({
        max_per_module <- as.integer(Sys.getenv("BRIDGE_CACHE_MAX_PER_MODULE", "10"))

        tryCatch({
            all_keys <- cache$list()
            if (length(all_keys) == 0L) return(invisible(NULL))

            # 1. Remove entries for tables no longer in the database
            db_tables <- tryCatch(
                DBI::dbGetQuery(con, "SELECT table_name FROM table_metadata")$table_name,
                error = function(e) character(0)
            )
            if (length(db_tables) > 0L) {
                stale <- vapply(all_keys, function(k) {
                    # Keys start with table_id (before first module-specific suffix)
                    # or use pipe-separated format (enrich|table|...)
                    if (grepl("^enrich\\|", k)) {
                        tbl <- strsplit(k, "\\|", fixed = FALSE)[[1]][2]
                    } else {
                        # Table name is the prefix before the column data
                        # Match against known table names
                        tbl <- NA_character_
                        for (db_tbl in db_tables) {
                            if (startsWith(k, paste0(db_tbl, "_"))) {
                                tbl <- db_tbl
                                break
                            }
                        }
                    }
                    !is.na(tbl) && !(tbl %in% db_tables)
                }, logical(1))
                for (k in all_keys[stale]) cache$del(k)
                all_keys <- all_keys[!stale]
            }

            # 2. Cap entries per table+module type (keep most recent N)
            # Group keys by table + module suffix pattern
            module_patterns <- c("_dep$", "_raw_heatmap_task$", "_dep_heatmap_task$",
                                 "_dep_volcano_task$", "_pca_v1$")
            for (pat in module_patterns) {
                matching <- grep(pat, all_keys, value = TRUE)
                if (length(matching) <= max_per_module) next
                # Group by table prefix (everything before the column data is table-specific)
                # Since we can't easily group, just cap total per pattern
                to_drop <- head(matching, length(matching) - max_per_module)
                for (k in to_drop) cache$del(k)
            }

            # Also cap enrichment entries
            enrich_keys <- grep("^enrich\\|", all_keys, value = TRUE)
            if (length(enrich_keys) > max_per_module * 3L) {
                to_drop <- head(enrich_keys, length(enrich_keys) - max_per_module * 3L)
                for (k in to_drop) cache$del(k)
            }

            # 3. Reclaim storage from orphaned hashes
            cache$gc()
        }, error = function(e) {
            message("BRIDGE: cache maintenance skipped: ", conditionMessage(e))
        })
    })

    rv <- reactiveValues(tables = list(), table_names = character(), data_cols = list(), datatype = character(), constrasts = list()) # variable that stores most of the important values for each table

    # Populate species choices from table_metadata
    local({
        metadata <- tryCatch(
            DBI::dbGetQuery(con, "SELECT table_name FROM table_metadata"),
            error = function(e) data.frame(table_name = character(0))
        )
        species_raw <- unique(sub("_.*", "", metadata$table_name))
        species_raw <- species_raw[nzchar(species_raw)]
        species_display <- tools::toTitleCase(species_raw)
        choices <- setNames(species_raw, species_display)
        shiny::updateSelectInput(session, "species", choices = choices)
    })

    valid_tables <- shiny::reactive({
        shiny::req(input$species)
        species <- tolower(input$species)

        # Query metadata table
        metadata <- DBI::dbGetQuery(con, "SELECT table_name FROM table_metadata")
        matches <- grep(paste0("^", species, "_"), metadata$table_name, value = TRUE)

        if (length(matches) == 0) {
            return(NULL)
        } else {
            return(matches)
        }
    })


    # Rendering of which table to choose in function of species selected
    output$table_selector <- shiny::renderUI({
        tables <- valid_tables()
        if (is.null(tables)) {
            return(shiny::helpText("No tables found for this species."))
        } else {
            shiny::selectInput("selected_table", "Select a data table:", choices = tables)
        }
    })

    # obtain desired metadata table from the species and table selected
    table_metadata <- shiny::reactive({
        shiny::req(input$selected_table)
        query <- sprintf("SELECT * FROM table_metadata WHERE table_name = '%s'", input$selected_table)
        metadata <- DBI::dbGetQuery(con, query)

        if (nrow(metadata) == 0) {
            return(NULL)
        }

        # Get all non-null fields in the metadata row
        column_types <- names(metadata)[!is.na(metadata[1, ]) & metadata[1, ] != ""]

        # Exclude 'table_name'
        column_types <- setdiff(column_types, "table_name")

        # Create a named list with each type and its corresponding columns
        meta_list <- lapply(column_types, function(ct) unlist(strsplit(metadata[[ct]], ",")))
        names(meta_list) <- column_types

        return(meta_list)
    })

    annotation_metadata <- shiny::reactive({
        shiny::req(input$species)
        annotation_table <- input$annotation_table
        query <- sprintf("SELECT * FROM '%s'", annotation_table)
        return(query)
    })

    # COLUMN SELECTOR INPUT

    # From the metadata table retrieved show the possible datapoints columns to choose
    output$column_selector <- shiny::renderUI({
        meta <- table_metadata()
        if (is.null(meta)) {
            return(shiny::helpText("No metadata available for this table."))
        }

        datapoint_cols <- unlist(strsplit(meta$datapoint_columns, ","))

        shinyWidgets::pickerInput(
            inputId = "datapoints_selected",
            label = "Select datapoint columns to load:",
            choices = c(datapoint_cols),
            selected = NULL,
            multiple = TRUE,
            options = shinyWidgets::pickerOptions(
                actionsBox = TRUE,
                liveSearch = TRUE,
                noneSelectedText = "Select columns",
                deselectAllText = "None",
                selectAllText = "Select All",
                dropupAuto = FALSE
            ),
            width = "100%"
        )
    })

    valid_anno <- shiny::reactive({
        shiny::req(input$species)
        species <- tolower(input$species)

        # Query metadata table
        metadata <- DBI::dbGetQuery(con, "SELECT table_name FROM annotation_metadata")
        matches <- grep(paste0("^", species, "_"), metadata$table_name, value = TRUE)

        if (length(matches) == 0) {
            return(NULL)
        } else {
            return(matches)
        }
    })


    output$annotation_selector <- shiny::renderUI({
        anno_files <- valid_anno()
        if (is.null(anno_files)) {
            return(shiny::helpText("No annotation files found for this species."))
        } else {
            shiny::selectInput("annotation_table", "Select an annotation file:", choices = anno_files)
        }
    })

    # LOAD DATA BUTTON SERVER

    # Server side of loading the data, selecting the desired columns
    shiny::observeEvent(input$load_data, {
        extract_valid_contrasts <- function(dep_obj) {
            rd_names <- character(0)

            if (methods::is(dep_obj, "DEGdata")) {
                tr <- dep_obj@test_result
                if (!is.null(tr)) rd_names <- colnames(tr)
            } else {
                rd <- SummarizedExperiment::rowData(dep_obj)
                if (!is.null(rd)) rd_names <- colnames(rd)
            }

            if (is.null(rd_names) || !length(rd_names)) return(character(0))

            sig_cols <- grep("_significant$", rd_names, value = TRUE)
            unique(sub("_significant$", "", sig_cols))
        }

        tryCatch({
            shiny::withProgress(message = "Loading and processing dataset", value = 0, {
            shiny::req(input$selected_table, input$datapoints_selected, input$annotation_table)

            shiny::incProgress(0.1, detail = "Reading metadata")
            meta <- table_metadata()
            annotation <- annotation_metadata()
            if (is.null(meta)) {
                return()
            }

            # Always include identifiers
            id_cols <- trimws(unlist(strsplit(meta$identifier_columns, ",")))
            datapoint_cols <- trimws(input$datapoints_selected)

            # Validate experimental design: need >= 2 conditions with >= 2 replicates each
            conditions <- sub("_[^_]+$", "", datapoint_cols)
            cond_counts <- table(conditions)
            if (length(cond_counts) < 2) {
                shiny::showNotification(
                    "At least 2 conditions are required for differential expression analysis. Please select columns from more conditions.",
                    type = "error", duration = NULL
                )
                return()
            }
            underpowered <- names(cond_counts)[cond_counts < 2]
            if (length(underpowered) > 0) {
                shiny::showNotification(
                    paste0("Each condition needs at least 2 replicates. The following have fewer: ",
                           paste(underpowered, collapse = ", "), ". Please select more replicates."),
                    type = "error", duration = NULL
                )
                return()
            }

            shiny::incProgress(0.2, detail = "Querying selected table")
            safe_cols <- sprintf('"%s"', unique(c(id_cols, datapoint_cols)))
            col_string <- paste(safe_cols, collapse = ", ")
            query <- sprintf("SELECT %s FROM %s", col_string, input$selected_table)
            new_data <- DBI::dbGetQuery(con, query)

            shiny::incProgress(0.15, detail = "Loading annotation")
            annotation_data <- DBI::dbGetQuery(con, annotation)
            table_id <- input$selected_table

            name_parts <- strsplit(table_id, "_", fixed = TRUE)[[1]]
            if (length(name_parts) < 2) {
                stop(sprintf("Could not infer datatype from table name '%s'. Expected '<species>_<datatype>_...'.", table_id))
            }
            datatype <- name_parts[2]

            is_reload <- table_id %in% rv$table_names

            # Always refresh table-specific state so new column selections take effect.
            per_table_slots <- c(
                "tables", "id_cols", "data_cols", "datatype",
                "species", "annotation", "dep_output", "contrasts", "constrasts",
                "current_dep_heatmap_key", "current_dep_volcano_key"
            )
            for (slot in per_table_slots) {
                if (is.list(rv[[slot]]) && !is.null(rv[[slot]][[table_id]])) {
                    rv[[slot]][[table_id]] <- NULL
                }
            }
            if (!is.null(rv$heatmap_state) && !is.null(rv$heatmap_state[[table_id]])) {
                rv$heatmap_state[[table_id]] <- NULL
            }

            shiny::incProgress(0.15, detail = "Updating app state")
            rv$tables[[table_id]] <- new_data
            rv$id_cols[[table_id]] <- id_cols
            rv$data_cols[[table_id]] <- input$datapoints_selected
            rv$datatype[[table_id]] <- datatype
            rv$species[[table_id]] <- input$species
            rv$annotation[[table_id]] <- annotation_data

            if (!is_reload) {
                rv$table_names <- c(rv$table_names, table_id)
            }

            cache_key <- paste(table_id, paste(rv$data_cols[[table_id]], collapse = "_"), "dep", sep = "_")

            shiny::incProgress(0.3, detail = "Preparing processed object")
            if (cache$exists(cache_key)) {
                dep_output <- cache$get(cache_key)
            } else {
                if (rv$datatype[[table_id]] == "proteomics") {
                    dep_output <- dep2_proteomics(rv$tables[[table_id]], table_id, rv)
                } else if (rv$datatype[[table_id]] == "phosphoproteomics") {
                    dep_output <- dep2_phosphoproteomics(rv$tables[[table_id]], table_id, rv)
                } else if (rv$datatype[[table_id]] == "rnaseq") {
                    dep_output <- dep2_rnaseq(rv$tables[[table_id]], table_id, rv)
                } else {
                    stop(sprintf("Unsupported datatype '%s' for table '%s'.", rv$datatype[[table_id]], table_id))
                }
                cache$set(cache_key, dep_output)
            }

            shiny::incProgress(0.1, detail = "Finalizing")
            valid_contrasts <- extract_valid_contrasts(dep_output)
            rv$dep_output[[table_id]] <- dep_output
            rv$contrasts[[table_id]] <- valid_contrasts

            })
        }, error = function(e) {
            showNotification(
                paste("Failed to load data:", conditionMessage(e)),
                type = "error",
                duration = NULL
            )
        })
    })

    # DELETE DATA SELECTOR

    shiny::observe({
        shiny::updateSelectInput(session, "remove", choices = rv$table_names, selected = NULL)
    })

    clear_table_state <- function(table_id) {
        if (is.null(table_id) || !nzchar(table_id)) return(invisible(NULL))

        remove_slot_entry <- function(slot_name, key) {
            obj <- rv[[slot_name]]
            if (is.null(obj)) return(invisible(NULL))

            # list-like slots
            if (is.list(obj)) {
                if (!is.null(obj[[key]])) obj[[key]] <- NULL
                rv[[slot_name]] <- obj
                return(invisible(NULL))
            }

            # named vectors (e.g., character)
            nms <- names(obj)
            if (!is.null(nms) && key %in% nms) {
                rv[[slot_name]] <- obj[nms != key]
                return(invisible(NULL))
            }

            invisible(NULL)
        }

        per_table_slots <- c(
            "tables", "id_cols", "data_cols", "datatype",
            "species", "annotation", "dep_output", "contrasts", "constrasts",
            "current_dep_heatmap_key", "current_dep_volcano_key"
        )

        for (slot in per_table_slots) {
            remove_slot_entry(slot, table_id)
        }

        if (!is.null(rv$heatmap_state) && !is.null(rv$heatmap_state[[table_id]])) {
            rv$heatmap_state[[table_id]] <- NULL
        }

        rv$table_names <- setdiff(rv$table_names, table_id)
        invisible(NULL)
    }

    # DELETE DATA BUTTON SERVER

    shiny::observeEvent(input$delete_data, {
        shiny::req(input$remove)

        if (input$remove %in% rv$table_names) {
            # Remove all per-table state
            clear_table_state(input$remove)

            # Optionally reset the dropdown selection
            shiny::updateSelectInput(session, "remove", choices = rv$table_names, selected = NULL)
        }
    })

    # LOADED TABLES INFO

    output$loaded_info <- shiny::renderUI({
        format_dims <- function(x) {
            if (is.null(x)) return("n/a")
            d <- tryCatch(dim(x), error = function(e) NULL)
            if (!is.null(d) && length(d) >= 2 && !any(is.na(d[1:2]))) {
                return(sprintf("%s x %s", format(d[1], big.mark = ","), format(d[2], big.mark = ",")))
            }
            "n/a"
        }

        format_dep_dims <- function(obj) {
            if (is.null(obj)) return("processing")

            # SummarizedExperiment and data.frame/matrix-like objects
            dims <- format_dims(obj)
            if (!identical(dims, "n/a")) return(dims)

            # Fallback for DEGdata-like objects
            tr <- tryCatch(obj@test_result, error = function(e) NULL)
            tr_dims <- format_dims(tr)
            if (!identical(tr_dims, "n/a")) return(tr_dims)

            "n/a"
        }

        if (length(rv$table_names) == 0) {
            return("No table loaded yet.")
        } else {
            # Create a list with each table name and download button
            htmltools::tagList(
                shiny::tags$ul(
                    lapply(rv$table_names, function(tbl_name) {
                        download_id <- paste0("download_", gsub("[^a-zA-Z0-9]", "_", tbl_name))
                        raw_dims <- format_dims(rv$tables[[tbl_name]])
                        dep_obj <- rv$dep_output[[tbl_name]]
                        dep_class <- if (is.null(dep_obj)) "pending" else class(dep_obj)[1]
                        dep_dims <- format_dep_dims(dep_obj)

                        shiny::tags$li(
                            shiny::tags$div(
                                shiny::tags$strong(tbl_name),
                                shiny::tags$br(),
                                shiny::tags$small(sprintf("Table: %s | RDS (%s): %s", raw_dims, dep_class, dep_dims)),
                                shiny::tags$br(),
                                shiny::downloadButton(download_id, label = "Download .rds.gz", class = "btn-sm")
                            )
                        )
                    })
                )
            )
        }
    })
    
    # Download handlers for each loaded table
    observeEvent(rv$table_names, {
        lapply(rv$table_names, function(tbl_name) {
            download_id <- paste0("download_", gsub("[^a-zA-Z0-9]", "_", tbl_name))
            
            output[[download_id]] <- downloadHandler(
                filename = function() {
                    paste0(tbl_name, ".rds.gz")
                },
                content = function(file) {
                    # Get the appropriate object based on datatype
                    obj <- rv$dep_output[[tbl_name]]
                    
                    if (is.null(obj)) {
                        showNotification(paste("Object for", tbl_name, "not yet processed."), type = "error")
                        return(NULL)
                    }
                    
                    # Save with gzip compression
                    saveRDS(obj, file = file, compress = "gzip")
                }
            )
        })
    })

    #### END OF DATA RETRIEVING SERVER ####


    #### MAIN UI FOR EACH TABLE #### (includes functions for heatmap ui)

    output$all_tables_ui <- renderUI({
        req(length(rv$tables) > 0, !is.null(rv$contrasts))
        tagList(lapply(rv$table_names, function(tbl_name) {
            shinydashboard::box(
                title = paste("Table:", tbl_name), width = 12, solidHeader = TRUE,
                status = "primary", style = "overflow-x: auto", collapsible = TRUE,
                shinydashboard::box(
                    title = "Table", width = 12, solidHeader = TRUE, status = "primary",
                    style = "overflow-x: auto", collapsible = TRUE,
                    h3(), DT::DTOutput(paste0("table_", tbl_name)), h3()
                ),
                shinydashboard::tabBox(
                    title = shinyWidgets::actionBttn("individual_help",
                        "Help",
                        color = "primary", icon = icon("question-circle"),
                        size = "sm", style = "bordered"
                    ),
                    width = 12, id = paste0("plotTabs_", tbl_name), selected = "Raw Heatmap",
                    tabPanel(
                        "Raw Heatmap",
                        fluidRow(RawHeatmapUI(paste0("RawHeatmap_", tbl_name), tbl_name))
                    ),
                    tabPanel(
                        "DE Heatmap",
                        fluidRow(DepHeatmapUI(paste0("DEHeatmap_", tbl_name), tbl_name))
                    ),
                    tabPanel(
                        "Volcano Plot",
                        fluidRow(VolcanoUI(paste0("Volcano_", tbl_name), tbl_name, rv$contrasts[[tbl_name]]))
                    ),
                    tabPanel(
                        "Gene Expression",
                        fluidRow(datapointsUI(paste0("datapoints_", tbl_name), tbl_name))
                        # remove datapointsTableUI if merged
                    ),
                    tabPanel(
                        "Enrichment analysis",
                        fluidRow(EnrichmentUI(paste0("Enrichment_", tbl_name), tbl_name))
                    ),
                    tabPanel(
                        "PCA",
                        fluidRow(pcaUI(paste0("pca_", tbl_name), tbl_name))
                    )
                )
            )
        }))
    })

    #### END OF MAIN UI FOR EACH TABLE ####


    # RENDER OF RAW DATA TABLES

    shiny::observe({
        lapply(rv$table_names, function(tbl_name) {
            dep_output <- rv$dep_output[[tbl_name]]
            df_raw <- rv$tables[[tbl_name]]
            if (is.null(df_raw) || is.null(dep_output)) {
                output[[paste0("table_", tbl_name)]] <- DT::renderDT({
                    DT::datatable(data.frame(Message = "Table is being removed or not loaded."),
                        options = list(dom = "t", ordering = FALSE, paging = FALSE, searching = FALSE),
                        rownames = FALSE
                    )
                })
                return(invisible(NULL))
            }
            if (identical(rv$datatype[[tbl_name]], "rnaseq")){
                df_dep <- DEP2::get_results(dep_output) %>%
                    dplyr::rename(Gene_ID = ID) %>%
                    dplyr::mutate(Gene_ID = gsub(".*_(ENS.*G[0-9]+)", "\\1", Gene_ID))
                    complete <- dplyr::left_join(df_raw, df_dep, by = "Gene_ID")
                     # %>% dplyr::select(-any_of(c("name")))
            } else if (identical(rv$datatype[[tbl_name]], "phosphoproteomics")) {
                df_dep <- DEP2::get_results(dep_output)                
                complete <- dplyr::left_join(df_raw %>% dplyr::mutate(
                        ID = paste0(Protein_ID, "_", pepG)), df_dep, by = "ID") 
                        #%>% dplyr::select(-any_of(c("name")))                        
            } else {
                df_dep <- DEP2::get_results(dep_output) %>%
                    dplyr::rename(Protein_ID = ID)
                # message("DF_DEP: ", paste0(colnames(df_dep), collapse = ","), "ID: ", head(df_dep$Protein_ID, 3), " Gene_Name: ", head(df_dep$Gene_Name, 3))
                complete <- dplyr::left_join(df_raw, df_dep, by = "Protein_ID")
                # %>% dplyr::select(-any_of(c("name")))
            }
            output[[paste0("table_", tbl_name)]] <- DT::renderDT({
                DT::datatable(complete %>% dplyr::select(where(~ !is.numeric(.)), where(is.numeric)), extensions = "Buttons", filter = "top", options = list(scrollX = TRUE, pageLength = 10, lengthMenu = c(5, 10, 25, 50, 100), dom = "Blfrtip", buttons = list(
                        list(extend = "copy", exportOptions = list(modifier = list(page = "all"))),
                        list(extend = "csv",  exportOptions = list(modifier = list(page = "all"))),
                        list(extend = "excel",exportOptions = list(modifier = list(page = "all"))),
                        list(extend = "pdf",  exportOptions = list(modifier = list(page = "all"))),
                        "print"
                    )))
            })
        })
    })

    # RawHeatmapServer(paste0("RawHeatmap_", tbl_name), rv) #Function in heatmap_server.R
    wired <- reactiveValues(
        datapoints = character(),
        raw      = character(),
        dep      = character(),
        pca      = character(),
        enrich   = character()
    )

    # Helper predicates
    has_dep <- function(tbl) !is.null(rv$dep_output[[tbl]])
    has_contr <- function(tbl) !is.null(rv$contrasts[[tbl]]) && length(rv$contrasts[[tbl]]) > 0

    observeEvent(list(rv$table_names, rv$dep_output, rv$contrasts), ignoreInit = FALSE, {
        tbls <- rv$table_names

        for (tbl in tbls) {
            # Always-ok modules (need raw table + data cols)
            if (!(tbl %in% wired$datapoints) && !is.null(rv$tables[[tbl]]) && !is.null(rv$data_cols[[tbl]])) {
                datapointsServer(paste0("datapoints_", tbl), rv, tbl)
                wired$datapoints <- union(wired$datapoints, tbl)
            }
            if (!(tbl %in% wired$raw) && !is.null(rv$tables[[tbl]])) {
                RawHeatmapServer(paste0("RawHeatmap_", tbl), rv, tbl, cache)
                wired$raw <- union(wired$raw, tbl)
            }

            # DEP / PCA depend on dep_output
            if (!(tbl %in% wired$dep) && has_dep(tbl)) {
                DepHeatmapServer(paste0("DEHeatmap_", tbl), rv, cache, tbl)
                VolcanoServer(paste0("Volcano_", tbl), rv, cache, tbl)
                wired$dep <- union(wired$dep, tbl)
            }
            if (!(tbl %in% wired$pca) && has_dep(tbl)) {
                pcaServer(paste0("pca_", tbl), rv, cache, tbl)
                wired$pca <- union(wired$pca, tbl)
            }

            # Enrichment also needs contrasts
            if (!(tbl %in% wired$enrich) && has_dep(tbl) && has_contr(tbl)) {
                EnrichmentServer(paste0("Enrichment_", tbl), rv, cache, tbl)
                wired$enrich <- union(wired$enrich, tbl)
            }
        }
    })

    # Optional cleanup when a table is removed
    observeEvent(input$delete_data, {
        wired$datapoints <- setdiff(wired$datapoints, input$remove)
        wired$raw <- setdiff(wired$raw, input$remove)
        wired$dep <- setdiff(wired$dep, input$remove)
        wired$pca <- setdiff(wired$pca, input$remove)
        wired$enrich <- setdiff(wired$enrich, input$remove)
    })

    ############################

    integration_ui(input, output, session, rv)

    ###########################

    help_buttons(input, output, session)
}

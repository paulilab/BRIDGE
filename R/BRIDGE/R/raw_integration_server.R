#' @export
raw_integration <- function(input, output, session, rv, combined_data) {
    shiny::observeEvent(input$integrate_data, {
        shiny::req(length(input$integration) > 0)
        selected_tables <- input$integration
        selected_types <- sapply(selected_tables, function(tbl) rv$datatype[[tbl]])

        # Extract selected columns
        selected_lists <- lapply(selected_tables, function(tbl) input[[paste0("cols_selected_int", tbl)]])

        # Validate selection
        if (any(sapply(selected_lists, is.null))) {
            combined_data(data.frame(Message = "Select columns from each table before integrating."))
            rv$integration_reference_cols <- NULL
            return()
        }

        # Store all unique data column names across all selected tables
        reference_data_names <- unique(unlist(selected_lists))
        tables_to_join <- list()

        for (i in seq_along(selected_tables)) {
            tbl <- selected_tables[[i]]
            selected_data_cols <- selected_lists[[i]]
            df <- rv$tables[[tbl]][, c(rv$id_cols[[tbl]], selected_data_cols), drop = FALSE]

            # Construct ID column
            if (rv$datatype[[tbl]] == "phosphoproteomics" && all(c("Gene_Name", "pepG") %in% names(df))) {
                df$unique_id <- paste(df$Gene_Name, df$pepG, sep = "_")
            } else if (rv$datatype[[tbl]] == "proteomics" && all(c("Gene_Name", "Gene_ID") %in% names(df))) {
                df$unique_id <- paste(df$Gene_Name, df$Gene_ID, sep = "_")
            } else {
                df$unique_id <- df$Gene_Name
            }

            # Ensure ID columns match across tables
            id_cols <- setdiff(c("unique_id", rv$id_cols[[tbl]]), selected_data_cols)
            all_id_names <- unique(unlist(lapply(selected_tables, function(x) c("unique_id", rv$id_cols[[x]]))))

            # Add missing ID columns with NA
            missing_ids <- setdiff(all_id_names, names(df))
            for (col in missing_ids) df[[col]] <- NA

            # Keep original data column names (don't force renaming to first table)
            # This allows integration of datasets with different numbers of columns
            
            # Reorder columns: ID columns first, then data columns, then source
            data_cols_order <- setdiff(names(df), c(all_id_names, "source"))
            df <- df[, c(all_id_names, data_cols_order), drop = FALSE]
            df$source <- tbl
            tables_to_join[[tbl]] <- df
        }

        # Combine all
        combined_df <- dplyr::bind_rows(tables_to_join) %>% dplyr::select(where(~ !is.numeric(.)), where(is.numeric))
        combined_data(combined_df)
        rv$integration_reference_cols <- reference_data_names

        # Update search box
        updateSelectizeInput(session, "search_gene_integration", choices = sort(unique(gsub("_.*", "", combined_df$unique_id))), server = TRUE)
    })
}

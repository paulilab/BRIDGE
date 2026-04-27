suppressPackageStartupMessages({
    library(shiny)
})

has_shinyfiles <- requireNamespace("shinyFiles", quietly = TRUE)

get_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
    }
    normalizePath(getwd())
}

project_root <- get_script_dir()
python_dir <- file.path(project_root, "Python")
raw_script_path <- file.path(python_dir, "db_adding.py")
annotation_script_path <- file.path(python_dir, "db_adding_annotation.py")

find_python <- function() {
    py3 <- Sys.which("python3")
    if (nzchar(py3)) return(py3)
    py <- Sys.which("python")
    if (nzchar(py)) return(py)
    ""
}

python_exec <- find_python()

pick_path <- function(text_path, upload_input) {
    if (!is.null(text_path) && nzchar(trimws(text_path))) {
        return(normalizePath(trimws(text_path), mustWork = FALSE))
    }
    if (!is.null(upload_input) && !is.null(upload_input$datapath)) {
        return(upload_input$datapath)
    }
    ""
}

detect_sep <- function(file_path) {
    if (grepl("\\.gz$", file_path, ignore.case = TRUE)) {
        con <- gzfile(file_path, "rt")
    } else {
        con <- file(file_path, "rt")
    }
    on.exit(close(con), add = TRUE)

    first_line <- readLines(con, n = 1, warn = FALSE)
    if (length(first_line) == 0) return(",")
    if (grepl("\\t", first_line)) "\t" else ","
}

read_columns_preview <- function(file_path) {
    sep <- detect_sep(file_path)
    df <- utils::read.table(
        file = file_path,
        sep = sep,
        header = TRUE,
        nrows = 1,
        quote = "\"",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    data.frame(index = seq_along(colnames(df)), column = colnames(df))
}

run_python_script <- function(script_path, args) {
    output <- tryCatch(
        {
            out <- system2(
                command = python_exec,
                args = c(script_path, args),
                stdout = TRUE,
                stderr = TRUE
            )
            list(status = attr(out, "status") %||% 0L, output = out)
        },
        error = function(e) {
            list(status = 1L, output = paste("Failed to run script:", conditionMessage(e)))
        }
    )
    output
}

`%||%` <- function(a, b) {
    if (!is.null(a)) a else b
}

ui <- fluidPage(
    titlePanel("BRIDGE Database Builder"),
    tags$p(
        "This app wraps the Python database scripts with an interactive UI. ",
        "Use table previews to get 1-based column indices for identifier/datapoint selections."
    ),
    tags$details(
        tags$summary(tags$b("Naming and file requirements (click to expand)")),
        tags$div(
            style = "margin-top: 10px;",
            tags$h4("Raw / omics input"),
            tags$ul(
                tags$li("Accepted file types: CSV/TSV (optionally .gz compressed)."),
                tags$li("Raw table name pattern: <species>_<datatype>_<optional info>_<id>."),
                tags$li("Example: zebrafish_proteomics_test_1."),
                tags$li("Identifier columns should include: Gene_Name, Gene_ID, Protein_ID."),
                tags$li("For phosphoproteomics, include identifier column pepG."),
                tags$li("Datapoint column names should end with replicate index, e.g. X6.hpf_1.")
            ),
            tags$h4("Annotation input"),
            tags$ul(
                tags$li("Accepted file types: CSV/TSV (optionally .gz compressed)."),
                tags$li("Annotation table name pattern: <species>_annotation_<version>."),
                tags$li("Example: zebrafish_annotation_GRCz11."),
                tags$li("Required columns: Gene_ID, Gene_Name, Chromosome, Gene_Start, Gene_End, Gene_Type, Strand."),
                tags$li("Use the same species prefix as your raw data tables.")
            ),
            tags$h4("Processed RDS (optional)"),
            tags$ul(
                tags$li("If you attach processed data, provide an .rds file with a SummarizedExperiment-compatible object."),
                tags$li("It should correspond to the same raw table/datapoint columns for correct cache matching.")
            )
        )
    ),
    tags$hr(),
    fluidRow(
        column(
            width = 8,
            textInput("db_path", "Database file path", value = file.path(getwd(), "user_database.db")),
            if (has_shinyfiles) {
                shinyFiles::shinySaveButton(
                    id = "browse_db_path",
                    label = "Browse DB path",
                    title = "Choose where to create/select the database file",
                    filetype = list(db = c("db"), sqlite = c("sqlite", "sqlite3"), all = c("*"))
                )
            } else {
                tags$small(
                    style = "color:#666;",
                    "Install 'shinyFiles' to enable filesystem browsing for DB path."
                )
            }
        ),
        column(
            width = 4,
            br(),
            actionButton("create_db", "Create empty DB file")
        )
    ),
    verbatimTextOutput("global_status"),
    tabsetPanel(
        tabPanel(
            "Raw / omics table",
            fluidRow(
                column(
                    width = 6,
                    fileInput("raw_csv_upload", "Select CSV/TSV(.gz) (optional upload)"),
                    textInput("raw_csv_path", "Or enter CSV/TSV path on disk", value = ""),
                    actionButton("preview_raw_cols", "Preview columns"),
                    tableOutput("raw_cols_table")
                ),
                column(
                    width = 6,
                    textInput("raw_table", "Database table name", value = ""),
                    textInput("id_cols", "Identifier columns (1-based, e.g. 1,2,3 or 1:3)", value = ""),
                    textInput("tp_cols", "Datapoint columns (1-based, e.g. 4:20)", value = ""),
                    checkboxInput("add_processed", "Also attach processed RDS", value = FALSE),
                    fileInput("rds_upload", "Select processed RDS (optional upload)"),
                    textInput("rds_path", "Or enter RDS path on disk", value = ""),
                    actionButton("run_raw", "Add raw table", class = "btn-primary")
                )
            )
        ),
        tabPanel(
            "Annotation table",
            fluidRow(
                column(
                    width = 6,
                    fileInput("ann_csv_upload", "Select annotation CSV/TSV(.gz) (optional upload)"),
                    textInput("ann_csv_path", "Or enter annotation file path on disk", value = ""),
                    textInput("ann_table", "Annotation table name", value = ""),
                    actionButton("run_annotation", "Add annotation table", class = "btn-primary")
                ),
                column(
                    width = 6,
                    tags$h4("Naming reminder"),
                    tags$ul(
                        tags$li("Raw table: <species>_<datatype>_<optional info>_<id>"),
                        tags$li("Annotation table: <species>_annotation_<version>")
                    )
                )
            )
        )
    ),
    tags$hr(),
    tags$h4("Script output"),
    verbatimTextOutput("script_log")
)

server <- function(input, output, session) {
    script_log <- reactiveVal("Ready.\n")

    append_log <- function(lines) {
        old <- script_log()
        new_text <- paste(c(old, lines), collapse = "\n")
        script_log(new_text)
    }

    output$global_status <- renderText({
        checks <- c(
            sprintf("Project root: %s", project_root),
            sprintf("Python executable: %s", ifelse(nzchar(python_exec), python_exec, "NOT FOUND")),
            sprintf("Raw script found: %s", file.exists(raw_script_path)),
            sprintf("Annotation script found: %s", file.exists(annotation_script_path))
        )
        paste(checks, collapse = "\n")
    })

    output$script_log <- renderText(script_log())

    if (has_shinyfiles) {
        roots <- c(home = path.expand("~"), wd = normalizePath(getwd()), root = "/")
        shinyFiles::shinyFileSave(input, "browse_db_path", roots = roots, session = session)

        observeEvent(input$browse_db_path, {
            save_info <- shinyFiles::parseSavePath(roots, input$browse_db_path)
            if (nrow(save_info) > 0) {
                selected_path <- as.character(save_info$datapath[1])
                updateTextInput(session, "db_path", value = selected_path)
                append_log(sprintf("Selected DB path via browser: %s", selected_path))
            }
        })
    }

    observeEvent(input$create_db, {
        db_path <- trimws(input$db_path)
        if (!nzchar(db_path)) {
            showNotification("Please provide a database path.", type = "error")
            return()
        }
        dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
        ok <- file.create(db_path)
        if (!ok && !file.exists(db_path)) {
            showNotification("Could not create database file.", type = "error")
            return()
        }
        append_log(sprintf("Created/verified DB file: %s", normalizePath(db_path, mustWork = FALSE)))
        showNotification("Database file ready.", type = "message")
    })

    observeEvent(input$preview_raw_cols, {
        csv_path <- pick_path(input$raw_csv_path, input$raw_csv_upload)
        if (!nzchar(csv_path) || !file.exists(csv_path)) {
            showNotification("Provide a valid raw CSV/TSV path first.", type = "error")
            return()
        }

        cols <- tryCatch(
            read_columns_preview(csv_path),
            error = function(e) {
                showNotification(paste("Could not read column preview:", conditionMessage(e)), type = "error")
                NULL
            }
        )
        output$raw_cols_table <- renderTable(cols)
    })

    observeEvent(input$run_raw, {
        if (!nzchar(python_exec)) {
            showNotification("Python executable not found in PATH.", type = "error")
            return()
        }
        if (!file.exists(raw_script_path)) {
            showNotification("Could not find Python/db_adding.py.", type = "error")
            return()
        }

        csv_path <- pick_path(input$raw_csv_path, input$raw_csv_upload)
        rds_path <- pick_path(input$rds_path, input$rds_upload)
        db_path <- trimws(input$db_path)
        table_name <- trimws(input$raw_table)
        id_cols <- trimws(input$id_cols)
        tp_cols <- trimws(input$tp_cols)

        if (!nzchar(csv_path) || !file.exists(csv_path)) {
            showNotification("Provide a valid raw CSV/TSV path.", type = "error")
            return()
        }
        if (!nzchar(db_path)) {
            showNotification("Provide a database path.", type = "error")
            return()
        }
        if (!nzchar(table_name)) {
            showNotification("Provide a table name.", type = "error")
            return()
        }
        if (!nzchar(id_cols) || !nzchar(tp_cols)) {
            showNotification("Provide both identifier and datapoint column selections.", type = "error")
            return()
        }
        if (isTRUE(input$add_processed) && (!nzchar(rds_path) || !file.exists(rds_path))) {
            showNotification("Processed RDS selected, but file path is missing/invalid.", type = "error")
            return()
        }

        args <- c(
            "--csv", csv_path,
            "--db", db_path,
            "--table", table_name,
            "--id-cols", id_cols,
            "--tp-cols", tp_cols
        )

        if (isTRUE(input$add_processed)) {
            args <- c(args, "--processed", "--rds", rds_path)
        }

        append_log(sprintf("\n[RAW] Running: %s %s", python_exec, raw_script_path))
        res <- run_python_script(raw_script_path, args)
        append_log(res$output)

        if (identical(as.integer(res$status), 0L)) {
            showNotification("Raw table successfully added.", type = "message")
        } else {
            showNotification("Raw table upload failed. See script output.", type = "error")
        }
    })

    observeEvent(input$run_annotation, {
        if (!nzchar(python_exec)) {
            showNotification("Python executable not found in PATH.", type = "error")
            return()
        }
        if (!file.exists(annotation_script_path)) {
            showNotification("Could not find Python/db_adding_annotation.py.", type = "error")
            return()
        }

        csv_path <- pick_path(input$ann_csv_path, input$ann_csv_upload)
        db_path <- trimws(input$db_path)
        table_name <- trimws(input$ann_table)

        if (!nzchar(csv_path) || !file.exists(csv_path)) {
            showNotification("Provide a valid annotation CSV/TSV path.", type = "error")
            return()
        }
        if (!nzchar(db_path)) {
            showNotification("Provide a database path.", type = "error")
            return()
        }
        if (!nzchar(table_name)) {
            showNotification("Provide an annotation table name.", type = "error")
            return()
        }

        args <- c(
            "--csv", csv_path,
            "--db", db_path,
            "--table", table_name
        )

        append_log(sprintf("\n[ANNOTATION] Running: %s %s", python_exec, annotation_script_path))
        res <- run_python_script(annotation_script_path, args)
        append_log(res$output)

        if (identical(as.integer(res$status), 0L)) {
            showNotification("Annotation table successfully added.", type = "message")
        } else {
            showNotification("Annotation upload failed. See script output.", type = "error")
        }
    })
}

app <- shiny::shinyApp(ui = ui, server = server)

if (interactive() || Sys.getenv("RUN_STANDALONE_DB_BUILDER", "0") == "1") {
    shiny::runApp(app, port = as.integer(Sys.getenv("BRIDGE_DB_BUILDER_PORT", "3839")))
}

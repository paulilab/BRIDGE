if (!require(BRIDGE)) {
    library(devtools)
    devtools::document("R/BRIDGE")
    devtools::install("R/BRIDGE", keep_source = T, upgrade = FALSE)    
}
#print("Loading BRIDGE package")
library(BRIDGE)
#print("Loaded BRIDGE package")

suppressPackageStartupMessages({
    library(future.callr) # or future for multisession
    library(promises)
    library(tidyverse)
    library(shiny)
    library(shinybusy)
    library(DEP2)
    library(NbClust)
    library(SummarizedExperiment)
    library(ggplot2)
    library(plotly)
    library(ComplexHeatmap)
    library(matrixStats)
    library(gtools)
    library(tidySummarizedExperiment)
})
# library(pool)
# library(shinydashboard)
# library(shinyWidgets)
ht_opt$message <- FALSE

parse_cli_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    opts <- list()
    positional <- character(0)

    for (arg in args) {
        if (startsWith(arg, "--")) {
            kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
            key <- kv[[1]]
            value <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "true"
            opts[[key]] <- value
        } else {
            positional <- c(positional, arg)
        }
    }

    list(options = opts, positional = positional)
}

cli_args <- parse_cli_args()

safe_int_or_na <- function(x) {
    val <- suppressWarnings(as.integer(x))
    if (length(val) < 1L || is.na(val[[1]])) return(NA_integer_)
    val[[1]]
}

configure_future_backend <- function(cli_options = list()) {
    backend_cli <- cli_options[["backend"]]
    workers_cli <- safe_int_or_na(cli_options[["workers"]])
    max_workers_cli <- safe_int_or_na(cli_options[["max-workers"]])

    backend <- tolower(if (!is.null(backend_cli) && nzchar(backend_cli)) {
        backend_cli
    } else {
        Sys.getenv("BRIDGE_FUTURE_BACKEND", "auto")
    })

    workers_env <- safe_int_or_na(Sys.getenv("BRIDGE_FUTURE_WORKERS", NA))
    max_workers_env <- safe_int_or_na(Sys.getenv("BRIDGE_FUTURE_MAX_WORKERS", NA))
    max_workers <- if (!is.na(max_workers_cli) && max_workers_cli >= 1L) {
        max_workers_cli
    } else if (!is.na(max_workers_env) && max_workers_env >= 1L) {
        max_workers_env
    } else {
        2L
    }

    auto_workers <- min(max_workers, max(1L, future::availableCores() - 1L))
    workers_raw <- if (!is.na(workers_cli) && workers_cli >= 1L) {
        workers_cli
    } else if (!is.na(workers_env) && workers_env >= 1L) {
        workers_env
    } else {
        auto_workers
    }
    workers <- min(max_workers, max(1L, workers_raw))

    if (!is.na(workers_raw) && workers_raw > max_workers) {
        warning(sprintf(
            "Requested workers=%d exceeds hard cap %d; using %d.",
            workers_raw, max_workers, workers
        ))
    }

    if (identical(backend, "callr")) {
        future::plan(future.callr::callr, workers = workers)
        return(invisible(list(backend = "callr", workers = workers)))
    }

    if (identical(backend, "multisession")) {
        future::plan(future::multisession, workers = workers)
        return(invisible(list(backend = "multisession", workers = workers)))
    }

    # auto: prefer persistent workers on Linux for lower per-task startup overhead,
    # and fallback to callr for portability/stability.
    if (identical(backend, "auto")) {
        if (tolower(Sys.info()[["sysname"]]) == "linux") {
            future::plan(future::multisession, workers = workers)
            return(invisible(list(backend = "multisession", workers = workers)))
        }

        future::plan(future.callr::callr, workers = workers)
        return(invisible(list(backend = "callr", workers = workers)))
    }

    warning(
        sprintf(
            "Unknown BRIDGE_FUTURE_BACKEND='%s'. Falling back to auto.",
            backend
        )
    )

    if (tolower(Sys.info()[["sysname"]]) == "linux") {
        future::plan(future::multisession, workers = workers)
        return(invisible(list(backend = "multisession", workers = workers)))
    }

    future::plan(future.callr::callr, workers = workers)
    invisible(list(backend = "callr", workers = workers))
}

future_cfg <- configure_future_backend(cli_args$options)
message(sprintf("BRIDGE async backend: %s (%d worker%s)",
    future_cfg$backend,
    future_cfg$workers,
    ifelse(future_cfg$workers == 1L, "", "s")
))
set.seed(42)

#print("Getting CLI arguments")

# Determine database path
get_db_path <- function() {
    # 1. Try environment variable (for Docker/Shiny Server)
    db_path <- Sys.getenv("BRIDGE_DB_PATH", unset = NA)
    if (!is.na(db_path) && file.exists(db_path)) return(db_path)
    # 2. Try default Docker path
    if (file.exists("/srv/data/database.db")) return("/srv/data/database.db")
    # 3. If interactive, prompt user
    if (interactive()) return(file.choose())
    # 4. If command line, use first argument
    args <- cli_args$positional
    if (length(args) > 0 && file.exists(args[1])) return(args[1])
    stop("No valid database path found.")
}
db_path <- get_db_path()
#print(paste("db_path:", db_path))

get_port <- function() {
    # 1. Try environment variable (for Docker/Shiny Server)
    port_env <- Sys.getenv("BRIDGE_PORT", unset = NA)
    if (!is.na(port_env) && nzchar(port_env)) return(as.integer(port_env))
    # 2. Try command line argument
    args <- cli_args$positional
    if (length(args) > 1 && !is.na(as.integer(args[2]))) return(as.integer(args[2]))
    # 3. Default
    return(3838L)
}
port <- get_port()
#print(paste("port:", port))
#print("About to define ui and server")

# Define UI and server
ui <- BRIDGE::ui
#print("UI")
#print(exists("ui"))
#print(exists("ui", where=asNamespace("BRIDGE")))

server <- function(input, output, session) {
    BRIDGE::server_function(input, output, session, db_path)
}
#print("SERVER")
#print(exists("server"))

app <- shiny::shinyApp(ui = ui, server = server)

# Only run the app if not under Shiny Server
if (interactive() || (Sys.getenv("RUN_STANDALONE", "0") == "1")) {
    app <- shiny::runApp(
        shiny::shinyApp(ui = ui, server = server),
        port = port
    )
}

app
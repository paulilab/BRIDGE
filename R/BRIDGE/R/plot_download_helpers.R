plot_download_controls <- function(ns, id_prefix, default = "pdf") {
    shiny::tagList(
        shiny::div(
            style = "display:flex; gap:10px; align-items:end; flex-wrap:wrap; margin-bottom:8px;",
            shiny::selectInput(
                ns(paste0(id_prefix, "_download_format")),
                "Download format",
                choices = c("pdf", "svg"),
                selected = default,
                width = "150px"
            ),
            shiny::downloadButton(ns(paste0(id_prefix, "_download")), "Download plot")
        )
    )
}

register_plot_download <- function(input, output, session, id_prefix, filename_prefix, plot_fun, width = 10, height = 7) {
    output[[paste0(id_prefix, "_download")]] <- shiny::downloadHandler(
        filename = function() {
            fmt <- input[[paste0(id_prefix, "_download_format")]]
            if (is.null(fmt) || !(fmt %in% c("pdf", "svg"))) fmt <- "pdf"
            paste0(filename_prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", fmt)
        },
        content = function(file) {
            fmt <- input[[paste0(id_prefix, "_download_format")]]
            if (is.null(fmt) || !(fmt %in% c("pdf", "svg"))) fmt <- "pdf"

            if (identical(fmt, "svg")) {
                grDevices::svg(file, width = width, height = height)
            } else {
                grDevices::pdf(file, width = width, height = height, onefile = FALSE)
            }
            on.exit(grDevices::dev.off(), add = TRUE)

            p <- plot_fun()
            if (is.null(p)) {
                # plot_fun drew the plot in-place (e.g. draw + decorate)
            } else if (inherits(p, "Heatmap") || inherits(p, "HeatmapList")) {
                ComplexHeatmap::draw(p, merge_legend = TRUE, newpage = FALSE)
            } else {
                print(p)
            }
        }
    )
}

#' @export
pcaUI <- function(id, tbl_name) {
    ns <- NS(id)

    shinydashboard::box(
        title = "PCA", width = 12, solidHeader = TRUE, status = "info",
        fluidRow(
            column(
                3,
                shinyWidgets::actionBttn(ns("compute"),
                    span("Compute PCA", style = "color: white;"),
                    style = "simple", color = "primary", size = "sm"
                )
            ),
            column(
                3,
                shinyWidgets::switchInput(
                    ns("interactive"),
                    "Interactive plot",
                    value = FALSE,
                    onLabel = "YES",
                    offLabel = "NO",
                    width = "auto"
                )
            ),
            column(6, uiOutput(ns("plot_download_ui")))
        ),
        shinycssloaders::withSpinner(
            uiOutput(ns("plot_slot")),
            type = 8, color = "#2b8cbe", caption = "Loading..."
        ),
        h5(),
        shinydashboard::box(
            title = "PCA loadings", width = 12, solidHeader = TRUE,
            status = "info", style = "overflow-x: auto",
            collapsible = TRUE, collapsed = FALSE,
            # DT::DTOutput(ns("pcs"), height = "300px") |>
            #    shinycssloaders::withSpinner(
            #        type = 8, color = "#2b8cbe",
            #        caption = "Loading...", hide.ui = FALSE
            #    )
            uiOutput(ns("pcs_panel"))
        )
    )
}

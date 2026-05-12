#' @export
pcaUI <- function(id, tbl_name) {
    ns <- NS(id)

    shinydashboard::box(
        title = "PCA", width = 12, solidHeader = TRUE, status = "info",
        fluidRow(
            column(
                2,
                shinyWidgets::actionBttn(ns("compute"),
                    span("Compute PCA", style = "color: white;"),
                    style = "simple", color = "primary", size = "sm"
                )
            ),
            column(
                2,
                numericInput(
                    ns("top_n_features"),
                    "Top variable features",
                    value = 500,
                    min = 50,
                    max = 5000,
                    step = 50
                )
            ),
            column(4, uiOutput(ns("pc_axes_ui"))),
            column(4, uiOutput(ns("plot_download_ui")))
        ),
        plotOutput(ns("plot"), height = "480px"),
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

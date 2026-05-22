#' @export
ui <- shinydashboard::dashboardPage(
    skin = "blue",
    shinydashboard::dashboardHeader(title = "BRIDGE"),
    shinydashboard::dashboardSidebar(
        shinydashboard::sidebarMenu(
            shinydashboard::menuItem("Individual Exploration", tabName = "raw_data", icon = shiny::icon("table")),
            shinydashboard::menuItem("Integration", tabName = "Integration", icon = shiny::icon("th")),
            shinyWidgets::actionBttn("general_help", "Help", icon = shiny::icon("question-circle"), size = "sm", style = "bordered")
        ),
        shiny::tags$div(
            style = "padding: 10px 15px;",
            shiny::tags$a(
                href = "https://bridge-paulilab.readthedocs.io/en/latest/",
                target = "_blank",
                rel = "noopener noreferrer",
                shiny::icon("book"),
                shiny::span(" Documentation")
            )
        )
    ),
    shinydashboard::dashboardBody(
        # Custom font styling for the dashboard title
        shiny::tags$head(
            shiny::tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css?family=Inknut+Antiqua"),
            shiny::tags$style(shiny::HTML("
        .main-header .logo {
          font-family: 'Inknut Antiqua', sans-serif;
        }
        #debug_panel_wrapper {
          position: fixed;
          bottom: 10px;
          right: 10px;
          z-index: 9999;
          background: rgba(0, 0, 0, 0.85);
          color: #0f0;
          font-family: 'Courier New', monospace;
          font-size: 11px;
          padding: 8px 12px;
          border-radius: 6px;
          max-width: 320px;
          pointer-events: none;
        }
      "))
        ),
        shiny::uiOutput("debug_panel"),
        shiny::fluidRow(
            shiny::column(
                width = 9,
                shinydashboard::tabItems(
                    shinydashboard::tabItem(
                        tabName = "raw_data",
                        shiny::fluidRow(
                            shiny::uiOutput("all_tables_ui")
                        )
                    ),
                    shinydashboard::tabItem(
                        tabName = "Integration",
                        shiny::fluidRow(
                            shiny::uiOutput("integration_ui")
                        )
                    )
                )
            ),
            shiny::column(
                width = 3,
                data_ui_object
            )
        )
    )
)

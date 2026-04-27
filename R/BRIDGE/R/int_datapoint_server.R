#' @export
int_datapoints_server <- function(input, output, session, data_combined, reference_data_names) {
    get_data_combined <- function() {
        if (is.function(data_combined)) data_combined() else data_combined
    }

    get_reference_cols <- function() {
        if (is.function(reference_data_names)) reference_data_names() else reference_data_names
    }

    processed_data <- reactive({
        shiny::req(input$scale_integration)
        scale_input <- input$scale_integration
        raw_data <- get_data_combined()
        shiny::req(!is.null(raw_data), is.data.frame(raw_data), nrow(raw_data) > 0)
        shiny::validate(
            shiny::need(all(c("Gene_Name", "source", "unique_id") %in% colnames(raw_data)),
                "Integrate data first to enable the datapoints plot."
            )
        )
        if (scale_input == "Continous") {
            raw_data
        } else if (scale_input == "Log-scale") {
            log2_transform(raw_data)
        }
    })

    selected_data <- reactive({
        shiny::req(input$search_gene_integration)
        gene <- stringr::str_to_lower(trimws(input$search_gene_integration))
        data <- processed_data()
        data <- subset(data, stringr::str_to_lower(Gene_Name) %in% gene | stringr::str_to_lower(gsub("_.*", "", Gene_Name)) %in% gene)
        validate(need(nrow(data) > 0, "No data found for the entered gene(s)"))
        data
    })

    data_long <- reactive({
        data <- selected_data()
        ref_cols <- get_reference_cols()
        shiny::req(!is.null(ref_cols), length(ref_cols) > 0)
        shiny::validate(
            shiny::need(all(ref_cols %in% colnames(data)),
                "Selected integration columns are not available in the current integrated table."
            )
        )

        # Clean stage names
        unique_datapoints <- ref_cols %>%
            gsub("_[0-9]+$", "", .) %>%
            gsub("[_.-]+$", "", .) %>%
            unique() %>%
            str_sort(decreasing = FALSE, numeric = TRUE)         

        # Reshape
        data_long <- data %>%
            pivot_longer(cols = all_of(ref_cols), names_to = "Stage", values_to = "Expression") %>%
            mutate(
                StageGroup = gsub("[0-9]+$", "", Stage) %>% gsub("[_.-]+$", "", .),
                #StageGroup = factor(StageGroup, levels = mixedrank(unique(Stage)))
                StageGroup = factor(StageGroup, levels = unique_datapoints)
            )
        
        data_long
    })
    data_avg <- reactive({
        data_long() %>%
            group_by(source, StageGroup, unique_id) %>%
            summarize(MeanExpression = mean(Expression, na.rm = TRUE), .groups = "drop")
    })

    integration_datapoints_plot <- reactive({
        long <- data_long()
        avg <- data_avg()

        ggplot2::ggplot() +
            geom_point(data = long, aes(x = StageGroup, y = Expression, color = unique_id), size = 3, alpha = 0.5) +
            geom_line(data = avg, aes(x = StageGroup, y = MeanExpression, color = unique_id, group = unique_id), linewidth = 1.2) +
            geom_point(data = avg, aes(x = StageGroup, y = MeanExpression, color = unique_id), size = 4, shape = 17) +
            facet_wrap(~source, ncol = 2, scales = "free_y") +
            labs(x = "Stage", y = "Expression", color = "Gene") +
            theme_minimal() +
            theme(
                text = element_text(size = 14),
                axis.text.x = element_text(angle = 45, hjust = 1),
                strip.text = element_text(face = "bold")
            )
    })

    output$integration_datapoints_plot <- shiny::renderPlot({
        integration_datapoints_plot()
    })

    output$integration_datapoints_download_ui <- shiny::renderUI({
        plot_download_controls(session$ns, "integration_datapoints")
    })

    register_plot_download(
        input = input,
        output = output,
        session = session,
        id_prefix = "integration_datapoints",
        filename_prefix = "integration_datapoints",
        plot_fun = function() integration_datapoints_plot(),
        width = 11,
        height = 7
    )
}

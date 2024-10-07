library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(shinydashboard)

# Load data
joined_cholera <- read.csv("joined_cholera")

# UI
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        font-family: Georgia, serif;
        background-color: #f0f4f8; 
      }
      h1, h2, h3 {
        font-weight: bold;
        color: #2c3e50;  
      }
    "))),
  
  # Title with logo
  titlePanel(
    div(
      img(src = "logo.png", height = 100),
      HTML("<b>CholeraWatch</b>")
    )),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("country", "Select Country:", choices = unique(joined_cholera$Country)),
      sliderInput("yearRange", "Select Year Range:", 
                  min = min(joined_cholera$Year), 
                  max = max(joined_cholera$Year), 
                  value = c(min(joined_cholera$Year), max(joined_cholera$Year)),
                  step = 1,
                  sep = "")),  # To remove the comma separator from the years
    
    mainPanel(
      fluidRow(
        column(4, infoBoxOutput("totalCasesBox")),
        column(4, infoBoxOutput("totalDeathsBox")),
        column(4, infoBoxOutput("cfrBox"))
      ),
      br(),  # Line break to add space between the info boxes and the plot
      
      div(class = "mainPanel-spacing",  # to increase space between info boxes and plot
          plotOutput("choleraPlot")),
      
      # New UI output for the country summary
      uiOutput("countrySummary")  # This line will display the summary
    ))
)

# Server
server <- function(input, output) {
  
  # Filtered dataset based on input
  filtered_data <- reactive({
    joined_cholera %>%
      filter(Country == input$country, Year >= input$yearRange[1], Year <= input$yearRange[2])
  })
  
  # InfoBox for Total Cases
  output$totalCasesBox <- renderInfoBox({
    data <- filtered_data()
    total_cases <- sum(data$Number_cases, na.rm = TRUE)
    
    infoBox(
      "Total Cases", total_cases, icon = icon("bar-chart"),
      color = "blue", fill = TRUE
    )
  })
  
  # InfoBox for Total Deaths
  output$totalDeathsBox <- renderInfoBox({
    data <- filtered_data()
    total_deaths <- sum(data$Number_deaths, na.rm = TRUE)
    
    infoBox(
      "Total Deaths", total_deaths, icon = icon("times-circle"),
      color = "orange", fill = TRUE
    )
  })
  
  # InfoBox for Case-Fatality Rate
  output$cfrBox <- renderInfoBox({
    data <- filtered_data()
    total_cases <- sum(data$Number_cases, na.rm = TRUE)
    total_deaths <- sum(data$Number_deaths, na.rm = TRUE)
    cfr <- ifelse(total_cases > 0, (total_deaths / total_cases) * 100, 0)
    
    infoBox(
      "Case-Fatality Rate", paste0(round(cfr, 2), "%"), icon = icon("percent"),
      color = "green", fill = TRUE
    )
  })
  
  # Reshape the dataset for faceting
  output$choleraPlot <- renderPlot({
    data <- filtered_data()
    
    validate(
      need(nrow(data) > 0, "No data available for this selection")
    )
    
    # Prepare the data for faceting
    long_data <- data %>%
      pivot_longer(cols = c(Number_cases, Number_deaths, Case_fatality), 
                   names_to = "Metric", values_to = "Value") %>%
      mutate(Metric = recode(Metric,
                             "Number_cases" = "Total Cases",
                             "Number_deaths" = "Total Fatalities",
                             "Case_fatality" = "Case-Fatality Rate (%)"))
    
    # Plot with faceting
    long_data$Metric <- factor(long_data$Metric, 
                               levels = c("Total Cases", "Total Fatalities", "Case-Fatality Rate (%)"))
    
    ggplot(long_data, aes(x = Year, y = Value, group = Metric, color = Metric)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      facet_grid(Metric ~ ., scales = "free_y", switch = "y") +  # Facet for each metric
      labs(title = paste("Trends in Cholera Cases, Fatalities, and Case-Fatality Rate in", input$country),
           x = "Year", y = "Count / Percentage") +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),   # Bold title, center-aligned
        strip.text = element_text(face = "bold"),                # Bold facet labels
        axis.title.x = element_text(face = "bold"),              # Bold x-axis label
        axis.title.y = element_text(face = "bold"),              # Bold y-axis label
        legend.position = "none",                                # Remove legend
        strip.background = element_rect(fill = "gray90"),        # Background for facet labels
        panel.spacing = unit(1, "lines")                         # Increase space between panels
      ) +
      scale_color_manual(values = c("Total Cases" = "darkblue", 
                                    "Total Fatalities" = "darkorange", 
                                    "Case-Fatality Rate (%)" = "darkgreen")) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10), 
                         labels = scales::number_format(accuracy = 1))   # No decimal points for years
  })
  
  # Generate the country summary
  output$countrySummary <- renderUI({
    data <- filtered_data()
    
    # Check if data is available before performing calculations
    if (nrow(data) == 0) {
      return(HTML("<b>No data available for this selection.</b>"))
    }
    
    # Summary calculations
    max_cases <- max(data$Number_cases, na.rm = TRUE)
    min_cases <- min(data$Number_cases, na.rm = TRUE)
    year_max_cases <- data$Year[which.max(data$Number_cases)]
    year_min_cases <- data$Year[which.min(data$Number_cases)]
    
    max_fatalities <- max(data$Number_deaths, na.rm = TRUE)
    min_fatalities <- min(data$Number_deaths, na.rm = TRUE)
    year_max_fatalities <- data$Year[which.max(data$Number_deaths)]
    year_min_fatalities <- data$Year[which.min(data$Number_deaths)]
    
    max_cfr <- max(data$Case_fatality, na.rm = TRUE)
    min_cfr <- min(data$Case_fatality, na.rm = TRUE)
    year_max_cfr <- data$Year[which.max(data$Case_fatality)]
    year_min_cfr <- data$Year[which.min(data$Case_fatality)]
    
    # Display the summary in an HTML format
    HTML(paste0(
      "<b>Summary:</b><br><br>",
      "1. Number of Cases:<br>",
      "Highest: ", max_cases, " (Year: ", year_max_cases, ")<br>",
      "Lowest: ", min_cases, " (Year: ", year_min_cases, ")<br><br>",
      "2. Number of Fatalities:<br>",
      "Highest: ", max_fatalities, " (Year: ", year_max_fatalities, ")<br>",
      "Lowest: ", min_fatalities, " (Year: ", year_min_fatalities, ")<br><br>",
      "3. Case-Fatality Rate (CFR):<br>",
      "Highest: ", max_cfr, "% (Year: ", year_max_cfr, ")<br>",
      "Lowest: ", min_cfr, "% (Year: ", year_min_cfr, ")"
    ))
  })
}

# Run the App
shinyApp(ui = ui, server = server)

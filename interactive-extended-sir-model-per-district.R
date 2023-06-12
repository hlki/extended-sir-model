# Additional libraries
library(shiny)

# UI
ui <- fluidPage(
  titlePanel("Interactive SIR model"),
  sidebarLayout(
    sidebarPanel(
      numericInput("T", "Transmissibility", value = 0.02, min = 0, max = 1, step = 0.01),
      numericInput("DI", "Duration of Infectiousness", value = 10, min = 1, max = 100, step = 1),
      numericInput("PNH", "% Needing Hospitalization", value = 0.138, min = 0, max = 1, step = 0.01),
      numericInput("PNIC", "% Needing ICU Care", value = 0.047, min = 0, max = 1, step = 0.01),
      numericInput("MR", "Mortality rate", value = 0.0326, min = 0, max = 1, step = 0.01),
      numericInput("CR", "Contact rate", value = 5.0, min = 0, max = 100, step = 0.1),
      numericInput("VR", "Vaccination rate", value = 0.00102, min = 0, max = 1, step = 0.0001),
      sliderInput("starting_time", "Starting time for CR increase", value = 550, min = 0, max = 900, step = 5),
      sliderInput("time_interval", "Time interval for CR increase", value = 25, min = 0, max = 900, step = 5),
      sliderInput("CR_increase", "Contact rate increase per time interval", value = 2.2, min = 0, max = 10, step = 0.1)
    ),
    mainPanel(
      plotOutput("sirPlot"),
      plotOutput("comparisonPlot"),
      plotOutput("healthcareImpactPlot")
    )
  )
)

# Server logic
server <- function(input, output) {
  
  cr_function <- reactive({
    function(time) {
      base_CR <- input$CR
      periods_passed <- max(0, (time - input$starting_time) %/% input$time_interval)
      CR_increase <- periods_passed * input$CR_increase
      new_CR <- base_CR + CR_increase
      return(new_CR)
    }
  })
  
  sir_model <- reactive({
    function(time, state, parameters) {
      with(as.list(c(state, parameters)), {
        CR <- cr_function()(time)  # use the time-dependent function here
        N <- S + I + R + V
        ECR <- T * CR
        RR <- 1 / DI
        
        if(time >= 410) {
          dS <- -ECR * S * I / N - VR * S
          dV <- VR * S
        } else {
          dS <- -ECR * S * I / N
          dV <- 0
        }
        
        dI <- ECR * S * I / N - RR * I
        dR <- RR * I
        return(list(c(dS, dI, dR, dV)))
      })
    }
  })
  
  # Output for SIR plot
  output$sirPlot <- renderPlot({
    input_data <- list(T = input$T, DI = input$DI, PNH = input$PNH, PNIC = input$PNIC, MR = input$MR, CR = input$CR, VR = input$VR) # converted to list
    initial_state <- c(S = 1152400, I = 11, R = 1, V = 0)
    time_steps <- seq_len(900)
    time_series <- ode(y = initial_state, times = time_steps, func = sir_model(), parms = input_data)
    
    time_series_df <- as.data.frame(time_series)
    time_series_df$NH <- time_series_df$I * input_data$PNH
    time_series_df$NIC <- time_series_df$I * input_data$PNIC
    time_series_df$D <- time_series_df$I * input_data$MR
    
    ggplot(time_series_df, aes(x = time)) +
      geom_line(aes(y = S, colour = "S")) +
      geom_line(aes(y = I, colour = "I")) +
      geom_line(aes(y = R, colour = "R")) +
      geom_line(aes(y = V, colour = "V")) +
      labs(x = "Time", y = "Count", colour = "Compartment", title = "Susceptible-Infected-Recovered-Vaccinated Model") +
      theme_minimal()
  })
  
  # Output for comparison plot
  output$comparisonPlot <- renderPlot({
    input_data <- list(T = input$T, DI = input$DI, PNH = input$PNH, PNIC = input$PNIC, MR = input$MR, CR = input$CR, VR = input$VR) # converted to list
    initial_state <- c(S = 1152400, I = 11, R = 1, V = 0)
    time_steps <- seq_len(900)
    time_series <- ode(y = initial_state, times = time_steps, func = sir_model(), parms = input_data)
    
    time_series_df <- as.data.frame(time_series)
    time_series_df$NH <- time_series_df$I * input_data$PNH
    time_series_df$NIC <- time_series_df$I * input_data$PNIC
    time_series_df$D <- time_series_df$I * input_data$MR
    
    plot(NA, xlim = range(time_steps), ylim = range(c(time_series_df$I, observed_data$total_infected)), xlab = "Time", ylab = "Infected")
    lines(time_series_df$time, time_series_df$I, lwd = 2)
    points(observed_data$time, observed_data$total_infected, pch = 19, col = "red")
    legend("topleft", legend = c("Model predictions", "Observed data"), lty = c(1, NA), pch = c(NA, 19), col = c("black", "red"))
  })
  
  # Output for healthcare impact plot
  output$healthcareImpactPlot <- renderPlot({
    input_data <- list(T = input$T, DI = input$DI, PNH = input$PNH, PNIC = input$PNIC, MR = input$MR, CR = input$CR, VR = input$VR) # converted to list
    initial_state <- c(S = 1152400, I = 11, R = 1, V = 0)
    time_steps <- seq_len(900)
    time_series <- ode(y = initial_state, times = time_steps, func = sir_model(), parms = input_data)
    
    time_series_df <- as.data.frame(time_series)
    time_series_df$NH <- time_series_df$I * input_data$PNH
    time_series_df$NIC <- time_series_df$I * input_data$PNIC
    time_series_df$D <- time_series_df$I * input_data$MR
    
    time_series_melt <- melt(time_series_df, id = "time")
    
    subset <- time_series_melt %>% filter(variable %in% c("NH", "NIC", "D"))
    
    ggplot(subset, aes(x = time, y = value, colour = variable)) +
      geom_line() +
      scale_color_manual(values = c("NH" = "green", "NIC" = "blue", "D" = "red")) +
      labs(x = "Time", y = "Count", colour = "Compartment", title = "Impact on healthcare system and mortality") +
      theme_minimal()
  })
}

# Run the app
shinyApp(ui = ui, server = server)

# Function to check, install and load packages
load_package <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    library(package_name, character.only = TRUE)
  }
}

# List of necessary packages
packages <- c("deSolve", "ggplot2", "reshape2", "tidyverse", "shiny")

# Apply the function to each package
lapply(packages, load_package)


# UI
ui <- fluidPage(
  titlePanel("Інтерактивна SIRV модель"),
  sidebarLayout(
    sidebarPanel(
      numericInput("a", "Імовірність зараження", value = 0.02, min = 0, max = 1, step = 0.01),
      numericInput("d", "Тривалість інфекційності", value = 10, min = 1, max = 100, step = 1),
      numericInput("pNH", "Особи які потребують госпіталізації (у %)", value = 0.112, min = 0, max = 1, step = 0.01),
      numericInput("pNICU", "Особи які потребують реанімаційної допомоги (у %)", value = 0.047, min = 0, max = 1, step = 0.01),
      numericInput("m", "Коефіцієнт смертності (на 1 мільйон)", value = 0.026, min = 0, max = 1, step = 0.01),
      numericInput("c", "Середня кількість контактів", value = 6.7, min = 0, max = 100, step = 0.1),
      numericInput("nu", "Рівень вакцинації", value = 0.00102, min = 0, max = 1, step = 0.0001),
      sliderInput("time_steps", "Час", value = 900, min = 10, max = 1095, step = 5),
      checkboxInput("use_dynamic_c", "Використовувати динамічне збільшення контактів", value = FALSE),
      conditionalPanel(
        condition = "input.use_dynamic_c == true",
        sliderInput("starting_time", "Момент часу для збільшення кількості контактів", value = 165, min = 0, max = 900, step = 5),
        sliderInput("time_interval", "Часовий інтервал для збільшення кількості контактів", value = 120, min = 0, max = 900, step = 5),
        sliderInput("c_increase", "Зростання частоти контактів за інтервал часу", value = 2, min = 0, max = 10, step = 0.1)
      ),
    ),
    mainPanel(
      plotOutput("sirPlot"),
      plotOutput("healthcareImpactPlot"),
      plotOutput("infectedComparisonPlot"),
      plotOutput("vaccinatedComparisonPlot"),
      plotOutput("recoveredComparisonPlot"),
      plotOutput("deathComparisonPlot")
    )
  )
)

# Load the observed data
observed_data <- read.csv("observed_data_per_district.csv")

# Server logic
server <- function(input, output, session) {
  observe({
    updateSliderInput(session, "starting_time", max = input$time_steps)
    updateSliderInput(session, "time_interval", max = input$time_steps)
  })
  
  cr_function <- reactive({
    function(time) {
      base_c <- input$c
      periods_passed <- max(0, (time - input$starting_time) %/% input$time_interval)
      c_increase <- periods_passed * input$c_increase
      new_c <- base_c + c_increase
      return(new_c)
    }
  })
  
  sir_model <- reactive({
    function(time, state, parameters) {
      with(as.list(c(state, parameters)), {
        # Calculate dynamic or static contact rate
        c_value <- if(input$use_dynamic_c) cr_function()(time) else input$c
        
        # Calculate model parameters
        N <- S + I + R + V
        beta <- a * c_value
        gamma <- 1 / d
        nu_value <- if(time >= 410) nu else 0
        
        # Print parameters for each time step (for debugging or review)
        formatted_output <- sprintf("t: %0.2f, S: %0.2f, I: %0.2f, R: %0.2f, V: %0.2f, beta: %0.4f, gamma: %0.4f, nu: %0.4f, NH: %0.2f, NICU: %0.2f, D: %0.2f",
                                    time, S, I, R, V, beta, gamma, nu_value, I * pNH, I * pNICU, I * m)
        print(formatted_output)
        
        # Differential equations
        dS <- -beta * S * I / N - nu_value * S
        dI <- beta * S * I / N - gamma * I
        dR <- gamma * I
        dV <- nu_value * S
        
        return(list(c(dS, dI, dR, dV)))
      })
    }
  })
  
  # Output for SIR plot
  output$sirPlot <- renderPlot({
    input_data <- list(a = input$a, d = input$d, pNH = input$pNH, pNICU = input$pNICU, m = input$m, c = input$c, nu = input$nu) # converted to list
    initial_state <- c(S = 1152400, I = 11, R = 1, V = 0)
    time_steps <- seq_len(input$time_steps)
    time_series <- ode(y = initial_state, times = time_steps, func = sir_model(), parms = input_data)
    
    time_series_df <- as.data.frame(time_series)
    time_series_df$NH <- time_series_df$I * input_data$pNH
    time_series_df$NICU <- time_series_df$I * input_data$pNICU
    time_series_df$D <- time_series_df$I * input_data$m
    
    ggplot(time_series_df, aes(x = time)) +
      geom_line(aes(y = S, colour = "S"), size = 1.5) + 
      geom_line(aes(y = I, colour = "I"), size = 1.5) +
      geom_line(aes(y = R, colour = "R"), size = 1.5) +
      geom_line(aes(y = V, colour = "V"), size = 1.5) +
      labs(x = "Час", y = "Кількість", colour = "", title = "Результати SIRV моделі") +
      theme_minimal() +
      theme(text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 14, face = "bold"),
            plot.title = element_text(size = 14, face = "bold"))
  })
  
  # Output for healthcare impact plot
  output$healthcareImpactPlot <- renderPlot({
    input_data <- list(a = input$a, d = input$d, pNH = input$pNH, pNICU = input$pNICU, m = input$m, c = input$c, nu = input$nu) # converted to list
    initial_state <- c(S = 1152400, I = 11, R = 1, V = 0)
    time_steps <- seq_len(input$time_steps)
    time_series <- ode(y = initial_state, times = time_steps, func = sir_model(), parms = input_data)

    time_series_df <- as.data.frame(time_series)
    time_series_df$NH <- time_series_df$I * input_data$pNH
    time_series_df$NICU <- time_series_df$I * input_data$pNICU
    time_series_df$D <- time_series_df$I * input_data$m

    time_series_melt <- melt(time_series_df, id = "time")

    subset <- time_series_melt %>% filter(variable %in% c("NH", "NICU", "D"))

    ggplot(subset, aes(x = time, y = value, colour = variable)) +
      geom_line(size = 1.5) +
      scale_color_manual(values = c("NH" = "green", "NICU" = "blue", "D" = "red")) +
      labs(x = "Час", y = "Кількість", colour = "", title = "Кількість хворих що потребують медичної допомоги та померлих від хвороби") +
      theme_minimal() +
      theme(text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 14, face = "bold"),
            plot.title = element_text(size = 14, face = "bold"))
  })

  # Output for infected comparison plot
  output$infectedComparisonPlot <- renderPlot({
    input_data <- list(a = input$a, d = input$d, pNH = input$pNH, pNICU = input$pNICU, m = input$m, c = input$c, nu = input$nu)
    initial_state <- c(S = 1152400, I = 11, R = 1, V = 0)
    time_steps <- seq_len(input$time_steps)
    time_series <- ode(y = initial_state, times = time_steps, func = sir_model(), parms = input_data)

    time_series_df <- as.data.frame(time_series)
    time_series_df$NH <- time_series_df$I * input_data$pNH
    time_series_df$NICU <- time_series_df$I * input_data$pNICU
    time_series_df$D <- time_series_df$I * input_data$m

    par(cex.lab = 1.5, font.lab = 1)
    plot(NA, xlim = range(time_steps), ylim = range(c(time_series_df$I, observed_data$total_infected)), xlab = "Час", ylab = "Інфіковані")
    lines(time_series_df$time, time_series_df$I, lwd = 2)
    points(observed_data$time, observed_data$total_infected, pch = 19, col = "red")
    legend("topleft", legend = c("Прогноз моделі", "Дані спостережень"), lty = c(1, NA), pch = c(NA, 19), col = c("black", "red"))

    title("Порівняння прогнозованої кількості інфікованих з спостереженнями")
  })


  # Output for vaccinated comparison plot
  output$vaccinatedComparisonPlot <- renderPlot({

    # Handle NA values
    observed_data <- observed_data[!is.na(observed_data$fully_vaccinated), ]

    input_data <- list(a = input$a, d = input$d, pNH = input$pNH, pNICU = input$pNICU, m = input$m, c = input$c, nu = input$nu)
    initial_state <- c(S = 1152400, I = 11, R = 1, V = 0)
    time_steps <- seq_len(input$time_steps)
    time_series <- ode(y = initial_state, times = time_steps, func = sir_model(), parms = input_data)

    time_series_df <- as.data.frame(time_series)
    time_series_df$NH <- time_series_df$I * input_data$pNH
    time_series_df$NICU <- time_series_df$I * input_data$pNICU
    time_series_df$D <- time_series_df$I * input_data$m

    par(cex.lab = 1.5, font.lab = 1)
    plot(NA, xlim = range(time_steps), ylim = range(c(time_series_df$V, observed_data$fully_vaccinated)), xlab = "Час", ylab = "Вакциновані")
    lines(time_series_df$time, time_series_df$V, lwd = 2)
    points(observed_data$time, observed_data$fully_vaccinated, pch = 19, col = "red")
    legend("topleft", legend = c("Прогноз моделі", "Дані спостережень"), lty = c(1, NA), pch = c(NA, 19), col = c("black", "red"))

    title("Порівняння прогнозованої кількості вакцинованих з спостереженнями")
  })

  # Output for recovered comparison plot
  output$recoveredComparisonPlot <- renderPlot({

    # Handle NA values
    observed_data <- observed_data[!is.na(observed_data$total_recovered), ]

    input_data <- list(a = input$a, d = input$d, pNH = input$pNH, pNICU = input$pNICU, m = input$m, c = input$c, nu = input$nu)
    initial_state <- c(S = 1152400, I = 11, R = 1, V = 0)
    time_steps <- seq_len(input$time_steps)
    time_series <- ode(y = initial_state, times = time_steps, func = sir_model(), parms = input_data)

    time_series_df <- as.data.frame(time_series)
    time_series_df$NH <- time_series_df$I * input_data$pNH
    time_series_df$NICU <- time_series_df$I * input_data$pNICU
    time_series_df$D <- time_series_df$I * input_data$m

    par(cex.lab = 1.5, font.lab = 1)
    plot(NA, xlim = range(time_steps), ylim = range(c(time_series_df$R, observed_data$total_recovered)), xlab = "Час", ylab = "Одужавші")
    lines(time_series_df$time, time_series_df$R, lwd = 2)
    points(observed_data$time, observed_data$total_recovered, pch = 19, col = "red")
    legend("topleft", legend = c("Прогноз моделі", "Дані спостережень"), lty = c(1, NA), pch = c(NA, 19), col = c("black", "red"))

    title("Порівняння прогнозованої кількості одужавших з спостереженнями")
  })

  # Output for death comparison plot
  output$deathComparisonPlot <- renderPlot({

    # Handle NA values
    observed_data <- observed_data[!is.na(observed_data$total_death), ]

    input_data <- list(a = input$a, d = input$d, pNH = input$pNH, pNICU = input$pNICU, m = input$m, c = input$c, nu = input$nu)
    initial_state <- c(S = 1152400, I = 11, R = 1, V = 0)
    time_steps <- seq_len(input$time_steps)
    time_series <- ode(y = initial_state, times = time_steps, func = sir_model(), parms = input_data)

    time_series_df <- as.data.frame(time_series)
    time_series_df$NH <- time_series_df$I * input_data$pNH
    time_series_df$NICU <- time_series_df$I * input_data$pNICU
    time_series_df$D <- time_series_df$I * input_data$m

    par(cex.lab = 1.5, font.lab = 1)
    plot(NA, xlim = range(time_steps), ylim = range(c(time_series_df$D, observed_data$total_death)), xlab = "Час", ylab = "Померлі")
    lines(time_series_df$time, time_series_df$D, lwd = 2)
    points(observed_data$time, observed_data$total_death, pch = 19, col = "red")
    legend("topleft", legend = c("Прогноз моделі", "Дані спостережень"), lty = c(1, NA), pch = c(NA, 19), col = c("black", "red"))

    title("Порівняння прогнозованої кількості померлих від хвороби з спостереженнями")
  })
}

# Run the app
shinyApp(ui = ui, server = server)


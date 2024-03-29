# Required packages
required_packages <- c("tidyverse", "reshape2", "deSolve")

for(pkg in required_packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

library(tidyverse)
library(reshape2)
library(deSolve)

# Input data
input_data <- tibble(
  # Transmissibility 
  T = 0.02,
  # Duration of Infectiousness
  DI = 10,
  # % Needing Hospitalization
  PNH = 0.138,
  # % Needing ICU Care
  PNIC = 0.047,
  # Mortality rate
  MR = 0.0326,
  # Contact rate
  CR = 7.4,
  # Vaccination rate
  VR = 0.00102
)

# Initial state
initial_state <- c(S = 1152400, I = 11, R = 1, V = 0)

# Load the observed data
observed_data <- read.csv("observed_data_per_district.csv")

# Define a function for CR that varies with time
cr_function <- function(time) {
  base_CR <- input_data$CR
  
  # Calculate the number of 20-day periods that have passed, minus one
  #periods_passed <- (time - 1) %/% 20
  # Calculate the number of 20-day periods that have passed since the 100th time step
  periods_passed <- max(0, (time - 301) %/% 100)
  
  # Calculate the increase in the contact rate
  CR_increase <- periods_passed * 3.0
  
  # Add the increase to the base contact rate
  new_CR <- base_CR + CR_increase
  
  return(new_CR)
}


# Define the SIR model
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    CR <- cr_function(time)  # use the time-dependent function here
    N <- S + I + R + V
    ECR <- T * CR
    RR <- 1 / DI
    
    # Check if the current time is after the 410th time step
    if(time >= 410) {
      dS <- -ECR * S * I / N - VR * S
      dV <- VR * S
    } else {
      # Before 410 time step, no vaccinations are done
      dS <- -ECR * S * I / N
      dV <- 0
    }
    
    dI <- ECR * S * I / N - RR * I
    dR <- RR * I
    return(list(c(dS, dI, dR, dV)))
  })
}

# Time steps
time_steps <- seq_len(900)

# Solve the differential equations
time_series <- ode(y = initial_state, times = time_steps, func = sir_model, parms = input_data)

# Convert the output to a data frame and calculate NH, NIC and D
time_series_df <- as.data.frame(time_series)
time_series_df$NH <- time_series_df$I * input_data$PNH
time_series_df$NIC <- time_series_df$I * input_data$PNIC
time_series_df$D <- time_series_df$I * input_data$MR

# Reshape the time series data for plotting
time_series_long <- melt(time_series_df, id = "time")

# Create the SIR plot
ggplot(time_series_long, aes(x = time, y = value, colour = variable)) +
  geom_line() +
  scale_color_manual(values = c("S" = "green", "I" = "red", "R" = "blue", "V" = "pink")) +
  labs(x = "Time", y = "Count", colour = "Compartment", title = "Susceptible-Infected-Recovered-Vaccinated Model") +
  theme_minimal()

# Subset the dataframe for the second plot
subset <- time_series_long %>% filter(variable %in% c("NH", "NIC", "D"))

# Create the impact on healthcare system and mortality plot
ggplot(subset, aes(x = time, y = value, colour = variable)) +
  geom_line() +
  scale_color_manual(values = c("NH" = "green", "NIC" = "blue", "D" = "red")) +
  labs(x = "Time", y = "Count", colour = "Compartment", title = "Impact on healthcare system and mortality") +
  theme_minimal()

# Run the model
model_output <- ode(y = initial_state, times = time_steps, func = sir_model, parms = input_data)
model_output_df <- as.data.frame(model_output)

# Check if there are non-finite values in the 'I' column of observed_data
print(any(!is.finite(observed_data$total_infected)))

# If the above line prints "TRUE", then there are non-finite values. 
# Remove these rows with the following code:

observed_data <- observed_data[is.finite(observed_data$total_infected), ]

# Then, proceed with the plotting code
plot(NA, xlim = range(time_steps), ylim = range(c(model_output_df$I, observed_data$total_infected)), xlab = "Time", ylab = "Infected")
lines(model_output_df$time, model_output_df$I, lwd = 2)
points(observed_data$time, observed_data$total_infected, pch = 19, col = "red")
legend("topleft", legend = c("Model predictions", "Observed data"), lty = c(1, NA), pch = c(NA, 19), col = c("black", "red"))
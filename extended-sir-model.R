# Required packages
# install.packages(c("tidyverse", "reshape2", "deSolve"))
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
  T = 0.02,
  DI = 10,
  PNH = 0.138,
  PNIC = 0.047,
  MR = 0.0326,
  CR = 13.4,
  vaccination_rate = 0.00102
)

# Initial state
initial_state <- c(S = 39701744, I = 25789, R = 34987, V = 15142679)

# Load the observed data
observed_data <- read.csv("observed_data.csv")

# Define a function for CR that varies with time
cr_function <- function(time) {
  if (time >= 1 && time <= 188) {
    return(6.7)  # Modify this function according to your needs
  }
  else {
    return(input_data$CR)
  }
}

# Define the SIR model
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    CR <- cr_function(time)  # use the time-dependent function here
    N <- S + I + R + V
    ECR <- T * CR
    RR <- 1 / DI
    dS <- -ECR * S * I / N - vaccination_rate * S
    dI <- ECR * S * I / N - RR * I
    dR <- RR * I
    dV <- vaccination_rate * S
    return(list(c(dS, dI, dR, dV)))
  })
}

# Time steps
time_steps <- seq_len(365)

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

# Define the cost function
cost_function <- function(params) {
  input_data$CR <- params[1]
  input_data$vaccination_rate <- params[2]
  model_output <- ode(y = initial_state, times = time_steps, func = sir_model, parms = input_data)
  model_output_df <- as.data.frame(model_output)
  observed_data <- observed_data
  sum_of_squared_errors <- sum((observed_data$S - model_output_df$S)^2 
                               + (observed_data$I - model_output_df$I)^2
                               + (observed_data$R - model_output_df$R)^2
                               + (observed_data$V - model_output_df$V)^2
  )
  return(sum_of_squared_errors)
}

# Optimize the model
optimization_result <- optim(c(input_data$CR, input_data$vaccination_rate), cost_function, method = "BFGS")

# Print the best-fit parameters
print(optimization_result$par)

# Run the model
model_output <- ode(y = initial_state, times = time_steps, func = sir_model, parms = input_data)
model_output_df <- as.data.frame(model_output)

# Check if there are non-finite values in the 'I' column of observed_data
print(any(!is.finite(observed_data$I)))

# If the above line prints "TRUE", then there are non-finite values. 
# Remove these rows with the following code:

observed_data <- observed_data[is.finite(observed_data$I), ]

# Then, proceed with the plotting code
plot(NA, xlim = range(time_steps), ylim = range(c(model_output_df$I, observed_data$I)), xlab = "Time", ylab = "Infected")
lines(model_output_df$time, model_output_df$I, lwd = 2)
points(observed_data$time, observed_data$I, pch = 19, col = "red")
legend("topleft", legend = c("Model predictions", "Observed data"), lty = c(1, NA), pch = c(NA, 19), col = c("black", "red"))
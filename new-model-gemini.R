# Multiwave SIR Model in R

library(deSolve)
library(changepoint)

# SIR Model Equations
SIR <- function(t, state, params) {
  with(as.list(c(state, params)), {
    S <- state[1]
    I <- state[2]
    N <- params$N
    beta <- params$beta
    gamma <- params$gamma
    dS <- -beta * (S * I) / N
    dI <- beta * (S * I) / N - gamma * I
    dR <- gamma * I
    return(c(dS, dI, dR))
  })
}

# Change Point Detection
change_point <- function(infection_data) {
  # Use cpts function to identify change points in infection rate
  changepoint_model <- cpts(infection_data, method = "PELT")
  change_points <- changepoint_model$segmentMeans$where
  return(change_points)
}

# Dynamic Parameter Adjustment (Example)
update_parameters <- function(current_time, change_points, params, intervention_data) {
  # Check if current time corresponds to a change point
  if (current_time %in% change_points) {
    # Update beta based on intervention (replace with your logic)
    if (intervention_data$lockdown[current_time]) {
      params$beta <- params$beta * 0.5 # Reduce transmission rate during lockdown
    }
  }
  return(params)
}

# Simulation Loop
SIR_model <- function(initial_state, params, time_points, intervention_data = NULL) {
  # Define empty lists to store SIR values
  susceptible <- infected <- recovered <- vector(length = length(time_points))
  
  # Start simulation loop
  change_points_detected <- change_point(infection_data$infections)
  for (i in 1:length(time_points)) {
    current_time <- time_points[i]
    # Update parameters if a change point is detected
    params <- update_parameters(current_time, change_points_detected, params, intervention_data)
    # Solve SIR model for current time step
    model_out <- solve.ode(y = initial_state, times = current_time, func = SIR, params = params)
    # Store SIR values for this time step
    susceptible[i] <- model_out[nrow(model_out), 1]
    infected[i] <- model_out[nrow(model_out), 2]
    recovered[i] <- model_out[nrow(model_out), 3]
    # Update initial state for next iteration based on model output
    initial_state <- c(susceptible[i], infected[i], recovered[i])
  }
  
  # Return simulation results
  return(list(time = time_points, susceptible = susceptible, infected = infected, recovered = recovered))
}

# Example Usage (replace with your data and desired interventions)
# Initial conditions and parameters
initial_state <- c(S0 = 1000, I0 = 1, R0 = 0)
params <- c(N = 10000, beta = 0.2, gamma = 0.1)

# Time points for simulation
time_points <- seq(from = 0, to = 50, by = 1)

# Intervention data (example - lockdown during specific period)
intervention_data <- list(lockdown = c(rep(FALSE, 20), rep(TRUE, 15), rep(FALSE, 15)))

# Run simulation
simulation_results <- SIR_model(initial_state, params, time_points, intervention_data)

# Visualize results (using ggplot2 - replace with your preferred package)
library(ggplot2)
ggplot(data = simulation_results, aes(x = time)) +
  geom_line(aes(y = susceptible), color = "blue", linetype = "solid", size = 1) +
  geom_line(aes(y = infected), color = "red", linetype = "solid", size = 1) +
  geom_line(aes(y = recovered), color = "green", linetype = "solid", size = 1) +
  labs(title = "SIR Model with Intervention (Lockdown)",
       x = "Time", y = "Population")

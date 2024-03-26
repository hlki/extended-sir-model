library(deSolve)
library(changepoint)

# SIR model equations
sir_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I
    return(list(c(dS, dI, dR)))
  })
}

# Initial conditions
initial_conditions <- c(S = 9900000, I = 100000, R = 0)  # Assuming a population of 10 million

# Model parameters
parameters <- c(beta = 0.5, gamma = 0.2)

# Change-point detection
detect_changepoint <- function(beta_values) {
  cpt <- cpt.mean(beta_values)
  changepoints <- cpt@cpts
  return(changepoints)
}

# Dynamic parameter adjustment
update_parameters <- function(changepoints, current_time, beta, gamma) {
  if (length(changepoints) > 0 && current_time > changepoints[1]) {
    # Update parameters based on the detected change point
    beta <- beta - 0.1  # Assuming a reduction of 0.1 in transmission rate
    gamma <- gamma + 0.05  # Assuming an increase of 0.05 in recovery rate
    changepoints <- changepoints[-1]  # Remove the handled change point
  }
  return(list(beta = beta, gamma = gamma, changepoints = changepoints))
}

# Simulation loop
times <- seq(0, 100, by = 0.1)
output <- list(time = times[1], S = initial_conditions[1], I = initial_conditions[2], R = initial_conditions[3])
changepoints <- numeric(0)  # Initialize an empty vector for change points
beta_values <- numeric(0)   # Initialize an empty vector to store beta values

for (t in times[-1]) {
  current_state <- c(output$S[length(output$S)], output$I[length(output$I)], output$R[length(output$R)])
  updated_params <- update_parameters(changepoints, t, parameters$beta, parameters$gamma)
  changepoints <- updated_params$changepoints
  parameters <- c(beta = updated_params$beta, gamma = updated_params$gamma)
  
  new_state <- ode(y = current_state, times = c(t - 0.1, t), func = sir_model, parms = parameters)
  output$time <- c(output$time, t)
  output$S <- c(output$S, new_state[2, 1])
  output$I <- c(output$I, new_state[2, 2])
  output$R <- c(output$R, new_state[2, 3])
  
  # Monitor beta values and detect change points
  beta_values <- c(beta_values, parameters$beta)
  changepoints <- c(changepoints, detect_changepoint(beta_values))
}

# Visualization
plot(output$time, output$S, type = "l", col = "green", xlab = "Time", ylab = "Number of Individuals", ylim = c(0, 10000000))
lines(output$time, output$I, col = "red")
lines(output$time, output$R, col = "blue")
legend("topright", c("Susceptible", "Infected", "Recovered"), col = c("green", "red", "blue"), lty = 1)
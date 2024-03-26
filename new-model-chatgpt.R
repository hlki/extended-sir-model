# Load necessary packages
if (!requireNamespace("deSolve", quietly = TRUE)) install.packages("deSolve")
if (!requireNamespace("changepoint", quietly = TRUE)) install.packages("changepoint")
library(deSolve)
library(changepoint)

# Define the SIR model
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

# Initial conditions
initial_state <- c(S = 418581000, I = 16, R = 0)  # 1% of the population initially infected
parameters <- c(beta = 0.35, gamma = 0.1)     # Initial values for beta and gamma
N <- 1                                        # Total population normalized to 1
times <- seq(0, 500, by = 1)                  # Time span for the simulation

# Simulate the initial wave
solution <- ode(y = initial_state, times = times, func = sir_model, parms = parameters)

# Hypothetical infection rate data showing changes
infection_rates <- c(rep(0.35, 80), rep(0.15, 80))  # Assume a significant drop in beta after 80 days

# Detect change points in the hypothetical infection rate data
cpt <- cpt.mean(infection_rates)

# Adjust parameters based on detected change point
if (length(cpt@cpts) > 0) {
  parameters$beta <- 0.15  # Adjusting beta based on change point detection
}

# Re-run the simulation with updated parameters from the point of intervention
solution_post_intervention <- ode(y = solution[81, -1], times = seq(81, 160, by = 1), func = sir_model, parms = parameters)

# Combine solutions for plotting
combined_solution <- rbind(solution[1:80, ], solution_post_intervention)

# Plot the results
plot(combined_solution, main = "Enhanced SIR Model Simulation with Intervention")
lines(solution_post_intervention, col = "red")
legend("right", legend = c("Before Intervention", "After Intervention"), col = c("black", "red"), lty = 1)

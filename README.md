# extended-sir-model
Extended SIR Model with Vaccination in R
This repository contains an R script for simulating the spread of an infectious disease using an extended SIR (Susceptible-Infected-Recovered) model that also incorporates vaccination.

The model is extended from the traditional SIR model to take into account the following factors:

Contact Rate (CR): The rate at which an infected individual comes into contact with a susceptible individual. This can be adjusted over time to simulate varying conditions such as lockdowns or social distancing measures.

Vaccination: The model considers the vaccinated population and the rate of vaccination.

Hospitalization and ICU Care: The model calculates the number of individuals needing hospitalization and ICU care based on the percentage of infected individuals who need such care.

Mortality Rate: The model also calculates the number of deaths based on the mortality rate among the infected individuals.

The script uses the deSolve package in R to solve the differential equations that govern the SIR model. The parameters of the model can be optimized to fit observed data, using a least-squares optimization technique. The output of the model is visualized using ggplot2 package.

The data used for this model should be in a CSV format, containing the number of Infected individuals for each time step. The script includes code for reading this data from a CSV file and comparing it with the output of the model.

Please note that this is a simple model meant for illustrative purposes, and does not take into account all the factors that might influence the spread of an actual infectious disease in a real-world population.
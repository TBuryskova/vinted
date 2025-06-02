library(ggplot2)
library(rootSolve)
library(tidyverse)
library(stats)

# Parameters
sigma_pl <- 3
gamma <- 0.6
lambda <- 1
epsilon <- 3



# Expected utilities
E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
E_CZ <- -exp(gamma^2 / 2)
epsilon <- abs(E_PL+E_CZ)/2


# Grid of prior beliefs
p_grid <- seq(0.01, 0.99, length.out = 100)
# ── Helper functions ────────────────────────────────────────────────────────
expected_util <- function(r, gamma, sigma_pl, epsilon) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  epsilon <- abs(E_PL+E_CZ)/2
  
  pmax((1 - r) * E_CZ + r * E_PL, -epsilon)
}


objective <- function(a, p, gamma, sigma_pl, epsilon, lambda) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  epsilon <- abs(E_PL+E_CZ)/2
  

  # Expected utility
  EU <- p * (a * E_PL + (1 - a) * (-epsilon)) + 
    (1 - p) * ((1-a)* E_CZ + (1 - (1-a)) * (-epsilon))
  
  # Joint probabilities
  f_PL1 <- p * a
  f_PL0 <- p * (1 - a)
  f_CZ1 <- (1 - p) * (1-a)
  f_CZ0 <- (1 - p) * (1 - (1-a))
  
  f_y1 <- f_PL1 + f_CZ1
  f_y0 <- f_PL0 + f_CZ0
  
  # Mutual Information
  I <- 0
  if (f_PL1 > 0) I <- I + f_PL1 * log(f_PL1 / (p * f_y1))
  if (f_PL0 > 0) I <- I + f_PL0 * log(f_PL0 / (p * f_y0))
  if (f_CZ1 > 0) I <- I + f_CZ1 * log(f_CZ1 / ((1 - p) * f_y1))
  if (f_CZ0 > 0) I <- I + f_CZ0 * log(f_CZ0 / ((1 - p) * f_y0))
  
  return(EU - lambda * I)
}
H <- function(x) {
  -x*log(x)-(1-x)*log(1-x)
}



solve_ab <- function(p, gamma, epsilon, lambda, sigma_pl) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  
  res <- optimize(
    maximum = TRUE,
    f = objective,
    lower = 0,
    upper = 1,
    p = p,
    gamma = gamma,
    sigma_pl = sigma_pl,
    epsilon = epsilon,
    lambda = lambda
  )
  
  return(res)  # res$maximum is optimal a, res$objective is optimal value
}



V_RI <- vapply(
  p_grid,
  function(p) solve_ab(p = p, gamma=gamma, epsilon = epsilon, lambda = lambda, sigma_pl=sigma_pl)$objective,
  numeric(1)
  )





# ── Value functions ─────────────────────────────
V_no_info   <- expected_util(p_grid, gamma, sigma_pl, epsilon)
V_full_info <- -p_grid * epsilon + (1 - p_grid) * E_CZ
V_perfect   <- pmax(
  -lambda * H(p_grid) + (-p_grid * epsilon + (1 - p_grid) * E_CZ),
  V_no_info
)


# Assemble dataframe
df <- data.frame(
  p = rep(p_grid, 4),
  value = c(V_no_info, V_full_info, V_perfect, V_RI),
  type = rep(c("No info",
               "Full info",
               "Perfect signal",
               "RI optimal"), each = length(p_grid)))

# Plot
ggplot(df, aes(x = p, y = value, color = type), alpha=0.1) +
  geom_line(linewidth = 1) +
  labs(
    title = "Expected Utility under Different Information Structures",
    x = "Prior p (probability of PL)",
    y = "Value",
    color = "Scenario"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

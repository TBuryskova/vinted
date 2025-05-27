library(ggplot2)
library(rootSolve)

# Parameters
sigma_pl <- 3
gamma <- 0.6
epsilon <- 2
lambda <- 0.1

# Expected utilities
E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
E_CZ <- -exp(gamma^2 / 2)
delta <- E_PL - E_CZ

# Grid of prior beliefs
p_grid <- seq(0.01, 0.99, length.out = 200)

# Function to solve for r*(p)
solve_r_star <- function(p) {
  foc <- function(r) {
    if (r <= 0 || r >= 1) return(1e6)
    lhs <- delta * (r - p)
    rhs <- lambda * log((r * (1 - p)) / ((1 - r) * p))
    return(lhs - rhs)
  }
  result <- tryCatch(
    uniroot(foc, lower = 1e-6, upper = 1 - 1e-6)$root,
    error = function(e) NA
  )
  return(result)
}

r_star <- sapply(p_grid, solve_r_star)

# Calculate all utilities
V_red <- p_grid * epsilon + (1 - p_grid) * E_CZ
V_red <- pmax(V_red, -epsilon)

V_blue <- -p_grid * epsilon + (1 - p_grid) * (E_CZ)

H_p <- -p_grid * log(p_grid) - (1 - p_grid) * log(1 - p_grid)
V_perfect <- -lambda * H_p + (-p_grid * epsilon + (1 - p_grid) * E_CZ)
V_perfect <- pmax(V_perfect, V_red)

info_cost <- lambda * (log((r_star^(r_star) / p_grid^p_grid)) / log((1-r_star)^(1-r_star) / (1-p_grid)^(1-p_grid))) 
expected_util <- r_star * E_PL + (1 - r_star) * E_CZ
expected_util <- pmax(expected_util, -epsilon)
V_ri <- expected_util - info_cost

# Assemble dataframe
df <- data.frame(
  p = rep(p_grid, 4),
  value = c(V_red, V_blue, V_perfect, V_ri),
  type = rep(c("No info",
               "Fixed penalty",
               "Perfect signal",
               "RI optimal"), each = length(p_grid))
)

# Plot
ggplot(df, aes(x = p, y = value/10, color = type, alpha=0.5)) +
  geom_line(size = 1) +
  labs(
    title = "Expected Utility under Different Information Structures",
    x = "Prior p (probability of PL)",
    y = "Value",
    color = "Scenario"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

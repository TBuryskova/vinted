library(ggplot2)
library(rootSolve)
library(tidyverse)
# Parameters
sigma_pl <- 6
gamma <- 0.6
lambda <- 0.1

# Expected utilities
E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
E_CZ <- -exp(gamma^2 / 2)
delta <- E_PL - E_CZ
epsilon <- -(E_PL+E_CZ)/2



# Grid of prior beliefs
p_grid <- seq(0.01, 0.99, length.out = 100)

expected_util <- function(r) {
  # stage-2 utility after observing posterior belief r
  r*pmax(r * E_PL,-epsilon) + (1-r)*pmax((1 - r) * E_CZ, -epsilon)
}

info_cost <- function(r, p) {
  if (r <= 0 || r >= 1 || p <= 0 || p >= 1) {
    return(Inf)  # high cost = strongly penalized = avoided
  }
  # per-posterior Shannon cost used earlier
  -lambda * (log((r^(r) / p^p)) -log((1-r)^(1-r) / (1-p)^(1-p))) }

V_RI <- function(r, p) {
expected_util(r) -info_cost(r, p)}

## ---- solve r*(p) by direct maximisation -----------------------------------
solve_r_star <- function(p) {
  # maximise V_RI over r in (0,1); 1e-6 guards against log(0)
  opt <- optimize(
    f        = V_RI,
    interval = c(1e-6, 1 - 1e-6),
    p        = p,
    maximum  = TRUE
  )
  return(opt$maximum)   # r* that maximises V_RI
}

## ---- compute r*(p) and values on a grid -----------------------------------
r_star  <- sapply(p_grid, solve_r_star)

# Value functions
# Value functions
V_no_info   <- expected_util(p_grid)
V_full_info <- -p_grid *  epsilon + (1 - p_grid) * (E_CZ)           # blue
H_p         <- -p_grid * log(p_grid) - (1 - p_grid) * log(1 - p_grid)

V_perfect   <- pmax(-lambda * info_cost(1,p_grid) + (-p_grid * epsilon + (1 - p_grid) * E_CZ), V_no_info) # purple
V_ri        <- mapply(function(r, p) V_RI(r, p), r_star, p_grid)     # green


# Assemble dataframe
df <- data.frame(
  p = rep(p_grid, 4),
  value = c(V_no_info, V_full_info, V_perfect, V_ri),
  type = rep(c("No info",
               "Full info",
               "Perfect signal",
               "RI optimal"), each = length(p_grid))
) %>% mutate(  drop_value = pmax(pmin(value,0),E_PL-2))

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

library(ggplot2)
library(rootSolve)
library(tidyverse)
# Parameters
sigma_pl <- 2
gamma <- 0.6
lambda <- 1
epsilon <- 3


# Expected utilities
E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
E_CZ <- -exp(gamma^2 / 2)
delta <- E_PL - E_CZ



# Grid of prior beliefs
p_grid <- seq(0.01, 0.99, length.out = 100)
# ── Helper functions ────────────────────────────────────────────────────────
expected_util <- function(r, E_PL, E_CZ, epsilon) {
  pmax((1 - r) * E_CZ + r * E_PL, -epsilon)
}

H <- function(x) {
  result <- -x * log(x) - (1 - x) * log(1 - x)
  result[is.nan(result) | is.infinite(result)] <- 0
  result
}

info_cost <- function(r, p, lambda){
  result <-  case_when(lambda * (H(p) - H(r))>0 ~ lambda * (H(p) - H(r)),
                       TRUE ~ NA)
  result
  
  }
V_RI <- function(r, p, E_PL, E_CZ, epsilon, lambda) {
  
  
  expected_util(r, E_PL, E_CZ, epsilon) - info_cost(r, p, lambda)
}

# ── Root-finding for r* ─────────────────────────────────────────────────────
solve_r_star <- function(p, E_PL, E_CZ, epsilon, lambda) {
  optimize(
    f        = V_RI,
    interval = c(1e-6, 1 - 1e-6),
    maximum  = TRUE,
    p        = p,
    E_PL     = E_PL,
    E_CZ     = E_CZ,
    epsilon  = epsilon,
    lambda   = lambda
  )$maximum
}

r_star <- vapply(
  p_grid,
  function(p) solve_r_star(p, E_PL, E_CZ, epsilon, lambda),
  numeric(1)
)

# ── Value functions (now using the new helpers) ─────────────────────────────
V_no_info   <- expected_util(p_grid, E_PL, E_CZ, epsilon)
V_full_info <- -p_grid * epsilon + (1 - p_grid) * E_CZ
V_perfect   <- pmax(
  -lambda * H(p_grid) + (-p_grid * epsilon + (1 - p_grid) * E_CZ),
  V_no_info
)
V_ri <- mapply(
  function(r, p) V_RI(r, p, E_PL, E_CZ, epsilon, lambda),
  r_star, p_grid
)

# Assemble dataframe
df <- data.frame(
  p = rep(p_grid, 4),
  value = c(V_no_info, V_full_info, V_perfect, V_ri),
  r_star = r_star,
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

ggplot(df, aes(x = p, y = r_star), alpha=0.1) +
  geom_line(linewidth = 1) +
  labs(
    title = "R star",
    x = "Prior p (probability of PL)",
    y = "r"
  ) +
  theme_minimal() +
  ylim(0,1) +
  theme(legend.position = "bottom")


# ── Grid over r and p ──
r_vals <- seq(0.01, 0.99, length.out = 100)
p_vals <- seq(0.01, 0.99, length.out = 100)
lambda_val <- 1

# ── Compute values on grid ──
df <- expand.grid(r = r_vals, p = p_vals) %>%
  mutate(cost = info_cost(r, p, lambda_val), U=expected_util(r, E_PL, E_CZ, epsilon))

# ── Plot ──
library(ggplot2)

ggplot(df, aes(x = p, y = r, fill = cost)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    title = paste0("Information Cost: λ = ", lambda_val),
    x = "Prior belief p",
    y = "Posterior belief r",
    fill = "Cost"
  ) +
  theme_minimal()

ggplot(df, aes(x = p, y = r, fill = U)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    title = paste0("U: λ = ", lambda_val),
    x = "Prior belief p",
    y = "Exp. util",
    fill = "Cost"
  ) +
  theme_minimal()


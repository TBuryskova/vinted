# Shiny app: expected-utility curves and optimal posteriors under RI
# Sliders control σ_PL, γ, λ, ε and both graphs recompute instantly.
# -------------------------------------------------------------------------

library(shiny)
library(tidyverse)
library(ggplot2)

# ---------- Helper math ----------------------------------------------------
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







compute_curves <- function(sigma_pl, gamma, lambda, epsilon, n_grid = 100) {
  # Expected utilities
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  
  p_grid <- seq(0.01, 0.99, length.out = n_grid)
  
  
  # Value curves
  V_no_info   <- expected_util(p_grid, epsilon=epsilon, gamma=gamma,  sigma_pl=sigma_pl  )
  V_full_info <- -p_grid * epsilon + (1 - p_grid) * E_CZ
  V_perfect   <- pmax(-lambda * H(p_grid) + (-p_grid * epsilon + (1 - p_grid) * E_CZ), V_no_info)
  V_RI <- vapply(
    p_grid,
    function(p) solve_ab(p = p, gamma=gamma, epsilon = epsilon, lambda = lambda, sigma_pl=sigma_pl)$objective,
    numeric(1)
  )
  a_star <- vapply(
    p_grid,
    function(p) solve_ab(p = p, gamma=gamma, epsilon = epsilon, lambda = lambda, sigma_pl=sigma_pl)$maximum,
    numeric(1)
  )
  
  list(
    df = data.frame(
      p     = rep(p_grid, 4),
      value = c(V_no_info, V_full_info, V_perfect, V_RI),
      type  = rep(c("No info", "Full info", "Perfect signal", "RI optimal"), each = n_grid),
      a_star= rep(a_star, 4)
    ),
    E_PL = E_PL,
    E_CZ = E_CZ
  )
}

# ---------- UI --------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Expected Utility under Different Information Structures"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("sigma_pl", "σ_PL", min = 1, max = 5, value = 3, step = 0.1),
      sliderInput("gamma",     "γ",     min = 0.1, max = 2, value = 0.6, step = 0.1),
      sliderInput("lambda",    "λ",     min = 0,   max = 5, value = 1, step = 0.1),
      uiOutput("epsilon_slider")
    ),
    mainPanel(
      plotOutput("curvePlot", height = "500px"),
      tags$br(),
      plotOutput("aStarPlot", height = "350px")
    )
  )
)

# ---------- Server ----------------------------------------------------------
server <- function(input, output, session) {
  base_curves <- reactive({
    compute_curves(input$sigma_pl, input$gamma, input$lambda, epsilon = 0)
  })
  
  output$epsilon_slider <- renderUI({
    E_vals <- base_curves()
    
    # Defensive check
    if (is.null(E_vals$E_PL) || is.null(E_vals$E_CZ) ||
        !is.finite(E_vals$E_PL) || !is.finite(E_vals$E_CZ)) {
      return(sliderInput("epsilon", "ε", min = 0, max = 1, value = 0.5))
    }
    
    sliderInput("epsilon", "ε",
                min = -E_vals$E_CZ,
                max = -E_vals$E_PL,
                value = (-E_vals$E_CZ - E_vals$E_PL) / 2,
                step = 0.1)
  })
  
  curves <- reactive({
    req(input$epsilon)  # wait until epsilon exists
    compute_curves(input$sigma_pl, input$gamma, input$lambda, input$epsilon)
  })
  
  output$curvePlot <- renderPlot({
    ggplot(curves()$df, aes(x = p, y = value, colour = type)) +
      geom_line(linewidth = 1, alpha = 0.8) +
      labs(
        x = "Prior p (probability PL)",
        y = "Value",
        colour = "Scenario"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  }, res = 96)
  
  output$aStarPlot <- renderPlot({
    p_vals <- seq(0.01, 0.99, length.out = 100)
    a_star_vals <- rep(curves()$df$a_star[1:100], 1)  # Extract only one set
    
    df_astar <- data.frame(p = p_vals, a_star = a_star_vals)
    
    ggplot(df_astar, aes(x = p, y = a_star)) +
      geom_line(linewidth = 1, colour = "steelblue", alpha = 0.9) +
      labs(
        x = "Prior p (probability PL)",
        y = expression(a["*"](p)),
        title = "Optimal a*(p) under Rational Inattention"
      ) +
      ylim(0, 1) +
      theme_minimal()
  }, res = 96)
}

shinyApp(ui, server)

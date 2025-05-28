## app.R
# Shiny app: expected‑utility curves under different information structures
# Sliders control σ_PL, γ, λ, ε and the graph recomputes instantly.
# -------------------------------------------------------------------------

library(shiny)
library(ggplot2)

# ---------- Helper math ----------------------------------------------------
expected_util <- function(r, E_PL, E_CZ, epsilon) {
  pmax((1 - r) * E_CZ + r * E_PL, -epsilon)
}

H <- function(x) {
  -x * log(x) - (1 - x) * log(1 - x)
}

info_cost <- function(r, p, lambda) {
  lambda * (H(p) - H(r))
}

V_RI_scalar <- function(r, p, E_PL, E_CZ, epsilon, lambda) {
  expected_util(r, E_PL, E_CZ, epsilon) - info_cost(r, p, lambda)
}

solve_r_star <- function(p, E_PL, E_CZ, epsilon, lambda) {
  optimize(
    f        = V_RI_scalar,
    interval = c(1e-6, 1 - 1e-6),
    maximum  = TRUE,
    p        = p,
    E_PL     = E_PL,
    E_CZ     = E_CZ,
    epsilon  = epsilon,
    lambda   = lambda
  )$maximum
}

compute_curves <- function(sigma_pl, gamma, lambda, epsilon, n_grid = 100) {
  # Expected utilities
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  
  p_grid <- seq(0.01, 0.99, length.out = n_grid)
  
  # Optimal posterior r*(p)
  r_star <- vapply(
    p_grid,
    solve_r_star,
    numeric(1),
    E_PL = E_PL,
    E_CZ = E_CZ,
    epsilon = epsilon,
    lambda = lambda
  )
  
  # Value curves
  V_no_info   <- expected_util(p_grid, E_PL, E_CZ, epsilon)
  V_full_info <- -p_grid * epsilon + (1 - p_grid) * E_CZ
  V_perfect   <- pmax(-lambda * H(p_grid) + (-p_grid * epsilon + (1 - p_grid) * E_CZ), V_no_info)
  V_ri        <- mapply(
    V_RI_scalar,
    r = r_star,
    p = p_grid,
    MoreArgs = list(E_PL = E_PL, E_CZ = E_CZ, epsilon = epsilon, lambda = lambda)
  )
  
  list(
    df = data.frame(
      p     = rep(p_grid, 4),
      value = c(V_no_info, V_full_info, V_perfect, V_ri),
      type  = rep(c("No info", "Full info", "Perfect signal", "RI optimal"), each = n_grid)
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
      plotOutput("curvePlot", height = "550px")
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
}
shinyApp(ui, server)

            
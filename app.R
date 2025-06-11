
library(shiny)
library(tidyverse)
library(ggplot2)

# ---------- Helper math ----------------------------------------------------
p_grid <- seq(0.02, 0.98, length.out = 100)
# ── Helper functions ────────────────────────────────────────────────
expected_util <- function(r, gamma, sigma_pl, epsilon) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  epsilon <- abs(E_PL+E_CZ)/2
  
  pmax((1 - r) * E_CZ + r * E_PL, -epsilon)
}


objective <- function(params, p, gamma, sigma_pl, lambda) {
  a <- params[1]  # P(y = 1 | PL)
  b <- params[2]  # P(y = 1 | CZ)
  
  # Expected utilities in each state
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  epsilon <- abs(E_PL + E_CZ) / 2
  
  # Expected utility
  EU <- p * (a * E_PL + (1 - a) * (-epsilon)) + 
    (1 - p) * (b * E_CZ + (1 - b) * (-epsilon))
  
  # Joint probabilities
  f_PL1 <- p * a
  f_PL0 <- p * (1 - a)
  f_CZ1 <- (1 - p) * b
  f_CZ0 <- (1 - p) * (1 - b)
  
  # Marginal probabilities
  f_y1 <- f_PL1 + f_CZ1
  f_y0 <- f_PL0 + f_CZ0
  
  # Mutual information
  I <- 0
  if (f_PL1 > 0) I <- I + f_PL1 * log(f_PL1 / (p * f_y1))
  if (f_PL0 > 0) I <- I + f_PL0 * log(f_PL0 / (p * f_y0))
  if (f_CZ1 > 0) I <- I + f_CZ1 * log(f_CZ1 / ((1 - p) * f_y1))
  if (f_CZ0 > 0) I <- I + f_CZ0 * log(f_CZ0 / ((1 - p) * f_y0))
  
  # Return minus net utility
  return(-EU+lambda * I)
}

H <- function(x) {
  -x*log(x)-(1-x)*log(1-x)
}

solve_ab <- function(p, gamma, epsilon, lambda, sigma_pl) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  
  res <- optim(
    par = c(0.5, 0.5),  
    fn = objective,
    p = p,
    gamma = gamma,
    method= "L-BFGS-B",
    sigma_pl = sigma_pl,
    lambda = lambda,
    lower = c(0, 0),
    upper = c(1, 1)
  )
  
  return(res) 
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
  
  V_RI <- mapply(
    function(p) {-(solve_ab(p = p, gamma=gamma, epsilon = epsilon, lambda = lambda, sigma_pl=sigma_pl)$value)},
    p_grid
  )
  
  a_star_RI <- vapply(
    p_grid,
    function(p) solve_ab(p = p, gamma=gamma, epsilon = epsilon, lambda = lambda, sigma_pl=sigma_pl)$par[1],
    numeric(1)
  )
  
  b_star_RI<- vapply(
    p_grid,
    function(p) (solve_ab(p = p, gamma=gamma, epsilon = epsilon, lambda = lambda, sigma_pl=sigma_pl)$par[2]),
    numeric(1)
  )
  
  buy_RI <- mapply(
    function(p, a, b) {p*a+(1-p)*b},
    p_grid, a_star_RI, b_star_RI
  )
  
  
  buy_no_info <- vapply(
    p_grid,
    function(p) (p*E_PL+(1-p)*E_CZ>-epsilon),
    numeric(1)
    
  )
  
  buy_full_info<- vapply(
    p_grid,
    function(p) 1-p,
    numeric(1)
  )
  
  
  
  buy_perfect <-  mapply(
    function(p, v_no) {
      signal_aq <- as.numeric(-lambda * H(p) + (-p * epsilon + (1 - p) * E_CZ) >= v_no)
      signal_aq * (1-p) + (1-signal_aq)*(p*E_PL+(1-p)*E_CZ>-epsilon)

    },
    p_grid, V_no_info
  )
  
  info_no_info <- vapply(
    p_grid,
    function(p) 0,
    numeric(1)
  )
  
  
  info_full_info <- vapply(
    p_grid,
    function(p) H(p),
    numeric(1)
  )
  
  info_perfect <-mapply(
    function(p, v_no) {
      signal_aq <- as.numeric(-lambda * H(p) + (-p * epsilon + (1 - p) * E_CZ) >= v_no)
      signal_aq*H(p) 
    },
    p_grid, V_no_info
  )
  
  
  info_RI <- mapply(
    function(p, a,b) {
      # Joint probabilities
      f_PL1 <- p * a
      f_PL0 <- p * (1 - a)
      f_CZ1 <- (1 - p) * b
      f_CZ0 <- (1 - p) * (1 - b)
      
      f_y1 <- f_PL1 + f_CZ1
      f_y0 <- f_PL0 + f_CZ0
      
      # Mutual Information
      I <- 0
      if (f_PL1 > 0) I <- I + f_PL1 * log(f_PL1 / (p * f_y1))
      if (f_PL0 > 0) I <- I + f_PL0 * log(f_PL0 / (p * f_y0))
      if (f_CZ1 > 0) I <- I + f_CZ1 * log(f_CZ1 / ((1 - p) * f_y1))
      if (f_CZ0 > 0) I <- I + f_CZ0 * log(f_CZ0 / ((1 - p) * f_y0))
      return(I)
    },
    p_grid, a_star_RI, b_star_RI
  )  
  
  
  
  # Behavioral discrimination measure
  behavioral_discrimination_RI <- mapply(
    function(a,b) {b-a},
    a_star_RI, b_star_RI
  )  
  
  behavioral_discrimination_no_info <- vapply(
    p_grid,
    function(p) 0,
    numeric(1)
  )
  
  
  behavioral_discrimination_perfect <- mapply(
    function(p, v_no) {
      signal_aq <- as.numeric(-lambda * H(p) + (-p * epsilon + (1 - p) * E_CZ) >= v_no)
      signal_aq * 1
    },
    p_grid, V_no_info
  )
  
  behavioral_discrimination_full_info <- vapply(
    p_grid,
    function(p) 1,
    numeric(1)
  )
  
  
  
  list(
    df = data.frame(
      p     = rep(p_grid, 4),
      value = c(V_no_info, V_full_info, V_perfect, V_RI),
      type  = rep(c("No info", "Full info", "Perfect signal", "RI optimal"), each = 100),
      behavioral_discrimination = c(behavioral_discrimination_no_info,behavioral_discrimination_full_info,behavioral_discrimination_perfect,behavioral_discrimination_RI),
      info = c(info_no_info,info_full_info,info_perfect,info_RI),
      buy= c(buy_no_info,buy_full_info,buy_perfect,buy_RI)
    ),
    E_PL = E_PL,
    E_CZ = E_CZ,
    star=data.frame(
      p=p_grid,
      a=a_star_RI,
      b=b_star_RI
    )
  )
}

# ---------- UI --------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Outcomes under Different Information Structures"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("sigma_pl", "σ_PL", min = 1, max = 5, value = 3, step = 0.1),
      sliderInput("gamma",     "γ",     min = 0.1, max = 2, value = 0.6, step = 0.1),
      sliderInput("lambda",    "λ",     min = 0,   max = 5, value = 1, step = 0.1),
      uiOutput("epsilon_slider")
      
    ),
    mainPanel(
      fluidRow(
        column(6, plotOutput("curvePlot", height = "400px")),
        column(6, plotOutput("buyPlot", height = "400px"))
      ),
      fluidRow(
        column(6, plotOutput("discrimPlot", height = "400px")),
        column(6, plotOutput("info_acquired", height = "400px"))
      ),
    fluidRow(
      column(6, plotOutput("astarplot", height = "400px")),
      column(6, plotOutput("bstarplot", height = "400px"))
    )
    )
  )
)


# ---------- Server ----------------------------------------------------------
server <- function(input, output, session) {
  base_curves <- reactive({
    compute_curves(input$sigma_pl, input$gamma, input$lambda, epsilon = 1)
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
    req(input$epsilon)
    compute_curves(input$sigma_pl, input$gamma, input$lambda, input$epsilon)
  })
  
  output$curvePlot <- renderPlot({
    ggplot(curves()$df, aes(x = p, y = value, colour = type)) +
      geom_line(linewidth = 1, alpha = 0.8) +
      labs(
        x = "Prior p (probability PL)",
        y = "Value",
        colour = "Scenario",
        title = "Expected payoff"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$buyPlot <- renderPlot({
    ggplot(curves()$df, aes(x = p, y = buy, colour = type)) +
      geom_line(linewidth = 1, alpha = 0.8) +
      labs(
        x = "Prior p (probability PL)",
        y = "Value",
        colour = "Scenario",
        title = "Probability of buying"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  
  output$discrimPlot <- renderPlot({
    ggplot(curves()$df, aes(x = p, y = behavioral_discrimination, colour = type)) +
      geom_line(linewidth = 1, alpha = 0.9) +
      labs(
        x = "Prior p (probability PL)",
        y = "Behavioral Discrimination",
        colour = "Scenario",
        title = "Behavioral Discrimination"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$info_acquired <- renderPlot({
    ggplot(curves()$df, aes(x = p, y = info, colour = type)) +
      geom_line(linewidth = 1, alpha = 0.9) +
      labs(
        x = "Prior p (probability PL)",
        y = "Mutual entropy",
        colour = "Scenario",
        title = "Amount of info acquired"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  
  
  output$astarplot <- renderPlot({
    ggplot(curves()$star, aes(x = p, y = a)) +
      geom_line(linewidth = 1, alpha = 0.9) +
      labs(
        x = "Prior p (probability PL)",
        y = "a star",
        title = "Optimal a=P(buy|PL)"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  
  
  output$bstarplot <- renderPlot({
    ggplot(curves()$star, aes(x = p, y = b)) +
      geom_line(linewidth = 1, alpha = 0.9) +
      labs(
        x = "Prior p (probability PL)",
        y = "a star",
        title = "Optimal b=P(buy|CZ)"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
}

shinyApp(ui, server)

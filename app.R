library(shiny)
library(tidyverse)
library(ggplot2)

# ---------- Helper functions ------------------------------------------------
H <- function(x) {
  -x * log(x) - (1 - x) * log(1 - x)
}

expected_util <- function(r, gamma, sigma_pl, epsilon) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  pmax((1 - r) * E_CZ + r * E_PL, -epsilon)
}

objective <- function(params, p, gamma, sigma_pl, lambda, epsilon) {
  a <- params[1]
  b <- params[2]
  
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  
  EU <- p * (a * E_PL + (1 - a) * (-epsilon)) + 
    (1 - p) * (b * E_CZ + (1 - b) * (-epsilon))
  
  f_PL1 <- p * a
  f_PL0 <- p * (1 - a)
  f_CZ1 <- (1 - p) * b
  f_CZ0 <- (1 - p) * (1 - b)
  
  f_y1 <- f_PL1 + f_CZ1
  f_y0 <- f_PL0 + f_CZ0
  
  I <- 0
  if (f_PL1 > 0) I <- I + f_PL1 * log(f_PL1 / (p * f_y1))
  if (f_PL0 > 0) I <- I + f_PL0 * log(f_PL0 / (p * f_y0))
  if (f_CZ1 > 0) I <- I + f_CZ1 * log(f_CZ1 / ((1 - p) * f_y1))
  if (f_CZ0 > 0) I <- I + f_CZ0 * log(f_CZ0 / ((1 - p) * f_y0))
  
  return(-EU + lambda * I)
}

solve_ab <- function(p, gamma, sigma_pl, lambda, epsilon) {
  res <- optim(
    par = c(0.5, 0.5),
    fn = objective,
    p = p,
    gamma = gamma,
    sigma_pl = sigma_pl,
    lambda = lambda,
    epsilon = epsilon,
    method = "L-BFGS-B",
    lower = c(0, 0),
    upper = c(1, 1)
  )
  return(res)
}

compute_curves <- function(sigma_pl, gamma, lambda, epsilon, n_grid = 1000) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  p_grid <- seq(0.001, 0.999, length.out = n_grid)

  
  V_no_info <- expected_util(p_grid, gamma = gamma, sigma_pl = sigma_pl, epsilon = epsilon)
  V_full_info <- -p_grid * epsilon + (1 - p_grid) * E_CZ
  V_perfect <- pmax(-lambda * H(p_grid) + V_full_info, V_no_info)
  
  V_RI <- mapply(
    function(p) {-(solve_ab(p, gamma, sigma_pl, lambda, epsilon)$value)},
    p_grid
  )
  
  a_star_RI <- vapply(
    p_grid, function(p) solve_ab(p, gamma, sigma_pl, lambda, epsilon)$par[1], numeric(1)
  )
  
  b_star_RI <- vapply(
    p_grid, function(p) solve_ab(p, gamma, sigma_pl, lambda, epsilon)$par[2], numeric(1)
  )
  
  buy_RI <- p_grid * a_star_RI + (1 - p_grid) * b_star_RI
  buy_no_info <- (p_grid * E_PL + (1 - p_grid) * E_CZ > -epsilon)
  buy_full_info <- 1 - p_grid
  buy_perfect <- mapply(
    function(p, v_no) {
      signal_aq <- as.numeric(-lambda * H(p) + (-p * epsilon + (1 - p) * E_CZ) >= v_no)
      signal_aq * (1 - p) + (1 - signal_aq) * (p * E_PL + (1 - p) * E_CZ > -epsilon)
    },
    p_grid, V_no_info
  )
  
  info_RI <- mapply(
    function(p, a, b) {
      f_PL1 <- p * a
      f_PL0 <- p * (1 - a)
      f_CZ1 <- (1 - p) * b
      f_CZ0 <- (1 - p) * (1 - b)
      f_y1 <- f_PL1 + f_CZ1
      f_y0 <- f_PL0 + f_CZ0
      I <- 0
      if (f_PL1 > 0) I <- I + f_PL1 * log(f_PL1 / (p * f_y1))
      if (f_PL0 > 0) I <- I + f_PL0 * log(f_PL0 / (p * f_y0))
      if (f_CZ1 > 0) I <- I + f_CZ1 * log(f_CZ1 / ((1 - p) * f_y1))
      if (f_CZ0 > 0) I <- I + f_CZ0 * log(f_CZ0 / ((1 - p) * f_y0))
      return(I)
    },
    p_grid, a_star_RI, b_star_RI
  )
  
  info_full_info <- H(p_grid)
  info_no_info <- rep(0, length(p_grid))
  info_perfect <- mapply(
    function(p, v_no) {
      signal_aq <- as.numeric(-lambda * H(p) + (-p * epsilon + (1 - p) * E_CZ) >= v_no)
      signal_aq * H(p)
    },
    p_grid, V_no_info
  )
  
  behavioral_discrimination_RI <- b_star_RI - a_star_RI
  behavioral_discrimination_full_info <- rep(1, length(p_grid))
  behavioral_discrimination_no_info <- rep(0, length(p_grid))
  behavioral_discrimination_perfect <- mapply(
    function(p, v_no) {
      signal_aq <- as.numeric(-lambda * H(p) + (-p * epsilon + (1 - p) * E_CZ) >= v_no)
      signal_aq * 1
    },
    p_grid, V_no_info
  )
  
  list(
    df = data.frame(
      p = rep(p_grid, 4),
      value = c(V_no_info, V_full_info, V_perfect, V_RI),
      type = rep(c("No info", "Full info", "Perfect signal", "RI optimal"), each = n_grid),
      behavioral_discrimination = c(behavioral_discrimination_no_info,
                                    behavioral_discrimination_full_info,
                                    behavioral_discrimination_perfect,
                                    behavioral_discrimination_RI),
      info = c(info_no_info, info_full_info, info_perfect, info_RI),
      buy = c(buy_no_info, buy_full_info, buy_perfect, buy_RI)
    ),
    E_PL = E_PL,
    E_CZ = E_CZ,
    star = data.frame(p = p_grid, a = a_star_RI, b = b_star_RI)
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
        column(6, plotOutput("curvePlot")),
        column(6, plotOutput("buyPlot"))
      ),
      fluidRow(
        column(6, plotOutput("discrimPlot")),
        column(6, plotOutput("info_acquired"))
      ),
      fluidRow(
        column(6, plotOutput("astarplot")),
        column(6, plotOutput("bstarplot"))
      )
    )
  )
)

# ---------- Server ----------------------------------------------------------
server <- function(input, output, session) {
  params <- reactive({
    gamma <- input$gamma
    sigma_pl <- input$sigma_pl
    E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
    E_CZ <- -exp(gamma^2 / 2)
    epsilon_default <- abs(E_PL + E_CZ) / 2
    list(gamma = gamma, sigma_pl = sigma_pl, E_PL = E_PL, E_CZ = E_CZ, epsilon_default = epsilon_default)
  })
  
  output$epsilon_slider <- renderUI({
    p <- params()
    sliderInput("epsilon", "ε",
                min = round(-p$E_CZ, 1),
                max = round(-p$E_PL, 1),
                value = round(p$epsilon_default, 1),
                step = 0.1)
  })
  
  curves <- reactive({
    p <- params()
    req(input$epsilon)
    compute_curves(
      sigma_pl = p$sigma_pl,
      gamma = p$gamma,
      lambda = input$lambda,
      epsilon = input$epsilon
    )
  })
  
  output$curvePlot <- renderPlot({
    ggplot(curves()$df, aes(x = p, y = value, colour = type)) +
      geom_line(linewidth = 1, alpha = 0.8) +
      labs(x = "Prior p (probability PL)", y = "Value", colour = "Scenario", title = "Expected payoff") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$buyPlot <- renderPlot({
    ggplot(curves()$df, aes(x = p, y = buy, colour = type)) +
      geom_line(linewidth = 1, alpha = 0.8) +
      labs(x = "Prior p (probability PL)", y = "Probability of Buying", colour = "Scenario", title="Probability of buying") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$discrimPlot <- renderPlot({
    ggplot(curves()$df, aes(x = p, y = behavioral_discrimination, colour = type)) +
      geom_line(linewidth = 1) +
      labs(x = "Prior p", y = "Discrimination", title = "Behavioral Discrimination", colour = "Scenario") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$info_acquired <- renderPlot({
    ggplot(curves()$df, aes(x = p, y = info, colour = type)) +
      geom_line(linewidth = 1) +
      labs(x = "Prior p", y = "Mutual Information", title = "Information Acquired", colour = "Scenario") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$astarplot <- renderPlot({
    ggplot(curves()$star, aes(x = p, y = a)) +
      geom_line(linewidth = 1) +
      labs(x = "Prior p", y = "a*", title = "Optimal a = P(buy | PL)") +
      theme_minimal()
  })
  
  output$bstarplot <- renderPlot({
    ggplot(curves()$star, aes(x = p, y = b)) +
      geom_line(linewidth = 1) +
      labs(x = "Prior p", y = "b*", title = "Optimal b = P(buy | CZ)") +
      theme_minimal()
  })
}

shinyApp(ui, server)

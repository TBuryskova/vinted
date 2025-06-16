library(shiny)
library(tidyverse)
library(ggplot2)

# ---------- Helper functions ------------------------------------------------
safe_log <- function(x) log(pmax(x, 1e-10))

H <- function(x) {
  -x * safe_log(x) - (1 - x) * safe_log(1 - x)
}

expected_util <- function(p, gamma, sigma_pl, epsilon) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  pmax((1 - p) * E_CZ + p * E_PL, -epsilon)
}

objective <- function(params, p, gamma, sigma_pl, lambda, epsilon) {
  a <- params[1]; b <- params[2]
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  EU <- p * (a * E_PL + (1 - a) * (-epsilon)) +
    (1 - p) * (b * E_CZ + (1 - b) * (-epsilon))
  
  f_PL1 <- p * a; f_PL0 <- p * (1 - a)
  f_CZ1 <- (1 - p) * b; f_CZ0 <- (1 - p) * (1 - b)
  f_y1 <- f_PL1 + f_CZ1; f_y0 <- f_PL0 + f_CZ0
  
  I <- 0
  if (f_PL1 > 0) I <- I + f_PL1 * safe_log(f_PL1 / (p * f_y1))
  if (f_PL0 > 0) I <- I + f_PL0 * safe_log(f_PL0 / (p * f_y0))
  if (f_CZ1 > 0) I <- I + f_CZ1 * safe_log(f_CZ1 / ((1 - p) * f_y1))
  if (f_CZ0 > 0) I <- I + f_CZ0 * safe_log(f_CZ0 / ((1 - p) * f_y0))
  
  -EU + lambda * I
}

solve_ab <- function(p, gamma, sigma_pl, lambda, epsilon) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  
  if (p < 1e-3 || p > 1 - 1e-3) {
    a_star <- as.numeric(E_PL > -epsilon)
    b_star <- as.numeric(E_CZ > -epsilon)
    EU <- p * (a_star * E_PL + (1 - a_star) * (-epsilon)) +
      (1 - p) * (b_star * E_CZ + (1 - b_star) * (-epsilon))
    return(list(value = -EU, par = c(a_star, b_star)))
  }
  
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

compute_curves_by_lambda <- function(p, gamma, sigma_pl, epsilon, lambda_grid) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  
  V_no_info <- rep(expected_util(p, gamma, sigma_pl, epsilon), length(lambda_grid))
  V_full_info <- rep(-p * epsilon + (1 - p) * E_CZ, length(lambda_grid))
  H_p <- H(p)
  V_perfect <- pmax(-lambda_grid * H_p + rep(V_full_info[1], length(lambda_grid)),
                    rep(V_no_info[1], length(lambda_grid)))  
  solve_res <- lapply(lambda_grid, function(lambda)
    solve_ab(p, gamma, sigma_pl, lambda, epsilon)
  )
  a_star <- sapply(solve_res, function(res) res$par[1])
  b_star <- sapply(solve_res, function(res) res$par[2])
  V_RI <- sapply(solve_res, function(res) -res$value)
  
  buy_no_info <- rep(as.numeric(p * E_PL + (1 - p) * E_CZ > -epsilon), length(lambda_grid))
  buy_full_info <- rep(1 - p, length(lambda_grid))
  buy_perfect <- mapply(function(lambda, v_no) {
    signal <- as.numeric(-lambda * H(p) + (-p * epsilon + (1 - p) * E_CZ) >= v_no)
    signal * (1 - p) + (1 - signal) * as.numeric(p * E_PL + (1 - p) * E_CZ > -epsilon)
  }, lambda_grid, V_no_info)
  
  buy_RI <- p * a_star + (1 - p) * b_star
  
  info_perfect <- sapply(lambda_grid, function(lambda) {
    signal <- as.numeric(-lambda * H(p) + (-p * epsilon + (1 - p) * E_CZ) >= expected_util(p, gamma, sigma_pl, epsilon))
    signal * H(p)
  })
  
  info_RI <- mapply(function(a, b) {
    f_PL1 <- p * a; f_PL0 <- p * (1 - a)
    f_CZ1 <- (1 - p) * b; f_CZ0 <- (1 - p) * (1 - b)
    f_y1 <- f_PL1 + f_CZ1; f_y0 <- f_PL0 + f_CZ0
    I <- 0
    if (f_PL1 > 0) I <- I + f_PL1 * safe_log(f_PL1 / (p * f_y1))
    if (f_PL0 > 0) I <- I + f_PL0 * safe_log(f_PL0 / (p * f_y0))
    if (f_CZ1 > 0) I <- I + f_CZ1 * safe_log(f_CZ1 / ((1 - p) * f_y1))
    if (f_CZ0 > 0) I <- I + f_CZ0 * safe_log(f_CZ0 / ((1 - p) * f_y0))
    return(I)
  }, a_star, b_star)
  
  data.frame(
    lambda = rep(lambda_grid, 4),
    value = c(V_no_info, V_full_info, V_perfect, V_RI),
    type = rep(c("No info", "Full info", "Perfect signal", "RI optimal"), each = length(lambda_grid)),
    info = c(rep(0, length(lambda_grid)),
             rep(H(p), length(lambda_grid)),
             info_perfect,
             info_RI),
    buy = c(buy_no_info, buy_full_info, buy_perfect, buy_RI),
    discrim = c(rep(0, length(lambda_grid)),
                rep(1, length(lambda_grid)),
                rep(1, length(lambda_grid)),
                b_star - a_star)
  )
  
}

# ---------- UI --------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Outcomes vs. λ under Fixed Prior"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("p", "Prior p (probability PL)", min = 0.01, max = 0.99, value = 0.5, step = 0.01),
      sliderInput("sigma_pl", "σ_PL", min = 1, max = 5, value = 3, step = 0.1),
      sliderInput("gamma", "γ", min = 0.1, max = 2, value = 0.6, step = 0.1),
      sliderInput("epsilon", "ε", min = 0.1, max = 5, value = 1, step = 0.1)
    ),
    mainPanel(
      plotOutput("utilityPlot"),
      plotOutput("buyPlot"),
      plotOutput("discrimPlot"),
      plotOutput("infoPlot")
    )
  )
)

# ---------- Server ----------------------------------------------------------
server <- function(input, output, session) {
  lambda_grid <- seq(0, 5, length.out = 200)
  
  curves <- reactive({
    compute_curves_by_lambda(
      p = input$p,
      gamma = input$gamma,
      sigma_pl = input$sigma_pl,
      epsilon = input$epsilon,
      lambda_grid = lambda_grid
    )
  })
  
  output$utilityPlot <- renderPlot({
    ggplot(curves(), aes(x = lambda, y = value, color = type)) +
      geom_line(linewidth = 1) +
      labs(title = "Expected Utility vs. λ", x = "λ (info cost)", y = "Expected Utility") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$buyPlot <- renderPlot({
    ggplot(curves(), aes(x = lambda, y = buy, color = type)) +
      geom_line(linewidth = 1) +
      labs(title = "Probability of Buying vs. λ", x = "λ (info cost)", y = "P(buy)") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$discrimPlot <- renderPlot({
    ggplot(curves(), aes(x = lambda, y = discrim, color = type)) +
      geom_line(linewidth = 1) +
      labs(title = "Behavioral Discrimination vs. λ", x = "λ (info cost)", y = "P(buy|CZ) - P(buy|PL)") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$infoPlot <- renderPlot({
    ggplot(curves(), aes(x = lambda, y = info, color = type)) +
      geom_line(linewidth = 1) +
      labs(title = "Information Acquired vs. λ", x = "λ (info cost)", y = "Mutual Information") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
}

shinyApp(ui, server)

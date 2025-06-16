# ========== App 3: Varying sigma_pl (for fixed p, gamma, lambda, epsilon) ==========

library(shiny)
library(tidyverse)
library(ggplot2)

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

compute_curves_by_sigma <- function(p, gamma, lambda, epsilon, sigma_grid) {
  H_p <- H(p)
  
  results <- lapply(sigma_grid, function(sigma_pl) {
    E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
    E_CZ <- -exp(gamma^2 / 2)
    V_no <- expected_util(p, gamma, sigma_pl, epsilon)
    V_full <- -p * epsilon + (1 - p) * E_CZ
    V_perf <- max(-lambda * H_p + V_full, V_no)
    
    res <- solve_ab(p, gamma, sigma_pl, lambda, epsilon)
    a <- res$par[1]; b <- res$par[2]; V_ri <- -res$value
    buy_ri <- p * a + (1 - p) * b
    discrim <- b - a
    
    info <- {
      f_PL1 <- p * a; f_PL0 <- p * (1 - a)
      f_CZ1 <- (1 - p) * b; f_CZ0 <- (1 - p) * (1 - b)
      f_y1 <- f_PL1 + f_CZ1; f_y0 <- f_PL0 + f_CZ0
      I <- 0
      if (f_PL1 > 0) I <- I + f_PL1 * safe_log(f_PL1 / (p * f_y1))
      if (f_PL0 > 0) I <- I + f_PL0 * safe_log(f_PL0 / (p * f_y0))
      if (f_CZ1 > 0) I <- I + f_CZ1 * safe_log(f_CZ1 / ((1 - p) * f_y1))
      if (f_CZ0 > 0) I <- I + f_CZ0 * safe_log(f_CZ0 / ((1 - p) * f_y0))
      I
    }
    
    tibble(
      sigma_pl = sigma_pl,
      type = c("No info", "Full info", "Perfect signal", "RI optimal"),
      value = c(V_no, V_full, V_perf, V_ri),
      buy = c(as.numeric(p * E_PL + (1 - p) * E_CZ > -epsilon),
              1 - p,
              ifelse(-lambda * H_p + V_full >= V_no, 1 - p, as.numeric(p * E_PL + (1 - p) * E_CZ > -epsilon)),
              buy_ri),
      discrim = c(0, 1, 1, discrim),
      info = c(0, H_p,
               ifelse(-lambda * H_p + V_full >= V_no, H_p, 0),
               info)
    )
  })
  
  bind_rows(results)
}

# ========== Shiny App UI & Server ==========

ui <- fluidPage(
  titlePanel("Outcomes vs. σ_PL (PL Quality Variance)"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("p", "Prior p (prob PL)", 0.01, 0.99, 0.5, 0.01),
      sliderInput("gamma", "γ (risk aversion)", 0.1, 2, 0.6, 0.1),
      sliderInput("lambda", "λ (info cost)", 0, 5, 1, 0.1),
      sliderInput("epsilon", "ε (outside option)", 0.1, 5, 1, 0.1)
    ),
    mainPanel(
      plotOutput("utilPlot"),
      plotOutput("buyPlot"),
      plotOutput("discrimPlot"),
      plotOutput("infoPlot")
    )
  )
)

server <- function(input, output, session) {
  sigma_grid <- seq(1, 5, length.out = 100)
  
  curves <- reactive({
    compute_curves_by_sigma(
      p = input$p,
      gamma = input$gamma,
      lambda = input$lambda,
      epsilon = input$epsilon,
      sigma_grid = sigma_grid
    )
  })
  
  output$utilPlot <- renderPlot({
    ggplot(curves(), aes(x = sigma_pl, y = value, color = type)) +
      geom_line(linewidth = 1) +
      labs(title = "Expected Utility vs σ_PL", x = "σ_PL", y = "Expected Utility") +
      theme_minimal() + theme(legend.position = "bottom")
  })
  
  output$buyPlot <- renderPlot({
    ggplot(curves(), aes(x = sigma_pl, y = buy, color = type)) +
      geom_line(linewidth = 1) +
      labs(title = "Buying Probability vs σ_PL", x = "σ_PL", y = "P(buy)") +
      theme_minimal() + theme(legend.position = "bottom")
  })
  
  output$discrimPlot <- renderPlot({
    ggplot(curves(), aes(x = sigma_pl, y = discrim, color = type)) +
      geom_line(linewidth = 1) +
      labs(title = "Discrimination vs σ_PL", x = "σ_PL", y = "P(buy|CZ) - P(buy|PL)") +
      theme_minimal() + theme(legend.position = "bottom")
  })
  
  output$infoPlot <- renderPlot({
    ggplot(curves(), aes(x = sigma_pl, y = info, color = type)) +
      geom_line(linewidth = 1) +
      labs(title = "Information Acquired vs σ_PL", x = "σ_PL", y = "Mutual Information") +
      theme_minimal() + theme(legend.position = "bottom")
  })
}

shinyApp(ui, server)

library(shiny)
library(tidyverse)
library(ggplot2)

# ---------- Helper functions ------------------------------------------------
H <- function(x) {
  -x * log(x) - (1 - x) * log(1 - x)
}
safe_log <- function(x) log(pmax(x, 1e-10))

expected_util <- function(r, E_PL, E_CZ, epsilon) {
  pmax((1 - r) * E_CZ + r * E_PL, epsilon)
}

objective <- function(params, p, E_PL, E_CZ, lambda, epsilon) {
  a <- params[1]
  b <- params[2]
  
  EU <- p * (a * E_PL + (1 - a) * (epsilon)) +
    (1 - p) * (b * E_CZ + (1 - b) * (epsilon))
  
  f_PL1 <- p * a
  f_PL0 <- p * (1 - a)
  f_CZ1 <- (1 - p) * b
  f_CZ0 <- (1 - p) * (1 - b)
  
  f_y1 <- f_PL1 + f_CZ1
  f_y0 <- f_PL0 + f_CZ0
  
  I <- 0
  if (f_PL1 > 0) I <- I + f_PL1 * safe_log(f_PL1 / (p * f_y1))
  if (f_PL0 > 0) I <- I + f_PL0 * safe_log(f_PL0 / (p * f_y0))
  if (f_CZ1 > 0) I <- I + f_CZ1 * safe_log(f_CZ1 / ((1 - p) * f_y1))
  if (f_CZ0 > 0) I <- I + f_CZ0 * safe_log(f_CZ0 / ((1 - p) * f_y0))
  
  return(-EU + lambda * I)
}

solve_ab <- function(p, E_PL, E_CZ, lambda, epsilon) {
  res <- optim(
    par = c(0.5, 0.5),
    fn = objective,
    p = p,
    E_PL = E_PL,
    E_CZ = E_CZ,
    lambda = lambda,
    epsilon = epsilon,
    method = "L-BFGS-B",
    lower = c(0, 0),
    upper = c(1, 1)
  )
  return(res)
}

compute_curves <- function(E_PL, E_CZ, lambda, p, n_grid = 1000) {
  eps_grid <- seq(0.001, 0.999, length.out = n_grid)
  
  V_no_info <- (eps_grid<p*E_PL+(1-p)*E_CZ)*(p*E_PL+(1-p)*E_CZ)+(eps_grid>p*E_PL+(1-p)*E_CZ)*eps_grid
  V_full_info <-  (eps_grid<E_PL)*(p*E_PL+(1-p)*E_CZ)+(eps_grid<E_CZ &eps_grid>E_PL)*(p*eps_grid+(1-p)*E_CZ)+(eps_grid>E_CZ)*eps_grid
  V_perfect <- pmax(-lambda * H(p) + V_full_info, V_no_info)
  
  solve_res <- lapply(eps_grid, function(eps) solve_ab(p, E_PL, E_CZ, lambda, eps))
  V_RI <- sapply(solve_res, function(res) -res$value)
  a_star_RI <- sapply(solve_res, function(res) res$par[1])
  b_star_RI <- sapply(solve_res, function(res) res$par[2])
  
  buy_RI <- p * a_star_RI + (1 - p) * b_star_RI
  buy_no_info <- 1*(eps_grid<p*E_PL+(1-p)*E_CZ)+0     
  buy_full_info <- (eps_grid<E_PL)*1+(eps_grid<E_CZ &eps_grid>E_PL)*( (1-p))
  buy_perfect <- as.numeric(-lambda * H(p) + (p * eps_grid + (1 - p) * E_CZ) >= 
                              p * E_PL + (1 - p) * E_CZ) * (1 - p) +
    as.numeric(-lambda * H(p) + (p * eps_grid + (1 - p) * E_CZ) < 
                 p * E_PL + (1 - p) * E_CZ) * (p * E_PL + (1 - p) * E_CZ > eps_grid)
  
  info_RI <- mapply(
    function(a, b, eps) {
      f_PL1 <- p * a
      f_PL0 <- p * (1 - a)
      f_CZ1 <- (1 - p) * b
      f_CZ0 <- (1 - p) * (1 - b)
      f_y1 <- f_PL1 + f_CZ1
      f_y0 <- f_PL0 + f_CZ0
      I <- 0
      if (f_PL1 > 0) I <- I + f_PL1 * safe_log(f_PL1 / (p * f_y1))
      if (f_PL0 > 0) I <- I + f_PL0 * safe_log(f_PL0 / (p * f_y0))
      if (f_CZ1 > 0) I <- I + f_CZ1 * safe_log(f_CZ1 / ((1 - p) * f_y1))
      if (f_CZ0 > 0) I <- I + f_CZ0 * safe_log(f_CZ0 / ((1 - p) * f_y0))
      return(I)
    }, a_star_RI, b_star_RI, eps_grid
  )
  info_full  <- H(p)*(eps_grid<E_CZ &eps_grid>E_PL)
  info_perfect <- as.numeric(-lambda * H(p) + V_full_info >= V_no_info) * H(p)
  behavioral_discrimination_RI <- b_star_RI - a_star_RI
  behavioral_discrimination_perfect <- as.numeric(-lambda * H(p) + V_full_info >= V_no_info)
  
  list(
    df = data.frame(
      epsilon = rep(eps_grid, 4),
      value = c(V_no_info, V_full_info, V_perfect, V_RI),
      type = rep(c("No info", "Full info", "Perfect signal", "RI optimal"), each = n_grid),
      behavioral_discrimination = c(rep(0, n_grid), rep(1, n_grid), behavioral_discrimination_perfect, behavioral_discrimination_RI),
      info = c(rep(0, n_grid), info_full, info_perfect, info_RI),
      buy = c(buy_no_info, buy_full_info, buy_perfect, buy_RI)
    ),
    star = data.frame(epsilon = eps_grid, a = a_star_RI, b = b_star_RI)
  )
}

# ---------- UI --------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Outcomes vs Epsilon under Different Info Structures"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("lambda", "λ", min = 0, max = 5, value = 1, step = 0.1),
      sliderInput("E_CZ", "E_CZ", min = 0, max = 1, value = 0.6, step = 0.01),
      sliderInput("E_PL", "E_PL", min = 0, max = 1, value = 0.3, step = 0.01),
      sliderInput("p", "p", min = 0, max = 1, value = 0.2, step = 0.01)
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
  curves <- reactive({
    compute_curves(E_PL = input$E_PL, E_CZ = input$E_CZ, lambda = input$lambda, p = input$p)
  })
  
  output$curvePlot <- renderPlot({
    ggplot(curves()$df, aes(x = epsilon, y = value, colour = type)) +
      geom_line(linewidth = 1) +
      labs(x = expression(epsilon), y = "Value", title = "Expected payoff") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$buyPlot <- renderPlot({
    ggplot(curves()$df, aes(x = epsilon, y = buy, colour = type)) +
      geom_line(linewidth = 1) +
      labs(x = expression(epsilon), y = "Probability of Buying", title = "Probability of buying") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$discrimPlot <- renderPlot({
    ggplot(curves()$df, aes(x = epsilon, y = behavioral_discrimination, colour = type)) +
      geom_line(linewidth = 1) +
      labs(x = expression(epsilon), y = "Discrimination", title = "Behavioral Discrimination") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$info_acquired <- renderPlot({
    ggplot(curves()$df, aes(x = epsilon, y = info, colour = type)) +
      geom_line(linewidth = 1) +
      labs(x = expression(epsilon), y = "Mutual Information", title = "Information Acquired") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$astarplot <- renderPlot({
    ggplot(curves()$star, aes(x = epsilon, y = a)) +
      geom_line(linewidth = 1) +
      labs(x = expression(epsilon), y = "a*", title = "Optimal a = P(buy | PL)") +
      theme_minimal()
  })
  
  output$bstarplot <- renderPlot({
    ggplot(curves()$star, aes(x = epsilon, y = b)) +
      geom_line(linewidth = 1) +
      labs(x = expression(epsilon), y = "b*", title = "Optimal b = P(buy | CZ)") +
      theme_minimal()
  })
}

shinyApp(ui, server)

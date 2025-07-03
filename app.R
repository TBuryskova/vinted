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
  V_perfect <- pmax(-lambda * rep(H(p), n_grid) + V_full_info, V_no_info)
  
  solve_res <- lapply(eps_grid, function(eps) solve_ab(p, E_PL, E_CZ, lambda, eps))
  V_RI <- sapply(solve_res, function(res) -res$value)
  a_star_RI <- sapply(solve_res, function(res) res$par[1])
  b_star_RI <- sapply(solve_res, function(res) res$par[2])
  
  buy_RI <- p * a_star_RI + (1 - p) * b_star_RI
  buy_no_info <- 1*(eps_grid<p*E_PL+(1-p)*E_CZ)+0     
  buy_full_info <- (eps_grid<E_PL)*1+(eps_grid<E_CZ &eps_grid>E_PL)*( (1-p))
  buy_perfect <- (V_full_info-lambda*rep(H(p), n_grid)>V_no_info)*buy_full_info+(V_full_info-lambda*rep(H(p), n_grid)<V_no_info)*buy_no_info
  
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
  info_full  <- rep(H(p),n_grid)
  info_perfect <- (V_full_info-lambda*rep(H(p),n_grid)>V_no_info)*rep(H(p),n_grid)+(V_full_info-lambda*rep(H(p),n_grid)<V_no_info)*0
  
  
  behavioral_discrimination_RI <- b_star_RI - a_star_RI
  behavioral_discrimination_perfect <- (V_full_info-lambda*rep(H(p),n_grid)>V_no_info)*1+(V_full_info-lambda*rep(H(p),n_grid)<V_no_info)*0
  
  
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

area_under_curve <- function(x, y) {
  ord <- order(x)
 sum(diff(x[ord]) * (head(y[ord], -1) + tail(y[ord], -1)) / 2)
}
# ---------- UI --------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Outcomes vs Epsilon under Different Info Structures"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("lambda", "λ", min = 0, max = 1, value = 0.2, step = 0.05),
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
  areas_buy <- reactive({
    df <- curves()$df
    df %>%
      group_by(type) %>%
      summarize(area = area_under_curve(epsilon, buy), .groups = "drop")
  })
  
  areas_discrim <- reactive({
    df <- curves()$df
    df %>%
      group_by(type) %>%
      summarize(area = area_under_curve(epsilon, behavioral_discrimination), .groups = "drop")
  })
  
  areas_util <- reactive({
    df <- curves()$df
    df %>%
      group_by(type) %>%
      summarize(area = area_under_curve(epsilon, value), .groups = "drop")
  })
  
  areas_info <- reactive({
    df <- curves()$df
    df %>%
      group_by(type) %>%
      summarize(area = area_under_curve(epsilon, info), .groups = "drop")
  })
  
  curves <- reactive({
    compute_curves(E_PL = input$E_PL, E_CZ = input$E_CZ, lambda = input$lambda, p = input$p)
  })
  
  output$curvePlot <- renderPlot({
    df <- curves()$df
    areas <- areas_util()
    ggplot(curves()$df, aes(x = epsilon, y = value, colour = type)) +
      geom_line(linewidth = 1) +
      geom_text(
        data = areas,
        aes(x = 0.1, y = 0.95 - 0.04 * as.numeric(factor(type)), 
            label = paste0("∫ = ", round(area,4)), colour = type),
       show.legend = FALSE
      ) +
      labs(x = expression(epsilon), y = "Value", title = "Expected payoff") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$buyPlot <- renderPlot({
    df <- curves()$df
    areas <- areas_buy()
    
    ggplot(df, aes(x = epsilon, y = buy, colour = type)) +
      geom_line(linewidth = 1) +
      labs(x = expression(epsilon), y = "Probability of Buying", title = "Probability of buying") +
      geom_text(
        data = areas,
        aes(x = 0.05, y = 0.55 - 0.05 * as.numeric(factor(type)), 
            label = paste0("∫ = ", round(area, 4)), colour = type),
        hjust = 0, show.legend = FALSE
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$discrimPlot <- renderPlot({
    df <- curves()$df
    areas <- areas_discrim()
    
    ggplot(df, aes(x = epsilon, y = behavioral_discrimination, colour = type)) +
      geom_line(linewidth = 1) +
      labs(x = expression(epsilon), y = "Discrimination", title = "Behavioral Discrimination") +
      geom_text(
        data = areas,
        aes(x = 0.95, y = 0.95 - 0.07 * as.numeric(factor(type)), 
            label = paste0("∫ = ", round(area, 4)), colour = type),
        hjust = 1, show.legend = FALSE
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  
  output$info_acquired <- renderPlot({
    df <- curves()$df
    areas <- areas_info()
    
    ggplot(curves()$df, aes(x = epsilon, y = info, colour = type)) +
      geom_line(linewidth = 1) +
      labs(x = expression(epsilon), y = "Mutual Information", title = "Information Acquired") +
      theme_minimal() +
      theme(legend.position = "bottom") +
      geom_text(
        data = areas,
        aes(x = 0.4, y = 0.3- 0.04 * as.numeric(factor(type)), 
            label = paste0("∫ = ", round(area, 4)), colour = type),
        hjust = 1, show.legend = FALSE
      ) 
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

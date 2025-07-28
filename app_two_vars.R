library(shiny)
library(ggplot2)
library(dplyr)
library(scales)
library(metR) # Although metR might not be strictly needed for line plots, keep for consistency if other features are added back.

ui <- fluidPage(
  titlePanel("Firm-User Interaction Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("p", "p (probability of L):", min = 0.01, max = 0.99, value = 0.2, step = 0.1),
      sliderInput("epsilon", "ε", min = 0, max = 5, value = 1.5, step = 0.5),
      sliderInput("U_L", "U_L:", min = 0, max = 4, value = 1, step = 1),
      sliderInput("U_H", "U_H:", min = 0, max = 5, value = 2, step = 1),
      
      # Slider pro pevnou hodnotu w
      sliderInput("w_value", "w (fixed value):", min = 0.01, max = 5, value = 0.5, step = 0.01),
      
      actionButton("plot_all", "Generate Plots")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("User Utility (vs lambda)", plotOutput("plot_utility_surface")),
        tabPanel("Firm Profit (vs lambda)", plotOutput("plot_profit_surface")),
        tabPanel("Adoption Rate (vs lambda)", plotOutput("plot_adoption_rate_surface"))
      )
    )
  )
)

server <- function(input, output) {
  
  solve_ab <- function(lambda, p, epsilon, U_L, U_H, w, tol = 1e-8, max_iter = 1000) {
    U_L_bar <- U_L - (1 + w)
    U_H_bar <- U_H - (1 + w)
    U_0_bar <- epsilon
    
    e_L <- exp(U_L_bar / lambda)
    e_H <- exp(U_H_bar / lambda)
    e_0 <- exp(U_0_bar / lambda)
    
    a <- 0.5
    b <- 0.5
    
    for (i in 1:max_iter) {
      a_old <- a
      b_old <- b
      
      P <- p * a + (1 - p) * b
      Q <- 1 - P
      
      if (is.na(P) || is.na(Q) || P <= 0 || Q <= 0) {
        return(list(a = NA, b = NA, converged = FALSE))
      }
      
      a_new <- (e_L + e_0 * P / Q)^(-1) * e_L
      b_new <- (e_H + e_0 * P / Q)^(-1) * e_H
      
      a_new <- min(max(a_new, 0), 1)
      b_new <- min(max(b_new, 0), 1)
      
      if (abs(a_new - a_old) < tol && abs(b_new - b_old) < tol) {
        return(list(a = a_new, b = b_new, converged = TRUE))
      }
      
      a <- a_new
      b <- b_new
    }
    return(list(a = a, b = b, converged = FALSE))
  }
  
  compute_expected_utility <- function(lambda, p, epsilon, U_L, U_H, w, a, b) {
    U_L_bar <- U_L - (1 + w)
    U_H_bar <- U_H - (1 + w)
    
    EU <- p * (a * U_L_bar + (1 - a) * epsilon) +
      (1 - p) * (b * U_H_bar + (1 - b) * epsilon)
    
    info_cost <- compute_mutual_info(lambda, p, a, b)
    
    return(EU - info_cost)
  }
  
  compute_mutual_info <- function(lambda, p, a, b) {
    eps <- 1e-10
    a <- min(max(a, eps), 1 - eps)
    b <- min(max(b, eps), 1 - eps)
    
    P1 <- p * a + (1 - p) * b
    P0 <- 1 - P1
    
    if (P1 <= 0 || P0 <= 0) return(NA)
    
    mi <- - lambda * (
      p * a * log(a / P1) +
        p * (1 - a) * log((1 - a) / P0) +
        (1 - p) * b * log(b / P1) +
        (1 - p) * (1 - b) * log((1 - b) / P0)
    )
    
    return(mi)
  }
  
  grid_data_reactive <- eventReactive(input$plot_all, {
    fixed_lambda_range <- c(0.01, 100) 
    
    lambda_vals <- seq(fixed_lambda_range[1], fixed_lambda_range[2], length.out = 100) # Hustší mřížka
    # w_vals je nyní pevná hodnota ze slideru
    w_fixed <- input$w_value             
    
    grid <- expand.grid(lambda = lambda_vals, w = w_fixed) # w je pevné
    
    withProgress(message = 'Computing grid...', value = 0, {
      result_grid <- grid %>%
        rowwise() %>%
        mutate(
          res = list(tryCatch(
            solve_ab(lambda, input$p, input$epsilon, input$U_L, input$U_H, w),
            error = function(e) list(a = NA, b = NA, converged = FALSE)
          )),
          a = res$a,
          b = res$b
        ) %>%
        mutate(
          adoption_rate = ifelse(is.na(a) | is.na(b), NA, input$p * a + (1 - input$p) * b),
          info = ifelse(is.na(a) | is.na(b), NA, compute_mutual_info(lambda, input$p, a, b)),
          utility = ifelse(is.na(a) | is.na(b), NA, compute_expected_utility(lambda, input$p, input$epsilon, input$U_L, input$U_H, w, a, b)),
          profit = ifelse(is.na(adoption_rate), NA, adoption_rate * w) 
        ) %>%
        ungroup() %>%
        filter(!is.na(utility) & !is.na(profit) & !is.na(adoption_rate)) 
      
      incProgress(1/1)
      return(result_grid)
    })
  })
  
  # Odstraněna funkce gradients_data_reactive, protože už nepotřebujeme 2D derivace pro 1D ploty
  
  # Odstraněn output$equilibriumPlot
  
  output$plot_utility_surface <- renderPlot({
    grid_for_plot <- grid_data_reactive()
    req(grid_for_plot)
    
    ggplot(grid_for_plot, aes(x = lambda, y = utility)) + # Y-axis is now utility value
      geom_line(color = "blue", linewidth = 1) + # Changed to line plot
      labs(title = paste0("User Utility (w = ", input$w_value, ")"), 
           x = expression(lambda), y = "Utility") +
      theme_minimal()
  })
  
  output$plot_profit_surface <- renderPlot({
    grid_for_plot <- grid_data_reactive()
    req(grid_for_plot)
    
    ggplot(grid_for_plot, aes(x = lambda, y = profit)) + # Y-axis is now profit value
      geom_line(color = "darkgreen", linewidth = 1) + # Changed to line plot
      labs(title = paste0("Firm Profit (w = ", input$w_value, ")"), 
           x = expression(lambda), y = "Profit") +
      theme_minimal()
  })
  
  output$plot_adoption_rate_surface <- renderPlot({
    grid_for_plot <- grid_data_reactive()
    req(grid_for_plot)
    
    ggplot(grid_for_plot, aes(x = lambda, y = adoption_rate)) + # Y-axis is now adoption_rate value
      geom_line(color = "purple", linewidth = 1) + # Changed to line plot
      labs(title = paste0("User Adoption Rate (w = ", input$w_value, ")"), 
           x = expression(lambda), y = "Adoption Rate") +
      theme_minimal()
  })
}

shinyApp(ui = ui, server = server)
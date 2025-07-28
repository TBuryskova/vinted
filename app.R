library(shiny)
library(ggplot2)
library(reshape2)

ui <- fluidPage(
  titlePanel("Fixed-Point Solver for a and b"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("p", "p (probability of L):", value = 0.7, min = 0.01, max = 0.99),
      sliderInput("epsilon", "ε:", value = 1.5, min = 0, max = 5),
      sliderInput("U_L", "U_L:", value = 1, min = 0, max = 5),
      sliderInput("U_H", "U_H:", value = 2, min = 0, max = 5),
      sliderInput("w", "w :", value = 0.2, min = 0.01, max = 0.99),
      actionButton("plot", "Plot All")
    ),
    
    mainPanel(
      plotOutput("heatmap"),
      plotOutput("utilityPlot"),
      plotOutput("platformPlot") , # Third plot
      plotOutput("infoPlot")  # Third plot
      
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
      
      if (is.na(P) || is.na(Q) || P <= 0 || Q <= 0) break
      
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
    
    mi <- - lambda * (
      p * a * log(a / P1) +
        p * (1 - a) * log((1 - a) / P0) +
        (1 - p) * b * log(b / P1) +
        (1 - p) * (1 - b) * log((1 - b) / P0)
    )
    
    return(mi)
  }
  
  
  
  observeEvent(input$plot, {
    lambda_vals <- seq(0, 1, length.out = 100)

    grid <- expand.grid(lambda = lambda_vals)
    
    P_buy <- numeric(nrow(grid))
    Utility <- numeric(nrow(grid))
    Info <- numeric(nrow(grid))
    
    for (i in 1:nrow(grid)) {
      lambda <- grid$lambda[i]
      w <- input$w
      
      res <- tryCatch(
        solve_ab(lambda, input$p, input$epsilon, input$U_L, input$U_H, w),
        error = function(e) list(a = NA, b = NA)
      )
      
      a <- res$a
      b <- res$b
      
      if (!is.na(a) && !is.na(b)) {
        P_buy[i] <- input$p * a + (1 - input$p) * b
        Info[i] <- compute_mutual_info(lambda, input$p, a, b)
        Utility[i] <- compute_expected_utility(lambda, input$p, input$epsilon, input$U_L, input$U_H, w, a, b)
      } else {
        P_buy[i] <- NA
        Info[i] <- NA
        Utility[i] <- NA
      }
    }
    
    grid$result <- P_buy
    grid$utility <- Utility
    grid$info <- Info


    output$heatmap <- renderPlot({
      ggplot(grid, aes(x = lambda,  y = result)) +
        geom_line() +
        labs(
          title = paste("P(buy)"),
          x = expression(lambda),
        ) +
        theme_minimal()
    })
    
    
    output$infoPlot <- renderPlot({
      ggplot(grid, aes(x = lambda,  y = -info) )+
        geom_line() +
        labs(
          title = paste("Info"),
          x = expression(lambda),
        ) +
        theme_minimal()
    })
    
    output$utilityPlot <- renderPlot({
      ggplot(grid, aes(x = lambda,  y = utility)) +
        geom_line() +
        labs(
          title = "Expected Utility minus Info Cost",
          x = expression(lambda),
        ) +
        theme_minimal()
    })
    
    output$platformPlot <- renderPlot({
      ggplot(grid, aes(x = lambda, y = result*w-info)) +
        geom_line() +
        labs(
          title = "Platform profit",
          x = expression(lambda),
        ) +
        theme_minimal()
    })
  })
}

shinyApp(ui = ui, server = server)

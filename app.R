library(shiny)
library(ggplot2)
library(dplyr)

ui <- fluidPage(
  titlePanel("Fixed-Point Solver for a and b"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("p", "p (probability of L):", value = 0.5, min = 0.01, max = 0.99),
      sliderInput("epsilon", "ε:", value = 5, min = 0, max = 10),
      sliderInput("U_L", "U_L:", value = 0, min = 0, max = 5),
      sliderInput("U_H", "U_H:", value =10, min = 0, max = 10)
      # All buttons removed as requested
    ),
    
    mainPanel(
      # Summary Table is now defined here, outside the tabsetPanel, but within mainPanel
      h3("Summary of Optimal Points"),
      tableOutput("summaryTable"), 
      hr(), # Horizontal rule for separation
      
      # Tabbed interface for plots
      tabsetPanel(
        id = "mainTabs", # Give an ID to the tabsetPanel for potential future use
        
        tabPanel("Firm First", # Tab 1 for Firm First analysis
                 h4("Firm First Analysis: Optimal λ and w for given a,b"),
                 plotOutput("utilityPlot"), 
                 plotOutput("platformPlot") 
        ),
        
        tabPanel("User First", # Tab 2 for User First analysis
                 h4("User First Analysis: Firm's Optimal Response to a and b"),
                 plotOutput("profitPlot2"), 
                 plotOutput("utilityPlot2")  
        )
      ) # End of tabsetPanel
    )
  )
)

server <- function(input, output) {
  # --- Core Functions (unchanged) ---
  welfare_objective <- function(x, p, epsilon, U_L, U_H) {
    lambda <- x[1]
    w <- x[2]
    a <- x[3]
    b <- x[4]
    P_buy <- a * p + (1 - p) * b
    info <- compute_mutual_info(p, a, b)
    utility <- compute_expected_utility(lambda, p, epsilon, U_L, U_H, w, a, b)
    profit <- P_buy * w + info 
    welfare <- -(utility + profit) 
    return(welfare)
  }
  
  solve_ab <- function(lambda, p, epsilon, U_L, U_H, w, tol = 1e-8, max_iter = 1000) {
    U_L_bar <- U_L - w
    U_H_bar <- U_H - w
    U_0_bar <- epsilon
    if (lambda < 1e-6) { 
      return(list(a = NaN, b = NaN, converged = FALSE))
    }
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
      denom_a <- e_L + e_0 * P / Q
      denom_b <- e_H + e_0 * P / Q
      if (denom_a <= 0 || denom_b <= 0) break
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
    U_L_bar <- U_L - w
    U_H_bar <- U_H - w
    EU <- p * (a * U_L_bar + (1 - a) * epsilon) +
      (1 - p) * (b * U_H_bar + (1 - b) * epsilon)
    mi <- compute_mutual_info(p, a, b)
    return(EU - lambda * mi)
  }
  
  compute_mutual_info <- function(p, a, b) {
    eps <- 1e-10 
    a <- min(max(a, eps), 1 - eps)
    b <- min(max(b, eps), 1 - eps)
    P1 <- p * a + (1 - p) * b
    P0 <- 1 - P1
    if (P1 <= 0 || P0 <= 0) {
      return(0)
    }
    mi <- (
      p * a * log(a / P1) +
        p * (1 - a) * log((1 - a) / P0) +
        (1 - p) * b * log(b / P1) +
        (1 - p) * (1 - b) * log((1 - b) / P0)
    )
    return(mi)
  }
  
  firm_profit_objective <- function(x, p, epsilon, U_L, U_H, a, b) {
    w <- x[1]
    lambda <- x[2]
    if (lambda < 0.01) lambda <- 0.01 
    P_buy_val <- a * p + (1 - p) * b
    info_val <- compute_mutual_info(p, a, b)
    profit <- w * P_buy_val + lambda * info_val
    return(-profit) 
  }
  
  # --- Function to find Pareto optimal points (generalized) ---
  find_pareto_optimal <- function(data, utility_col_name, profit_col_name, x_col_name, y_col_name) {
    # Select relevant columns and rename them for generic processing
    temp_data <- data %>%
      select(x_val = !!sym(x_col_name), y_val = !!sym(y_col_name),
             utility = !!sym(utility_col_name), profit = !!sym(profit_col_name)) %>%
      filter(!is.na(utility) & !is.na(profit))
    
    if (nrow(temp_data) == 0) {
      # Return an empty data frame with the correct column names for plotting
      return(data.frame(x_val = numeric(0), y_val = numeric(0), utility = numeric(0), profit = numeric(0)))
    }
    
    # Sort data by utility ascending, then by profit descending for tie-breaking
    sorted_data_for_pareto <- temp_data %>%
      arrange(utility, desc(profit))
    
    pareto_indices <- c()
    current_max_profit <- -Inf
    
    # Iterate from the point with the highest utility (end of sorted_data_for_pareto)
    for (i in nrow(sorted_data_for_pareto):1) {
      if (sorted_data_for_pareto$profit[i] >= current_max_profit) {
        pareto_indices <- c(pareto_indices, i)
        current_max_profit <- sorted_data_for_pareto$profit[i]
      }
    }
    
    # Select the Pareto optimal points and re-sort them by utility ascending
    pareto_optimal_points <- sorted_data_for_pareto[pareto_indices, ] %>%
      arrange(utility) 
    
    # Rename columns back to original names for plotting
    pareto_optimal_points_final <- pareto_optimal_points %>%
      rename(!!x_col_name := x_val, !!y_col_name := y_val) %>%
      select(!!sym(x_col_name), !!sym(y_col_name), utility, profit) 
    
    return(pareto_optimal_points_final)
  }
  
  # --- Reactive expressions for computations ---
  
  # Reactive for Social Planner Optimization
  reactive_planner_opt <- reactive({
    req(input$p, input$epsilon, input$U_L, input$U_H) 
    
    opt_result <- optim(
      par = c(1, 2, 0.5, 0.5), 
      fn = welfare_objective,
      method = "L-BFGS-B",
      lower = c(0.02, 0.01, 0.01, 0.01),
      upper = c(10, input$U_H - input$epsilon, 0.99, 0.99),
      p = input$p,
      epsilon = input$epsilon,
      U_L = input$U_L,
      U_H = input$U_H
    )
    
    opt_lambda <- opt_result$par[1]
    opt_w <- opt_result$par[2]
    opt_a <- opt_result$par[3]
    opt_b <- opt_result$par[4]
    
    opt_P_buy <- input$p * opt_a + (1 - input$p) * opt_b
    opt_info <- compute_mutual_info(input$p, opt_a, opt_b)
    opt_utility <- compute_expected_utility(opt_lambda, input$p, input$epsilon, input$U_L, input$U_H, opt_w, opt_a, opt_b)
    opt_profit <- opt_P_buy * opt_w + opt_info 
    opt_welfare <- opt_utility + opt_profit 
    
    data.frame(
      Criterion = "Social Planner Optimum",
      Lambda = round(opt_lambda, 3),
      w = round(opt_w, 3),
      a = round(opt_a, 3), 
      b = round(opt_b, 3), 
      Utility = round(opt_utility, 3),
      Profit = round(opt_profit, 3),
      Welfare = round(opt_welfare, 3)
    )
  })
  
  # Reactive for Grid1 data (Firm First analysis)
  reactive_grid1_data <- reactive({
    req(input$p, input$epsilon, input$U_L, input$U_H)
    
    lambda_vals <- seq(0.02, 10, length.out = 100)
    w_vals <- seq(0.01, input$U_H - input$epsilon, length.out = 100)
    grid <- expand.grid(lambda = lambda_vals, w = w_vals)
    
    P_buy <- numeric(nrow(grid))
    Utility <- numeric(nrow(grid))
    Info <- numeric(nrow(grid))
    Profit <- numeric(nrow(grid))
    Welfare <- numeric(nrow(grid))
    a_vals <- numeric(nrow(grid)) 
    b_vals <- numeric(nrow(grid)) 
    
    for (i in 1:nrow(grid)) {
      lambda <- grid$lambda[i]
      w <- grid$w[i]
      
      res <- tryCatch(
        solve_ab(lambda, input$p, input$epsilon, input$U_L, input$U_H, w),
        error = function(e) list(a = NA, b = NA) 
      )
      
      a <- res$a
      b <- res$b
      
      if (!is.na(a) && !is.na(b)) {
        P_buy[i] <- input$p * a + (1 - input$p) * b
        Info[i] <- compute_mutual_info(input$p, a, b)
        Utility[i] <- compute_expected_utility(lambda, input$p, input$epsilon, input$U_L, input$U_H, w, a, b)
        Profit[i] <- P_buy[i] * w + Info[i] 
        Welfare[i] <- Utility[i] + Profit[i]
        a_vals[i] <- a 
        b_vals[i] <- b 
      } else {
        P_buy[i] <- NA
        Info[i] <- NA
        Utility[i] <- NA
        Profit[i] <- NA
        Welfare[i] <- NA
        a_vals[i] <- NA 
        b_vals[i] <- NA 
      }
    }
    
    grid$result <- P_buy
    grid$utility <- Utility
    grid$info <- Info
    grid$profit <- Profit
    grid$welfare <- Welfare
    grid$a <- a_vals 
    grid$b <- b_vals 
    
    return(grid)
  })
  
  # Reactive for Pareto Optimal points from Grid1
  reactive_pareto_points_grid1 <- reactive({
    req(reactive_grid1_data()) 
    find_pareto_optimal(reactive_grid1_data(), 
                        utility_col_name = "utility", 
                        profit_col_name = "profit", 
                        x_col_name = "lambda", 
                        y_col_name = "w")
  })
  
  # Reactive for Grid2 data (User First analysis)
  reactive_grid2_data <- reactive({
    req(input$p, input$epsilon, input$U_L, input$U_H)
    
    a_vals_grid2 <- seq(0.01, 0.99, length.out = 50) 
    b_vals_grid2 <- seq(0.01, 0.99, length.out = 50)
    grid2 <- expand.grid(a = a_vals_grid2, b = b_vals_grid2)
    
    max_firm_profit <- numeric(nrow(grid2))
    user_utility_at_firm_opt <- numeric(nrow(grid2))
    optimal_w_firm_grid2 <- numeric(nrow(grid2))
    optimal_lambda_firm_grid2 <- numeric(nrow(grid2))
    
    for (i in 1:nrow(grid2)) {
      current_a <- grid2$a[i]
      current_b <- grid2$b[i]
      
      firm_opt_result <- tryCatch(
        optim(
          par = c( (input$U_H - input$epsilon) / 2, 1), 
          fn = firm_profit_objective,
          method = "L-BFGS-B",
          lower = c(0.01, 0.02), 
          upper = c(input$U_H - input$epsilon, 10), 
          p = input$p,
          epsilon = input$epsilon,
          U_L = input$U_L,
          U_H = input$U_H,
          a = current_a,
          b = current_b
        ),
        error = function(e) {
          warning(paste("Optimization failed for a=", current_a, "b=", current_b, ":", e$message))
          list(value = NA, par = c(NA, NA)) 
        }
      )
      
      optimal_w_firm_grid2[i] <- firm_opt_result$par[1]
      optimal_lambda_firm_grid2[i] <- firm_opt_result$par[2] 
      
      max_firm_profit[i] <- -firm_opt_result$value
      
      user_utility_at_firm_opt[i] <- compute_expected_utility(
        optimal_lambda_firm_grid2[i], 
        input$p, 
        input$epsilon, 
        input$U_L, 
        input$U_H, 
        optimal_w_firm_grid2[i], 
        current_a, 
        current_b
      )
    }
    
    grid2$max_firm_profit <- max_firm_profit
    grid2$user_utility_at_firm_opt <- user_utility_at_firm_opt
    grid2$optimal_w_firm <- optimal_w_firm_grid2
    grid2$optimal_lambda_firm <- optimal_lambda_firm_grid2
    
    return(grid2)
  })
  
  # Reactive for Pareto Optimal points from Grid2
  reactive_pareto_points_grid2 <- reactive({
    req(reactive_grid2_data()) 
    find_pareto_optimal(reactive_grid2_data(), 
                        utility_col_name = "user_utility_at_firm_opt", 
                        profit_col_name = "max_firm_profit", 
                        x_col_name = "a", 
                        y_col_name = "b")
  })
  
  # Reactive for Nash Equilibrium calculation
  reactive_nash_equilibrium_data <- reactive({
    req(input$p, input$epsilon, input$U_L, input$U_H) 
    
    w_curr <- (input$U_H - input$epsilon) / 2
    lambda_curr <- 1.0
    
    initial_ab_res <- solve_ab(lambda_curr, input$p, input$epsilon, input$U_L, input$U_H, w_curr)
    a_curr <- initial_ab_res$a
    b_curr <- initial_ab_res$b
    
    if (is.na(a_curr) || is.na(b_curr)) {
      a_curr <- 0.5
      b_curr <- 0.5
      warning("Initial solve_ab for Nash failed, using default a=0.5, b=0.5.")
    }
    
    max_iter_nash <- 200
    tol_nash <- 1e-6
    
    converged <- FALSE
    
    for (i in 1:max_iter_nash) {
      w_prev <- w_curr
      lambda_prev <- lambda_curr
      a_prev <- a_curr
      b_prev <- b_curr
      
      res_ab <- tryCatch(
        solve_ab(lambda_curr, input$p, input$epsilon, input$U_L, input$U_H, w_curr),
        error = function(e) {
          warning(paste("solve_ab failed in Nash iteration", i, ":", e$message))
          list(a = NA, b = NA)
        }
      )
      
      a_next <- res_ab$a
      b_next <- res_ab$b
      
      if (is.na(a_next) || is.na(b_next)) {
        warning(paste("a or b became NA in Nash iteration", i, ". Using previous values."))
        a_next <- a_curr 
        b_next <- b_curr
      }
      
      firm_opt_res <- tryCatch(
        optim(
          par = c(w_curr, lambda_curr), 
          fn = firm_profit_objective,
          method = "L-BFGS-B",
          lower = c(0.01, 0.02), 
          upper = c(input$U_H - input$epsilon, 10), 
          p = input$p,
          epsilon = input$epsilon,
          U_L = input$U_L,
          U_H = input$U_H,
          a = a_next, 
          b = b_next
        ),
        error = function(e) {
          warning(paste("optim failed in Nash iteration", i, ":", e$message))
          list(value = NA, par = c(NA, NA))
        }
      )
      
      w_next <- firm_opt_res$par[1]
      lambda_next <- firm_opt_res$par[2]
      
      if (is.na(w_next) || is.na(lambda_next)) {
        warning(paste("w or lambda became NA in Nash iteration", i, ". Using previous values."))
        w_next <- w_curr 
        lambda_next <- lambda_curr
      }
      
      w_curr <- w_next
      lambda_curr <- lambda_next
      a_curr <- a_next
      b_curr <- b_next
      
      diff_w <- abs(w_curr - w_prev)
      diff_lambda <- abs(lambda_curr - lambda_prev)
      diff_a <- abs(a_curr - a_prev)
      diff_b <- abs(b_curr - b_prev)
      
      if (max(diff_w, diff_lambda, diff_a, diff_b, na.rm = TRUE) < tol_nash) {
        converged <- TRUE
        message(paste("Nash Equilibrium converged in", i, "iterations."))
        break
      }
    }
    
    if (!converged) {
      warning("Nash Equilibrium did not converge within max iterations.")
    }
    
    nash_P_buy <- input$p * a_curr + (1 - input$p) * b_curr
    nash_info <- compute_mutual_info(input$p, a_curr, b_curr)
    nash_utility <- compute_expected_utility(lambda_curr, input$p, input$epsilon, input$U_L, input$U_H, w_curr, a_curr, b_curr)
    nash_profit <- nash_P_buy * w_curr + nash_info
    nash_welfare <- nash_utility + nash_profit
    
    data.frame(
      Criterion = "Nash Equilibrium",
      Lambda = round(lambda_curr, 3),
      w = round(w_curr, 3),
      a = round(a_curr, 3), 
      b = round(b_curr, 3), 
      Utility = round(nash_utility, 3),
      Profit = round(nash_profit, 3),
      Welfare = round(nash_welfare, 3)
    )
  }) # End of reactive_nash_equilibrium_data
  
  # --- Plot Outputs (now depend directly on reactive data) ---
  
  output$utilityPlot <- renderPlot({
    req(reactive_grid1_data(), reactive_pareto_points_grid1())
    grid_data <- reactive_grid1_data()
    pareto_points <- reactive_pareto_points_grid1()
    max_point <- grid_data[which.max(grid_data$utility), ] 
    
    ggplot(grid_data, aes(x = lambda, y = w)) +
      geom_tile(aes(fill = utility)) +
      geom_contour(aes(z = utility), color = "white") +
      geom_point(data = max_point, aes(x = lambda, y = w), color = "red", size = 3) + 
      # Add Pareto optimal points
      geom_point(data = pareto_points, aes(x = lambda, y = w), 
                 shape = 4, color = "black", size = 4, stroke = 1.5) + # Black crosses
      labs(title = "Expected Utility - Info Cost", x = expression(lambda), y = "w") +
      scale_fill_viridis_c() +
      theme_minimal()
  },
  res = 96 
  )
  
  output$platformPlot <- renderPlot({
    req(reactive_grid1_data(), reactive_pareto_points_grid1())
    grid_data <- reactive_grid1_data()
    pareto_points <- reactive_pareto_points_grid1()
    max_point <- grid_data[which.max(grid_data$profit), ] 
    
    ggplot(grid_data, aes(x = lambda, y = w )) +
      geom_tile(aes(fill = profit)) +
      geom_contour(aes(z = profit), color = "white") +
      geom_point(data = max_point, aes(x = lambda, y = w), color = "red", size = 3) + 
      # Add Pareto optimal points
      geom_point(data = pareto_points, aes(x = lambda, y = w), 
                 shape = 4, color = "black", size = 4, stroke = 1.5) + # Black crosses
      labs(title = "Platform Profit", x = expression(lambda), y = "w") +
      scale_fill_viridis_c() +
      theme_minimal()
  },
  res = 96
  )
  
  output$profitPlot2 <- renderPlot({
    req(reactive_grid2_data(), reactive_pareto_points_grid2()) 
    grid2_data <- reactive_grid2_data()
    pareto_points_grid2 <- reactive_pareto_points_grid2() 
    max_point <- grid2_data[which.max(grid2_data$max_firm_profit), ] 
    
    ggplot(grid2_data, aes(x = a, y = b)) +
      geom_tile(aes(fill = max_firm_profit)) +
      geom_contour(aes(z = max_firm_profit), color = "white") +
      geom_point(data = max_point, aes(x = a, y = b), color = "red", size = 3) + 
      # Add Pareto optimal points for the second tab
      geom_point(data = pareto_points_grid2, aes(x = a, y = b), 
                 shape = 4, color = "black", size = 4, stroke = 1.5) + # Black crosses
      labs(title = "Firm's Max Profit (as function of a, b)", x = "a", y = "b") +
      scale_fill_viridis_c(option = "plasma", na.value = "grey") + 
      theme_minimal()
  },
  res = 96
  )
  
  output$utilityPlot2 <- renderPlot({
    req(reactive_grid2_data(), reactive_pareto_points_grid2()) 
    grid2_data <- reactive_grid2_data()
    pareto_points_grid2 <- reactive_pareto_points_grid2() 
    max_point <- grid2_data[which.max(grid2_data$user_utility_at_firm_opt), ] 
    
    ggplot(grid2_data, aes(x = a, y = b)) +
      geom_tile(aes(fill = user_utility_at_firm_opt)) +
      geom_point(data = max_point, aes(x = a, y = b), color = "red", size = 3) + 
      geom_contour(aes(z = user_utility_at_firm_opt), color = "white") +
      # Add Pareto optimal points for the second tab
      geom_point(data = pareto_points_grid2, aes(x = a, y = b), 
                 shape = 4, color = "black", size = 4, stroke = 1.5) + # Black crosses
      labs(title = "User Utility at Firm's Max Profit (as function of a, b)", x = "a", y = "b") +
      scale_fill_viridis_c(option = "magma", na.value = "grey") + 
      theme_minimal()
  },
  res = 96
  )
  
  # --- Helper function to safely extract a value or return NA ---
  # This prevents errors if a data frame is empty (e.g., from which.max returning integer(0))
  safe_extract <- function(df, col_name) {
    if (is.null(df) || nrow(df) == 0 || is.null(df[[col_name]])) {
      return(NA_real_) # Return NA_real_ for numeric NA
    }
    # Ensure the extracted value is treated as a scalar, even if it's a 1-element vector
    return(df[[col_name]][1]) 
  }
  
  # --- Single Definition for Summary Table ---
  output$summaryTable <- renderTable({
    # Ensure all necessary reactive data is available before rendering
    req(reactive_planner_opt(), reactive_grid1_data(), 
        reactive_grid2_data(), reactive_nash_equilibrium_data())
    
    planner_opt_data <- reactive_planner_opt()
    grid1_data <- reactive_grid1_data() 
    grid2_data <- reactive_grid2_data() 
    nash_data <- reactive_nash_equilibrium_data()
    
    # Re-find max points from the reactive grid data for the table
    # Using tryCatch with safe_extract to handle cases where which.max might return integer(0)
    max_profit_grid1_data <- tryCatch({
      grid1_data[which.max(grid1_data$profit), ]
    }, error = function(e) data.frame()) # Return empty dataframe on error
    
    max_user_utility_grid2_data <- tryCatch({
      grid2_data[which.max(grid2_data$user_utility_at_firm_opt), ]
    }, error = function(e) data.frame()) # Return empty dataframe on error
    
    base_table <- data.frame(
      Criterion = c("Max Profit (Grid1)", "Max User Utility (Grid2)"), 
      Lambda = c(
        round(safe_extract(max_profit_grid1_data, "lambda"), 3), 
        round(safe_extract(max_user_utility_grid2_data, "optimal_lambda_firm"), 3) 
      ),
      w = c(
        round(safe_extract(max_profit_grid1_data, "w"), 3), 
        round(safe_extract(max_user_utility_grid2_data, "optimal_w_firm"), 3) 
      ),
      a = c(
        round(safe_extract(max_profit_grid1_data, "a"), 3), 
        round(safe_extract(max_user_utility_grid2_data, "a"), 3) 
      ),
      b = c(
        round(safe_extract(max_profit_grid1_data, "b"), 3), 
        round(safe_extract(max_user_utility_grid2_data, "b"), 3) 
      ),
      Utility = c(
        round(safe_extract(max_profit_grid1_data, "utility"), 3), 
        round(safe_extract(max_user_utility_grid2_data, "user_utility_at_firm_opt"), 3) 
      ),
      Profit = c(
        round(safe_extract(max_profit_grid1_data, "profit"), 3), 
        round(safe_extract(max_user_utility_grid2_data, "max_firm_profit"), 3) 
      ),
      Welfare = c(
        round(safe_extract(max_profit_grid1_data, "welfare"), 3), 
        round(safe_extract(max_user_utility_grid2_data, "user_utility_at_firm_opt") + safe_extract(max_user_utility_grid2_data, "max_firm_profit"), 3) 
      )
    )
    
    # Combine all rows
    final_table <- rbind(base_table, nash_data, planner_opt_data)
    
    return(final_table)
  },
  digits = 3, 
  striped = TRUE, 
  bordered = TRUE 
  )
}

shinyApp(ui = ui, server = server)
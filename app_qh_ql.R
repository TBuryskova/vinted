library(shiny)
library(ggplot2)
library(dplyr)
library(viridis)
library(purrr)

# --- Mutual Information (analytical) ---
mutual_info_analytic <- function(r) {
  if (r <= 0 || r >= 1) return(0)
  1 + r * log2(r) + (1 - r) * log2(1 - r)
}

# --- Compute cutoff (new for log utility) ---
cutoff_pi <- function(w, VH, VL, eps) {
  # Robustness check to prevent NaN
  if (VH <= w || VL <= w || eps <= 0) {
    return(NA) # Indicate invalid parameters
  }
  
  if (VH == VL) {
    return(ifelse(log(VH - w) > log(eps), 0, 1))
  }
  
  pi_star <- log(eps / (VL - w)) / log((VH - w) / (VL - w))
  return(pi_star)
}

# --- Closed-form r* and V* (new for log utility) ---
alytic_solution <- function(VH, VL, w, eps, lambda) {
  # Robustness check
  if (VH <= w || VL <= w || eps <= 0) {
    return(list(r_star = NA, V_star = -Inf))
  }
  
  gain <- 0.5 * (log(VH - w) - log(VL - w))
  
  # The consumer's unconstrained optimal attention level
  r_star_raw <- 2^(gain / lambda) / (1 + 2^(gain / lambda))
  
  # The cutoff value for r, based on the firm's price
  pi_star <- cutoff_pi(w, VH, VL, eps)
  
  # Check if pi_star is a valid number
  if(is.na(pi_star) || VH <= w || VL <= w || eps <= 0) {
    return(list(r_star = NA, V_star = -Inf))
  }
  
  if (pi_star > 0.5) {
    r_cutoff <- (pi_star - 0.5) / 0.5
  } else {
    r_cutoff <- 1
  }
  
  if (r_star_raw > r_cutoff) {
    # Unconstrained solution is valid
    r_star <- r_star_raw
    EU <- 0.5 * (log(VL - w) + log(eps)) + 0.5 * r_star * (log(VH - w) - log(VL - w))
    MI <- mutual_info_analytic(r_star)
    V <- EU - lambda * MI
  } else {
    # Consumer chooses not to acquire information, r = 0.5
    r_star <- 0.5
    EU <- log(eps)
    MI <- 0
    V <- EU
  }
  
  list(r_star = r_star, V_star = V)
}

# --- Buy probability (new for log utility) ---
optimal_q <- function(r, VH, VL, w, eps, p) {
  # Robustness check
  if (VH <= w || VL <= w || eps <= 0) {
    return(NA)
  }
  
  if (VH == VL) {
    return(ifelse(log(VH - w) > log(eps), 1, 0))
  }
  
  pi_star <- cutoff_pi(w, VH, VL, eps)
  if (is.na(pi_star)) {
    return(NA)
  }
  
  alpha <- p * r + (1 - p) * (1 - r)
  
  # Handle potential division by zero if alpha or (1-alpha) is zero
  pi_h <- ifelse(alpha > 0, (p * r) / alpha, 0)
  pi_l <- ifelse(alpha < 1, (p * (1 - r)) / (1 - alpha), 0)
  
  q_h <- ifelse(pi_h > pi_star, 1, 0)
  q_l <- ifelse(pi_l > pi_star, 1, 0)
  q <- alpha * q_h + (1 - alpha) * q_l
  q
}

# --- Firm's objective (updated for log utility) ---
firm_profit <- function(r, w, lambda, VH, VL, eps, p) {
  q <- optimal_q(r, VH, VL, w, eps, p)
  mi <- mutual_info_analytic(r)
  
  # Robustness check for invalid parameters from `optimal_q`
  if (is.na(q)) {
    return(-Inf)
  }
  
  # Corrected: lambda * mi is added to the profit as requested.
  profit <- w * q + lambda * mi
  
  if (is.na(profit)) {
    return(-Inf) # Penalize invalid parameters
  }
  profit
}

# --- UI ---
ui <- fluidPage(
  titlePanel("Rational Inattention: Consumer and Firm Perspective"),
  sidebarLayout(
    sidebarPanel(
      h4("Model Parameters"),
      sliderInput("vh_slider", "Upper Utility (VH)", min = 5, max = 20, value = 10, step = 1),
      sliderInput("vl_slider", "Lower Utility (VL)", min = 0, max = 10, value = 5, step = 1),
      sliderInput("eps_slider", "Default Utility (ε)", min = 0.01, max = 10, value = 1, step = 0.5),
      hr(),
      uiOutput("paramWarning"), # New UI element for warning
      hr(),
      h4("Optimal Strategy Values"),
      htmlOutput("optimalValuesOutput")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Consumer",
                 fluidRow(
                   column(4, plotOutput("rStarPlot")),
                   column(4, plotOutput("VStarPlot")),
                   column(4, plotOutput("qStarPlot"))
                 )
        ),
        tabPanel("Firm (given r)",
                 mainPanel(
                   plotOutput("firmOptPlot"),
                   plotOutput("firmProfitPlot")
                 )
        ),
        tabPanel("Firm's Optimization",
                 fluidRow(
                   column(12, plotOutput("firmOptimizationPlot", height = "600px"))
                 )
        ),
        tabPanel("Game Equilibrium",
                 fluidRow(
                   actionButton("runEquilibrium", "Run Simulation", class = "btn-primary"),
                   hr(),
                   column(12, tableOutput("equilibriumTable"))
                 )
        )
      )
    )
  )
)

# --- Server ---
server <- function(input, output) {
  
  # New: Reactive value to check parameter validity
  param_validity <- reactive({
    # We must ensure w is always less than VL
    input$vh_slider > input$vl_slider && input$vl_slider > 0.5
  })
  
  output$paramWarning <- renderUI({
    if (!param_validity()) {
      div(
        class = "alert alert-danger",
        "Warning: For the log utility to be valid, VH must be > VL, and VL must be > 0.5. Please adjust sliders.",
        style = "color: red; font-weight: bold;"
      )
    }
  })
  
  reactive_data <- reactive({
    req(param_validity()) # Only proceed if parameters are valid
    
    VH <- input$vh_slider
    VL <- input$vl_slider
    eps <- input$eps_slider
    p <- 0.5 # p is still fixed
    
    # --- Firm's Optimization Tab Calculations ---
    # Updated ranges for w and lambda
    firm_lambda_vals <- seq(0.01, 10, length.out = 50) 
    firm_w_vals <- seq(0.01, VL - 0.1, length.out = 50) 
    
    firm_optimization_grid <- expand.grid(lambda = firm_lambda_vals, w = firm_w_vals)
    
    firm_optimization_grid$consumer_r_star <- mapply(
      function(w, lambda) alytic_solution(VH, VL, w, eps, lambda)$r_star,
      firm_optimization_grid$w,
      firm_optimization_grid$lambda
    )
    
    firm_optimization_grid$firm_profit <- mapply(
      function(r, w, lambda) firm_profit(r, w, lambda, VH, VL, eps, p),
      firm_optimization_grid$consumer_r_star,
      firm_optimization_grid$w,
      firm_optimization_grid$lambda
    )
    
    optimal_strategy <- firm_optimization_grid[which.max(firm_optimization_grid$firm_profit), ]
    
    optimal_r <- optimal_strategy$consumer_r_star
    optimal_w <- optimal_strategy$w
    optimal_lambda <- optimal_strategy$lambda
    
    pbuy_at_opt <- optimal_q(optimal_r, VH, VL, optimal_w, eps, p)
    consumer_utility_at_opt <- alytic_solution(VH, VL, optimal_w, eps, optimal_lambda)$V_star
    social_welfare_at_opt <- optimal_strategy$firm_profit + consumer_utility_at_opt
    
    # --- Consumer Tab Calculations ---
    # Updated ranges for w and lambda
    lambda_vals <- seq(0.01, 10, length.out = 100)
    w_vals <- seq(0.01, VL - 0.1, length.out = 30) 
    grid <- expand.grid(lambda = lambda_vals, w = w_vals)
    
    sol <- mapply(function(l, w) alytic_solution(VH, VL, w, eps, l), grid$lambda, grid$w, SIMPLIFY = FALSE)
    grid$r_star <- sapply(sol, function(x) x$r_star)
    grid$V_star <- sapply(sol, function(x) x$V_star)
    grid$q_star <- mapply(function(r, w) optimal_q(r, VH, VL, w, eps, p), grid$r_star, grid$w)
    
    # --- Firm Tab Calculations ---
    r_vals <- seq(0.51, 0.99, by = 0.01)
    df_firm_opt <- lapply(r_vals, function(r) {
      grid_firm <- expand.grid(lambda = lambda_vals, w = w_vals) 
      grid_firm$profit <- mapply(function(w, lambda) firm_profit(r, w, lambda, VH, VL, eps, p), grid_firm$w, grid_firm$lambda)
      opt <- grid_firm[which.max(grid_firm$profit), ]
      data.frame(r = r, profit = opt$profit)
    }) %>% bind_rows()
    
    results_firm_profit <- lapply(r_vals, function(r) {
      grid_firm <- expand.grid(lambda = lambda_vals, w = w_vals) 
      grid_firm$profit <- mapply(function(w, lambda) firm_profit(r, w, lambda, VH, VL, eps, p), grid_firm$w, grid_firm$lambda)
      opt <- grid_firm[which.max(grid_firm$profit), ]
      data.frame(r = r, lambda_opt = opt$lambda, w_opt = opt$w, profit = opt$profit)
    })
    df_firm_profit <- bind_rows(results_firm_profit)
    
    list(
      grid = grid,
      firm_optimization_grid = firm_optimization_grid,
      optimal_strategy = optimal_strategy,
      optimal_r = optimal_r,
      optimal_w = optimal_w,
      optimal_lambda = optimal_lambda,
      pbuy_at_opt = pbuy_at_opt,
      consumer_utility_at_opt = consumer_utility_at_opt,
      social_welfare_at_opt = social_welfare_at_opt,
      df_firm_opt = df_firm_opt,
      df_firm_profit = df_firm_profit
    )
  })
  
  # --- Game Equilibrium Reactive Value ---
  equilibrium_data <- reactiveVal(NULL)
  
  # --- New: Function to find firm's optimal w, lambda for a given r ---
  find_firm_optimal_given_r <- function(r_current, VH, VL, eps, p) {
    firm_lambda_vals <- seq(0.01, 10, length.out = 50)
    firm_w_vals <- seq(0.01, VH - 0.1, length.out = 50) 
    
    firm_optimization_grid <- expand.grid(lambda = firm_lambda_vals, w = firm_w_vals)
    
    firm_optimization_grid$firm_profit <- mapply(
      function(w, lambda) firm_profit(r_current, w, lambda, VH, VL, eps, p),
      firm_optimization_grid$w,
      firm_optimization_grid$lambda
    )
    
    optimal_strategy <- firm_optimization_grid[which.max(firm_optimization_grid$firm_profit), ]
    return(list(w = optimal_strategy$w, lambda = optimal_strategy$lambda, profit = optimal_strategy$firm_profit))
  }
  
  # --- New: Observe button click to run simulation ---
  observeEvent(input$runEquilibrium, {
    req(param_validity()) # Only run if parameters are valid
    
    # Initialize variables for the game
    VH <- input$vh_slider
    VL <- input$vl_slider
    eps <- input$eps_slider
    p <- 0.5
    
    w_prev <- -1
    lambda_prev <- -1
    r_prev <- 0.51 # Starting r for the consumer
    
    # Tolerance for convergence
    tol <- 1e-4
    
    # Data frame to store results
    history <- data.frame(
      iteration = integer(),
      w = numeric(),
      lambda = numeric(),
      r = numeric(),
      profit = numeric()
    )
    
    # Loop for the game
    for (i in 1:100) { # Max 100 iterations to prevent infinite loops
      # Step 1: Firm chooses w, lambda given current r
      firm_optimal <- find_firm_optimal_given_r(r_prev, VH, VL, eps, p)
      w_new <- firm_optimal$w
      lambda_new <- firm_optimal$lambda
      
      # Step 2: User responds with optimal r
      consumer_optimal <- alytic_solution(VH, VL, w_new, eps, lambda_new)
      r_new <- consumer_optimal$r_star
      
      # Store this iteration's results
      new_row <- data.frame(
        iteration = i,
        w = w_new,
        lambda = lambda_new,
        r = r_new,
        profit = firm_optimal$profit
      )
      history <- rbind(history, new_row)
      
      # Check for convergence
      if (abs(w_new - w_prev) < tol && abs(lambda_new - lambda_prev) < tol && abs(r_new - r_prev) < tol) {
        break
      }
      
      # Update values for the next iteration
      w_prev <- w_new
      lambda_prev <- lambda_new
      r_prev <- r_new
    }
    
    # Store the results in the reactive value
    equilibrium_data(history)
  })
  
  # --- UI Outputs ---
  output$optimalValuesOutput <- renderUI({
    req(param_validity())
    data <- reactive_data()
    div(
      p(HTML(paste0("<b>Optimal Price (w):</b> ", round(data$optimal_w, 4)))),
      p(HTML(paste0("<b>Optimal Inattention Cost (λ):</b> ", round(data$optimal_lambda, 4)))),
      p(HTML(paste0("<b>Consumer Attention (r):</b> ", round(data$optimal_r, 4)))),
      p(HTML(paste0("<b>Firm Profit:</b> ", round(data$optimal_strategy$firm_profit, 4)))),
      p(HTML(paste0("<b>Consumer Utility (V):</b> ", round(data$consumer_utility_at_opt, 4)))),
      p(HTML(paste0("<b>Buy Probability (q):</b> ", round(data$pbuy_at_opt, 4)))),
      p(HTML(paste0("<b>Social Welfare:</b> ", round(data$social_welfare_at_opt, 4))))
    )
  })
  
  # --- Consumer Tab Plots ---
  output$rStarPlot <- renderPlot({
    req(param_validity())
    data <- reactive_data()
    ggplot(data$grid, aes(x = lambda, y = w, fill = r_star)) +
      geom_tile() +
      scale_fill_viridis_c(name = "r*") +
      labs(title = "Optimal Precision r*", x = "λ", y = "w") +
      theme_minimal()
  })
  
  output$VStarPlot <- renderPlot({
    req(param_validity())
    data <- reactive_data()
    ggplot(data$grid, aes(x = lambda, y = w, fill = V_star)) +
      geom_tile() +
      scale_fill_viridis_c(name = "V*") +
      labs(title = "Value Function V*", x = "λ", y = "w") +
      theme_minimal()
  })
  
  output$qStarPlot <- renderPlot({
    req(param_validity())
    data <- reactive_data()
    ggplot(data$grid, aes(x = lambda, y = w, fill = q_star)) +
      geom_tile() +
      scale_fill_viridis_c(name = "P(buy)") +
      labs(title = "Buy Probability q*", x = "λ", y = "w") +
      theme_minimal()
  })
  
  # --- Firm Tab Plots ---
  output$firmOptPlot <- renderPlot({
    req(param_validity())
    data <- reactive_data()
    ggplot(data$df_firm_opt, aes(x = r, y = profit)) +
      geom_line(color = "black", size = 1) +
      labs(title = "Firm's Maximum Profit vs Consumer Attention r",
           x = "r (Consumer Attention)", y = "Max Profit") +
      theme_minimal()
  })
  
  output$firmProfitPlot <- renderPlot({
    req(param_validity())
    data <- reactive_data()
    ggplot(data$df_firm_profit, aes(x = r)) +
      geom_line(aes(y = lambda_opt, color = "λ*"), size = 1) +
      geom_line(aes(y = w_opt, color = "w*"), size = 1) +
      scale_color_manual(values = c("λ*" = "blue", "w*" = "green", "Profit" = "black")) +
      labs(title = "Firm’s Optimal Strategy vs Exogenous r",
           y = "Optimal Values", x = "r (Consumer Attention)", color = "") +
      theme_minimal()
  })
  
  # --- Firm's Optimization Tab Plot ---
  output$firmOptimizationPlot <- renderPlot({
    req(param_validity())
    data <- reactive_data()
    p <- ggplot(data$firm_optimization_grid, aes(x = w, y = lambda, fill = firm_profit)) +
      geom_tile() +
      scale_fill_viridis_c(name = "Firm Profit") +
      labs(title = "Firm's Optimal Strategy", x = "Price (w)", y = "Cost of Inattention (λ)") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "right"
      )
    
    p + geom_point(
      data = data$optimal_strategy,
      aes(x = w, y = lambda),
      color = "red",
      shape = 4, # 'x' shape
      size = 8,
      stroke = 2,
      inherit.aes = FALSE
    )
  })
  
  # --- New: Output for the equilibrium table ---
  output$equilibriumTable <- renderTable({
    req(param_validity())
    equilibrium_data()
  })
}

shinyApp(ui, server)

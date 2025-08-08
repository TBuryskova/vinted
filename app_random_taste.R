# app_random_taste_firm.R
library(shiny)
library(ggplot2)
library(dplyr)
library(purrr) # Added to load the map_df function
library(tidyr) # Added for pivot_longer

# --- Model primitives ---

# Logistic function
Lambda <- function(x) 1 / (1 + exp(-x))

# Mutual information in bits for binary symmetric channel
I_bits <- function(r) {
  ifelse(r == 0.5, 0,
         ifelse(r == 1, 1,
                1 + (1 - r) * (log(1 - r) / log(2)) + r * (log(r) / log(2))))
}

# Expected utility given signal precision r
EU_r <- function(r, UH, UL, c, beta, R) {
  # Net payoffs from buying
  uH_net <- UH - c
  uL_net <- UL - c
  
  # Expected payoff after high signal
  EU_high <- R + Lambda(beta * (r * uH_net + (1 - r) * uL_net - R)) *
    ((r * uH_net + (1 - r) * uL_net) - R)
  
  # Expected payoff after low signal
  EU_low <- R + Lambda(beta * ((1 - r) * uH_net + r * uL_net - R)) *
    (((1 - r) * uH_net + r * uL_net) - R)
  
  # Average over prior p=0.5
  0.5 * EU_high + 0.5 * EU_low
}

# Consumer's objective: EU - kappa * MI
U_r <- function(r, UH, UL, c, beta, R, kappa) {
  EU_r(r, UH, UL, c, beta, R) - kappa * I_bits(r)
}

# Ex-ante purchase probability
P_buy <- function(r, UH, UL, c, beta, R) {
  uH_net <- UH - c
  uL_net <- UL - c
  
  p_buy_high <- Lambda(beta * (r * uH_net + (1 - r) * uL_net - R))
  p_buy_low  <- Lambda(beta * ((1 - r) * uH_net + r * uL_net - R))
  
  0.5 * p_buy_high + 0.5 * p_buy_low
}

# Firm's profit function
Firm_Profit <- function(r, c, kappa, UH, UL, R, beta) {
  # P_buy and I_bits depend on the consumer's optimal r, which is passed here
  profit <- P_buy(r, UH, UL, c, beta, R) * c + kappa * I_bits(r)
  return(profit)
}

# Social Welfare function
Social_Welfare <- function(r, c, kappa, UH, UL, R, beta) {
  firm_profit <- Firm_Profit(r, c, kappa, UH, UL, R, beta)
  consumer_utility <- U_r(r, UH, UL, c, beta, R, kappa)
  firm_profit + consumer_utility
}


# --- UI ---
ui <- fluidPage(
  titlePanel("Rational Inattention – Firm Profit Maximization"),
  tabsetPanel(
    tabPanel("Firm Acts First",
             sidebarLayout(
               sidebarPanel(
                 h4("Consumer's Parameters"),
                 sliderInput("UH", "Utility High (UH):", min = 0, max = 10, value = 5, step = 0.5),
                 sliderInput("UL", "Utility Low (UL):", min = 0, max = 10, value = 0, step = 0.5),
                 sliderInput("R", "Reservation payoff (R):", min = 0, max = 10, value = 2.5, step = 0.5),
                 hr(),
                 h4("Firm's Choice Parameters"),
                 sliderInput("c_range", "Cost Range (c):", min = 0, max = 10, value = c(0, 5), step = 0.5),
                 sliderInput("k_range", "Kappa Range (k):", min = 0, max = 10, value = c(0, 5), step = 0.5),
                 sliderInput("steps", "Grid steps:", min = 10, max = 100, value = 25, step = 1)
               ),
               mainPanel(
                 h3("Firm's Profit Heatmap"),
                 plotOutput("heat_profit", height = "400px"),
                 hr(),
                 h3("Optimal Solution"),
                 textOutput("optimal_solution")
               )
             )
    ),
    tabPanel("Backward Induction",
             sidebarLayout(
               sidebarPanel(
                 h4("Consumer's Parameters"),
                 sliderInput("UH_cf", "Utility High (UH):", min = 0, max = 10, value = 5, step = 0.5),
                 sliderInput("UL_cf", "Utility Low (UL):", min = 0, max = 10, value = 0, step = 0.5),
                 sliderInput("R_cf", "Reservation payoff (R):", min = 0, max = 10, value = 0.5, step = 0.5),
                 hr(),
                 h4("Firm's Best Response Parameters"),
                 sliderInput("r_range_cf", "Consumer's r range:", min = 0.5, max = 1, value = c(0.5, 1), step = 0.01),
                 sliderInput("c_range_cf", "Firm's c range:", min = 1, max = 10, value = c(1, 5), step = 0.5),
                 sliderInput("k_range_cf", "Firm's kappa range:", min = 0, max = 10, value = c(1, 5), step = 0.5),
                 sliderInput("steps_cf", "Grid steps:", min = 10, max = 50, value = 25, step = 1)
               ),
               mainPanel(
                 h3("Firm's Best Response Functions"),
                 plotOutput("best_response_plot", height = "400px"),
                 hr(),
                 h3("Equilibrium Solution (Backward Induction)"),
                 textOutput("optimal_response_cf")
               )
             )
    ),
    tabPanel("Iterative Best Responses",
             sidebarLayout(
               sidebarPanel(
                 h4("Simulation Parameters"),
                 p("Initial Cost (c): 2"),
                 p("Initial Kappa (κ): 2"),
                 p("Number of Iterations: 10")
               ),
               mainPanel(
                 h3("Best Response Dynamics"),
                 plotOutput("iterative_plot", height = "400px"),
                 hr(),
                 h3("Equilibrium Values"),
                 textOutput("equilibrium_values")
               )
             )
    ),
    tabPanel("Social Planner",
             sidebarLayout(
               sidebarPanel(
                 h4("Social Planner Parameters"),
                 sliderInput("UH_sp", "Utility High (UH):", min = 0, max = 10, value = 5, step = 0.5),
                 sliderInput("UL_sp", "Utility Low (UL):", min = 0, max = 10, value = 0, step = 0.5),
                 sliderInput("R_sp", "Reservation payoff (R):", min = 0, max = 10, value = 0.5, step = 0.5),
                 hr(),
                 h4("Optimization Ranges"),
                 sliderInput("c_range_sp", "Cost Range (c):", min = 0, max = 10, value = c(0, 5), step = 0.5),
                 sliderInput("k_range_sp", "Kappa Range (k):", min = 0, max = 10, value = c(0, 5), step = 0.5),
                 sliderInput("r_range_sp", "r Range (r):", min = 0.5, max = 1, value = c(0.5, 1), step = 0.01),
                 sliderInput("steps_sp", "Grid steps:", min = 10, max = 50, value = 25, step = 1)
               ),
               mainPanel(
                 h3("Social Planner's Optimal Solution"),
                 textOutput("social_planner_solution"),
                 hr(),
                 h3("Social Welfare at Optimum"),
                 textOutput("social_welfare_value")
               )
             )
    ),
    tabPanel("Summary of All Scenarios",
             mainPanel(
               h3("Comparison of Optimal Outcomes"),
               tableOutput("all_scenarios_summary")
             )
    )
  )
)

# --- Server ---
server <- function(input, output, session) {
  beta <- 2
  
  # --- Firm Acts First Tab Logic ---
  firm_sweep_data <- reactive({
    c_vals <- seq(input$c_range[1], input$c_range[2], length.out = input$steps)
    k_vals <- seq(input$k_range[1], input$k_range[2], length.out = input$steps)
    
    expand.grid(c = c_vals, kappa = k_vals) %>%
      rowwise() %>%
      mutate(
        sol_consumer = list(optimize(function(r) -U_r(r, input$UH, input$UL, c, beta, input$R, kappa),
                                     interval = c(0.5, 1), maximum = FALSE)),
        r_star = sol_consumer$minimum,
        profit = Firm_Profit(r_star, c, kappa, input$UH, input$UL, input$R, beta)
      ) %>%
      select(-sol_consumer) %>%
      ungroup()
  })
  
  firm_optimal_solution <- reactive({
    data <- firm_sweep_data()
    data[which.max(data$profit),]
  })
  
  output$heat_profit <- renderPlot({
    ggplot(firm_sweep_data(), aes(x = c, y = kappa, fill = profit)) +
      geom_tile() +
      scale_fill_viridis_c(name = "Firm Profit") +
      labs(title = "Firm's Profit as a Function of c and kappa", x = "Cost c", y = "kappa") +
      theme_minimal()
  })
  
  output$optimal_solution <- renderText({
    sol <- firm_optimal_solution()
    paste0("The firm's optimal choices are:\n",
           "  - Optimal Cost (c*): ", round(sol$c, 2), "\n",
           "  - Optimal Kappa (κ*): ", round(sol$kappa, 2), "\n",
           "The consumer's optimal response is:\n",
           "  - Optimal Signal Precision (r*): ", round(sol$r_star, 2), "\n",
           "Leading to a maximum firm profit of: ", round(sol$profit, 2))
  })
  
  # --- Backward Induction Tab Logic ---
  backward_induction_data <- reactive({
    r_vals <- seq(input$r_range_cf[1], input$r_range_cf[2], length.out = input$steps_cf)
    c_vals <- seq(input$c_range_cf[1], input$c_range_cf[2], length.out = input$steps_cf)
    k_vals <- seq(input$k_range_cf[1], input$k_range_cf[2], length.out = input$steps_cf)
    
    # This function finds the firm's best response (c*, kappa*) for a given r
    firm_best_response <- function(r) {
      grid <- expand.grid(c = c_vals, kappa = k_vals)
      grid$profit <- Firm_Profit(r, grid$c, grid$kappa, input$UH_cf, input$UL_cf, input$R_cf, beta)
      best_firm_response <- grid[which.max(grid$profit), ]
      return(best_firm_response)
    }
    
    # Now, the consumer chooses the r that maximizes their utility,
    # anticipating the firm's best response
    map_df(r_vals, function(r) {
      firm_sol <- firm_best_response(r)
      c_star <- firm_sol$c
      k_star <- firm_sol$kappa
      
      consumer_utility <- U_r(r, input$UH_cf, input$UL_cf, c_star, beta, input$R_cf, k_star)
      
      data.frame(r = r,
                 c_star = c_star,
                 k_star = k_star,
                 firm_profit = firm_sol$profit,
                 consumer_utility = consumer_utility)
    })
  })
  
  consumer_first_optimal_solution <- reactive({
    data <- backward_induction_data()
    data[which.max(data$consumer_utility),]
  })
  
  output$best_response_plot <- renderPlot({
    data <- backward_induction_data()
    
    # Reshape data for plotting both c_star and k_star on the same y-axis
    plot_data <- data %>%
      pivot_longer(cols = c(c_star, k_star),
                   names_to = "parameter",
                   values_to = "value")
    
    ggplot(plot_data, aes(x = r, y = value, color = parameter)) +
      geom_line(size = 1) +
      labs(title = "Firm's Best Response to Consumer's Information Choice",
           x = "Consumer's Signal Precision (r)",
           y = "Firm's Optimal Choice",
           color = "Parameter") +
      scale_color_manual(values = c("c_star" = "blue", "k_star" = "red"),
                         labels = c("Optimal Cost (c*)", "Optimal Kappa (κ*)")) +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$optimal_response_cf <- renderText({
    sol <- consumer_first_optimal_solution()
    paste0("The user's optimal choice is:\n",
           "  - Optimal Signal Precision (r*): ", round(sol$r, 2), "\n",
           "The firm's optimal response is:\n",
           "  - Optimal Cost (c*): ", round(sol$c_star, 2), "\n",
           "  - Optimal Kappa (κ*): ", round(sol$k_star, 2), "\n",
           "Leading to a maximum user utility of: ", round(sol$consumer_utility, 2))
  })
  
  # --- Iterative Best Responses Tab Logic ---
  iterative_results <- reactive({
    
    # Fixed initial values and iterations
    c_current <- 2
    k_current <- 2
    iterations <- 10
    
    # Create empty data frame to store results
    history <- data.frame(iteration = 0, c = c_current, kappa = k_current, r = NA)
    
    for (i in 1:iterations) {
      # Consumer's Best Response to current firm parameters
      r_new_opt <- optimize(function(r) -U_r(r, 5, 0, c_current, beta, 0.5, k_current),
                            interval = c(0.5, 1), maximum = FALSE)
      r_new <- r_new_opt$minimum
      
      # Firm's Best Response to new consumer parameter
      firm_profit_grid <- expand.grid(c = seq(0, 10, length.out = 100),
                                      kappa = seq(0, 10, length.out = 100))
      firm_profit_grid$profit <- Firm_Profit(r_new, firm_profit_grid$c, firm_profit_grid$kappa, 5, 0, 0.5, beta)
      
      best_firm_response <- firm_profit_grid[which.max(firm_profit_grid$profit), ]
      c_new <- best_firm_response$c
      k_new <- best_firm_response$kappa
      
      # Add new values to history
      history <- rbind(history, data.frame(iteration = i, c = c_new, kappa = k_new, r = r_new))
      
      # Update current values for next iteration
      c_current <- c_new
      k_current <- k_new
    }
    
    history
  })
  
  output$iterative_plot <- renderPlot({
    plot_data <- iterative_results() %>%
      pivot_longer(cols = c(c, kappa, r),
                   names_to = "parameter",
                   values_to = "value")
    
    ggplot(plot_data, aes(x = iteration, y = value, color = parameter)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      labs(title = "Best Response Dynamics Over Iterations",
           x = "Iteration",
           y = "Parameter Value",
           color = "Parameter") +
      scale_color_manual(values = c("c" = "blue", "kappa" = "red", "r" = "green4")) +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$equilibrium_values <- renderText({
    final_state <- tail(iterative_results(), 1)
    paste0("The simulation has reached the following approximate equilibrium after ", final_state$iteration, " iterations:\n",
           "  - Firm's Optimal Cost (c): ", round(final_state$c, 2), "\n",
           "  - Firm's Optimal Kappa (κ): ", round(final_state$kappa, 2), "\n",
           "  - Consumer's Optimal r: ", round(final_state$r, 2))
  })
  
  # --- Social Planner Tab Logic ---
  social_planner_results <- reactive({
    c_vals <- seq(input$c_range_sp[1], input$c_range_sp[2], length.out = input$steps_sp)
    k_vals <- seq(input$k_range_sp[1], input$k_range_sp[2], length.out = input$steps_sp)
    r_vals <- seq(input$r_range_sp[1], input$r_range_sp[2], length.out = input$steps_sp)
    
    grid <- expand.grid(c = c_vals, kappa = k_vals, r = r_vals)
    
    grid$welfare <- Social_Welfare(grid$r, grid$c, grid$kappa, input$UH_sp, input$UL_sp, input$R_sp, beta)
    
    best_solution <- grid[which.max(grid$welfare), ]
    
    best_solution
  })
  
  output$social_planner_solution <- renderText({
    sol <- social_planner_results()
    paste0("The social planner's optimal choices are:\n",
           "  - Optimal Cost (c*): ", round(sol$c, 2), "\n",
           "  - Optimal Kappa (κ*): ", round(sol$kappa, 2), "\n",
           "  - Optimal Signal Precision (r*): ", round(sol$r, 2))
  })
  
  output$social_welfare_value <- renderText({
    sol <- social_planner_results()
    paste0("The maximum social welfare is: ", round(sol$welfare, 2))
  })
  
  # --- Summary of All Scenarios ---
  firm_first_summary <- reactive({
    sol <- firm_optimal_solution()
    r_star <- sol$r_star
    c_star <- sol$c
    k_star <- sol$kappa
    
    list(
      r = r_star,
      c = c_star,
      kappa = k_star,
      p_buy = P_buy(r_star, input$UH, input$UL, c_star, beta, input$R),
      profit = sol$profit,
      utility = U_r(r_star, input$UH, input$UL, c_star, beta, input$R, k_star),
      welfare = Social_Welfare(r_star, c_star, k_star, input$UH, input$UL, input$R, beta)
    )
  })
  
  backward_induction_summary <- reactive({
    sol <- consumer_first_optimal_solution()
    r_star <- sol$r
    c_star <- sol$c_star
    k_star <- sol$k_star
    
    list(
      r = r_star,
      c = c_star,
      kappa = k_star,
      p_buy = P_buy(r_star, input$UH_cf, input$UL_cf, c_star, beta, input$R_cf),
      profit = Firm_Profit(r_star, c_star, k_star, input$UH_cf, input$UL_cf, input$R_cf, beta),
      utility = sol$consumer_utility,
      welfare = Social_Welfare(r_star, c_star, k_star, input$UH_cf, input$UL_cf, input$R_cf, beta)
    )
  })
  
  iterative_summary <- reactive({
    final_state <- tail(iterative_results(), 1)
    r_val <- final_state$r
    c_val <- final_state$c
    k_val <- final_state$kappa
    
    list(
      r = r_val,
      c = c_val,
      kappa = k_val,
      p_buy = P_buy(r_val, 5, 0, c_val, beta, 0.5),
      profit = Firm_Profit(r_val, c_val, k_val, 5, 0, 0.5, beta),
      utility = U_r(r_val, 5, 0, c_val, beta, 0.5, k_val),
      welfare = Social_Welfare(r_val, c_val, k_val, 5, 0, 0.5, beta)
    )
  })
  
  social_planner_summary <- reactive({
    sol <- social_planner_results()
    r_star <- sol$r
    c_star <- sol$c
    k_star <- sol$kappa
    
    list(
      r = r_star,
      c = c_star,
      kappa = k_star,
      p_buy = P_buy(r_star, input$UH_sp, input$UL_sp, c_star, beta, input$R_sp),
      profit = Firm_Profit(r_star, c_star, k_star, input$UH_sp, input$UL_sp, input$R_sp, beta),
      utility = U_r(r_star, input$UH_sp, input$UL_sp, c_star, beta, input$R_sp, k_star),
      welfare = sol$welfare
    )
  })
  
  output$all_scenarios_summary <- renderTable({
    firm_sol <- firm_first_summary()
    bi_sol <- backward_induction_summary()
    iter_sol <- iterative_summary()
    sp_sol <- social_planner_summary()
    
    data.frame(
      Scenario = c("Firm Acts First", "Backward Induction", "Iterative Best Responses", "Social Planner"),
      `r*` = c(firm_sol$r, bi_sol$r, iter_sol$r, sp_sol$r),
      `c*` = c(firm_sol$c, bi_sol$c, iter_sol$c, sp_sol$c),
      `κ*` = c(firm_sol$kappa, bi_sol$kappa, iter_sol$kappa, sp_sol$kappa),
      `P(buy)` = c(firm_sol$p_buy, bi_sol$p_buy, iter_sol$p_buy, sp_sol$p_buy),
      `Profit` = c(firm_sol$profit, bi_sol$profit, iter_sol$profit, sp_sol$profit),
      `Utility` = c(firm_sol$utility, bi_sol$utility, iter_sol$utility, sp_sol$utility),
      `Welfare` = c(firm_sol$welfare, bi_sol$welfare, iter_sol$welfare, sp_sol$welfare)
    ) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
  }, striped = TRUE, bordered = TRUE, rownames = FALSE, digits = 2, sanitize.text.function = function(x) x)
}

shinyApp(ui, server)

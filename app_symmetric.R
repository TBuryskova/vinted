library(shiny)
library(ggplot2)
library(dplyr)
library(DT) # For interactive tables

# Helper functions

# Function to calculate a_star (from equation 10)
a_star_val <- function(lambda, w, UL, epsilon) {
  if (lambda <= 0) return(NA) # Lambda must be positive
  val <- (UL - w - epsilon) / (2 * lambda)
  
  # Clip a_star between 0 and 1
  if (is.infinite(val) || is.na(val)) return(NA)
  return(pmax(0, pmin(1, 0.5 + val)))
}

# Function to calculate b_star (from equation 11)
b_star_val <- function(lambda, w, UH, epsilon) {
  if (lambda <= 0) return(NA) # Lambda must be positive
  val <- (UH - w - epsilon) / (2 * lambda)
  
  # Clip b_star between 0 and 1
  if (is.infinite(val) || is.na(val)) return(NA)
  return(pmax(0, pmin(1, 0.5 + val)))
}

# NEW: Robust and consistent mutual information calculation
calculate_mutual_information <- function(a, b) {
  # Clip a, b to ensure they are within [0, 1] for robustness
  a <- pmax(0, pmin(1, a))
  b <- pmax(0, pmin(1, b))
  
  s <- a + b
  
  # Adjust s and 2-s for denominators if they are exactly zero
  # Using a small positive epsilon to prevent division by zero or log(0)
  s_denom <- ifelse(s == 0, 1e-9, s) 
  two_minus_s_denom <- ifelse( (2-s) == 0, 1e-9, 2-s) 
  
  # Helper for individual log terms: prob * log(ratio)
  log_term_safe <- function(prob, numerator, denominator) {
    if (prob == 0) return(0) # If probability is zero, the term is zero
    if (denominator <= 0) return(NA) # Should ideally not happen with denom adjustments
    
    ratio <- numerator / denominator
    # If ratio is non-positive, clip it to a small positive value before log
    if (ratio <= 0) ratio <- 1e-9 
    
    prob * log(ratio)
  }
  
  mi_term_a1 <- log_term_safe(a, 2*a, s_denom)
  mi_term_a2 <- log_term_safe(1-a, 2*(1-a), two_minus_s_denom)
  mi_term_b1 <- log_term_safe(b, 2*b, s_denom)
  mi_term_b2 <- log_term_safe(1-b, 2*(1-b), two_minus_s_denom)
  
  # If any term is NA, the whole MI is problematic
  if (any(is.na(c(mi_term_a1, mi_term_a2, mi_term_b1, mi_term_b2)))) return(NA)
  
  mi_val <- 0.5 * (mi_term_a1 + mi_term_a2 + mi_term_b1 + mi_term_b2)
  
  if (is.infinite(mi_val) || is.na(mi_val)) return(NA)
  return(mi_val)
}

# Updated: Function to calculate mutual information I(s, \hat{s}) (equation 9)
# Now uses the common calculate_mutual_information helper
Optimal_Mutual_Information <- function(lambda, w, UL, UH, epsilon) {
  a_s <- a_star_val(lambda, w, UL, epsilon)
  b_s <- b_star_val(lambda, w, UH, epsilon)
  
  if (is.na(a_s) || is.na(b_s)) return(NA)
  
  return(calculate_mutual_information(a_s, b_s))
}

# Updated: User Utility Function (equation 8 with a_star, b_star substituted)
# Now uses the common calculate_mutual_information helper
User_Utility_Function <- function(lambda, w, UL, UH, epsilon) {
  a_s <- a_star_val(lambda, w, UL, epsilon)
  b_s <- b_star_val(lambda, w, UH, epsilon)
  
  if (is.na(a_s) || is.na(b_s)) return(NA) # Propagate NA if a_s or b_s are invalid
  
  term_a <- 0.5 * (a_s * (UL - w) + (1 - a_s) * epsilon)
  term_b <- 0.5 * (b_s * (UH - w) + (1 - b_s) * epsilon)
  
  mutual_info <- calculate_mutual_information(a_s, b_s) # Use common MI function
  
  if (is.na(mutual_info)) return(NA)
  
  user_u <- term_a + term_b - lambda * mutual_info
  
  if (is.infinite(user_u) || is.na(user_u)) return(NA)
  
  return(user_u)
}

# Unconditional Probability of Buying (equation 12)
Unconditional_Prob_Buying <- function(lambda, w, UL, UH, epsilon) {
  a_s <- a_star_val(lambda, w, UL, epsilon)
  b_s <- b_star_val(lambda, w, UH, epsilon)
  
  if (is.na(a_s) || is.na(b_s)) return(NA)
  
  prob <- 0.5 * a_s + 0.5 * b_s
  
  if (is.infinite(prob) || is.na(prob)) return(NA)
  
  return(prob)
}

# Firm Profit Function (equation 13)
Firm_Profit_Function <- function(lambda, w, UL, UH, epsilon) {
  prob_buying <- Unconditional_Prob_Buying(lambda, w, UL, UH, epsilon)
  
  if (is.na(prob_buying)) return(NA)
  
  profit <- prob_buying * w
  
  if (is.infinite(profit) || is.na(profit)) return(NA)
  
  return(profit)
}

# Updated: User Utility Function given direct a,b inputs
# Now uses the common calculate_mutual_information helper
User_Utility_Given_ab_lambda_w <- function(a, b, lambda, w, UL, UH, epsilon) {
  # Validate a, b are within [0, 1] with a small tolerance for floating point
  a <- pmax(0, pmin(1, a)) # Clip a,b to be within [0,1]
  b <- pmax(0, pmin(1, b))
  
  if (lambda <= 0 || is.infinite(lambda) || is.na(lambda)) return(NA)
  
  term1_part <- 0.5 * (a * (UL - w) + (1 - a) * epsilon) +
    0.5 * (b * (UH - w) + (1 - b) * epsilon)
  
  mutual_info_calculated <- calculate_mutual_information(a, b) # Use common MI function
  
  if (is.na(mutual_info_calculated)) return(NA)
  
  result <- term1_part - lambda * mutual_info_calculated
  
  if (is.na(result) || is.infinite(result)) return(NA)
  return(result)
}

# NEW: Social Welfare Objective Function for optim
Social_Welfare_Objective <- function(params, UL, UH, epsilon) {
  a <- params[1]
  b <- params[2]
  lambda <- params[3]
  w <- params[4]
  
  # Return Inf if lambda is non-positive (optim minimizes, so Inf is bad)
  if (lambda <= 0) return(Inf) 
  
  # Firm Profit based on direct a,b (not a_star from lambda,w)
  # a and b are decision variables for the social planner
  prob_buying_sp <- 0.5 * a + 0.5 * b
  firm_profit_val <- prob_buying_sp * w
  
  # User Utility using the robust function for direct a,b, lambda, w
  user_utility_val <- User_Utility_Given_ab_lambda_w(a, b, lambda, w, UL, UH, epsilon)
  
  # If any calculation results in NA or Inf, return Inf for minimization
  if (is.na(firm_profit_val) || is.infinite(firm_profit_val) ||
      is.na(user_utility_val) || is.infinite(user_utility_val)) {
    return(Inf)
  }
  
  welfare <- firm_profit_val + user_utility_val
  
  # optim minimizes, so return negative welfare for maximization
  return(-welfare)
}


# --- Shiny UI ---
ui <- fluidPage(
  titlePanel("Optimal Point Plotter for Firm Profit and Scenario Comparison"),
  sidebarLayout(
    sidebarPanel(
      h3("Model Parameters"),
      sliderInput("UL", "U_L:", min = 0, max = 10, value = 0, step = 0.1),
      sliderInput("UH", "U_H:", min = 0, max = 10, value = 10, step = 0.1),
      sliderInput("epsilon", "epsilon:", min =0, max = 10, value = 5, step = 0.1),
      hr(),
      numericInput("grid_res", "General Grid Resolution (points per axis):", value = 50, min = 10, max = 200)
    ),
    mainPanel(
      tabsetPanel(id = "main_tabs", # Added an id to reference the tabsetPanel
                  tabPanel("Firm First Plot", 
                           fluidRow(
                             column(6, sliderInput("lambda_range", "Lambda (λ) Range (min, max):", min = 0.01, max = 10, value = c(0.01, 10), step = 0.1)),
                             column(6, sliderInput("w_range", "W Range (min, max):", min = 0, max = 10, value = c(0, 10), step = 0.1))
                           ),
                           plotOutput("profit_plot", height = "600px")
                  ),
                  tabPanel("User First Plot", 
                           plotOutput("user_utility_plot", height = "600px")
                  ), 
                  # NEW: Social Planner Tab
                  tabPanel("Social Planner", 
                           actionButton("calculate_social_planner", "Calculate Social Planner Optimum"),
                           hr(),
                           DTOutput("social_planner_table") # Table to show optimization results
                  ),
                  tabPanel("Scenario Comparison", DTOutput("optimal_values_table")),
                  tabPanel("About", uiOutput("about_text"))
      )
    )
  )
)

# --- Shiny Server ---
server <- function(input, output) {
  
  # Reactive expression to generate the grid and calculate Objective Function values (Firm First)
  grid_data <- reactive({
    req(input$lambda_range, input$w_range, input$grid_res)
    
    lambda_vals <- seq(input$lambda_range[1], input$lambda_range[2], length.out = input$grid_res)
    w_vals <- seq(input$w_range[1], input$w_range[2], length.out = input$grid_res)
    
    plot_df <- expand.grid(lambda = lambda_vals, w = w_vals)
    
    plot_df %>%
      rowwise() %>%
      mutate(
        firm_profit = Firm_Profit_Function(lambda, w, input$UL, input$UH, input$epsilon)
      ) %>%
      ungroup()
  })
  
  # Reactive expression to find the global maximum of the firm profit on the grid (Firm First)
  firm_first_optimal_point <- reactive({
    df <- grid_data() %>% drop_na(firm_profit)
    if (nrow(df) == 0) return(NULL)
    
    max_row <- df %>% slice_max(firm_profit, n = 1, with_ties = FALSE)
    
    if (nrow(max_row) > 0) {
      optimal_lambda <- max_row$lambda
      optimal_w <- max_row$w
      
      max_row$user_utility <- User_Utility_Function(optimal_lambda, optimal_w, input$UL, input$UH, input$epsilon)
      max_row$uncond_prob_buying <- Unconditional_Prob_Buying(optimal_lambda, optimal_w, input$UL, input$UH, input$epsilon)
      max_row$optimal_mutual_info <- Optimal_Mutual_Information(optimal_lambda, optimal_w, input$UL, input$UH, input$epsilon)
      
      max_row$a_star_opt <- a_star_val(optimal_lambda, optimal_w, input$UL, input$epsilon)
      max_row$b_star_opt <- b_star_val(optimal_lambda, optimal_w, input$UH, input$epsilon)
    }
    max_row
  })
  
  # Reactive expression for User First scenario results (fixed lambda_max, w_max)
  user_first_results <- reactive({
    # We need to ensure lambda_range and w_range are available before using them
    req(input$lambda_range, input$w_range) 
    
    optimal_lambda_uf <- input$lambda_range[2] # Max lambda from UI slider
    optimal_w_uf <- input$w_range[2]       # Max w from UI slider
    
    a_star_uf <- a_star_val(optimal_lambda_uf, optimal_w_uf, input$UL, input$epsilon)
    b_star_uf <- b_star_val(optimal_lambda_uf, optimal_w_uf, input$UH, input$epsilon)
    
    firm_profit_uf <- (0.5 * a_star_uf + 0.5 * b_star_uf) * optimal_w_uf # Firm profit at these a_star, b_star, w
    user_utility_uf <- User_Utility_Function(optimal_lambda_uf, optimal_w_uf, input$UL, input$UH, input$epsilon)
    uncond_prob_buying_uf <- Unconditional_Prob_Buying(optimal_lambda_uf, optimal_w_uf, input$UL, input$UH, input$epsilon)
    optimal_mutual_info_uf <- Optimal_Mutual_Information(optimal_lambda_uf, optimal_w_uf, input$UL, input$UH, input$epsilon)
    
    data.frame(
      Scenario = "User First",
      Optimal_Lambda = optimal_lambda_uf,
      Optimal_W = optimal_w_uf,
      Firm_Profit = firm_profit_uf,
      User_Utility = user_utility_uf,
      Welfare = firm_profit_uf + user_utility_uf,
      Uncond_Prob_Buying = uncond_prob_buying_uf,
      Optimal_Mutual_Info = optimal_mutual_info_uf,
      Optimal_a_star = a_star_uf,
      Optimal_b_star = b_star_uf,
      stringsAsFactors = FALSE
    )
  })
  
  # Reactive expression for User First plot data (user utility surface in a-b space)
  user_first_plot_data <- reactive({
    # Use a higher resolution for user plot for better visual accuracy
    user_plot_grid_res <- 100 # Increased resolution
    req(input$lambda_range, input$w_range) # Ensure these are available
    
    a_vals <- seq(0, 1, length.out = user_plot_grid_res)
    b_vals <- seq(0, 1, length.out = user_plot_grid_res)
    
    # Use the lambda and w from the user first scenario (max values from sliders)
    fixed_lambda <- input$lambda_range[2]
    fixed_w <- input$w_range[2]
    
    plot_df_user <- expand.grid(a = a_vals, b = b_vals)
    
    plot_df_user %>%
      rowwise() %>%
      mutate(
        user_utility_val = User_Utility_Given_ab_lambda_w(a, b, fixed_lambda, fixed_w, input$UL, input$UH, input$epsilon)
      ) %>%
      ungroup()
  })
  
  # NEW: Reactive for Social Planner results (triggered by button)
  social_planner_results <- eventReactive(input$calculate_social_planner, {
    req(input$UL, input$UH, input$epsilon, input$lambda_range, input$w_range)
    
    # Define bounds for optimization
    lower_bounds <- c(a = 0, b = 0, lambda = input$lambda_range[1], w = input$w_range[1])
    upper_bounds <- c(a = 1, b = 1, lambda = input$lambda_range[2], w = input$w_range[2])
    
    # Initial parameters for optim: midpoints of ranges (or can try other strategies)
    initial_params <- c(
      a = 0.5, 
      b = 0.5, 
      lambda = (input$lambda_range[1] + input$lambda_range[2]) / 2,
      w = (input$w_range[1] + input$w_range[2]) / 2
    )
    
    # Add progress bar for optimization
    withProgress(message = 'Calculating Social Planner Optimum...', value = 0.5, {
      optim_result <- optim(
        par = initial_params,
        fn = Social_Welfare_Objective,
        method = "L-BFGS-B",
        lower = lower_bounds,
        upper = upper_bounds,
        UL = input$UL,
        UH = input$UH,
        epsilon = input$epsilon,
        control = list(factr = 1e7) # For convergence tolerance; higher is more tolerant
      )
    })
    
    # Extract results
    optimal_a_sp <- optim_result$par[1]
    optimal_b_sp <- optim_result$par[2]
    optimal_lambda_sp <- optim_result$par[3]
    optimal_w_sp <- optim_result$par[4]
    
    # Calculate metrics at the found optimal point for consistency
    sp_firm_profit <- (0.5 * optimal_a_sp + 0.5 * optimal_b_sp) * optimal_w_sp
    sp_user_utility <- User_Utility_Given_ab_lambda_w(optimal_a_sp, optimal_b_sp, optimal_lambda_sp, optimal_w_sp, input$UL, input$UH, input$epsilon)
    sp_welfare <- sp_firm_profit + sp_user_utility # This should match -optim_result$value
    sp_uncond_prob_buying <- 0.5 * optimal_a_sp + 0.5 * optimal_b_sp
    sp_mutual_info <- calculate_mutual_information(optimal_a_sp, optimal_b_sp)
    
    data.frame(
      Scenario = "Social Planner",
      Optimal_Lambda = optimal_lambda_sp,
      Optimal_W = optimal_w_sp,
      Firm_Profit = sp_firm_profit,
      User_Utility = sp_user_utility,
      Welfare = sp_welfare,
      Uncond_Prob_Buying = sp_uncond_prob_buying,
      Optimal_Mutual_Info = sp_mutual_info,
      Optimal_a_star = optimal_a_sp,
      Optimal_b_star = optimal_b_sp,
      stringsAsFactors = FALSE
    )
  })
  
  # Render the plot for Firm First
  output$profit_plot <- renderPlot({
    df_clean <- grid_data() %>% drop_na(firm_profit)
    max_point <- firm_first_optimal_point()
    
    p <- ggplot(df_clean, aes(x = lambda, y = w)) +
      geom_contour_filled(aes(z = firm_profit), alpha = 0.7) +
      scale_fill_viridis_d(option = "plasma", name = "Firm Profit") +
      labs(title = "Firm First Optimal Point",
           subtitle = paste0("U_L = ", input$UL, ", U_H = ", input$UH, ", epsilon = ", input$epsilon),
           x = "Lambda (λ)",
           y = "W") +
      theme_minimal() +
      theme(legend.position = "bottom", # Changed legend position
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) +
      coord_cartesian(xlim = input$lambda_range, ylim = input$w_range)
    
    if (!is.null(max_point) && nrow(max_point) > 0) {
      p <- p + 
        annotate("point", x = max_point$lambda, y = max_point$w,
                 color = "red", shape = 8, size = 5, stroke = 2) +
        annotate("text", x = max_point$lambda, y = max_point$w,
                 label = paste0("(λ=", round(max_point$lambda, 2), ", w=", round(max_point$w, 2), ")"),
                 vjust = -1.5, hjust = 0.5, color = "red", size = 5)
    }
    p
  })
  
  # Render the plot for User First (Utility in a-b space)
  output$user_utility_plot <- renderPlot({
    df_user_clean <- user_first_plot_data() %>% drop_na(user_utility_val)
    uf_opt_point <- user_first_results()
    
    p <- ggplot(df_user_clean, aes(x = a, y = b)) +
      geom_contour_filled(aes(z = user_utility_val), alpha = 0.7) +
      scale_fill_viridis_d(option = "cividis", name = "User Utility") + # Changed color scale for distinction
      labs(title = "User First Optimal Point (User Utility in a-b space)",
           subtitle = paste0("Fixed λ = ", round(uf_opt_point$Optimal_Lambda, 2), ", Fixed w = ", round(uf_opt_point$Optimal_W, 2),
                             ", U_L = ", input$UL, ", U_H = ", input$UH, ", epsilon = ", input$epsilon),
           x = "Probability of Buying (a)",
           y = "Probability of Buying (b)") +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) # a and b are probabilities [0,1]
    
    # Mark the optimal a*, b* for the User First scenario
    if (!is.null(uf_opt_point) && nrow(uf_opt_point) > 0 && 
        !is.na(uf_opt_point$Optimal_a_star) && !is.na(uf_opt_point$Optimal_b_star)) {
      
      opt_a <- uf_opt_point$Optimal_a_star
      opt_b <- uf_opt_point$Optimal_b_star
      
      p <- p + 
        annotate("point", x = opt_a, y = opt_b,
                 color = "red", shape = 8, size = 5, stroke = 2) +
        annotate("text", x = opt_a, y = opt_b,
                 label = paste0("(a*=", round(opt_a, 2), ", b*=", round(opt_b, 2), ")"),
                 vjust = -1.5, hjust = 0.5, color = "red", size = 5)
    }
    p
  })
  
  # NEW: Render the Social Planner optimal values output as a DT table
  output$social_planner_table <- renderDT({
    sp_df <- social_planner_results() # This reactive will be NULL until button is clicked
    if (is.null(sp_df) || nrow(sp_df) == 0) {
      return(datatable(data.frame(Note = "Click 'Calculate Social Planner Optimum' to run the optimization."), 
                       options = list(dom = 't', paging = FALSE, searching = FALSE)))
    }
    datatable(sp_df %>% mutate(across(where(is.numeric), ~ round(., 4))), 
              options = list(dom = 't', paging = FALSE, searching = FALSE),
              caption = paste0("Social Planner Optimal Outcome (U_L = ", input$UL, ", U_H = ", input$UH, ", epsilon = ", input$epsilon, ")")
    )
  })
  
  # Render the overall optimal values output as a DT table
  output$optimal_values_table <- renderDT({
    firm_first_df <- firm_first_optimal_point()
    if (is.null(firm_first_df) || nrow(firm_first_df) == 0) {
      firm_first_df_for_combine <- data.frame(
        Scenario = "Firm First", Optimal_Lambda = NA, Optimal_W = NA,
        Firm_Profit = NA, User_Utility = NA, Welfare = NA,
        Uncond_Prob_Buying = NA, Optimal_Mutual_Info = NA,
        Optimal_a_star = NA, Optimal_b_star = NA, stringsAsFactors = FALSE
      )
    } else {
      firm_first_df_for_combine <- firm_first_df %>%
        mutate(Scenario = "Firm First", Welfare = user_utility + firm_profit) %>%
        select(Scenario, Optimal_Lambda = lambda, Optimal_W = w,
               Firm_Profit = firm_profit, User_Utility = user_utility, Welfare,
               Uncond_Prob_Buying = uncond_prob_buying, Optimal_Mutual_Info = optimal_mutual_info,
               Optimal_a_star = a_star_opt, Optimal_b_star = b_star_opt)
    }
    
    user_first_df_for_combine <- user_first_results() %>%
      select(Scenario, Optimal_Lambda, Optimal_W,
             Firm_Profit, User_Utility, Welfare,
             Uncond_Prob_Buying, Optimal_Mutual_Info,
             Optimal_a_star, Optimal_b_star) # Ensure these columns exist for consistency
    
    # Get social planner results if they've been calculated
    social_planner_df_raw <- social_planner_results() 
    
    if (is.null(social_planner_df_raw) || nrow(social_planner_df_raw) == 0) {
      # Placeholder for social planner if not calculated yet
      social_planner_df_for_combine <- data.frame(
        Scenario = "Social Planner", Optimal_Lambda = NA, Optimal_W = NA,
        Firm_Profit = NA, User_Utility = NA, Welfare = NA,
        Uncond_Prob_Buying = NA, Optimal_Mutual_Info = NA,
        Optimal_a_star = NA, Optimal_b_star = NA, stringsAsFactors = FALSE
      )
    } else {
      social_planner_df_for_combine <- social_planner_df_raw %>%
        select(Scenario, Optimal_Lambda, Optimal_W,
               Firm_Profit, User_Utility, Welfare,
               Uncond_Prob_Buying, Optimal_Mutual_Info,
               Optimal_a_star, Optimal_b_star)
    }
    
    combined_df <- bind_rows(firm_first_df_for_combine, user_first_df_for_combine, social_planner_df_for_combine) %>%
      mutate(across(where(is.numeric), ~ round(., 4)))
    
    datatable(combined_df, options = list(dom = 't', paging = FALSE, searching = FALSE),
              caption = paste0("Comparison of Optimal Outcomes (U_L = ", input$UL, ", U_H = ", input$UH, ", epsilon = ", input$epsilon, ")")
    )
  })
  
  # About section content
  output$about_text <- renderUI({
    HTML("
      <h4>About This Application</h4>
      <p>This R Shiny app visualizes different optimization scenarios for a firm and its users.</p>
      <p>The app calculates and displays:</p>
      <ol>
        <li>The firm's profit (objective function) value across a user-defined grid of $\\lambda$ and $w$ values for the <b>'Firm First'</b> scenario (Platform as leader).</li>
        <li>The user's utility surface for the <b>'User First'</b> scenario (User as leader), assuming the platform's best response is to set $\\lambda$ and $w$ to their maximum allowed values (based on the UI sliders).</li>
        <li>The optimal outcome for a <b>'Social Planner'</b>, who aims to maximize the sum of Firm Profit and User Utility simultaneously across $a, b, \\lambda,$ and $w$.</li>
        <li>A comparison of key metrics across all three scenarios in a table.</li>
      </ol>
      <h5>How to Use:</h5>
      <ul>
        <li>Adjust the core parameters $U_L$, $U_H$, and $\\epsilon$ using the sliders on the left.</li>
        <li>In the 'Firm First Plot' tab, set the desired range for $\\lambda$ and $w$. These ranges define the search space for the 'Firm First' global maximum and also set the upper bounds for $\\lambda$ and $w$ in the 'User First' and 'Social Planner' scenarios.</li>
        <li>Increase the 'General Grid Resolution' for more detailed surface plots (Firm First, User First). This will increase computation time.</li>
        <li>For the 'Social Planner' scenario, click the 'Calculate Social Planner Optimum' button to run the multi-parameter optimization.</li>
      </ul>
      <h5>Plot Interpretations:</h5>
      <ul>
        <li><b>'Firm First Plot':</b> Shows the firm's profit as a function of $\\lambda$ and $w$. Darker colors indicate higher profit. The <span style='color:red;'>&#10038;</span> red star marks the optimal $(\\lambda, w)$ for the firm.</li>
        <li><b>'User First Plot':</b> Shows the user's utility as a function of conditional buying probabilities $a$ and $b$, given fixed $\\lambda$ and $w$. Darker colors indicate higher utility. The <span style='color:red;'>&#10038;</span> red star marks the optimal $(a^*, b^*)$ for the user under these conditions.</li>
      </ul>
      <h5>Scenario Comparison Tab:</h5>
      <p>This tab displays a table summarizing the optimal outcomes for all three scenarios, including optimal parameters, profits, utilities, welfare, and buying probabilities.</p>
      <h5>Important Notes:</h5>
      <ul>
        <li>The 'Global Maximum on Grid' for 'Firm First' is the highest point found within the discrete grid. It might not be the exact global maximum if the grid resolution is too coarse or if the true maximum lies outside the defined lambda/w range.</li>
        <li>The 'User First' scenario's results are based on the specific assumption that the platform's best response is to always set $\\lambda$ and $w$ to their maximum feasible values, as described in the provided theoretical model.</li>
        <li>The 'Social Planner' optimization uses a numerical method ('L-BFGS-B'). While generally robust for bounded problems, it may find a local optimum rather than the global one, especially for complex, non-convex functions. The initial parameter guess is set to the mid-point of the parameter ranges.</li>
        <li>Numerical instabilities (e.g., very large/small numbers, division by zero, or invalid logarithm arguments) might occur for certain parameter combinations or at the edges of the plotting/optimization range. Points where calculations result in non-finite values (NA, Inf, NaN) are automatically excluded or penalized.</li>
      </ul>
    ")
  })
}

shinyApp(ui = ui, server = server)
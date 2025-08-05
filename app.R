library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# Define wide, fixed bounds for lambda and w for the entire app
# THESE ARE NOW GLOBAL AND ACCESSIBLE TO ALL FUNCTIONS
lambda_min_global <- 1e-9 # Changed from 0 to a small positive value to avoid issues when lambda=0
lambda_max_global <- 10  # Set to a sufficiently large value
w_min_global <- 0   # Allow for negative prices (subsidiies)
w_max_global <- 10   # Set to a sufficiently large value

# Helper functions

# Corrected and numerically stable function to calculate a_star
a_star_val <- function(lambda, w, UL, epsilon) {
  # Handle case where lambda is zero or very close to it
  if (lambda <= 1e-9) { 
    if (UL - w > epsilon) return(1)
    else return(0)
  }
  
  # Numerically stable calculation using log-sum-exp trick
  exponent <- (epsilon - (UL - w)) / lambda
  result <- 1 / (1 + exp(exponent))
  
  # Ensure the result is within [0, 1] for robustness
  return(pmax(0, pmin(1, result)))
}

# Corrected and numerically stable function to calculate b_star
b_star_val <- function(lambda, w, UH, epsilon) {
  # Handle case where lambda is zero or very close to it
  if (lambda <= 1e-9) { 
    if (UH - w > epsilon) return(1)
    else return(0)
  }
  
  # Numerically stable calculation using log-sum-exp trick
  exponent <- (epsilon - (UH - w)) / lambda
  result <- 1 / (1 + exp(exponent))
  
  # Ensure the result is within [0, 1] for robustness
  return(pmax(0, pmin(1, result)))
}
# Robust and consistent mutual information calculation
calculate_mutual_information <- function(a, b) {
  # Clip a, b to ensure they are within [0, 1] for robustness
  a <- pmax(0, pmin(1, a))
  b <- pmax(0, pmin(1, b))
  
  s <- a + b
  
  # Adjust s and 2-s for denominators if they are exactly zero
  # Using a small positive epsilon to prevent division by zero or log(0)
  s_denom <- ifelse(s == 0, 1e-9, s) # Use 1e-9 for non-zero denominator
  two_minus_s_denom <- ifelse( (2-s) == 0, 1e-9, 2-s) # Use 1e-9 for non-zero denominator
  
  # Helper for individual log terms: prob * log(ratio)
  log_term_safe <- function(prob, numerator, denominator) {
    # CRITICAL FIX: If probability is zero, the term is always zero.
    # This prevents 0 * log(0) or 0 * log(negative) leading to NaN/Inf.
    if (prob == 0) return(0)
    
    # Defensive check: if denominator is problematic (though s_denom/two_minus_s_denom should prevent this)
    if (denominator <= 0) return(NA)
    
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

# Mutual information I(s, \hat{s}) (equation 9)
Optimal_Mutual_Information <- function(lambda, w, UL, UH, epsilon) {
  a_s <- a_star_val(lambda, w, UL, epsilon)
  b_s <- b_star_val(lambda, w, UH, epsilon)
  
  if (is.na(a_s) || is.na(b_s)) return(NA)
  
  return(calculate_mutual_information(a_s, b_s))
}

# User Utility Function (equation 8 with a_star, b_star substituted)
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
# This function calculates firm profit given lambda and w, assuming user's best response
# Assuming calculate_mutual_information() is available and correct
Firm_Profit_Function <- function(lambda, w, UL, UH, epsilon) {
  a_s <- a_star_val(lambda, w, UL, epsilon)
  b_s <- b_star_val(lambda, w, UH, epsilon)
  
  if (is.na(a_s) || is.na(b_s)) return(NA)
  
  # Calculate the unconditional probability of buying
  prob_buying <- 0.5 * a_s + 0.5 * b_s
  
  # Calculate the mutual information term
  mutual_info <- calculate_mutual_information(a_s, b_s)
  
  if (is.na(prob_buying) || is.na(mutual_info)) return(NA)
  
  # Correct profit calculation: Revenue from transactions + Revenue from attention
  profit <- (prob_buying * w) + (lambda * mutual_info)
  
  if (is.infinite(profit) || is.na(profit)) return(NA)
  
  return(profit)
}

# User Utility Function given direct a,b inputs (used for social planner and user first plot)
User_Utility_Given_ab_lambda_w <- function(a, b, lambda, w, UL, UH, epsilon) {
  # Clip a, b to be within [0,1]
  a <- pmax(0, pmin(1, a)) 
  b <- pmax(0, pmin(1, b))
  
  if (is.infinite(lambda) || is.na(lambda)) return(NA)
  
  term1_part <- 0.5 * (a * (UL - w) + (1 - a) * epsilon) +
    0.5 * (b * (UH - w) + (1 - b) * epsilon)
  
  mutual_info_calculated <- calculate_mutual_information(a, b) 
  
  if (is.na(mutual_info_calculated)) return(NA)
  
  result <- term1_part - lambda * mutual_info_calculated
  
  if (is.na(result) || is.infinite(result)) return(NA)
  return(result)
}

# Firm's Best Response (lambda, w) to user's choices (a, b)
Firm_Best_Response_Lambda_W <- function(a, b, UL, UH, epsilon) {
  # Clip a, b to ensure they are within [0, 1] for robust calculation of MI
  a <- pmax(0, pmin(1, a))
  b <- pmax(0, pmin(1, b))
  
  C1_w_coeff <- 0.5 * (a + b)
  C2_lambda_coeff <- calculate_mutual_information(a, b)
  
  if (is.na(C2_lambda_coeff)) {
    return(list(lambda = NA, w = NA, profit = NA, explanation = "Mutual Information calculation resulted in NA."))
  }
  
  optimal_w <- NA
  optimal_lambda <- NA
  w_explanation <- ""
  lambda_explanation <- ""
  
  if (C1_w_coeff > 0) {
    optimal_w <- w_max_global
    w_explanation <- paste0("Coefficient of w (0.5*(a+b)) is positive (", round(C1_w_coeff, 3), "), so w is set to w_max_global.")
  } else if (C1_w_coeff < 0) {
    optimal_w <- w_min_global
    w_explanation <- paste0("Coefficient of w (0.5*(a+b)) is negative (", round(C1_w_coeff, 3), "), so w is set to w_min_global.")
  } else { # C1_w_coeff is 0
    optimal_w <- w_min_global # Arbitrarily pick min if indifferent
    w_explanation <- paste0("Coefficient of w (0.5*(a+b)) is zero, so w does not affect profit. Arbitrarily chosen as w_min_global (", w_min_global, ").")
  }
  
  if (C2_lambda_coeff > 0) {
    optimal_lambda <- lambda_max_global
    lambda_explanation <- paste0("Coefficient of lambda (MI) is positive (", round(C2_lambda_coeff, 3), "), so lambda is set to lambda_max_global.")
  } else if (C2_lambda_coeff < 0) {
    optimal_lambda <- lambda_min_global
    lambda_explanation <- paste0("Coefficient of lambda (MI) is negative (", round(C2_lambda_coeff, 3), "), so lambda is set to lambda_min_global.")
  } else { # C2_lambda_coeff is 0
    optimal_lambda <- lambda_min_global # Arbitrarily pick min if indifferent
    lambda_explanation <- paste0("Coefficient of lambda (MI) is zero, so lambda does not affect profit. Arbitrarily chosen as lambda_min_global (", lambda_min_global, ").")
  }
  
  profit <- C1_w_coeff * optimal_w + optimal_lambda * C2_lambda_coeff
  
  return(list(
    lambda = optimal_lambda,
    w = optimal_w,
    profit = profit,
    explanation = paste(w_explanation, lambda_explanation, sep = "\n")
  ))
}


# Social Welfare Objective Function for optim (First-Best)
Social_Welfare_Objective <- function(params, UL, UH, epsilon) {
  a <- params[1]
  b <- params[2]
  lambda <- params[3]
  w <- params[4]
  
  if (lambda < lambda_min_global) return(Inf)
  
  a_bounded <- pmax(0, pmin(1, a))
  b_bounded <- pmax(0, pmin(1, b))
  
  mi_val_sp <- calculate_mutual_information(a_bounded, b_bounded)
  if (is.na(mi_val_sp)) return(Inf)
  
  firm_profit_val <- 0.5 * (a_bounded + b_bounded) * w + lambda * mi_val_sp
  user_utility_val <- User_Utility_Given_ab_lambda_w(a_bounded, b_bounded, lambda, w, UL, UH, epsilon)
  
  if (is.na(firm_profit_val) || is.infinite(firm_profit_val) ||
      is.na(user_utility_val) || is.infinite(user_utility_val)) {
    return(Inf)
  }
  
  welfare <- firm_profit_val + user_utility_val
  
  return(-welfare) # optim minimizes
}

# Objective function to minimize the distance to a mutual best response (Nash Equilibrium)
Nash_Equilibrium_Objective <- function(params, UL, UH, epsilon) {
  a <- params[1]
  b <- params[2]
  lambda <- params[3]
  w <- params[4]
  
  a_br <- a_star_val(lambda, w, UL, epsilon)
  b_br <- b_star_val(lambda, w, UH, epsilon)
  
  # Call firm best response function and extract values
  firm_br_result <- Firm_Best_Response_Lambda_W(a, b, UL, UH, epsilon)
  lambda_br <- firm_br_result$lambda
  w_br <- firm_br_result$w
  
  if (is.na(a_br) || is.na(b_br) || is.na(lambda_br) || is.na(w_br)) {
    return(Inf)
  }
  
  distance <- (a - a_br)^2 + (b - b_br)^2 +
    (lambda - lambda_br)^2 + (w - w_br)^2
  
  return(distance)
}

# --- NEW OBJECTIVE FOR SECOND-BEST SOCIAL PLANNER ---
# Objective function to maximize welfare by controlling only lambda and w.
# User's behavior (a,b) is a best response to lambda and w.
Second_Best_Welfare_Objective <- function(params, UL, UH, epsilon) {
  lambda <- params[1]
  w <- params[2]
  
  # Constraints check: must be within global bounds
  if (lambda < lambda_min_global || w < w_min_global || lambda > lambda_max_global || w > w_max_global) {
    return(Inf)
  }
  
  # Calculate user's best responses to the given lambda and w
  a_s <- a_star_val(lambda, w, UL, epsilon)
  b_s <- b_star_val(lambda, w, UH, epsilon)
  
  if (is.na(a_s) || is.na(b_s)) {
    return(Inf) # If user's best response is not well-defined, this is a bad choice for the planner
  }
  
  # Calculate firm's profit given these parameters
  firm_profit <- Firm_Profit_Function(lambda, w, UL, UH, epsilon)
  # Calculate user's utility given these parameters
  user_utility <- User_Utility_Function(lambda, w, UL, UH, epsilon)
  
  if (is.na(firm_profit) || is.infinite(firm_profit) ||
      is.na(user_utility) || is.infinite(user_utility)) {
    return(Inf)
  }
  
  welfare <- firm_profit + user_utility
  
  return(-welfare) # optim minimizes, so return negative welfare for maximization
}


# --- Shiny UI ---
ui <- fluidPage(
  titlePanel("Outcomes for binary signal"),
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
      tabsetPanel(id = "main_tabs",
                  tabPanel("Firm First Plot",
                           plotOutput("profit_plot", height = "600px")
                  ),
                  tabPanel("User First Plot",
                           plotOutput("user_utility_plot", height = "600px")
                  ),
                  tabPanel("Scenario Comparison",
                           tableOutput("optimal_values_table"),
                           hr(),
                           fluidRow(
                             column(12,
                                    h4("Select Parameters to Plot:"),
                                    checkboxGroupInput("params_to_plot", label = NULL,
                                                       choices = c(
                                                         "Optimal_Lambda", "Optimal_W", "Firm_Profit",
                                                         "User_Utility", "Welfare", "Uncond_Prob_Buying",
                                                         "Optimal_Mutual_Info", "Optimal_a_star", "Optimal_b_star"
                                                       ),
                                                       inline = TRUE),
                                    uiOutput("param_plots") # This will be dynamic plot outputs
                             )
                           )
                  ),
                  tabPanel("About", uiOutput("about_text"))
      )
    )
  )
)

# --- Shiny Server ---
server <- function(input, output) {
  
  # Reactive expression to generate the grid and calculate Firm Profit values (Firm First)
  # Firm as Leader: Firm chooses lambda and w to maximize profit, user responds with a_star, b_star
  grid_data <- reactive({
    req(input$grid_res) # Only depends on grid_res now
    
    lambda_vals <- seq(lambda_min_global, lambda_max_global, length.out = input$grid_res)
    w_vals <- seq(w_min_global, w_max_global, length.out = input$grid_res)
    
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
  
  # Reactive expression for User First plot data (user utility surface in a-b space)
  # User as Leader: User chooses a,b. Firm responds with optimal lambda, w. User anticipates this.
  user_first_plot_data <- reactive({
    req(input$UL, input$UH, input$epsilon)
    
    user_plot_grid_res <- input$grid_res # Use general grid res for consistency
    a_vals <- seq(0, 1, length.out = user_plot_grid_res)
    b_vals <- seq(0, 1, length.out = user_plot_grid_res)
    
    plot_df_user <- expand.grid(a = a_vals, b = b_vals)
    
    results_df <- plot_df_user %>%
      rowwise() %>%
      mutate(
        # Call Firm_Best_Response_Lambda_W without bounds arguments (it uses global fixed bounds)
        firm_br = list(Firm_Best_Response_Lambda_W(a, b, input$UL, input$UH, input$epsilon)),
        responsive_lambda = firm_br$lambda,
        responsive_w = firm_br$w,
        firm_profit_br = firm_br$profit,
        user_utility_val = User_Utility_Given_ab_lambda_w(a, b, responsive_lambda, responsive_w, input$UL, input$UH, input$epsilon)
      ) %>%
      ungroup() %>%
      select(-firm_br) # Remove the list column
    
    return(results_df)
  })
  
  # Reactive expression to find the optimal point for the User First scenario
  user_first_optimal_point <- reactive({
    df <- user_first_plot_data() %>% drop_na(user_utility_val)
    if (nrow(df) == 0) return(NULL)
    
    max_row <- df %>% slice_max(user_utility_val, n = 1, with_ties = FALSE)
    
    if (nrow(max_row) > 0) {
      optimal_a_uf <- max_row$a
      optimal_b_uf <- max_row$b
      optimal_lambda_uf <- max_row$responsive_lambda
      optimal_w_uf <- max_row$responsive_w
      
      max_row$Scenario <- "User First"
      max_row$Optimal_Lambda <- optimal_lambda_uf
      max_row$Optimal_W <- optimal_w_uf
      max_row$Firm_Profit <- max_row$firm_profit_br # Firm profit at this BR
      max_row$User_Utility <- max_row$user_utility_val # User utility at their optimal a,b
      max_row$Welfare <- max_row$Firm_Profit + max_row$User_Utility
      max_row$Uncond_Prob_Buying <- 0.5 * optimal_a_uf + 0.5 * optimal_b_uf
      max_row$Optimal_Mutual_Info <- calculate_mutual_information(optimal_a_uf, optimal_b_uf)
      max_row$Optimal_a_star <- optimal_a_uf
      max_row$Optimal_b_star <- optimal_b_uf
      
      # Select relevant columns for consistency with other scenarios
      max_row <- max_row %>%
        select(Scenario, Optimal_Lambda, Optimal_W, Firm_Profit, User_Utility, Welfare,
               Uncond_Prob_Buying, Optimal_Mutual_Info, Optimal_a_star, Optimal_b_star)
    }
    max_row
  })
  
  # Reactive for Social Planner (First-Best) results
  social_planner_results <- reactive({
    req(input$UL, input$UH, input$epsilon)
    
    lower_bounds <- c(a = 0, b = 0, lambda = lambda_min_global, w = w_min_global)
    upper_bounds <- c(a = 1, b = 1, lambda = lambda_max_global, w = w_max_global)
    
    initial_params <- c(
      a = 0.5,
      b = 0.5,
      lambda = (lambda_min_global + lambda_max_global) / 2,
      w = (w_min_global + w_max_global) / 2
    )
    
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
        control = list(factr = 1e7)
      )
    })
    
    optimal_a_sp <- optim_result$par[1]
    optimal_b_sp <- optim_result$par[2]
    optimal_lambda_sp <- optim_result$par[3]
    optimal_w_sp <- optim_result$par[4]
    
    sp_firm_profit <- (0.5 * optimal_a_sp + 0.5 * optimal_b_sp) * optimal_w_sp + optimal_lambda_sp * calculate_mutual_information(optimal_a_sp, optimal_b_sp)
    sp_user_utility <- User_Utility_Given_ab_lambda_w(optimal_a_sp, optimal_b_sp, optimal_lambda_sp, optimal_w_sp, input$UL, input$UH, input$epsilon)
    sp_welfare <- sp_firm_profit + sp_user_utility
    sp_uncond_prob_buying <- 0.5 * optimal_a_sp + 0.5 * optimal_b_sp
    sp_mutual_info <- calculate_mutual_information(optimal_a_sp, optimal_b_sp)
    
    data.frame(
      Scenario = "Social Planner (First-Best)",
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
  
  # Reactive for the Nash Equilibrium results
  nash_equilibrium_results <- reactive({
    req(input$UL, input$UH, input$epsilon)
    
    lower_bounds <- c(a = 0, b = 0, lambda = lambda_min_global, w = w_min_global)
    upper_bounds <- c(a = 1, b = 1, lambda = lambda_max_global, w = w_max_global)
    
    initial_params <- c(
      a = 0.5,
      b = 0.5,
      lambda = (lambda_min_global + lambda_max_global) / 2,
      w = (w_min_global + w_max_global) / 2
    )
    
    withProgress(message = 'Calculating Nash Equilibrium...', value = 0.5, {
      optim_result <- optim(
        par = initial_params,
        fn = Nash_Equilibrium_Objective,
        method = "L-BFGS-B",
        lower = lower_bounds,
        upper = upper_bounds,
        UL = input$UL,
        UH = input$UH,
        epsilon = input$epsilon,
        control = list(factr = 1e7)
      )
    })
    
    optimal_a_ne <- optim_result$par[1]
    optimal_b_ne <- optim_result$par[2]
    optimal_lambda_ne <- optim_result$par[3]
    optimal_w_ne <- optim_result$par[4]
    
    ne_firm_profit <- Firm_Profit_Function(optimal_lambda_ne, optimal_w_ne, input$UL, input$UH, input$epsilon)
    ne_user_utility <- User_Utility_Function(optimal_lambda_ne, optimal_w_ne, input$UL, input$UH, input$epsilon)
    ne_welfare <- ne_firm_profit + ne_user_utility
    ne_uncond_prob_buying <- Unconditional_Prob_Buying(optimal_lambda_ne, optimal_w_ne, input$UL, input$UH, input$epsilon)
    ne_mutual_info <- Optimal_Mutual_Information(optimal_lambda_ne, optimal_w_ne, input$UL, input$UH, input$epsilon)
    
    data.frame(
      Scenario = "Nash Equilibrium",
      Optimal_Lambda = optimal_lambda_ne,
      Optimal_W = optimal_w_ne,
      Firm_Profit = ne_firm_profit,
      User_Utility = ne_user_utility,
      Welfare = ne_welfare,
      Uncond_Prob_Buying = ne_uncond_prob_buying,
      Optimal_Mutual_Info = ne_mutual_info,
      Optimal_a_star = optimal_a_ne,
      Optimal_b_star = optimal_b_ne,
      stringsAsFactors = FALSE
    )
  })
  
  # --- NEW REACTIVE FOR SECOND-BEST SOCIAL PLANNER ---
  second_best_social_planner_results <- reactive({
    req(input$UL, input$UH, input$epsilon)
    
    # Define bounds for optimization. Only lambda and w are controlled.
    lower_bounds <- c(lambda = lambda_min_global, w = w_min_global)
    upper_bounds <- c(lambda = lambda_max_global, w = w_max_global)
    
    # Initial guess
    initial_params <- c(
      lambda = (lambda_min_global + lambda_max_global) / 2,
      w = (w_min_global + w_max_global) / 2
    )
    
    withProgress(message = 'Calculating Second-Best Social Planner Optimum...', value = 0.5, {
      optim_result <- optim(
        par = initial_params,
        fn = Second_Best_Welfare_Objective,
        method = "L-BFGS-B",
        lower = lower_bounds,
        upper = upper_bounds,
        UL = input$UL,
        UH = input$UH,
        epsilon = input$epsilon,
        control = list(factr = 1e7)
      )
    })
    
    # Extract results
    optimal_lambda_sb <- optim_result$par[1]
    optimal_w_sb <- optim_result$par[2]
    
    # Calculate the resulting user choices using the *correct* function calls
    # a_star requires UL, b_star requires UH
    optimal_a_sb <- a_star_val(optimal_lambda_sb, optimal_w_sb, input$UL, input$epsilon)
    optimal_b_sb <- b_star_val(optimal_lambda_sb, optimal_w_sb, input$UH, input$epsilon)
    
    # Calculate all metrics at the found optimal point
    sb_firm_profit <- Firm_Profit_Function(optimal_lambda_sb, optimal_w_sb, input$UL, input$UH, input$epsilon)
    sb_user_utility <- User_Utility_Function(optimal_lambda_sb, optimal_w_sb, input$UL, input$UH, input$epsilon)
    sb_welfare <- sb_firm_profit + sb_user_utility
    sb_uncond_prob_buying <- Unconditional_Prob_Buying(optimal_lambda_sb, optimal_w_sb, input$UL, input$UH, input$epsilon)
    sb_mutual_info <- Optimal_Mutual_Information(optimal_lambda_sb, optimal_w_sb, input$UL, input$UH, input$epsilon)
    
    data.frame(
      Scenario = "Social Planner (Second-Best)",
      Optimal_Lambda = optimal_lambda_sb,
      Optimal_W = optimal_w_sb,
      Firm_Profit = sb_firm_profit,
      User_Utility = sb_user_utility,
      Welfare = sb_welfare,
      Uncond_Prob_Buying = sb_uncond_prob_buying,
      Optimal_Mutual_Info = sb_mutual_info,
      Optimal_a_star = optimal_a_sb,
      Optimal_b_star = optimal_b_sb,
      stringsAsFactors = FALSE
    )
  })
  
  # Reactive for the combined data table and plots
  optimal_values_data <- reactive({
    
    # Firm First
    firm_first_optimal <- firm_first_optimal_point()
    firm_first_df <- if (is.null(firm_first_optimal) || nrow(firm_first_optimal) == 0) {
      data.frame(Scenario = "Firm First", Optimal_Lambda = NA, Optimal_W = NA,
                 Firm_Profit = NA, User_Utility = NA, Welfare = NA,
                 Uncond_Prob_Buying = NA, Optimal_Mutual_Info = NA,
                 Optimal_a_star = NA, Optimal_b_star = NA, stringsAsFactors = FALSE)
    } else {
      data.frame(
        Scenario = "Firm First",
        Optimal_Lambda = firm_first_optimal$lambda,
        Optimal_W = firm_first_optimal$w,
        Firm_Profit = firm_first_optimal$firm_profit,
        User_Utility = firm_first_optimal$user_utility,
        Welfare = firm_first_optimal$firm_profit + firm_first_optimal$user_utility,
        Uncond_Prob_Buying = firm_first_optimal$uncond_prob_buying,
        Optimal_Mutual_Info = firm_first_optimal$optimal_mutual_info,
        Optimal_a_star = firm_first_optimal$a_star_opt,
        Optimal_b_star = firm_first_optimal$b_star_opt,
        stringsAsFactors = FALSE
      )
    }
    
    # User First
    user_first_optimal <- user_first_optimal_point()
    user_first_df <- if (is.null(user_first_optimal) || nrow(user_first_optimal) == 0) {
      data.frame(Scenario = "User First", Optimal_Lambda = NA, Optimal_W = NA,
                 Firm_Profit = NA, User_Utility = NA, Welfare = NA,
                 Uncond_Prob_Buying = NA, Optimal_Mutual_Info = NA,
                 Optimal_a_star = NA, Optimal_b_star = NA, stringsAsFactors = FALSE)
    } else {
      user_first_optimal
    }
    
    # Social Planner (First-Best)
    social_planner_optimal <- social_planner_results()
    social_planner_df <- if (is.null(social_planner_optimal) || nrow(social_planner_optimal) == 0) {
      data.frame(Scenario = "Social Planner (First-Best)", Optimal_Lambda = NA, Optimal_W = NA,
                 Firm_Profit = NA, User_Utility = NA, Welfare = NA,
                 Uncond_Prob_Buying = NA, Optimal_Mutual_Info = NA,
                 Optimal_a_star = NA, Optimal_b_star = NA, stringsAsFactors = FALSE)
    } else {
      social_planner_optimal
    }
    
    # Nash Equilibrium
    nash_eq_optimal <- nash_equilibrium_results()
    nash_eq_df <- if (is.null(nash_eq_optimal) || nrow(nash_eq_optimal) == 0) {
      data.frame(Scenario = "Nash Equilibrium", Optimal_Lambda = NA, Optimal_W = NA,
                 Firm_Profit = NA, User_Utility = NA, Welfare = NA,
                 Uncond_Prob_Buying = NA, Optimal_Mutual_Info = NA,
                 Optimal_a_star = NA, Optimal_b_star = NA, stringsAsFactors = FALSE)
    } else {
      nash_eq_optimal
    }
    
    # Second-Best Social Planner
    sb_social_planner_optimal <- second_best_social_planner_results()
    sb_social_planner_df <- if (is.null(sb_social_planner_optimal) || nrow(sb_social_planner_optimal) == 0) {
      data.frame(Scenario = "Social Planner (Second-Best)", Optimal_Lambda = NA, Optimal_W = NA,
                 Firm_Profit = NA, User_Utility = NA, Welfare = NA,
                 Uncond_Prob_Buying = NA, Optimal_Mutual_Info = NA,
                 Optimal_a_star = NA, Optimal_b_star = NA, stringsAsFactors = FALSE)
    } else {
      sb_social_planner_optimal
    }
    
    # Combine all five data frames
    combined_df <- bind_rows(firm_first_df, user_first_df, social_planner_df, nash_eq_df, sb_social_planner_df) %>%
      mutate(across(where(is.numeric), ~ round(., 2)))
    
    return(combined_df)
  })
  
  # Render the plot for Firm First
  output$profit_plot <- renderPlot({
    df_clean <- grid_data() %>% drop_na(firm_profit)
    max_point <- firm_first_optimal_point()
    
    p <- ggplot(df_clean, aes(x = lambda, y = w)) +
      geom_contour_filled(aes(z = firm_profit), alpha = 0.7) +
      labs(title = "Firm First Optimal Point (Firm as Leader, User as Follower)",
           subtitle = paste0("U_L = ", input$UL, ", U_H = ", input$UH, ", epsilon = ", input$epsilon),
           x = "Lambda (λ)",
           y = "W") +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) +
      coord_cartesian(xlim = c(lambda_min_global, lambda_max_global), ylim = c(w_min_global, w_max_global))
    
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
    uf_opt_point <- user_first_optimal_point()
    
    p <- ggplot(df_user_clean, aes(x = a, y = b)) +
      geom_contour_filled(aes(z = user_utility_val), alpha = 0.7) +
      labs(title = "User First Optimal Point (User as Leader, Firm as Follower)",
           subtitle = paste0("Firm chooses λ, w as best response to User's a,b. U_L = ", input$UL, ", U_H = ", input$UH, ", epsilon = ", input$epsilon),
           x = "User's Chosen a",
           y = "User's Chosen b") +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
    
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
  
  # Render the overall optimal values output as a DT table
  output$optimal_values_table <- renderTable({
    (optimal_values_data())
    
  })
  
  # Dynamic UI for parameter plots
  output$param_plots <- renderUI({
    req(input$params_to_plot) # Ensure at least one checkbox is ticked
    plot_list <- lapply(input$params_to_plot, function(param_name) {
      plot_id <- paste0("plot_", param_name)
      plotOutput(plot_id, height = "300px") # Adjust height as needed
    })
    do.call(tagList, plot_list)
  })
  
  # Define renderPlot for each possible parameter
  all_plot_params <- c(
    "Optimal_Lambda", "Optimal_W", "Firm_Profit", "User_Utility", "Welfare",
    "Uncond_Prob_Buying", "Optimal_Mutual_Info", "Optimal_a_star", "Optimal_b_star"
  )
  
  for (param in all_plot_params) {
    local({
      my_param <- param # Localize the parameter name for each iteration
      output_id <- paste0("plot_", my_param)
      
      output[[output_id]] <- renderPlot({
        req(my_param %in% input$params_to_plot) # Only render if checkbox is ticked
        current_data <- optimal_values_data()
        req(current_data)
        if (!(my_param %in% colnames(current_data))) {
          return(NULL) # Or show a message if column is missing
        }
        
        # Get appropriate y-axis label
        y_label <- switch(my_param,
                          "Optimal_Lambda" = "Optimal Lambda (λ)",
                          "Optimal_W" = "Optimal W",
                          "Firm_Profit" = "Firm Profit",
                          "User_Utility" = "User Utility",
                          "Welfare" = "Total Welfare",
                          "Uncond_Prob_Buying" = "Unconditional Probability of Buying",
                          "Optimal_Mutual_Info" = "Optimal Mutual Information",
                          "Optimal_a_star" = "User's Optimal a*",
                          "Optimal_b_star" = "User's Optimal b*",
                          gsub("_", " ", my_param))
        
        ggplot(current_data, aes_string(x = "Scenario", y = my_param, fill = "Scenario")) +
          geom_col(position = position_dodge(width = 0.9), color = "black") +
          labs(title = paste0("Comparison of ", y_label),
               y = y_label,
               x = "") +
          theme_minimal() +
          theme(legend.position = "none",
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                axis.title.y = element_text(size = 12)) +
          geom_text(aes_string(label = my_param), vjust = -0.5, size = 3)
      })
    })
  }
  
  # About section content
  output$about_text <- renderUI({
    HTML("
      <h4>About This Application</h4>
      <p>This R Shiny app visualizes different optimization scenarios for a firm and its users.</p>
      <p>The app calculates and displays:</p>
      <ol>
        <li>The firm's profit (objective function) value across a predefined grid of $\\lambda$ and $w$ values for the <b>'Firm First'</b> scenario (Platform as leader). In this scenario, the firm chooses $\\lambda$ and $w$ to maximize its profit, anticipating the user's optimal buying probabilities ($a^*$ and $b^*$).</li>
        <li>The user's utility surface for the <b>'User First'</b> scenario (User as leader). Here, the user chooses their desired conditional buying probabilities ($a$ and $b$). The firm, acting as a follower, observes these choices and optimally sets its $\\lambda$ and $w$ to maximize its profit, such that the user's actual best responses match the chosen $a$ and $b$. The user anticipates this firm's best response when making their initial choice of $a$ and $b$.</li>
        <li>The optimal outcome for a <b>'Social Planner (First-Best)'</b>, who aims to maximize the sum of Firm Profit and User Utility simultaneously across $a, b, \\lambda,$ and $w$. This is a theoretical benchmark.</li>
        <li>The optimal outcome for a <b>'Social Planner (Second-Best)'</b>, who maximizes total welfare by choosing only $\\lambda$ and $w$ and anticipating the user's best response behavior.</li>
        <li>The <b>'Nash Equilibrium'</b>, where both the firm's choices and the user's choices are mutual best responses to each other.</li>
        <li>A comparison of key metrics across all scenarios in a table.</li>
      </ol>
      <h5>How to Use:</h5>
      <ul>
        <li>Adjust the core parameters $U_L$, $U_H$, and $\\epsilon$ using the sliders on the left.</li>
        <li>The ranges for $\\lambda$ (0 to 10) and $w$ (0 to 10) are fixed internally and define the search space for all optimizations.</li>
        <li>Increase the 'General Grid Resolution' for more detailed surface plots. This will increase computation time.</li>
      </ul>
      <h5>Important Notes:</h5>
      <ul>
        <li>The 'Social Planner (First-Best)', 'Nash Equilibrium', and 'Social Planner (Second-Best)' scenarios all rely on numerical optimization ('L-BFGS-B'). While generally robust, it may find a local optimum rather than the global one, especially for complex functions.</li>
        <li>Numerical instabilities (e.g., very large/small numbers, division by zero, or invalid logarithm arguments) might occur for certain parameter combinations. Points where calculations result in non-finite values (NA, Inf, NaN) are automatically excluded or penalized.</li>
      </ul>
    ")
  })
}

shinyApp(ui = ui, server = server)
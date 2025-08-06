library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT) # For a nice-looking summary table

# User-defined parameters (global for simplicity in this app)
UH <- 10      # utility of buying when X = H
UL <- 0       # utility of buying when X = L
epsilon <- 8  # utility of leaving
p <- 0.5      # prior
eps <- 10^(-3) # Small epsilon for numerical stability
lambda_range <- c(0,0.3)
w_range <- c(0,0.1)


# Function: probability of buying (P_buy)
P_buy <- function(r, w) {
  # Posterior when signal = H (High)
  prob_Y_H <- r * p + (1 - r) * (1 - p)
  posterior_H_given_Y_H <- (r * p) / prob_Y_H
  
  # Expected utility from action H given signal H (Principal buys)
  EU_buy_given_signal_H <- posterior_H_given_Y_H * (UH - w) + (1 - posterior_H_given_Y_H) * (UL - w)
  
  # Posterior when signal = L (Low)
  prob_Y_L <- (1 - r) * p + r * (1 - p)
  posterior_H_given_Y_L <- ((1 - r) * p) / prob_Y_L
  
  # Expected utility from action H given signal L (Principal buys)
  EU_buy_given_signal_L <- posterior_H_given_Y_L * (UH - w) + (1 - posterior_H_given_Y_L) * (UL - w)
  
  # Determine if Principal buys based on signal H or L
  # If EU of buying is >= epsilon (utility of leaving), principal buys.
  # Otherwise, principal leaves.
  
  # P_buy when signal is H
  buy_if_signal_H <- ifelse(EU_buy_given_signal_H >= epsilon, 1, 0)
  
  # P_buy when signal is L
  buy_if_signal_L <- ifelse(EU_buy_given_signal_L >= epsilon, 1, 0)
  
  # Overall probability of buying, considering the probability of each signal
  # P(Buy) = P(Buy|Signal=H)P(Signal=H) + P(Buy|Signal=L)P(Signal=L)
  # P(Signal=H) = prob_Y_H
  # P(Signal=L) = prob_Y_L
  
  p_buy <- buy_if_signal_H * prob_Y_H + buy_if_signal_L * prob_Y_L
  
  return(p_buy)
}

# Function to calculate entropy of a Bernoulli variable (helper for mutual_info_bernoulli)
entropy <- function(q) {
  # Handle log(0) case for q=0 or q=1, where entropy is 0
  if (q == 0 || q == 1) {
    return(0)
  }
  return(-q * log2(q) - (1 - q) * log2(1 - q))
}

# Function to compute mutual information for Bernoulli variables
mutual_info_bernoulli <- function(p, r) {
  # Calculate the probability of Y=1 (signal H)
  p_y_H <- r * p + (1 - r) * (1 - p)
  
  # Calculate the entropy of Y
  h_y <- entropy(p_y_H)
  
  # Calculate the conditional entropy H(Y|X) which is the entropy of the channel noise
  # This is the entropy of a Bernoulli variable with probability r (precision)
  h_y_given_x <- entropy(r)
  
  # Calculate mutual information
  mutual_info <- h_y - h_y_given_x
  
  return(mutual_info)
}

# Total profit of the agent
total_profit <- function(r, lambda, w) {
  # If r=0.5, mutual information is 0, so agent's profit is w * P_buy(0.5, w).
  # If r != 0.5 and lambda is infinite, profit will be infinite.
  
  if (r == 0.5) { 
    profit <- w * P_buy(r, w) 
  } else {
    profit <- w * P_buy(r, w) + lambda * mutual_info_bernoulli(p, r)
  }
  
  return(profit)
}


# Total utility of the principal
total_utility <- function(r, lambda, w) {
  # Check for the condition UH - epsilon < w
  if (UH - epsilon < w) {
    EU_total <- epsilon # Principal always chooses to leave, gets epsilon
  } else {
    # Posterior probability when signal = H (High)
    prob_Y_H <- r * p + (1 - r) * (1 - p)
    posterior_H_given_Y_H <- ifelse(prob_Y_H == 0, 0, (r * p) / prob_Y_H)
    
    # Posterior probability when signal = L (Low)
    prob_Y_L <- (1 - r) * p + r * (1 - p)
    posterior_H_given_Y_L <- ifelse(prob_Y_L == 0, 0, ((1 - r) * p) / prob_Y_L)
    
    # Expected utility from action H (High) given signal H
    EU_H_given_signal_H <- posterior_H_given_Y_H * (UH - w) + (1 - posterior_H_given_Y_H) * (UL - w)
    
    # Expected utility from action L (Low) given signal H (Action L is to leave, utility epsilon)
    EU_L_given_signal_H <- epsilon
    
    # Expected utility from action H given signal L
    EU_H_given_signal_L <- posterior_H_given_Y_L * (UH - w) + (1 - posterior_H_given_Y_L) * (UL - w)
    
    # Expected utility from action L given signal L (Action L is to leave, utility epsilon)
    EU_L_given_signal_L <- epsilon
    
    # Optimal action given signal H
    EU_optimal_given_signal_H <- max(EU_H_given_signal_H, EU_L_given_signal_H)
    
    # Optimal action given signal L
    EU_optimal_given_signal_L <- max(EU_H_given_signal_L, EU_L_given_signal_L)
    
    # Expected utility overall, considering the probability of each signal
    EU_total <- (prob_Y_H * EU_optimal_given_signal_H) + 
      (prob_Y_L * EU_optimal_given_signal_L) - 
      lambda * mutual_info_bernoulli(p, r)
  }
  
  return(EU_total)
}

# Social Welfare Function
social_welfare <- function(r, lambda, w) {
  total_utility(r, lambda, w) + total_profit(r, lambda, w)
}



# UI
ui <- fluidPage(
  titlePanel("Optimal Signal Design"),
  sidebarLayout(
    sidebarPanel(
      h4("Grid and Parameter Controls"),
   
      numericInput("grid_size", "Grid Points per Axis",
                   value = 50, min = 10, max = 500, step = 10),
      helpText("Note: Increasing the grid size may slow down calculations.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary",
                 tableOutput("summary_table")
        ),
        tabPanel("Principal's Best Response",
                 plotOutput("r_star_plot", height = "400px"),
                 plotOutput("profit_plot", height = "400px")
        ),
        tabPanel("Agent's Best Response",
                 plotOutput("optimalLambdaWPlot", height = "400px"),
                 plotOutput("totalUtilityPlot", height = "400px")
        ),
        tabPanel("Nash Equilibrium",
                 verbatimTextOutput("nash_equilibrium_text"),
                 plotOutput("nash_plot", height = "600px")
        ),
        tabPanel("Social Optimum",
                 verbatimTextOutput("social_optimum_text"),
                 plotOutput("social_optimum_plot", height = "600px")
        )
      )
    )
  )
)

# SERVER
server <- function(input, output) {
  # Optimization for r* given lambda and w (Principal's best response)
  optimal_r <- function(lambda, w) {
    # If UH - epsilon < w, the principal's utility is fixed at epsilon,
    # regardless of r. In this case, the principal would choose r=0.5
    # to minimize the cost of information (as mutual_info(p, 0.5) = 0).
    if (UH - epsilon < w|lambda==0.3) {
      return(0.5) 
    } else {
      obj <- function(r_val) -total_utility(r_val, lambda, w)
      # Optimize over a range of r values.
      # The range is from 0.5 + eps to 1 - eps to avoid log(0) and ensure valid precision.
      # If the optimum is exactly 0.5, this range will find it if it's the boundary.
      opt <- optimize(obj, c(0.5 + eps, 1 - eps))
      return(opt$minimum)
    }
  }
  
  # Reactive grid parameters for PLOTTING RANGES AND COMPUTATIONS
  lambda_seq <-  reactive({
    seq(lambda_range[1], lambda_range[2], length.out = input$grid_size)
  })
  
  w_seq <-  reactive({
    seq(w_range[1], w_range[2], length.out = input$grid_size)
  })
  
  r_seq <- reactive({
    # Ensure r covers the full range [0.5, 1] for optimization, avoiding exact 0.5 or 1 for log
    seq(0.5 + eps, 1 - eps, length.out = input$grid_size)
  })
  
  # Principal's best response surface
  grid_data <- reactive({
    req(lambda_seq(), w_seq()) # Now using lambda_seq and w_seq for computations
    expand.grid(lambda = lambda_seq(), w = w_seq()) %>%
      distinct() %>%
      rowwise() %>%
      mutate(
        r_star = optimal_r(lambda, w),
        profit = total_profit(r_star, lambda, w)
      ) %>%
      ungroup()
  })
  
  output$r_star_plot <- renderPlot({
    df <- grid_data()
    ggplot(df, aes(x = lambda, y = w, z = r_star)) +
      geom_contour_filled() +
      labs(title = "Principal's Optimal Signal Precision r*",
           x = expression(lambda),
           y = "w",
           fill = "r*") +
      theme_minimal()
  })
  
  output$profit_plot <- renderPlot({
    df <- grid_data()
    ggplot(df, aes(x = lambda, y = w, z = profit)) +
      geom_contour_filled() +
      geom_point(data = df[which.max(df$profit), ], aes(x = lambda, y = w), color = "blue", size = 3) +
      labs(title = "Total Profit",
           x = expression(lambda),
           y = "w",
           fill = "Profit") +
      theme_minimal()
  })
  
  # Agent's best response curve (calculates over finite lambda/w ranges)
  optimal_values_r <- reactive({
    req(lambda_seq(), w_seq(), r_seq()) # Now using lambda_seq and w_seq for computations
    
    result <- lapply(r_seq(), function(r_val) {
      grid <- expand.grid(lambda = lambda_seq(), w = w_seq()) %>%
        distinct()
      grid$Pbuy <- mapply(P_buy, r = r_val, w = grid$w)
      grid$profit <- grid$w * grid$Pbuy + grid$lambda * mutual_info_bernoulli(p, r_val)
      
      # Find all rows that maximize profit for this r_val
      max_profit_for_r <- max(grid$profit)
      opt_points <- grid %>% filter(profit == max_profit_for_r)
      
      # Return all optimal (lambda, w) pairs for this r_val
      data.frame(r = r_val, lambda_star = opt_points$lambda, w_star = opt_points$w, profit_agent = max_profit_for_r)
    })
    
    df <- bind_rows(result)
    return(df)
  })
  
  output$optimalLambdaWPlot <- renderPlot({
    df <- optimal_values_r()
    # No filtering for finite values needed here as values are always finite
    ggplot(df, aes(x = r)) +
      geom_line(aes(y = lambda_star, color = "λ*"), size = 1) +
      geom_line(aes(y = w_star, color = "w*"), size = 1) +
      labs(
        title = "Agent's Optimal λ* and w* as functions of r",
        x = "r (signal accuracy)",
        y = "Optimal value",
        color = ""
      ) +
      scale_color_manual(values = c("λ*" = "blue", "w*" = "darkred")) +
      theme_minimal()
  })
  
  output$totalUtilityPlot <- renderPlot({
    df_opt <- optimal_values_r()
    df_utility <- df_opt %>%
      rowwise() %>%
      mutate(utility = total_utility(r, lambda_star, w_star)) %>%
      ungroup()
    
    max_utility_value <- max(df_utility$utility)
    max_utility_points <- df_utility %>% filter(utility == max_utility_value)
    
    ggplot(df_utility, aes(x = r, y = utility)) +
      geom_line(color = "darkgreen", size = 1.2) +
      geom_point(data = max_utility_points, aes(x = r, y = utility), color = "red", size = 4) +
      labs(
        title = "Principal's Utility on Agent's Best Response Curve",
        x = "r (signal accuracy)",
        y = "Total Utility"
      ) +
      theme_minimal()
  })
  
  # Nash Equilibrium calculation and plot
  nash_data <- reactive({
    req(r_seq())
    agent_br <- optimal_values_r() # Agent's best response, now with finite values
    
    nash_df <- agent_br %>%
      rowwise() %>%
      mutate(
        principal_best_r = optimal_r(lambda_star, w_star), # lambda_star and w_star are always finite
        is_nash = abs(r - principal_best_r) < 0.005 # Tolerance for Nash equilibrium
      ) %>%
      ungroup()
    
    nash_points <- nash_df %>% filter(is_nash)
    
    return(list(full_df = nash_df, nash = nash_points))
  })
  
  output$nash_equilibrium_text <- renderPrint({
    nash_points <- nash_data()$nash
    if (nrow(nash_points) > 0) {
      cat("Nash Equilibrium found at:\n")
      print(nash_points)
    } else {
      cat("No Nash Equilibrium found within the search grid. Try adjusting the grid range or size.")
    }
  })
  
  output$nash_plot <- renderPlot({
    data <- nash_data()
    # No filtering for finite values needed here as values are always finite
    principal_br_plot <- ggplot(data$full_df, aes(x = lambda_star, y = w_star, color = r)) +
      geom_point(aes(size = principal_best_r)) +
      scale_size_continuous(range = c(1, 5), name = "Principal's Best r") +
      labs(
        title = "Nash Equilibrium as an intersection of best responses",
        x = expression(lambda),
        y = "w",
        color = "Agent's Best r"
      ) +
      theme_minimal()
    
    if (nrow(data$nash) > 0) {
      # Add convex hull if there are 3 or more Nash points
      if (nrow(data$nash) >= 3) {
        hull_points <- data$nash[chull(data$nash$lambda_star, data$nash$w_star), ]
        principal_br_plot <- principal_br_plot +
          geom_polygon(data = hull_points, aes(x = lambda_star, y = w_star), 
                       fill = "black", alpha = 0.1, inherit.aes = FALSE) # Highlight area
      }
      
      # Plot individual Nash points
      principal_br_plot <- principal_br_plot +
        geom_point(data = data$nash, aes(x = lambda_star, y = w_star),
                   color = "black", shape = 4, size = 5, stroke = 1.5) # Slightly smaller points
    }
    
    principal_br_plot
  })
  
  # Social Optimum Calculation (calculates over finite lambda/w ranges)
  social_optimum <- reactive({
    req(r_seq(), lambda_seq(), w_seq()) # Now using lambda_seq and w_seq for computations
    
    # Create a 3D grid for all parameters
    full_grid <- expand.grid(
      r = r_seq(),
      lambda = lambda_seq(), # Using finite lambda range
      w = w_seq() # Using finite w range
    ) %>%
      distinct() # Remove duplicates
    
    # Calculate social welfare for each point in the grid
    full_grid$social_welfare <- mapply(social_welfare,
                                       r = full_grid$r,
                                       lambda = full_grid$lambda,
                                       w = full_grid$w)
    
    # Find all rows with the maximum social welfare
    max_welfare <- max(full_grid$social_welfare)
    optimum_points <- full_grid %>% filter(social_welfare == max_welfare)
    
    return(optimum_points)
  })
  
  # Agent First scenario calculation (Agent maximizes profit, Principal reacts with optimal r)
  agent_first_optimum <- reactive({
    df_opt <- optimal_values_r() # This now contains finite values
    max_profit_agent_value <- max(df_opt$profit_agent)
    max_profit_points <- df_opt %>%
      filter(profit_agent == max_profit_agent_value)
    
    return(max_profit_points)
  })
  
  output$social_optimum_text <- renderPrint({
    optimum <- social_optimum()
    cat("Social Optimum found at:\n")
    print(optimum)
  })
  
  # Summary table
  summary_data <- reactive({
    req(nash_data(), social_optimum(), agent_first_optimum())
    
    summary_rows <- list()
    
    # Principal's Max Profit (This is a single point from the grid_data)
    principal_opt_point <- grid_data()[which.max(grid_data()$profit), ]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Scenario = "Principal's Max Profit",
      r = principal_opt_point$r_star,
      w = principal_opt_point$w,
      lambda = principal_opt_point$lambda
    )
    
    # Agent First (can have multiple optimal points)
    agent_first_points <- agent_first_optimum()
    for (i in 1:nrow(agent_first_points)) {
      point <- agent_first_points[i, ]
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Scenario = paste0("Agent First (", i, ")"),
        r = point$r,
        w = point$w_star,
        lambda = point$lambda_star
      )
    }
    
    # Nash Equilibrium (can have multiple points)
    nash_points <- nash_data()$nash
    if (nrow(nash_points) > 0) {
      for (i in 1:nrow(nash_points)) {
        point <- nash_points[i, ]
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          Scenario = paste0("Nash Equilibrium (", i, ")"),
          r = point$r,
          w = point$w_star,
          lambda = point$lambda_star
        )
      }
    } else {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Scenario = "Nash Equilibrium",
        r = NA, w = NA, lambda = NA # Indicate no Nash if none found
      )
    }
    
    # Social Optimum (can have multiple optimal points)
    social_opt_points <- social_optimum()
    for (i in 1:nrow(social_opt_points)) {
      point <- social_opt_points[i, ]
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Scenario = paste0("Social Optimum (", i, ")"),
        r = point$r,
        w = point$w,
        lambda = point$lambda
      )
    }
    
    # Combine all rows into a single data frame
    summary_df <- bind_rows(summary_rows)
    
    # Calculate all metrics for each scenario
    summary_df_calculated <- summary_df %>%
      rowwise() %>%
      mutate(
        P_buy = P_buy(r, w),
        Utility_Principal = total_utility(r, lambda, w),
        Profit_Agent = total_profit(r, lambda, w),
        Total_Welfare = social_welfare(r, lambda, w)
      ) %>%
      ungroup() %>%
      mutate(
        across(c(r, w, lambda, P_buy, Utility_Principal, Profit_Agent, Total_Welfare), ~round(., 4))
      )
    
    return(summary_df_calculated)
  })
  
  output$summary_table <- renderTable({
    summary_data()
  })
}

shinyApp(ui = ui, server = server)

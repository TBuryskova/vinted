library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT) # For a nice-looking summary table

# User-defined parameters
a <- 10     # utility of buying when X = H
b <- 0    # utility of buying when X = L
c <- 8     # utility of leaving
p <- 0.5    # prior
eps <- 1e-8 # small value to prevent log from collapsing


# Function: probability of buying
P_buy <- function(r, w) {
  
  p_buy <- ifelse(r >= (1-p)*(c+w-a)/(-(1-p)*a+p*b+(c+w)*(1-2*p)), 1, 
                  ifelse(r >= ((c+w)*p-b*p)/((1-p)*a-p*b+(c+w)*(2*p-1)), 0.5,0))
        
  return(p_buy)
}

# Function: mutual information
mutual_info <- function(r) {
  r_clamped <- pmin(pmax(r, eps), 1 - eps)
  term1 <- ifelse(r_clamped > 0, r_clamped * log(r_clamped / p), 0)
  term2 <- ifelse(r_clamped < 1, (1 - r_clamped) * log((1 - r_clamped) / (1 - p)), 0)
  return(term1 + term2)
}

# Total profit of the agent
total_profit <- function(r, lambda, w) {
  w * P_buy(r, w) + lambda * mutual_info(r)
}

# Total utility of the principal
total_utility <- function(r, lambda, w) {
  EU_buy_s0 <- (1 - r) * a + r * b - w
  EU_buy_s1 <- r * a + (1-r) * b - w
  
  EU_s0 <- max(EU_buy_s0, c)
  EU_s1 <- max(EU_buy_s1, c)
  
  EU_total <- p * EU_s0 + (1 - p) * EU_s1 - lambda * mutual_info(r)
  
  return(EU_total)
}
# Social Welfare Function
social_welfare <- function(r, lambda, w) {
  total_utility(r, lambda, w) + total_profit(r, lambda, w)
}

# Optimization for r* given lambda and w
optimal_r <- function(lambda, w) {
  obj <- function(r) -total_utility(r, lambda, w)
  opt <- optimize(obj, c(0.5 + eps, 1 - eps))
  return(opt$minimum)
}

# UI
ui <- fluidPage(
  titlePanel("Optimal Signal Design"),
  sidebarLayout(
    sidebarPanel(
      h4("Grid and Parameter Controls"),
      sliderInput("lambda_range", "λ Range",
                  min = 0, max = 100, value = c(0.01, 1), step = 0.01),
      sliderInput("w_range", "w Range",
                  min = 0, max = 100, value = c(0.01, 5), step = 0.01),
      numericInput("grid_size", "Grid Points per Axis",
                   value = 50, min = 10, max = 500, step = 10),
      helpText("Note: Increasing the grid size may slow down calculations.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary",
                 DTOutput("summary_table")
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
# SERVER
server <- function(input, output) {
  
  # Reactive grid parameters
  lambda_seq <- reactive({
    seq(input$lambda_range[1], input$lambda_range[2], length.out = input$grid_size)
  })
  
  w_seq <- reactive({
    seq(input$w_range[1], input$w_range[2], length.out = input$grid_size)
  })
  
  r_seq <- reactive({
    seq(0.5000, 0.9999, length.out = input$grid_size)
  })
  
  # Principal's best response surface
  grid_data <- reactive({
    req(lambda_seq(), w_seq())
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
  
  # Agent's best response curve
  optimal_values_r <- reactive({
    req(lambda_seq(), w_seq(), r_seq())
    MI <- mutual_info
    
    result <- lapply(r_seq(), function(r) {
      grid <- expand.grid(lambda = lambda_seq(), w = w_seq()) %>%
        distinct()
      grid$Pbuy <- mapply(P_buy, r = r, w = grid$w)
      grid$profit <- grid$w * grid$Pbuy + grid$lambda * MI(r)
      opt <- grid[which.max(grid$profit), ]
      data.frame(r = r, lambda_star = opt$lambda, w_star = opt$w)
    })
    
    df <- bind_rows(result)
    return(df)
  })
  
  output$optimalLambdaWPlot <- renderPlot({
    df <- optimal_values_r()
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
    
    max_utility_point <- df_utility[which.max(df_utility$utility), ]
    
    ggplot(df_utility, aes(x = r, y = utility)) +
      geom_line(color = "darkgreen", size = 1.2) +
      geom_point(data = max_utility_point, aes(x = r, y = utility), color = "red", size = 4) +
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
    agent_br <- optimal_values_r()
    
    nash_df <- agent_br %>%
      rowwise() %>%
      mutate(
        principal_best_r = optimal_r(lambda_star, w_star),
        is_nash = abs(r - principal_best_r) < 0.005 # Tolerance
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
      principal_br_plot <- principal_br_plot +
        geom_point(data = data$nash, aes(x = lambda_star, y = w_star),
                   color = "black", shape = 4, size = 10, stroke = 2)
    }
    
    principal_br_plot
  })
  
  # Social Optimum Calculation
  social_optimum <- reactive({
    req(r_seq(), lambda_seq(), w_seq())
    
    # Create a 3D grid for all parameters
    full_grid <- expand.grid(
      r = r_seq(),
      lambda = lambda_seq(),
      w = w_seq()
    ) %>%
      distinct() # Add this line to remove duplicates
    
    # Calculate social welfare for each point in the grid
    full_grid$social_welfare <- mapply(social_welfare,
                                       r = full_grid$r,
                                       lambda = full_grid$lambda,
                                       w = full_grid$w)
    
    # Find the row with the maximum social welfare
    optimum_point <- full_grid[which.max(full_grid$social_welfare), ]
    
    return(optimum_point)
  })
  
  # Agent First scenario calculation
  agent_first_optimum <- reactive({
    df_opt <- optimal_values_r()
    max_profit_point <- df_opt %>%
      rowwise() %>%
      mutate(profit_agent = total_profit(r, lambda_star, w_star)) %>%
      ungroup() %>%
      filter(profit_agent == max(profit_agent)) %>%
      slice(1)
    
    return(max_profit_point)
  })
  
  output$social_optimum_text <- renderPrint({
    optimum <- social_optimum()
    cat("Social Optimum found at:\n")
    print(optimum)
  })
  
  # Summary table
  summary_data <- reactive({
    req(nash_data(), social_optimum(), agent_first_optimum())
    
    # Get data for each scenario
    principal_opt_point <- grid_data()[which.max(grid_data()$profit), ]
    nash_point <- nash_data()$nash
    social_opt_point <- social_optimum()
    agent_first_point <- agent_first_optimum()
    
    # Handle the case where no Nash Equilibrium is found
    if (nrow(nash_point) == 0) {
      nash_r <- NA
      nash_w <- NA
      nash_lambda <- NA
    } else {
      # If multiple Nash points are found, take the first one.
      nash_r <- nash_point$r[1]
      nash_w <- nash_point$w_star[1]
      nash_lambda <- nash_point$lambda_star[1]
    }
    
    # Extract Agent First values
    agent_r <- agent_first_point$r
    agent_w <- agent_first_point$w_star
    agent_lambda <- agent_first_point$lambda_star
    
    # Create a data frame for summary
    scenario <- c("Principal's Max Profit", "Agent First", "Nash Equilibrium", "Social Optimum")
    
    r_vals <- c(principal_opt_point$r_star, agent_r, nash_r, social_opt_point$r)
    w_vals <- c(principal_opt_point$w, agent_w, nash_w, social_opt_point$w)
    lambda_vals <- c(principal_opt_point$lambda, agent_lambda, nash_lambda, social_opt_point$lambda)
    
    # Calculate all metrics for each scenario
    pbuy_vals <- mapply(P_buy, r = r_vals, w = w_vals)
    utility_vals <- mapply(total_utility, r = r_vals, lambda = lambda_vals, w = w_vals)
    profit_vals <- mapply(total_profit, r = r_vals, lambda = lambda_vals, w = w_vals)
    welfare_vals <- mapply(social_welfare, r = r_vals, lambda = lambda_vals, w = w_vals)
    
    data.frame(
      Scenario = scenario,
      r = round(r_vals, 4),
      w = round(w_vals, 4),
      lambda = round(lambda_vals, 4),
      P_buy = round(pbuy_vals, 4),
      Utility_Principal = round(utility_vals, 4),
      Profit_Agent = round(profit_vals, 4),
      Total_Welfare = round(welfare_vals, 4)
    )
  })
  
  output$summary_table <- renderTable({
    summary_data()
  })
}

shinyApp(ui = ui, server = server)

shinyApp(ui = ui, server = server)


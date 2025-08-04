library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# Uživatelské parametry
a <- 2     # užitek z koupě když X = 0
b <- 8      # užitek z koupě když X = 1
c <- 5      # užitek z odchodu
p <- 0.5    # prior
eps <- 1e-8 # malá hodnota, aby se log nezhroutil

# Funkce: posterior belief po signálu
posterior <- function(r, s) {
  if (s == 1) {
    return((r * p) / (r * p + (1 - r) * (1 - p)))
  } else {
    return(((1 - r) * p) / ((1 - r) * p + r * (1 - p)))
  }
}

# Funkce: pravděpodobnost koupě
P_buy <- function(r, w) {
  p1 <- posterior(r, 1)
  p0 <- posterior(r, 0)
  
  u1 <- p1 * b + (1 - p1) * a - w
  u0 <- p0 * b + (1 - p0) * a - w
  
  p_buy <- ifelse(u1 >= c, r * p + (1 - r) * (1 - p), 0) +
    ifelse(u0 >= c, (1 - r) * p + r * (1 - p), 0)
  return(p_buy)
}

# Funkce: mutual information
mutual_info <- function(r) {
  r_clamped <- pmin(pmax(r, eps), 1 - eps)
  term1 <- ifelse(r_clamped > 0, r_clamped * log(r_clamped / p), 0)
  term2 <- ifelse(r_clamped < 1, (1 - r_clamped) * log((1 - r_clamped) / (1 - p)), 0)
  return(term1 + term2)
}

total_profit <- function(r, lambda, w) {
  w * P_buy(r, w) + lambda * mutual_info(r)
}

total_utility <- function(r, lambda, w) {
  EU_buy_s0 <- (1 - r) * a + r * b - w
  EU_buy_s1 <- r * a + (1-r) * b - w
  
  EU_s0 <- max(EU_buy_s0, c)
  EU_s1 <- max(EU_buy_s1, c)
  
  EU_total <- p * EU_s0 + (1 - p) * EU_s1 - lambda * mutual_info(r)
  
  return(EU_total)
}

# Optimalizace r* pro dané lambda a w
optimal_r <- function(lambda, w) {
  obj <- function(r) -total_utility(r, lambda, w)
  opt <- optimize(obj, c(eps, 1 - eps))
  return(opt$minimum)
}

# UI
ui <- fluidPage(
  titlePanel("Optimal Signal Design"),
  tabsetPanel(
    tabPanel("Principal's Best Response",
             plotOutput("r_star_plot", height = "600px"),
             plotOutput("profit_plot", height = "600px")
    ),
    tabPanel("Agent's Best Response",
             plotOutput("optimalLambdaWPlot", height = "400px"),
             plotOutput("totalUtilityPlot", height = "400px")
    ),
    tabPanel("Nash Equilibrium",
             verbatimTextOutput("nash_equilibrium_text"),
             plotOutput("nash_plot", height = "600px")
    )
  )
)

# SERVER
server <- function(input, output) {
  
  # Principal's best response surface
  grid_data <- reactive({
    lambda_seq <- seq(0.01, 1, length.out = 100)
    w_seq <- seq(0, 5, length.out = 100)
    
    expand.grid(lambda = lambda_seq, w = w_seq) %>%
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
    r_seq <- seq(0.5001, 0.999, length.out = 100)
    
    lambda_grid <- seq(0.01, 1, length.out = 30)
    w_grid <- seq(0.01, 5, length.out = 30)
    
    MI <- mutual_info
    
    result <- lapply(r_seq, function(r) {
      grid <- expand.grid(lambda = lambda_grid, w = w_grid)
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
    
    ggplot(df_utility, aes(x = r, y = utility)) +
      geom_line(color = "darkgreen", size = 1.2) +
      labs(
        title = "Principal's Utility on Agent's Best Response Curve",
        x = "r (signal accuracy)",
        y = "Total Utility"
      ) +
      theme_minimal()
  })
  
  # Nash Equilibrium calculation and plot
  nash_data <- reactive({
    agent_br <- optimal_values_r()
    
    # We are looking for a point on the agent's best response curve
    # where r is also the principal's best response to (lambda*, w*)
    nash_df <- agent_br %>%
      rowwise() %>%
      mutate(
        principal_best_r = optimal_r(lambda_star, w_star),
        # Check if r is approximately equal to the principal's best response
        is_nash = abs(r - principal_best_r) < 0.005 # Tolerance
      ) %>%
      ungroup()
    
    # Find the Nash equilibrium points
    nash_points <- nash_df %>% filter(is_nash)
    
    return(list(full_df = nash_df, nash = nash_points))
  })
  
  output$nash_equilibrium_text <- renderPrint({
    nash_points <- nash_data()$nash
    if (nrow(nash_points) > 0) {
      cat("Nash Equilibrium found at:\n")
      print(nash_points)
    } else {
      cat("No Nash Equilibrium found within the search grid.")
    }
  })
  
  output$nash_plot <- renderPlot({
    data <- nash_data()
    
    # Plot the principal's best response surface
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
    
    # Add the Nash point(s) to the plot
    if (nrow(data$nash) > 0) {
      principal_br_plot <- principal_br_plot +
        geom_point(data = data$nash, aes(x = lambda_star, y = w_star), 
                   color = "black", shape = 4, size = 10, stroke = 2)
    }
    
    principal_br_plot
  })
}

shinyApp(ui = ui, server = server)
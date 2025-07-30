# Install necessary packages if you haven't already
# install.packages(c("shiny", "ggplot2", "dplyr", "tidyr", "viridis"))

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis) # For better color scales in geom_contour_filled

# --- Helper Functions for Core Model Components ---

# A(lambda,w) = exp[ (U_L–w)/lambda ]
A <- function(lambda, w, UL) {
  # Handle potential division by zero or extreme values for lambda
  if (lambda <= 0 || is.infinite(lambda) || is.na(lambda)) return(NA)
  exp((UL - w) / lambda)
}

# B(lambda,w) = exp[ (U_H–w)/lambda ]
B <- function(lambda, w, UH) {
  if (lambda <= 0 || is.infinite(lambda) || is.na(lambda)) return(NA)
  exp((UH - w) / lambda)
}

# C(lambda) = exp[ epsilon/lambda ]
C <- function(lambda, epsilon) {
  if (lambda <= 0 || is.infinite(lambda) || is.na(lambda)) return(NA)
  exp(epsilon / lambda)
}

# D(lambda,w) = A + C
D <- function(lambda, w, UL, epsilon) {
  val_A <- A(lambda, w, UL)
  val_C <- C(lambda, epsilon)
  if (is.na(val_A) || is.na(val_C)) return(NA)
  val_A + val_C
}

# F(lambda,w) = B + C
F_val <- function(lambda, w, UH, epsilon) { # Renamed to F_val to avoid conflict with base R 'F'
  val_B <- B(lambda, w, UH)
  val_C <- C(lambda, epsilon)
  if (is.na(val_B) || is.na(val_C)) return(NA)
  val_B + val_C
}

# S(lambda,w) = A/D + B/F (This is a* + b*)
S_val <- function(lambda, w, UL, UH, epsilon) {
  val_A <- A(lambda, w, UL)
  val_D <- D(lambda, w, UL, epsilon)
  val_B <- B(lambda, w, UH)
  val_F <- F_val(lambda, w, UH, epsilon)
  
  # Check for division by zero or NA/Inf
  a_star_term <- if (val_D == 0 || is.na(val_D) || is.infinite(val_D)) NA else val_A / val_D
  b_star_term <- if (val_F == 0 || is.na(val_F) || is.infinite(val_F)) NA else val_B / val_F
  
  if (is.na(a_star_term) || is.na(b_star_term)) return(NA)
  a_star_term + b_star_term
}

# a*(lambda,w) = A/D
a_star_val <- function(lambda, w, UL, epsilon) {
  val_A <- A(lambda, w, UL)
  val_D <- D(lambda, w, UL, epsilon)
  if (val_D == 0 || is.na(val_D) || is.infinite(val_D)) NA else val_A / val_D
}

# b*(lambda,w) = B/F
b_star_val <- function(lambda, w, UH, epsilon) {
  val_B <- B(lambda, w, UH)
  val_F <- F_val(lambda, w, UH, epsilon)
  if (val_F == 0 || is.na(val_F) || is.infinite(val_F)) NA else val_B / val_F
}

# --- Firm Profit Function (Objective Function) ---
Firm_Profit_Function <- function(lambda, w, UL, UH, epsilon) {
  if (lambda <= 0 || is.infinite(lambda) || is.na(lambda)) return(NA)
  
  a_s <- a_star_val(lambda, w, UL, epsilon)
  b_s <- b_star_val(lambda, w, UH, epsilon)
  s <- S_val(lambda, w, UL, UH, epsilon) # s is a* + b*
  
  if (any(is.na(c(a_s, b_s, s))) || any(is.infinite(c(a_s, b_s, s)))) return(NA)
  
  # Term 1: 0.5 * (a* + b*) * w
  term1 <- 0.5 * (a_s + b_s) * w
  if (is.na(term1) || is.infinite(term1)) return(NA)
  
  # Term 2: lambda * 0.5 * [...] (Mutual Information part)
  safe_log <- function(x) {
    if (x <= 0 || is.na(x) || is.infinite(x)) NA else log(x)
  }
  
  log_part_a1 <- safe_log((2 * a_s) / s)
  log_part_a2 <- safe_log((2 * (1 - a_s)) / (2 - s))
  log_part_b1 <- safe_log((2 * b_s) / s)
  log_part_b2 <- safe_log((2 * (1 - b_s)) / (2 - s))
  
  if (any(is.na(c(log_part_a1, log_part_a2, log_part_b1, log_part_b2))) ||
      any(is.infinite(c(log_part_a1, log_part_a2, log_part_b1, log_part_b2)))) return(NA)
  
  term2_sum <- a_s * log_part_a1 + (1 - a_s) * log_part_a2 +
    b_s * log_part_b1 + (1 - b_s) * log_part_b2
  
  if (is.na(term2_sum) || is.infinite(term2_sum)) return(NA)
  
  term2 <- lambda * 0.5 * term2_sum
  if (is.na(term2) || is.infinite(term2)) return(NA)
  
  result <- term1 + term2
  if (is.na(result) || is.infinite(result)) return(NA)
  return(result)
}

# --- User Utility Function ---
User_Utility_Function <- function(lambda, w, UL, UH, epsilon) {
  if (lambda <= 0 || is.infinite(lambda) || is.na(lambda)) return(NA)
  a_s <- a_star_val(lambda, w, UL, epsilon)
  b_s <- b_star_val(lambda, w, UH, epsilon)
  if (is.na(a_s) || is.na(b_s)) return(NA)
  
  term_a_utility <- a_s * (UL - w) + (1 - a_s) * epsilon
  term_b_utility <- b_s * (UH - w) + (1 - b_s) * epsilon
  if (is.na(term_a_utility) || is.na(term_b_utility)) return(NA)
  
  0.5 * term_a_utility + 0.5 * term_b_utility
}

# --- Unconditional Probability of Buying ---
Unconditional_Prob_Buying <- function(lambda, w, UL, UH, epsilon) {
  a_s <- a_star_val(lambda, w, UL, epsilon)
  b_s <- b_star_val(lambda, w, UH, epsilon)
  if (is.na(a_s) || is.na(b_s)) return(NA)
  0.5 * (a_s + b_s)
}

# --- Optimal Mutual Information ---
Optimal_Mutual_Information <- function(lambda, w, UL, UH, epsilon) {
  if (lambda <= 0 || is.infinite(lambda) || is.na(lambda)) return(NA)
  
  a_s <- a_star_val(lambda, w, UL, epsilon)
  b_s <- b_star_val(lambda, w, UH, epsilon)
  s <- S_val(lambda, w, UL, UH, epsilon) # s is a* + b*
  
  if (any(is.na(c(a_s, b_s, s))) || any(is.infinite(c(a_s, b_s, s)))) return(NA)
  
  safe_log <- function(x) {
    if (x <= 0 || is.na(x) || is.infinite(x)) NA else log(x)
  }
  
  log_part_a1 <- safe_log((2 * a_s) / s)
  log_part_a2 <- safe_log((2 * (1 - a_s)) / (2 - s))
  log_part_b1 <- safe_log((2 * b_s) / s)
  log_part_b2 <- safe_log((2 * (1 - b_s)) / (2 - s))
  
  if (any(is.na(c(log_part_a1, log_part_a2, log_part_b1, log_part_b2))) ||
      any(is.infinite(c(log_part_a1, log_part_a2, log_part_b1, log_part_b2)))) return(NA)
  
  term_sum <- a_s * log_part_a1 + (1 - a_s) * log_part_a2 +
    b_s * log_part_b1 + (1 - b_s) * log_part_b2
  
  if (is.na(term_sum) || is.infinite(term_sum)) return(NA)
  
  lambda * 0.5 * term_sum
}


# --- Shiny UI ---
ui <- fluidPage(
  titlePanel("Optimal Point Plotter for Firm Profit"),
  sidebarLayout(
    sidebarPanel(
      h3("Model Parameters"),
      sliderInput("UL", "U_L:", min = 0, max = 10, value = 0, step = 0.1),
      sliderInput("UH", "U_H:", min = 0, max = 10, value = 10, step = 0.1),
      sliderInput("epsilon", "epsilon:", min =0, max = 10, value = 5, step = 0.1),
      hr(),
      h3("Plotting Range for (λ, w)"),
      sliderInput("lambda_range", "Lambda (λ) Range (min, max):", min = 0.01, max = 10, value = c(0.01, 10), step = 0.1),
      sliderInput("w_range", "W Range (min, max):", min = 0, max = 10, value = c(0, 10), step = 0.1),
      numericInput("grid_res", "Grid Resolution (points per axis):", value = 50, min = 10, max = 200)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Firm First Plot", plotOutput("profit_plot", height = "600px")),
        tabPanel("Optimal Values", htmlOutput("optimal_values_output")), # New output for optimal values
        tabPanel("About", htmlOutput("about_text"))
      )
    )
  )
)

# --- Shiny Server ---
server <- function(input, output) {
  
  # Reactive expression to generate the grid and calculate Objective Function values
  grid_data <- reactive({
    req(input$lambda_range, input$w_range, input$grid_res)
    
    lambda_vals <- seq(input$lambda_range[1], input$lambda_range[2], length.out = input$grid_res)
    w_vals <- seq(input$w_range[1], input$w_range[2], length.out = input$grid_res)
    
    # Create a data frame for all combinations of lambda and w
    plot_df <- expand.grid(lambda = lambda_vals, w = w_vals)
    
    # Calculate Firm Profit for each point in the grid
    plot_df %>%
      rowwise() %>%
      mutate(
        firm_profit = Firm_Profit_Function(lambda, w, input$UL, input$UH, input$epsilon)
      ) %>%
      ungroup()
  })
  
  # Reactive expression to find the global maximum of the firm profit on the grid
  global_max_point <- reactive({
    df <- grid_data() %>% drop_na(firm_profit) # Ensure only valid profit values are considered
    if (nrow(df) == 0) return(NULL)
    
    # Find the row with the maximum firm_profit
    max_row <- df %>% slice_max(firm_profit, n = 1, with_ties = FALSE)
    
    # Calculate additional metrics at this optimal point
    if (nrow(max_row) > 0) {
      optimal_lambda <- max_row$lambda
      optimal_w <- max_row$w
      
      max_row$user_utility <- User_Utility_Function(optimal_lambda, optimal_w, input$UL, input$UH, input$epsilon)
      max_row$uncond_prob_buying <- Unconditional_Prob_Buying(optimal_lambda, optimal_w, input$UL, input$UH, input$epsilon)
      max_row$optimal_mutual_info <- Optimal_Mutual_Information(optimal_lambda, optimal_w, input$UL, input$UH, input$epsilon)
      
      # Also calculate a* and b* at the optimal point
      max_row$a_star_opt <- a_star_val(optimal_lambda, optimal_w, input$UL, input$epsilon)
      max_row$b_star_opt <- b_star_val(optimal_lambda, optimal_w, input$UH, input$epsilon)
      
    }
    max_row
  })
  
  
  # Render the plot
  output$profit_plot <- renderPlot({
    df_clean <- grid_data() %>% drop_na(firm_profit) # Clean data for plotting objective function surface
    max_point <- global_max_point()
    
    p <- ggplot(df_clean, aes(x = lambda, y = w)) +
      # Plot the firm profit as a filled contour map
      geom_contour_filled(aes(z = firm_profit), alpha = 0.7) +
      # Changed to scale_fill_viridis_d for discrete fill values from geom_contour_filled
      scale_fill_viridis_d(option = "plasma", name = "Firm Profit") +
      labs(title = "Firm First and Global Maximum",
           subtitle = paste0("U_L = ", input$UL, ", U_H = ", input$UH, ", epsilon = ", input$epsilon),
           x = "Lambda (λ)",
           y = "W") +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) +
      # Ensure plot limits match slider ranges
      coord_cartesian(xlim = input$lambda_range, ylim = input$w_range)
    
    # Add the global maximum point if found
    if (!is.null(max_point) && nrow(max_point) > 0) {
      p <- p + geom_point(data = max_point,
                          aes(x = lambda, y = w),
                          color = "red", shape = 8, size = 5, stroke = 2) + # Red star shape for global max
        geom_text(data = max_point,
                  aes(x = lambda, y = w, label = paste0("(λ=", round(lambda, 2), ", w=", round(w, 2), ")")),
                  vjust = -1.5, hjust = 0.5, color = "red", size = 5)
    }
    
    p
  })
  
  # Render the optimal values output
  output$optimal_values_output <- renderUI({
    max_point <- global_max_point()
    if (is.null(max_point) || nrow(max_point) == 0) {
      return(HTML("<p>No optimal point found within the current grid and parameters.</p>"))
    }
    
    HTML(paste0(
      "<h4>Optimal Firm First:</h4>",
      "<p><b>Optimal &lambda;:</b> ", round(max_point$lambda, 4), "</p>",
      "<p><b>Optimal w:</b> ", round(max_point$w, 4), "</p>",
      "<p><b>Optimal Firm Profit:</b> ", round(max_point$firm_profit, 4), "</p>",
      "<p><b>User Utility:</b> ", round(max_point$user_utility, 4), "</p>",
      "<p><b>Optimal Firm Profit:</b> ", round(max_point$firm_profit, 4), "</p>",
      "<p><b>Welfare:</b> ", round(max_point$user_utility, 4)+round(max_point$firm_profit, 4), "</p>",
      "<p><b>Unconditional Probability of Buying (P*(Y=1)):</b> ", round(max_point$uncond_prob_buying, 4), "</p>",
      "<p><b>Optimal Mutual Information:</b> ", round(max_point$optimal_mutual_info, 4), "</p>",
      "<p><b>Optimal a*:</b> ", round(max_point$a_star_opt, 4), "</p>",
      "<p><b>Optimal b*:</b> ", round(max_point$b_star_opt, 4), "</p>"
    ))
  })
  
  
  # About section content
  output$about_text <- renderUI({
    HTML("
      <h4>About This Application</h4>
      <p>This R Shiny app visualizes the firm's profit function surface and identifies the global maximum point on the explored grid.</p>
      <p>The app works by:</p>
      <ol>
        <li>Calculating the firm's profit (objective function) value across a user-defined grid of $\\lambda$ and $w$ values.</li>
        <li>Plotting the profit function as a filled contour map, where different colors represent different profit values.</li>
        <li>Identifying the point on the entire grid that yields the highest profit value (Global Maximum on Grid).</li>
        <li>Computing and displaying several related metrics at this identified optimal point.</li>
      </ol>
      <h5>How to Use:</h5>
      <ul>
        <li>Adjust the values of $U_L$, $U_H$, and $\\epsilon$ using the sliders on the left.</li>
        <li>Set the desired range for $\\lambda$ and $w$ to focus your search for the global maximum.</li>
        <li>Increase the 'Grid Resolution' for a more detailed surface plot and potentially a more accurate global maximum identification, though this will increase computation time.</li>
      </ul>
      <h5>Plot Interpretation:</h5>
      <ul>
        <li><b>Filled Contours:</b> Represent the value of the firm's profit across the $(\\lambda, w)$ plane. Darker colors typically indicate higher values (depending on the color scale).</li>
        <li><span style='color:red;'>&#10038;</span> <b>Red Star:</b> Indicates the point on the grid where the firm's profit reaches its highest value (the Global Maximum on the Grid). A label next to it shows its coordinates.</li>
      </ul>
      <h5>Optimal Values Tab:</h5>
      <p>This tab displays the computed values for optimal $\\lambda$, optimal $w$, firm profit, user utility, unconditional probability of buying, optimal mutual information, and the optimal $a^*$ and $b^*$ at the identified global maximum point.</p>
      <h5>Important Notes:</h5>
      <ul>
        <li>The 'Global Maximum on Grid' is the highest point found within the discrete grid. It serves as a good visual confirmation but might not be the exact global maximum if the grid resolution is too coarse or if the true maximum lies outside the defined lambda/w range.</li>
        <li>Numerical instabilities (e.g., very large/small numbers, division by zero, or invalid logarithm arguments) might occur for certain parameter combinations or at the edges of the plotting range. Points where calculations result in non-finite values (NA, Inf, NaN) are automatically excluded.</li>
      </ul>
    ")
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

# ========== Combined Shiny App ==========
# This app integrates four separate Shiny apps into one.
# You can choose which variable to vary (p, lambda, sigma, or gamma)
# and the UI will adapt accordingly. Users can select which plots to display.

# --- 1. Load Libraries ---
library(shiny)
library(tidyverse)
library(ggplot2)

# --- 2. Core Mathematical & Helper Functions (Common to all apps) ---

#' A safe logarithm function that avoids log(0) errors.
safe_log <- function(x) {
  log(pmax(x, 1e-10))
}

#' Shannon Entropy function H(x).
H <- function(x) {
  # Use safe_log to prevent issues with p=0 or p=1 in calculations.
  -x * safe_log(x) - (1 - x) * safe_log(1 - x)
}

#' Calculates the expected utility of not acquiring information.
#' @param p Prior probability of product type PL.
#' @param gamma Risk aversion parameter.
#' @param sigma_pl Standard deviation of the PL product quality.
#' @param epsilon Utility of the outside option.
#' @return The maximum of the expected utility and the outside option.
expected_util <- function(p, gamma, sigma_pl, epsilon) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  pmax((1 - p) * E_CZ + p * E_PL, -epsilon)
}

#' The objective function to be minimized by optim().
#' This function calculates the negative expected utility minus the cost of information.
#' @param params A numeric vector with two values: a = P(buy|PL) and b = P(buy|CZ).
#' @return The value of the objective function.
objective <- function(params, p, gamma, sigma_pl, lambda, epsilon) {
  a <- params[1]
  b <- params[2]
  
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  
  EU <- p * (a * E_PL + (1 - a) * (-epsilon)) +
    (1 - p) * (b * E_CZ + (1 - b) * (-epsilon))
  
  f_PL1 <- p * a; f_PL0 <- p * (1 - a)
  f_CZ1 <- (1 - p) * b; f_CZ0 <- (1 - p) * (1 - b)
  f_y1 <- f_PL1 + f_CZ1; f_y0 <- f_PL0 + f_CZ0
  
  I <- 0
  if (f_PL1 > 0) I <- I + f_PL1 * safe_log(f_PL1 / (p * f_y1))
  if (f_PL0 > 0) I <- I + f_PL0 * safe_log(f_PL0 / (p * f_y0))
  if (f_CZ1 > 0) I <- I + f_CZ1 * safe_log(f_CZ1 / ((1 - p) * f_y1))
  if (f_CZ0 > 0) I <- I + f_CZ0 * safe_log(f_CZ0 / ((1 - p) * f_y0))
  
  -EU + lambda * I
}

#' Solves for the optimal a* and b* using optimization.
#' Handles edge cases for p near 0 or 1.
#' @return A list containing the result from optim() or the manually calculated values.
solve_ab <- function(p, gamma, sigma_pl, lambda, epsilon) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  
  # Handle edge cases where p is very close to 0 or 1, as optim can be unstable.
  if (p < 1e-3 || p > 1 - 1e-3) {
    a_star <- as.numeric(E_PL > -epsilon)
    b_star <- as.numeric(E_CZ > -epsilon)
    EU <- p * (a_star * E_PL + (1 - a_star) * (-epsilon)) +
      (1 - p) * (b_star * E_CZ + (1 - b_star) * (-epsilon))
    return(list(value = -EU, par = c(a_star, b_star)))
  }
  
  # Use L-BFGS-B for bounded optimization
  res <- optim(
    par = c(0.5, 0.5),
    fn = objective,
    p = p,
    gamma = gamma,
    sigma_pl = sigma_pl,
    lambda = lambda,
    epsilon = epsilon,
    method = "L-BFGS-B",
    lower = c(0, 0),
    upper = c(1, 1)
  )
  return(res)
}


# --- 3. Data Computation Functions (one for each variable) ---

# 3.1 Varying p
compute_curves_p <- function(sigma_pl, gamma, lambda, epsilon, n_grid = 100) {
  p_grid <- seq(0.001, 0.999, length.out = n_grid)
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  
  V_no_info <- expected_util(p_grid, gamma, sigma_pl, epsilon)
  V_full_info <- -p_grid * epsilon + (1 - p_grid) * E_CZ
  V_perfect <- pmax(-lambda * H(p_grid) + V_full_info, V_no_info)
  
  solve_res <- lapply(p_grid, function(p) solve_ab(p, gamma, sigma_pl, lambda, epsilon))
  V_RI <- sapply(solve_res, function(res) -res$value)
  a_star_RI <- sapply(solve_res, function(res) res$par[1])
  b_star_RI <- sapply(solve_res, function(res) res$par[2])
  
  buy_RI <- p_grid * a_star_RI + (1 - p_grid) * b_star_RI
  buy_no_info <- as.numeric(p_grid * E_PL + (1 - p_grid) * E_CZ > -epsilon)
  buy_full_info <- 1 - p_grid
  buy_perfect <- mapply(function(p, v_no) {
    signal_aq <- as.numeric(-lambda * H(p) + (-p * epsilon + (1 - p) * E_CZ) >= v_no)
    signal_aq * (1 - p) + (1 - signal_aq) * as.numeric(p * E_PL + (1 - p) * E_CZ > -epsilon)
  }, p_grid, V_no_info)
  
  info_RI <- mapply(function(p, a, b) {
    f_PL1 <- p * a; f_PL0 <- p * (1 - a); f_CZ1 <- (1 - p) * b; f_CZ0 <- (1 - p) * (1 - b)
    f_y1 <- f_PL1 + f_CZ1; f_y0 <- f_PL0 + f_CZ0
    I <- 0
    if (f_PL1 > 0) I <- I + f_PL1 * safe_log(f_PL1 / (p * f_y1)); if (f_PL0 > 0) I <- I + f_PL0 * safe_log(f_PL0 / (p * f_y0))
    if (f_CZ1 > 0) I <- I + f_CZ1 * safe_log(f_CZ1 / ((1 - p) * f_y1)); if (f_CZ0 > 0) I <- I + f_CZ0 * safe_log(f_CZ0 / ((1 - p) * f_y0))
    return(I)
  }, p_grid, a_star_RI, b_star_RI)
  
  info_perfect <- mapply(function(p, v_no) as.numeric(-lambda * H(p) + (-p * epsilon + (1 - p) * E_CZ) >= v_no) * H(p), p_grid, V_no_info)
  
  discrim_RI <- b_star_RI - a_star_RI
  discrim_perfect <- mapply(function(p, v_no) as.numeric(-lambda * H(p) + (-p * epsilon + (1 - p) * E_CZ) >= v_no), p_grid, V_no_info)
  
  df <- bind_rows(
    tibble(p = p_grid, value = V_no_info, type = "No info", discrim = 0, info = 0, buy = buy_no_info),
    tibble(p = p_grid, value = V_full_info, type = "Full info", discrim = 1, info = H(p_grid), buy = buy_full_info),
    tibble(p = p_grid, value = V_perfect, type = "Perfect signal", discrim = discrim_perfect, info = info_perfect, buy = buy_perfect),
    tibble(p = p_grid, value = V_RI, type = "RI optimal", discrim = discrim_RI, info = info_RI, buy = buy_RI)
  )
  
  star_df <- tibble(p = p_grid, a = a_star_RI, b = b_star_RI)
  
  list(df = df, star = star_df)
}

# 3.2 Varying lambda
compute_curves_lambda <- function(p, gamma, sigma_pl, epsilon, lambda_grid) {
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2); E_CZ <- -exp(gamma^2 / 2)
  H_p <- H(p)
  
  solve_res <- lapply(lambda_grid, function(lambda) solve_ab(p, gamma, sigma_pl, lambda, epsilon))
  a_star <- sapply(solve_res, function(res) res$par[1]); b_star <- sapply(solve_res, function(res) res$par[2])
  
  V_no <- rep(expected_util(p, gamma, sigma_pl, epsilon), length(lambda_grid))
  V_full <- rep(-p * epsilon + (1 - p) * E_CZ, length(lambda_grid))
  V_perf <- pmax(-lambda_grid * H_p + V_full, V_no)
  V_ri <- sapply(solve_res, function(res) -res$value)
  
  buy_no <- rep(as.numeric(p * E_PL + (1 - p) * E_CZ > -epsilon), length(lambda_grid))
  buy_full <- rep(1 - p, length(lambda_grid))
  buy_perf <- ifelse(-lambda_grid * H_p + V_full >= V_no, 1-p, buy_no)
  buy_ri <- p * a_star + (1 - p) * b_star
  
  discrim_ri <- b_star - a_star
  discrim_perf <- ifelse(-lambda_grid * H_p + V_full >= V_no, 1, 0)
  
  info_ri <- mapply(function(a, b) {
    f_PL1 <- p * a; f_PL0 <- p * (1 - a); f_CZ1 <- (1 - p) * b; f_CZ0 <- (1 - p) * (1 - b)
    f_y1 <- f_PL1 + f_CZ1; f_y0 <- f_PL0 + f_CZ0
    I <- 0
    if (f_PL1 > 0) I <- I + f_PL1 * safe_log(f_PL1 / (p * f_y1)); if (f_PL0 > 0) I <- I + f_PL0 * safe_log(f_PL0 / (p * f_y0))
    if (f_CZ1 > 0) I <- I + f_CZ1 * safe_log(f_CZ1 / ((1 - p) * f_y1)); if (f_CZ0 > 0) I <- I + f_CZ0 * safe_log(f_CZ0 / ((1 - p) * f_y0))
    I
  }, a_star, b_star)
  info_perf <- ifelse(-lambda_grid * H_p + V_full >= V_no, H_p, 0)
  
  bind_rows(
    tibble(lambda = lambda_grid, type = "No info", value = V_no, buy = buy_no, discrim = 0, info = 0),
    tibble(lambda = lambda_grid, type = "Full info", value = V_full, buy = buy_full, discrim = 1, info = H_p),
    tibble(lambda = lambda_grid, type = "Perfect signal", value = V_perf, buy = buy_perf, discrim = discrim_perf, info = info_perf),
    tibble(lambda = lambda_grid, type = "RI optimal", value = V_ri, buy = buy_ri, discrim = discrim_ri, info = info_ri)
  )
}

# 3.3 Varying sigma_pl
compute_curves_sigma <- function(p, gamma, lambda, epsilon, sigma_grid) {
  H_p <- H(p)
  
  results <- lapply(sigma_grid, function(sigma_pl) {
    E_PL <- -exp(gamma^2 * sigma_pl^2 / 2); E_CZ <- -exp(gamma^2 / 2)
    V_no <- expected_util(p, gamma, sigma_pl, epsilon)
    V_full <- -p * epsilon + (1 - p) * E_CZ
    V_perf <- max(-lambda * H_p + V_full, V_no)
    
    res <- solve_ab(p, gamma, sigma_pl, lambda, epsilon)
    a <- res$par[1]; b <- res$par[2]; V_ri <- -res$value
    buy_ri <- p * a + (1 - p) * b; discrim <- b - a
    
    info <- {
      f_PL1 <- p * a; f_PL0 <- p * (1 - a); f_CZ1 <- (1 - p) * b; f_CZ0 <- (1 - p) * (1 - b)
      f_y1 <- f_PL1 + f_CZ1; f_y0 <- f_PL0 + f_CZ0
      I <- 0; if (f_PL1 > 0) I <- I + f_PL1 * safe_log(f_PL1 / (p * f_y1)); if (f_PL0 > 0) I <- I + f_PL0 * safe_log(f_PL0 / (p * f_y0))
      if (f_CZ1 > 0) I <- I + f_CZ1 * safe_log(f_CZ1 / ((1 - p) * f_y1)); if (f_CZ0 > 0) I <- I + f_CZ0 * safe_log(f_CZ0 / ((1 - p) * f_y0))
      I
    }
    
    bind_rows(
      tibble(sigma_pl = sigma_pl, type = "No info", value = V_no, buy = as.numeric(p * E_PL + (1 - p) * E_CZ > -epsilon), discrim = 0, info = 0),
      tibble(sigma_pl = sigma_pl, type = "Full info", value = V_full, buy = 1 - p, discrim = 1, info = H_p),
      tibble(sigma_pl = sigma_pl, type = "Perfect signal", value = V_perf, buy = ifelse(-lambda * H_p + V_full >= V_no, 1 - p, as.numeric(p * E_PL + (1 - p) * E_CZ > -epsilon)), discrim = ifelse(-lambda * H_p + V_full >= V_no, 1, 0), info = ifelse(-lambda * H_p + V_full >= V_no, H_p, 0)),
      tibble(sigma_pl = sigma_pl, type = "RI optimal", value = V_ri, buy = buy_ri, discrim = discrim, info = info)
    )
  })
  bind_rows(results)
}

# 3.4 Varying gamma
compute_curves_gamma <- function(p, sigma_pl, lambda, epsilon, gamma_grid) {
  H_p <- H(p)
  
  results <- lapply(gamma_grid, function(gamma) {
    E_PL <- -exp(gamma^2 * sigma_pl^2 / 2); E_CZ <- -exp(gamma^2 / 2)
    V_no <- expected_util(p, gamma, sigma_pl, epsilon)
    V_full <- -p * epsilon + (1 - p) * E_CZ
    V_perf <- max(-lambda * H_p + V_full, V_no)
    
    res <- solve_ab(p, gamma, sigma_pl, lambda, epsilon)
    a <- res$par[1]; b <- res$par[2]; V_ri <- -res$value
    buy_ri <- p * a + (1 - p) * b; discrim <- b - a
    
    info <- {
      f_PL1 <- p * a; f_PL0 <- p * (1 - a); f_CZ1 <- (1 - p) * b; f_CZ0 <- (1 - p) * (1 - b)
      f_y1 <- f_PL1 + f_CZ1; f_y0 <- f_PL0 + f_CZ0
      I <- 0; if (f_PL1 > 0) I <- I + f_PL1 * safe_log(f_PL1 / (p * f_y1)); if (f_PL0 > 0) I <- I + f_PL0 * safe_log(f_PL0 / (p * f_y0))
      if (f_CZ1 > 0) I <- I + f_CZ1 * safe_log(f_CZ1 / ((1 - p) * f_y1)); if (f_CZ0 > 0) I <- I + f_CZ0 * safe_log(f_CZ0 / ((1 - p) * f_y0))
      I
    }
    
    bind_rows(
      tibble(gamma = gamma, type = "No info", value = V_no, buy = as.numeric(p * E_PL + (1 - p) * E_CZ > -epsilon), discrim = 0, info = 0),
      tibble(gamma = gamma, type = "Full info", value = V_full, buy = 1 - p, discrim = 1, info = H_p),
      tibble(gamma = gamma, type = "Perfect signal", value = V_perf, buy = ifelse(-lambda * H_p + V_full >= V_no, 1 - p, as.numeric(p * E_PL + (1 - p) * E_CZ > -epsilon)), discrim = ifelse(-lambda * H_p + V_full >= V_no, 1, 0), info = ifelse(-lambda * H_p + V_full >= V_no, H_p, 0)),
      tibble(gamma = gamma, type = "RI optimal", value = V_ri, buy = buy_ri, discrim = discrim, info = info)
    )
  })
  
  bind_rows(results)
}


# --- 4. UI Definition ---
ui <- fluidPage(
  titlePanel("Integrated Analysis of Information Acquisition"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Select X Variable"),
      selectInput("main_var", "Vary by:",
                  choices = list("Prior p" = "p",
                                 "Info Cost λ" = "lambda",
                                 "PL Quality Variance σ_PL" = "sigma",
                                 "Risk Aversion γ" = "gamma"),
                  selected = "p"),
      
      hr(),
      h4("2. Set Fixed Parameters"),
      
      conditionalPanel(
        condition = "input.main_var !== 'p'",
        sliderInput("p_in", "Prior p (prob PL)", 0.01, 0.99, 0.5, 0.01)
      ),
      conditionalPanel(
        condition = "input.main_var !== 'lambda'",
        sliderInput("lambda_in", "λ (info cost)", 0, 5, 1, 0.1)
      ),
      conditionalPanel(
        condition = "input.main_var !== 'sigma'",
        sliderInput("sigma_pl_in", "σ_PL (PL quality variance)", 1, 5, 2, 0.1)
      ),
      conditionalPanel(
        condition = "input.main_var !== 'gamma'",
        sliderInput("gamma_in", "γ (risk aversion)", 0.1, 2, 0.6, 0.1)
      ),
      
      uiOutput("epsilon_slider_ui"),
      
      hr(),
      h4("3. Select Plots to Display"),
      checkboxGroupInput("selected_plots", NULL,
                         choices = c("Expected Utility", 
                                     "Probability of Buying",
                                     "Behavioral Discrimination", 
                                     "Information Acquired"),
                         selected = c("Expected Utility", "Probability of Buying",
                                      "Behavioral Discrimination", "Information Acquired"))
    ),
    
    mainPanel(
      fluidRow(
        column(6, 
               # Use a conditional panel to show/hide this plot
               conditionalPanel(
                 condition = "input.selected_plots.includes('Expected Utility')",
                 plotOutput("utilPlot")
               )
        ),
        column(6, 
               conditionalPanel(
                 condition = "input.selected_plots.includes('Probability of Buying')",
                 plotOutput("buyPlot")
               )
        )
      ),
      fluidRow(
        column(6, 
               conditionalPanel(
                 condition = "input.selected_plots.includes('Behavioral Discrimination')",
                 plotOutput("discrimPlot")
               )
        ),
        column(6, 
               conditionalPanel(
                 condition = "input.selected_plots.includes('Information Acquired')",
                 plotOutput("infoPlot")
               )
        ))
  
  #     fluidRow(
  #       column(6,
  #              # This panel depends on BOTH the main variable AND the checkbox
  #              conditionalPanel(
  #                condition = FALSE,
  #                plotOutput("astarplot")
  #              )
  #       ),
  #       column(6,
  #              conditionalPanel(
  #                condition = FALSE,
  #                plotOutput("bstarplot")
  #              )
  #       )
  #     )
   )
   )
)


# --- 5. Server Logic ---
server <- function(input, output, session) {
  
  epsilon_params <- reactive({
    gamma <- if(input$main_var == "gamma") 0.6 else input$gamma_in
    sigma_pl <- if(input$main_var == "sigma") 2 else input$sigma_pl_in
    req(gamma, sigma_pl)
    
    E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
    E_CZ <- -exp(gamma^2 / 2)
    
    default_val <- abs(E_PL + E_CZ) / 2
    min_val <- min(-E_PL, -E_CZ)
    max_val <- max(-E_PL, -E_CZ)
    
    list(
      min = round(min_val, 1),
      max = round(max_val, 1) + 0.1,
      value = round(default_val, 1)
    )
  })
  
  output$epsilon_slider_ui <- renderUI({
    params <- epsilon_params()
    sliderInput("epsilon_in", "ε (outside option)",
                min = params$min,
                max = params$max,
                value = params$value,
                step = 0.1)
  })
  
  all_curves <- reactive({
    req(input$epsilon_in) 
    
    p_grid      <- seq(0.001, 0.999, length.out = 100)
    lambda_grid <- seq(0, 5, length.out = 100)
    sigma_grid  <- seq(1, 5, length.out = 100)
    gamma_grid  <- seq(0.1, 2, length.out = 100)
    
    switch(input$main_var,
           "p" = compute_curves_p(input$sigma_pl_in, input$gamma_in, input$lambda_in, input$epsilon_in, n_grid = 100),
           "lambda" = list(df = compute_curves_lambda(input$p_in, input$gamma_in, input$sigma_pl_in, input$epsilon_in, lambda_grid)),
           "sigma" = list(df = compute_curves_sigma(input$p_in, input$gamma_in, input$lambda_in, input$epsilon_in, sigma_grid)),
           "gamma" = list(df = compute_curves_gamma(input$p_in, input$sigma_pl_in, input$lambda_in, input$epsilon_in, gamma_grid))
    )
  })
  
  plot_aesthetics <- reactive({
    switch(input$main_var,
           "p"      = list(x_var = "p", x_lab = "Prior p (prob PL)"),
           "lambda" = list(x_var = "lambda", x_lab = "λ (info cost)"),
           "sigma"  = list(x_var = "sigma_pl", x_lab = "σ_PL (PL quality variance)"),
           "gamma"  = list(x_var = "gamma", x_lab = "γ (risk aversion)")
    )
  })
  
  render_standard_plot <- function(y_var, y_lab, title) {
    renderPlot({
      aes_list <- plot_aesthetics()
      df <- all_curves()$df
      ggplot(df, aes_string(x = aes_list$x_var, y = y_var, color = "type")) +
        geom_line(linewidth = 1.1, alpha = 0.8) +
        labs(title = title, x = aes_list$x_lab, y = y_lab, color = "Scenario") +
        theme_minimal(base_size = 14) +
        theme(legend.position = "bottom")
    })
  }
  
  output$utilPlot    <- render_standard_plot("value", "Expected Utility", "Expected Utility")
  output$buyPlot     <- render_standard_plot("buy", "P(buy)", "Probability of Buying")
  output$discrimPlot <- render_standard_plot("discrim", "P(buy|CZ) - P(buy|PL)", "Behavioral Discrimination")
  output$infoPlot    <- render_standard_plot("info", "Mutual Information", "Information Acquired")
  
  # The renderPlot functions are now always available, but the UI controls their visibility.
  output$astarplot <- renderPlot({
    req(input$main_var == "p")
    ggplot(all_curves()$star, aes(x = p, y = a)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Prior p", y = "a*", title = "Optimal a = P(buy | PL)") +
      theme_minimal(base_size = 14)
  })
  
  output$bstarplot <- renderPlot({
    req(input$main_var == "p")
    ggplot(all_curves()$star, aes(x = p, y = b)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Prior p", y = "b*", title = "Optimal b = P(buy | CZ)") +
      theme_minimal(base_size = 14)
  })
  
}

# --- 6. Run the App ---
shinyApp(ui, server)

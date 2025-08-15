library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- Helper Functions ---
pB <- function(w, lambda) {
  threshold <- lambda * log(cosh(1 / lambda))
  ifelse(abs(w) >= threshold,
         ifelse(w > 0, 0, 1),
         (cosh(1 / lambda) - exp(w / lambda)) /
           (2 * cosh(1 / lambda) - (exp(w / lambda) + exp(-w / lambda))))
}

u_fun <- function(w, lambda) {
  C  <- cosh(1 / lambda)
  X  <- exp(w / lambda)
  kc <- (X * C - 1) / (X * (C - X))
  threshold <- lambda * log(cosh(1 / lambda))
  
  EU_interior <- 0.5 * ((1 - w) / (1 + kc * exp(-(1 - w) / lambda)) -
                          (1 + w) / (1 + kc * exp((1 + w) / lambda)))
  
  ifelse(abs(w) >= threshold,
         ifelse(w > 0, 0, -w),
         EU_interior)
}

I_RI <- function(w, lambda) {
  C  <- cosh(1 / lambda)
  X  <- exp(w / lambda)
  threshold <- lambda * log(cosh(1 / lambda))
  
  pB_val <- (C - X) / (2 * C - (X + 1/X))
  
  clamp <- function(z) pmin(pmax(z, 1e-12), 1 - 1e-12)
  pB_val <- clamp(pB_val)
  
  qH <- 1 / (1 + ((1 - pB_val) / pB_val) * exp(-(1 - w) / lambda))
  qL <- 1 / (1 + ((1 - pB_val) / pB_val) * exp((1 + w) / lambda))
  
  qH <- clamp(qH)
  qL <- clamp(qL)
  
  I_interior <- 0.5 * (
    qH * log(qH / pB_val) + (1 - qH) * log((1 - qH) / (1 - pB_val)) +
      qL * log(qL / pB_val) + (1 - qL) * log((1 - qL) / (1 - pB_val))
  )
  
  ifelse(abs(w) >= threshold, 0, I_interior)
}

V_RI <- function(w, lambda) {
  u_fun(w, lambda) - lambda * I_RI(w, lambda)
}

# Conditional probabilities
pB_H <- function(w, lambda) {
  C  <- cosh(1 / lambda)
  X  <- exp(w / lambda)
  pB_val <- (C - X) / (2 * C - (X + 1/X))
  clamp <- function(z) pmin(pmax(z, 1e-12), 1 - 1e-12)
  pB_val <- clamp(pB_val)
  qH <- 1 / (1 + ((1 - pB_val) / pB_val) * exp(-(1 - w) / lambda))
  clamp(qH)
}

pB_L <- function(w, lambda) {
  C  <- cosh(1 / lambda)
  X  <- exp(w / lambda)
  pB_val <- (C - X) / (2 * C - (X + 1/X))
  clamp <- function(z) pmin(pmax(z, 1e-12), 1 - 1e-12)
  pB_val <- clamp(pB_val)
  qL <- 1 / (1 + ((1 - pB_val) / pB_val) * exp((1 + w) / lambda))
  clamp(qL)
}

# Mapping definition labels
def_labels <- c(
  "Fee + info",
  "User",
  "Social planner"
)

# Objective definitions
a_fun <- function(w, lambda, def) {
  if (def == "Fee + info") {
    return(I_RI(w, lambda) * lambda)
  } else if (def == "User") {
    return(V_RI(w, lambda))
  } else if (def == "Social planner") {
    return(V_RI(w, lambda) + pB(w, lambda) * w + I_RI(w, lambda) * lambda)
  }
}

# --- UI ---
ui <- fluidPage(
  titlePanel("Optimization over λ for Given w"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("w", "Choose w", min = 0, max = 1, value = 0.5, step = 0.1),
      sliderInput("lambda_range", "Lambda range", min = 0.01, max = 2, value = c(0.01, 1), step = 0.01),
      selectInput("profit_def", "Profit Definition",
                  choices = def_labels, selected = def_labels[1]),
      numericInput("res", "Resolution (points)", value = 200, min = 10, max = 500, step = 10)
    ),
    mainPanel(
      h4("Summary Table: Max over λ for all definitions"),
      tableOutput("summary_table"),
      plotOutput("plot_probs"),
      plotOutput("plot_profits")
    )
  )
)

# --- Server ---
server <- function(input, output, session) {
  
  lambda_data <- reactive({
    lambdas <- seq(input$lambda_range[1], input$lambda_range[2], length.out = input$res)
    w <- input$w
    
    profit_vals <- sapply(def_labels, function(def) sapply(lambdas, function(l) a_fun(w, l, def)))
    colnames(profit_vals) <- def_labels
    
    data.frame(
      lambda = lambdas,
      pB = pB(w, lambdas),
      pB_H = pB_H(w, lambdas),
      pB_L = pB_L(w, lambdas),
      profit_vals,
      check.names = FALSE  # keep spaces in names
    )
  })
  
  # Table: λ* for all definitions
  output$summary_table <- renderTable({
    df <- lambda_data()
    w <- input$w
    
    results_list <- lapply(def_labels, function(chosen_def) {
      vals <- df[[chosen_def]]
      if (length(vals) == 0) return(NULL) # Handle cases where a column might be missing
      i_max <- which.max(vals)
      l_opt <- df$lambda[i_max]
      
      data.frame(
        Definition = chosen_def,
        lambda_opt = l_opt,
        profit_max = a_fun(w, l_opt, chosen_def),
        pB = pB(w, l_opt),
        pB_H = pB_H(w, l_opt),
        pB_L = pB_L(w, l_opt),
        u = u_fun(w, l_opt),
        I_RI = I_RI(w, l_opt),
        V_RI = V_RI(w, l_opt),
        stringsAsFactors = FALSE
      )
    })
    
    # Combine the list of data frames into a single data frame
    do.call(rbind, results_list)
    
  }, digits = 4)
  
  output$plot_probs <- renderPlot({
    df <- lambda_data()
    ggplot(df, aes(x = lambda)) +
      geom_line(aes(y = pB, color = "pB (overall)"), linewidth = 1) +
      geom_line(aes(y = pB_H, color = "pB | H"), linewidth = 1) +
      geom_line(aes(y = pB_L, color = "pB | L"), linewidth = 1) +
      labs(y = "Probability", color = "Legend",
           title = paste("Buying Probabilities for w =", input$w)) +
      theme_minimal()
  })
  
  output$plot_profits <- renderPlot({
    df <- lambda_data()
    df_long <- pivot_longer(df, cols = all_of(def_labels), names_to = "Definition", values_to = "Profit")
    ggplot(df_long, aes(x = lambda, y = Profit, color = Definition)) +
      geom_line(linewidth = 1) +
      labs(y = "Profit", color = "Definition",
           title = paste("Profit for all definitions, w =", input$w)) +
      theme_minimal()
  })
}

shinyApp(ui, server)
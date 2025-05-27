## app.R
# Rational Inattention Plotter
# --------------------------------------------------
# Adjust the four sliders to see how the expected‑utility
# curves change under different parameter values.
# --------------------------------------------------

# ---- libraries ----
library(shiny)
library(ggplot2)

# ---- helper functions (pure R, no reactives here) ----
compute_curves <- function(sigma_pl, gamma, epsilon, lambda, n_grid = 100) {
  # Expected utilities
  E_PL <- -exp(gamma^2 * sigma_pl^2 / 2)
  E_CZ <- -exp(gamma^2 / 2)
  
  # Prior grid
  p_grid <- seq(0.01, 0.99, length.out = n_grid)
  
 
  expected_util <- function(r) {
    # stage-2 utility after observing posterior belief r
    r*pmax(r * E_PL,-epsilon) + (1-r)*pmax((1 - r) * E_CZ, -epsilon)
  }
  
  info_cost <- function(r, p) {
    if (r <= 0 || r >= 1 || p <= 0 || p >= 1) {
      return(Inf)  # high cost = strongly penalized = avoided
    }
    # per-posterior Shannon cost used earlier
    -lambda * (log((r^(r) / p^p)) -log((1-r)^(1-r) / (1-p)^(1-p))) }
  
  V_RI <- function(r, p) {
    expected_util(r) -info_cost(r, p)}
  
  # Solve r*(p) by maximise V_RI over r in (0,1)
  solve_r_star <- function(p) {
    opt <- optimize(V_RI, interval = c(1e-6, 1 - 1e-6), p = p, maximum = TRUE)
    opt$maximum
  }
  
  r_star <- vapply(p_grid, solve_r_star, numeric(1))
  
  V_no_info   <- pmax(p_grid * E_PL + (1 - p_grid) * E_CZ, -epsilon)  
  V_full_info <- -p_grid * epsilon + (1 - p_grid) * (E_CZ)            
  H_p         <- -p_grid * log(p_grid) - (1 - p_grid) * log(1 - p_grid)
  V_perfect   <- pmax(-lambda * H_p + (-p_grid * epsilon + (1 - p_grid) * E_CZ), V_no_info) 
  V_ri        <- mapply(function(r, p) V_RI(r, p), r_star, p_grid)     
  
  # Assemble tidy data frame
  df <- data.frame(
    p = rep(p_grid, 4),
    value = c(V_no_info, V_full_info, V_perfect, V_ri),
    type = rep(c("No info",
                 "Full info",
                 "Perfect signal",
                 "RI optimal"), each = length(p_grid))
  ) 

}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Rational Inattention – Expected Utility Curves"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("sigma_pl", "σ_PL (variance of PL quality)", min = 1, max = 5, value = 3, step = 0.1),
      sliderInput("gamma",     "γ (risk‑aversion)",            min = 0.1, max = 2, value = 0.6, step = 0.1),
      sliderInput("epsilon",   "ε (penalty for not buying)",   min = E_PL,   max = E_CZ, value = 3, step = 0.1),
      sliderInput("lambda",    "λ (info cost)",                min = 0,   max = 5, value = 3, step = 0.1)
    ),
    mainPanel(
      plotOutput("riPlot", height = "550px")
    )
  )
)

# ---- SERVER ----
server <- function(input, output, session) {
  
  curves <- reactive({
    compute_curves(input$sigma_pl, input$gamma, input$epsilon, input$lambda)
  })
  
  output$riPlot <- renderPlot({
    dat <- curves()
    df  <- dat$df
    ggplot(df, aes(x = p, y = value, colour = type)) +
      geom_line(linewidth = 1) +
      labs(x = "Prior p (probability PL)", y = "Value", colour = "Scenario") +
      ggtitle("Expected Utility under Different Information Structures") +
      theme_minimal() +
      theme(legend.position = "bottom")
  }, res = 96)
}

# ---- Run the app ----
shinyApp(ui = ui, server = server)

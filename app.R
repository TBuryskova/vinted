library(shiny)
library(ggplot2)
library(dplyr)

# --- Helper Functions ---
pB <- function(x, y) {
  threshold <- y * log(cosh(1 / y))
  ifelse(abs(x) >= threshold,
         ifelse(x > 0, 0, 1),
         (cosh(1 / y) - exp(x / y)) /
           (2 * cosh(1 / y) - (exp(x / y) + exp(-x / y))))
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
  
  # Clamp to avoid NaNs
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

# a(x,y) definitions
a_fun <- function(x, y, def) {
  if (def == "def1") {
    return(pB(x, y) * x)
  } else if (def == "def2") {
    return(I_RI(x, y) * y)
  } else if (def == "def3") {
    return(pB(x, y) * x + I_RI(x, y) * y)
  }    else if (def == "def4") {
    return(V_RI(x,y))
  }   else if (def == "def5") {
    return(V_RI(x,y)+pB(x, y) * x + I_RI(x, y) * y)
  }   
}

# --- UI ---
ui <- fluidPage(
  titlePanel("Different definitions of a(x,y)"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("xrange", "X range", min = 0, max = 2, value = c(0, 1), step = 0.05),
      sliderInput("yrange", "Y range", min = 0.01, max = 2, value = c(0.01, 1), step = 0.05),
      numericInput("res", "Resolution (grid points per axis)", value = 50, min = 10, max = 200, step = 5),
      selectInput("a_def", "Choose a(x,y) definition",
                  choices = list(
                    "pB * w" = "def1",
                    "lambda * I_RI" = "def2",
                    "pB * w + lambda * I_RI" = "def3",
                    "user" = "def4",
                    "social planner" = "def5"
                  ),
                  selected = "def1")
    ),
    mainPanel(
      tableOutput("summary_table"),
      plotOutput("plot_a"),
      verbatimTextOutput("opt_results"),
      h4("Summary Table for All Definitions"),
    )
  )
)

# --- Server ---
server <- function(input, output, session) {
  
  grid_data <- reactive({
    xs <- seq(input$xrange[1], input$xrange[2], length.out = input$res)
    ys <- seq(input$yrange[1], input$yrange[2], length.out = input$res)
    expand.grid(x = xs, y = ys) %>%
      mutate(a_val = mapply(a_fun, x, y, MoreArgs = list(def = input$a_def)))
  })
  
  output$plot_a <- renderPlot({
    ggplot(grid_data(), aes(x = x, y = y, fill = a_val)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma") +
      labs(fill = "a(x,y)", title = paste("Heatmap of a(x,y) -", input$a_def)) +
      theme_minimal()
  })
  
  output$opt_results <- renderPrint({
    df <- grid_data()
    opt_row <- df[which.max(df$a_val), ]
    w <- opt_row$x
    lambda <- opt_row$y
    a_val <- a_fun(w, lambda, input$a_def)
    
    list(
      chosen_definition = input$a_def,
      optimal_point = c(x = w, y = lambda),
      a_xy = a_val,
      pB_xy = pB(w, lambda),
      u_xy = u_fun(w, lambda),
      I_xy = I_RI(w, lambda),
      V_xy = V_RI(w, lambda)
    )
  })
  
  # Summary table for all definitions
  output$summary_table <- renderTable({
    defs <- c("def1", "def2", "def3","def4","def5")
    xs <- seq(input$xrange[1], input$xrange[2], length.out = input$res)
    ys <- seq(input$yrange[1], input$yrange[2], length.out = input$res)
    grid <- expand.grid(x = xs, y = ys)
    
    res_list <- lapply(defs, function(def) {
      grid$a_val <- mapply(a_fun, grid$x, grid$y, MoreArgs = list(def = def))
      opt_row <- grid[which.max(grid$a_val), ]
      w <- opt_row$x
      lambda <- opt_row$y
      data.frame(
        Definition = def,
        w_opt = w,
        lambda_opt = lambda,
        a_xy = a_fun(w, lambda, def),
        pB_xy = pB(w, lambda),
        u_xy = u_fun(w, lambda),
        I_xy = I_RI(w, lambda),
        V_xy = V_RI(w, lambda)
      )
    })
    
    do.call(rbind, res_list)
  }, digits = 4)
}

shinyApp(ui, server)

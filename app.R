library(shiny)
library(ggplot2)
library(dplyr)

# --- Helper Functions ---
pB <- function(x, y) {
  if(y==0) {
    if(x>=1) {
  return(0)} else {
    return((1-x)*0.5)
  }
  } else{
  threshold <- y * log(cosh(1 / y))
  pB <- ifelse(abs(x) >= threshold,
         ifelse(x > 0, 0, 1),
         (cosh(1 / y) - exp(x / y)) /
           (2 * cosh(1 / y) - (exp(x / y) + exp(-x / y))))
  return(pB)
  }
}

# pB_H function: calculates qH, the conditional probability of a positive outcome given the observer has a high observation.
pB_H <- function(w, lambda) {
  C  <- cosh(1 / lambda)
  X  <- exp(w / lambda)
  pB_val <- (C - X) / (2 * C - (X + 1/X))
  # Clamp to avoid numerical issues
  clamp <- function(z) pmin(pmax(z, 1e-12), 1 - 1e-12)
  pB_val <- clamp(pB_val)
  qH <- 1 / (1 + ((1 - pB_val) / pB_val) * exp(-(1 - w) / lambda))
  min(clamp(qH),1,na.rm=TRUE)
}

# pB_L function: calculates qL, the conditional probability of a positive outcome given the observer has a low observation.
pB_L <- function(w, lambda) {
  C  <- cosh(1 / lambda)
  X  <- exp(w / lambda)
  pB_val <- (C - X) / (2 * C - (X + 1/X))
  # Clamp to avoid numerical issues
  clamp <- function(z) pmin(pmax(z, 1e-12), 1 - 1e-12)
  pB_val <- clamp(pB_val)
  qL <- 1 / (1 + ((1 - pB_val) / pB_val) * exp((1 + w) / lambda))
  max(clamp(qL),0,na.rm=TRUE)
}

# u_fun: a utility function, now with correct handling of lambda = 0
u_fun <- function(w, lambda) {
  if (lambda == 0) {
    # Special case for lambda = 0
    return(0.5 * (1 - w))
  } else {
    # Calculations for lambda != 0
    C <- cosh(1 / lambda)
    X <- exp(w / lambda)
    kc <- (X * C - 1) / (X * (C - X))
    threshold <- lambda * log(cosh(1 / lambda))
    
    EU_interior <- 0.5 * ((1 - w) / (1 + kc * exp(-(1 - w) / lambda)) -
                            (1 + w) / (1 + kc * exp((1 + w) / lambda)))
    
    if (abs(w) >= threshold) {
      if (w > 0) {
        return(0)
      } else {
        return(-w)
      }
    } else {
      return(EU_interior)
    }
  }
}

# I_RI: a function that calculates information rent
I_RI <- function(w, lambda) {
  if (lambda == 0) {
    # Special case for lambda = 0
    return(0.5 * (1 - w))
  } else {
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
    
    I_interior <- ifelse(abs(w) >= threshold, 0, I_interior)
    return(I_interior)
  }
}

# V_RI: a function for user's value
V_RI <- function(w, lambda) {
  u_fun(w, lambda) - lambda * I_RI(w, lambda)
}

# a_fun: a function to choose which a(x,y) definition to use
a_fun <- function(x, y, def) {
  if (def == "def1") {
    if(y==0) {
      return(0.5*x)
    } else {
      return(pB(x, y) * x + I_RI(x, y) * y )
    }
  }  else if (def == "def4") {
    return(V_RI(x,y))
  } else if (def == "def5") {
    return(V_RI(x,y)+pB(x, y) * x + I_RI(x, y) * y)
  }
}

# --- UI ---
ui <- fluidPage(
  titlePanel("Different definitions of a(x,y)"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("xrange", "X range", min = 0, max = 2, value = c(0, 1), step = 0.05),
      sliderInput("yrange", "Y range", min = 0, max = 2, value = c(0, 1), step = 0.05),
      numericInput("res", "Resolution (grid points per axis)", value = 200, min = 10, max = 200, step = 5),
      selectInput("a_def", "Choose a(x,y) definition",
                  choices = list(
                    "profit" = "def1",
                    "user" = "def4",
                    "social planner" = "def5"
                  ),
                  selected = "def1")
    ),
    mainPanel(
      h4("Summary Table"),
      tableOutput("summary_table"),
      plotOutput("plot_a"),
      plotOutput("plot_buy"),
      plotOutput("plot_util"),
      
      # New plot for welfare heatmap
      plotOutput("plot_welfare")
    )
  )
)

# --- Server ---
server <- function(input, output, session) {
  
  grid_data <- reactive({
    xs <- seq(input$xrange[1], input$xrange[2], length.out = input$res)
    ys <- seq(input$yrange[1], input$yrange[2], length.out = input$res)
    
    # Calculate all necessary values in one go
    grid <- expand.grid(x = xs, y = ys) %>%
      mutate(
        a_val = mapply(a_fun, x, y, MoreArgs = list(def = input$a_def)),
        welfare = mapply(function(w, lambda) {
          V_RI(w, lambda) + pB(w, lambda) * w + I_RI(w, lambda) * lambda
        }, x, y),
        util = mapply(V_RI, x, y),
        pBuy = mapply(pB, x, y)
      )
    
    return(grid)
  })
  
  output$plot_a <- renderPlot({
    ggplot(grid_data(), aes(x = x, y = y, fill = a_val)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma") +
      labs(fill = "profit", x="w", y="lambda", title = paste("Profit as a function of w and lambda")) +
      theme_minimal()
  })
  
  output$plot_buy <- renderPlot({
    df=grid_data()
    
    ggplot(df, aes(x = x, y = y, fill = pBuy)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma") +
      labs(fill = "pBuy", x="w", y="lambda", title = paste("pBuy as a function of w and lambda")) +
      theme_minimal()
  })
  
  output$plot_util <- renderPlot({
    df=grid_data()
    
    ggplot(df, aes(x = x, y = y, fill = util)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma") +
      labs(fill = "util", x="w", y="lambda", title = paste("User utility as a function of w and lambda")) +
      theme_minimal()
  }) 
  
  # New plot output for welfare heatmap
  output$plot_welfare <- renderPlot({
    # Get the data from the reactive expression
    df <- grid_data()

    ggplot(df, aes(x = x, y = y, fill = welfare)) +
      geom_tile() +
      # Add a point to show the maximum welfare location
      geom_point(aes(x = 1, y = 0), color = "red", size = 5, shape = 4, inherit.aes = FALSE) +
      scale_fill_viridis_c(option = "viridis") + # Using a different color scale for visual distinction
      labs(fill = "Welfare", x="w", y="lambda", title = "Welfare as a function of w and lambda") +
      theme_minimal()
  })
  
  # Summary table for all definitions
  output$summary_table <- renderTable({
    defs <- c("def1", "def4", "def5")
    xs <- seq(input$xrange[1], input$xrange[2], length.out = input$res)
    ys <- seq(input$yrange[1], input$yrange[2], length.out = input$res)
    grid <- expand.grid(x = xs, y = ys)
    
    res_list <- lapply(defs, function(def) {
      grid$a_val <- mapply(a_fun, grid$x, grid$y, MoreArgs = list(def = def))
      opt_row <- grid[which.max(grid$a_val), ]
      w <- opt_row$x
      lambda <- opt_row$y
      
      # Calculate welfare for the optimal point
      welfare_val <- V_RI(w,lambda) + pB(w, lambda) * w + I_RI(w, lambda) * lambda
      
      data.frame(
        Definition = def,
        w_opt = w,
        lambda_opt = lambda,
        a_xy = a_fun(w, lambda, def),
        pB_H_xy = pB_H(w, lambda),
        pB_L_xy = pB_L(w, lambda),
        u_xy = u_fun(w, lambda),
        I_xy = I_RI(w, lambda),
        V_xy = V_RI(w, lambda),
        welfare = welfare_val
      )
    })
    
    do.call(rbind, res_list)
  }, digits = 4)
}

shinyApp(ui, server)

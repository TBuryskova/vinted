library(shiny)
library(ggplot2)
library(dplyr)

# --- Helper Functions (same as before) ---
pB <- function(x, y) {
  if(y==0) {
    if(x>=1) {
      return(0)
    } else {
      return((1-x)*0.5)
    }
  } else {
    threshold <- y * log(cosh(1 / y))
    pB <- ifelse(abs(x) >= threshold,
                 ifelse(x > 0, 0, 1),
                 (cosh(1 / y) - exp(x / y)) /
                   (2 * cosh(1 / y) - (exp(x / y) + exp(-x / y))))
    return(pB)
  }
}

pB_H <- function(w, lambda) {
  C  <- cosh(1 / lambda)
  X  <- exp(w / lambda)
  pB_val <- (C - X) / (2 * C - (X + 1/X))
  clamp <- function(z) pmin(pmax(z, 1e-12), 1 - 1e-12)
  pB_val <- clamp(pB_val)
  qH <- 1 / (1 + ((1 - pB_val) / pB_val) * exp(-(1 - w) / lambda))
  min(clamp(qH),1,na.rm=TRUE)
}

pB_L <- function(w, lambda) {
  C  <- cosh(1 / lambda)
  X  <- exp(w / lambda)
  pB_val <- (C - X) / (2 * C - (X + 1/X))
  clamp <- function(z) pmin(pmax(z, 1e-12), 1 - 1e-12)
  pB_val <- clamp(pB_val)
  qL <- 1 / (1 + ((1 - pB_val) / pB_val) * exp((1 + w) / lambda))
  max(clamp(qL),0,na.rm=TRUE)
}

u_fun <- function(w, lambda) {
  if (lambda == 0) {
    return(0.5 * (1 - w))
  } else {
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

I_RI <- function(w, lambda) {
  if (lambda == 0) {
    return(0.5 * (1 - w))
  } else {
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
    
    I_interior <- ifelse(abs(w) >= threshold, 0, I_interior)
    return(I_interior)
  }
}

V_RI <- function(w, lambda) {
  u_fun(w, lambda) - lambda * I_RI(w, lambda)
}

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
      plotOutput("schematic_plot"),
      plotOutput("plot_util"),
      plotOutput("plot_welfare")
    )
  )
)

# --- Server ---
server <- function(input, output, session) {
  
  grid_data <- reactive({
    xs <- seq(input$xrange[1], input$xrange[2], length.out = input$res)
    ys <- seq(input$yrange[1], input$yrange[2], length.out = input$res)
    
    grid <- expand.grid(x = xs, y = ys) %>%
      mutate(
        a_val = mapply(a_fun, x, y, MoreArgs = list(def = input$a_def)),
        welfare = mapply(function(w, lambda) {
          V_RI(w,lambda) + pB(w, lambda) * w + I_RI(w, lambda) * lambda
        }, x, y),
        util = mapply(V_RI, x, y),
        pBuy = mapply(pB, x, y)
      )
    return(grid)
  })
  
  # helper for bw contour plot
  bw_contour_plot <- function(df, zvar, fill_label, title) {
    ggplot(df, aes(x = x, y = y, z = !!sym(zvar))) +
      geom_contour_filled(bins = 15) +
      scale_fill_grey(start = 1, end = 0) +
      labs(fill = fill_label, x="w", y="lambda", title = title) +
      theme_minimal()
  }
  
  output$plot_a <- renderPlot({
    bw_contour_plot(grid_data(), "a_val", "profit", "Profit as a function of w and lambda")
  })
  
  output$plot_buy <- renderPlot({
    bw_contour_plot(grid_data(), "pBuy", "pBuy", "pBuy as a function of w and lambda")
  })
  
  output$plot_util <- renderPlot({
    bw_contour_plot(grid_data(), "util", "util", "User utility as a function of w and lambda")
  }) 
  
  output$plot_welfare <- renderPlot({
    bw_contour_plot(grid_data(), "welfare", "Welfare", "Welfare as a function of w and lambda")
  })
  
  output$schematic_plot <- renderPlot({
    df <- grid_data()
    ggplot(df, aes(x = x, y = y, z = pBuy)) +
      geom_contour_filled(breaks = seq(0, 0.5, by = 0.05)) +
      scale_fill_grey(start = 1, end = 0) +
      labs(title = "Schematic: pBuy contours", x = "w", y = "lambda", fill = "pBuy") +
      theme_minimal()
  })
  
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

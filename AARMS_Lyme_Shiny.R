#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Lyme Disease Model with Acaricide"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput(inputId = "a",
                        label ="Proportion of rodents with acaricide treatment",
                        min = 0.0,
                        max = 0.999,
                        value = 0.2)
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
           #plotOutput("plot")
           textOutput("selected_var")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    mod5 <- function(t, y , parms){
        #Pull state variables from y vector
        L = y[1]
        Ns = y[2]
        Ni = y[3]
        A = y[4]
        Hs = y[5]
        Hi = y[6]
        Hr = y[7]
        
        #Pull the required parameter values from the parms vector
        b_l=parms["b_l"]
        b_n=parms["b_n"]
        r=parms["r"]
        a =parms["a"]
        delta=parms["delta"]
        mu_l=parms["mu_l"]
        mu_n=parms["mu_n"]
        mu_a=parms["mu_a"]
        mu_h=parms["mu_h"]
        beta=parms["beta"]
        gamma=parms["gamma"]
        gamma_h=parms["gamma_h"]
        p_hl=parms["p_hl"]
        p_nh=parms["p_nh"]
        p_hn = parms["p_hn"]
        
        dL = beta*A*exp(-gamma*A) 
        - b_l*L*(1-r) - mu_l*L
        
        dNs = b_l*L*(1-r)*(1-a)*(1- 
                                   (1- (1-delta)^((Ni/((1-r)*(Hs+Hi+Hr)
                                   )
                                   )
                                   )
                                   )
        )*((Hs + (1-p_hl)*Hi + Hr)/(Hs+Hi+Hr)) 
        - b_n*Ns*(1-r) - mu_n*Ns
        
        dNi = b_l*L*(1-r)*(1-a)*p_hl*(Hi/(Hs+Hi+Hr)) + b_l*L*(1-r)*(1-a)*( 
          (1- (1-delta)^((Ni/((1-r)*(Hs+Hi+Hr)
          )
          )
          )
          )
        ) *((Hs + (1-p_hl)*Hi + Hr)/(Hs+Hi+Hr)) - b_n*Ni*(1-r) - mu_n*Ni
        
        dA = b_n*(Ns+Ni)*(1-r)*(1-a) - mu_a*A
        
        dHs = mu_h*(Hs+Hi+Hr) - b_n*p_nh*(1-r)*Ni*(Hs/(Hs+Hi+Hr)) - mu_h*Hs
        
        dHi = b_n*p_nh*(1-r)*Ni*(Hs/(Hs+Hi+Hr)) - gamma_h*Hi - mu_h*Hi
        
        dHr = gamma_h*Hi - mu_h*Hr
        
        res=c(dL, dNs, dNi, dA, dHs, dHi, dHr)
        list(res)
      }
    
    values <- reactive({
      req(input$a)
      times  = seq(0, 1000, by=1)
      parms= c(
        b_l=0.5,
        b_n=0.5,
        r=0,
        delta=0.7,
        mu_l=0.01,
        mu_n=0.002,
        mu_a=0.1,
        mu_h=0.001,
        beta=15,
        gamma=0.0005,
        gamma_h=0.1,
        p_hl=0.8,
        p_nh=0.9,
        p_hn = 0.8
      )
      start = c(L=6000, Ns=1000, Ni=1500, A=800, Hs=1000, Hi=2000, Hr=100)
      out <- ode(y = start, times = times, func = mod5, 
          parms = c(parms, a = input$a))
      out = as.data.frame(out)
    })
    
     output$selected_var <- renderText({ 
       paste("You have selected", out[2,2])
     })
    
    output$plot <- renderPlot({
        with(out, {
        plot(x = out$time, y = out$L, col = "yellow", ylab = "Number of individuals",
             xlab = "Time (days)", type = "l", xlim = c(0, 4000), ylim = c(-10,40000))
        #lines(x = out$time, y = out$L, col = "yellow")
        lines(x = out$time, y = out$Ns, col = "red")
        lines(x = out$time, y = out$A, col = "purple")
        lines(x = out$time, y = out$Hi, col = "magenta")
        lines(x = out$time, y = out$Ni, col = "pink")
        #lines(x = out$time, y = out$Ai, col = "orange")
        lines(x = out$time, y = out$Hr, col = "green")
        lines(x = out$time, y = out$Hs, col = "cornflowerblue")
        })
      })
}


# Run the application 
shinyApp(ui = ui, server = server)

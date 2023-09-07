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
                        value = 0.2),
            sliderInput(inputId = "b_l",
                        label ="Rate at which larval ticks attach on rodents",
                        min = 0.0,
                        max = 0.999,
                        value = 0.5),
            sliderInput(inputId = "b_n",
                        label ="Rate at which nymphal ticks bite rodents",
                        min = 0.0,
                        max = 0.999,
                        value = 0.5),
            sliderInput(inputId = "r",
                        label ="Proportion of rodents with repellent insecticide",
                        min = 0.0,
                        max = 0.999,
                        value = 0),
            sliderInput(inputId = "delta",
                        label ="Probability of co-feeding transmission",
                        min = 0.0,
                        max = 0.999,
                        value = 0.7),
            sliderInput(inputId = "mu_l",
                        label ="Death rate of larval ticks",
                        min = 0.0,
                        max = 0.999,
                        value = 0.01),
            sliderInput(inputId = "mu_n",
                        label ="Death rate of nymphal ticks",
                        min = 0.0,
                        max = 0.999,
                        value = 0.002),
            sliderInput(inputId = "mu_a",
                        label ="Death rate of adult ticks",
                        min = 0.0,
                        max = 0.999,
                        value = 0.1),
            sliderInput(inputId = "mu_h",
                        label ="Death rate of rodents",
                        min = 0.0,
                        max = 0.999,
                        value = 0.001),
            sliderInput(inputId = "beta",
                        label ="Maximum birth rate for larvae",
                        min = 0.0,
                        max = 50,
                        value = 15),
            sliderInput(inputId = "gamma",
                        label ="Intensity of density-dependence in larvae birth rate",
                        min = 0.0,
                        max = 0.999,
                        value = 0.0005),
            sliderInput(inputId = "gamma_h",
                        label ="Recovery rate of hosts",
                        min = 0.0,
                        max = 0.999,
                        value = 0.1),
            sliderInput(inputId = "p_hl",
                        label ="Probability of transmission from host to larva",
                        min = 0.0,
                        max = 0.999,
                        value = 0.8),
            sliderInput(inputId = "p_nh",
                        label ="Probability of transmission from nymph to host",
                        min = 0.0,
                        max = 0.999,
                        value = 0.9),
            sliderInput(inputId = "p_hn",
                        label ="Probability of transmission from host to nymph",
                        min = 0.0,
                        max = 0.999,
                        value = 0.8)
            
        ),

        # Show a plot of the generated distribution
        mainPanel(plotOutput("plot")
        )
           #textOutput("selected_var")
        
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    mod5 <- function(t, vars , parms){
      with(as.list(c(vars, parms)), { 
        #Pull state variables from y vector
        # L = vars[1]
        # Ns = vars[2]
        # Ni = vars[3]
        # A = vars[4]
        # Hs = vars[5]
        # Hi = vars[6]
        # Hr = vars[7]
        
        #Pull the required parameter values from the parms vector
        # b_l=parms["b_l"]
        # b_n=parms["b_n"]
        # r=parms["r"]
        # a =parms["a"]
        # delta=parms["delta"]
        # mu_l=parms["mu_l"]
        # mu_n=parms["mu_n"]
        # mu_a=parms["mu_a"]
        # mu_h=parms["mu_h"]
        # beta=parms["beta"]
        # gamma=parms["gamma"]
        # gamma_h=parms["gamma_h"]
        # p_hl=parms["p_hl"]
        # p_nh=parms["p_nh"]
        # p_hn = parms["p_hn"]
        
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
        return(list(res))
      })
      }
    
    values <- reactive({
      req(input$a, input$b_l, input$b_n, input$r, input$delta, input$mu_l,
          input$mu_n, input$mu_a, input$mu_h, input$beta, input$gamma, 
          input$gamma_h, input$p_hl, input$p_hn, input$p_nh)
      
      times  = seq(0, 5000, by=1)
      parms= c(a = input$a, b_l = input$b_l, b_n = input$b_n, r = input$r, 
               delta = input$delta, mu_l = input$mu_l,
               mu_n = input$mu_n, mu_a = input$mu_a, mu_h = input$mu_h, 
               beta = input$beta, gamma = input$gamma, 
               gamma_h = input$gamma_h, p_hl = input$p_hl, 
               p_hn = input$p_hn, p_nh = input$p_nh)
      start = c(L=6000, Ns=1000, Ni=1500, A=800, Hs=1000, Hi=2000, Hr=100)
      
      ode(y = start, times = times, func = mod5, 
          parms = parms)
      
    })
    
    
    output$plot <- renderPlot({
      out <- as.data.frame(values())
      
        with(out, {
        plot(x = time, y = L, col = "yellow", ylab = "Number of individuals",
             xlab = "Time (days)", type = "l", xlim = c(0, 400), ylim = c(-10,4000))
        lines(x = out$time, y = out$L, col = "yellow")
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

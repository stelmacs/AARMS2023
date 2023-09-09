# 
# Shiny app for SIR model with vaccinations
#
# Created by Claus Ekstrøm 2019
# @ClausEkstrom
#

library("shiny")
library("deSolve")
library("cowplot")
library("ggplot2")
library("tidyverse")
library("ggrepel")
library("shinydashboard")

## Create an SIR function
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
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
    return(list(c(dL, dNs, dNi, dA, dHs, dHi, dHr)))
  })
}


#
# Define UI 
#

ui <- dashboardPage(
  dashboardHeader(disable = TRUE),
  dashboardSidebar(
    # sliderInput("popsize",
    #             "Population size (millions):",
    #             min = 1, max = 300, value = 6
    # ),
    # sliderInput("connum",
    #             "Basic reproductive number (R0, # persons):",
    #             min = .5, max = 20, value = 5
    # ),
    # sliderInput("pinf",
    #             "# infected at outbreak:",
    #             min = 1, max = 50, value = 2
    # ),
    # sliderInput("pvac",
    #             "Proportion vaccinated / immune (%):",
    #             min = 0, max = 100, value = 75
    # ),
    # sliderInput("vaceff",
    #             "Vaccine effectiveness (%):",
    #             min = 0, max = 100, value = 85
    # ),
    # sliderInput("infper",
    #             "Infection period (days):",
    #             min = 1, max = 30, value = 7
    # ),
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
                value = 0.8),
    sliderInput("timeframe",
                "Time frame (days):",
                min = 1, max = 400, value = 200
    )
    
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
                              /* body */
                              .content-wrapper, .right-side {
                              background-color: #fffff8;
                              }                              
                              '))),
    
    #    mainPanel(
    fluidRow(plotOutput("distPlot")),
    br(),
    # fluidRow(
    #   # Dynamic valueBoxes
    #   valueBoxOutput("progressBox", width = 6),
    #   valueBoxOutput("approvalBox", width = 6),
    #   valueBoxOutput("BRRBox", width = 6),
    #   valueBoxOutput("HIBox", width = 6)
    # ),
    br(),
    br()
  )
)

#
# Define server 
#
server <- function(input, output) {
  # Create reactive input
  dataInput <- reactive({
    init       <-
      c(L=6000, Ns=1000, Ni=1500, A=800, Hs=1000, Hi=2000, Hr=100
      )
    parameters <- c(a = input$a, b_l = input$b_l, b_n = input$b_n, r = input$r, 
                    delta = input$delta, mu_l = input$mu_l,
                    mu_n = input$mu_n, mu_a = input$mu_a, mu_h = input$mu_h, 
                    beta = input$beta, gamma = input$gamma, 
                    gamma_h = input$gamma_h, p_hl = input$p_hl, 
                    p_hn = input$p_hn, p_nh = input$p_nh)

    ## Time frame
    times <- seq(0, input$timeframe, by = 1)
    
    ## Solve using ode (General Solver for Ordinary Differential Equations)
    out <- ode(
      y = init,
      times = times,
      func = sir,
      parms = parameters
    )   
    #    out
    as.data.frame(out)
  })
  
  output$distPlot <- renderPlot({
    out <-
      dataInput() %>%
      gather(key, value, -time) %>%
      mutate(
        id = row_number(),
        key2 = recode(
          key,   # dL, dNs, dNi, dA, dHs, dHi, dHr
          L = "Larvae",
          Ns = "Susceptible Nymphs",
          Ni = "Infected Nymphs",
          A = "Adult Ticks",
          Hs = "Susceptible Hosts",
          Hi = "Infected Hosts",
          Hr = "Treated Hosts"
        ),
        keyleft = recode(
          key,
          L = "Larvae",
          Ns = "Susceptible Nymphs",
          Ni = "Infected Nymphs",
          A = "Adult Ticks",
          Hs = "Susceptible Hosts",
          Hi = "Infected Hosts",
          Hr = "Treated Hosts"
        ),
        keyright = recode(
          key,
          L = "Larvae",
          Ns = "Susceptible Nymphs",
          Ni = "Infected Nymphs",
          A = "Adult Ticks",
          Hs = "Susceptible Hosts",
          Hi = "Infected Hosts",
          Hr = "Treated Hosts"
        )
      )
    
    ggplot(data = out,
           aes(
             x = time,
             y = value,
             group = key2,
             col = key2,
             label = key2,
             data_id = id
           )) + # ylim(0, 1) +
      ylab("Number of individuals") + xlab("Time (days)") +
      geom_line(size = 1) +
      geom_text_repel(
        data = subset(out, time == max(time)),
        aes(label = keyright),
        size = 2,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 1,
        direction = "y"
      ) +
      geom_text_repel(
        data = subset(out, time == min(time)),
        aes(label = keyleft),
        size = 6,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 0,
        direction = "y"
      ) +
      theme(legend.position = "none") +
      scale_colour_manual(values = c("yellow","red","purple","magenta","pink","green",
                                     "cornflowerblue")) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 4000)) +
      theme(
        rect=element_rect(size=0),
        legend.position="none",
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.background=element_rect(fill="transparent", colour=NA),
        legend.key = element_rect(fill = "transparent", colour = "transparent")
      )
    
  })
  
  # output$progressBox <- renderValueBox({
  #   valueBox(
  #     dataInput() %>% filter(time == max(time)) %>% select(R) %>% mutate(R = round(100 * R, 2)) %>% paste0("%"),
  #     "Proportion of full population that got the disease by end of time frame",
  #     icon = icon("thumbs-up", lib = "glyphicon"),
  #     color = "black"
  #   )
  # })
  
  # output$approvalBox <- renderValueBox({
  #   valueBox(
  #     paste0(round(
  #       100 * (dataInput() %>% filter(row_number() == n()) %>% mutate(res = (R + I) / (S + I + R)) %>% pull("res")), 2), "%"),
  #     "Proportion of susceptibles that will get the disease by end of time frame",
  #     icon = icon("thermometer-full"),
  #     color = "black"
  #   )
  # })
  
  # output$BRRBox <- renderValueBox({
  #   valueBox(
  #     paste0(round(input$connum *
  #                    (1 - input$pvac / 100 * input$vaceff / 100), 2), ""),
  #     "Effective R0 (for populationen at outbreak, when immunity is taken into account)",
  #     icon = icon("arrows-alt"),
  #     color = "red"
  #   )
  # })
  
  # output$HIBox <- renderValueBox({
  #   valueBox(
  #     paste0(round(100 * (1 - 1 / (input$connum)), 2), "%"),
  #     "Proportion of population that needs to be immune for herd immunity",
  #     icon = icon("medkit"),
  #     color = "blue"
  #   )
  # })  
}

# Run the application
shinyApp(ui = ui, server = server)


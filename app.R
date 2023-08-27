# Original R Code made by Carl Walters, Aug 2023
# modified to Shiny app by Villy Christensen, Aug-Sep 2023
#library(rsconnect)
#rsconnect::deployApp('path/to/your/app')

library(shiny)
library(shinyBS)
#library(shinycssloaders)

# declare global variables here (if really needed)
bls_txt = "Rule of thumb: fish swim with 1 bl/sec when searching, up to 2 bl/sec when in pursuit or actively searching, and 0.5 bl/sec while in feeding mode, so 1 bl/sec is a reasonable choice"
exit_txt = "Note that when exit probability > 0.05/day, we recommend you increase the number of model time steps (per day)"
turn_txt = "Fish tend to change direction after a prey hunt/capture, so how many times an hour does this species do that? Planktivores, often. Piscivores, less often"
home_txt = "Proportion of movement that is directed back toward the individual’s starting point or home range center"
nfish_txt = "More fish (or think of each as a school) = slower calculation. It may not change results noticeably, but try it ... " #Schools move together, so for those it's number of schools rather than fish"
lgt_txt = "Get this from your Ecospace basemap. Use the average of cell length and cell width"

Xend = array(0,1)
Yend = array(0,1)
xone = array(0,1)
yone= array(0,1)
Dmove= array(0,1)

# if updating this list, remember that it has to be done with names/dimensioning and in two places down at bottom
nm <- c('name',"body length","body length CV ","swim speed","swim speed CV" ,"hours active","hours CV" ,"turns", "home base",
        "cell length (km)",'model step' ,"N fish","Diffusivity (D; km^2/day)" ,"Mixing rate (yr-1)", "Ecospace dispersal (km/yr)", "time")
res <- matrix(ncol = length(nm), nrow = 0)
colnames(res) <- nm

# Define UI for application 
ui <- fluidPage(
  
  # Application title
  titlePanel("Ecospace dispersal IBM"),
  
  
  fluidRow(
    column(4,
           wellPanel(
             textInput('name','Species name',""),
             sliderInput("bl", "Body length (cm)", value = 60, min = 5, max = 200),
             numericInput("bls", "Swim speed (bl/sec)", value = 1, min = 0.001, max = 5),
             sliderInput("hours", "Hours active / day:  ", value = 18, min = 1, max = 24),   
             sliderInput("turns", "# turns / active hour", value = 20, min = 1, max = 120),
             numericInput("celllength", "Cell length (km)", value = 110, min = 0.001, max = 1000),
             numericInput("steps", "Model time steps (per day)", value = 1, min = 1, max = 100),
             numericInput("nfish", "Number of fish", value = 1000, min = 1, max = 10000),
             actionButton("batch", "Generate 20 runs"),
             downloadButton('downFile',"Save session output"),
             # downloadButton("downloadPlot", "Download plots"),
             
             bsPopover("bls", 'Swimming speed?', bls_txt, placement = "bottom", trigger = "hover",options = NULL),
             bsPopover("turns", 'Number of turns?', turn_txt, placement = "bottom", trigger = "hover",options = NULL),
             bsPopover("wthome", 'Home base?', home_txt, placement = "bottom", trigger = "hover",options = NULL),
             bsPopover("celllength", 'Ecospace cell length', lgt_txt, placement = "bottom", trigger = "hover",options = NULL),
             bsPopover("nfish", 'Number of fish (or shools)?', nfish_txt, placement = "bottom", trigger = "hover",options = NULL),
             
             br(),
             br(),
             textOutput("distturn"),
             textOutput("Pexit"),
             textOutput("Pexit_warning"),
             textOutput("diffus"),
             textOutput("annmix"),
             #textOutput("swimspeed"),
             #textOutput("mixecospace"),
             #textOutput("meanD"),
             textOutput("movedist"),
             #textOutput("meanDyr"),
             #br(),
             textOutput("time")
             
           )       
    ),
    
    column(4,
           wellPanel(
             numericInput("bl.cv", "Body length (CV)", value = 0.1, min = 0, max = 1),
             numericInput("bls.cv", "Swim speed (CV)", value = 0.1, min = 0, max = 1),
             numericInput("hours.cv", "Hours active (CV)", value = 0.1, min = 0, max = 1),
             sliderInput("t.sd", "Turn st.dev (deg.)", value = 40, min = 0, max = 180),   
             sliderInput("wthome", "Move towards home (prop):  ", value = 0, min = 0, max = 0.50),
             plotOutput("histPlot")
           )
    ),
    
    column(4,
           plotOutput("distPlot"),
           plotOutput("trajPlot")
    )
  )
)



# Define server logic 
server <- function(input, output) {
  
  do.calc = function()    {
    #read input variables to local variables (it's faster!)
    wthome = input$wthome
    celllength = input$celllength
    t.sd = input$t.sd * pi / 180   #convert from degrees to radians
    nfish = input$nfish
    hours = input$hours
    nday = input$turns * hours
    steps = input$steps
    bl.cv = input$bl.cv
    bls.cv= input$bls.cv
    hours.cv = input$hours.cv
    
    swimspeed = input$bl * input$bls * 60*60 * hours /100/1000 /steps  # swimspeed to set below in km/day
    dday = swimspeed / celllength         # distance moved per day in cell lengths
    Dstep1 = dday / nday                  # distance moved per movement step (cell lengths)
    ex=array(0,nfish)
    Dmove <<- array(0,nfish)
    Xend<<-array(0,nfish)
    Yend<<-array(0,nfish)
    xone<<-array(0,nday)
    yone<<-array(0,nday)
    for(ifish in 1:nfish){
      Dstep = Dstep1 * rnorm(1,1,bl.cv)  * rnorm(1,1,bls.cv)  * rnorm(1,1,hours.cv)     
      # body length cv multiplied on    # swim speed cv   # hours active cv
      th = runif(1,0,1) * 2 * pi
      X = runif(1,0,1)
      Y = runif(1,0,1)
      Xs = X
      Ys = Y
      xmhome=0
      ymhome=0
      for(i in 1:nday){
        if(ifish==1) {
          xone[i] <<- X-Xs
          yone[i] <<- Y-Ys
        }
        th = th + rnorm(1,0,t.sd)
        if(wthome>0) {
          X = X + (1-wthome)*Dstep * sin(th)+wthome*xmhome
          Y = Y + (1-wthome)*Dstep * cos(th)+wthome*ymhome
        } else {
          X = X + Dstep * sin(th)
          Y = Y + Dstep * cos(th)
        }
        if (X < 0){ex[ifish] = 1}
        if (X > 1){ex[ifish] = 1}
        if (Y < 0) {ex[ifish] = 1}
        if (Y > 1) {ex[ifish] = 1}
        d=sqrt((X-Xs)^2+(Y-Ys)^2)
        thy=acos((Ys-Y)/d)
        thx=asin((Xs-X)/d)
        #th=wthome*thhome+(1-wthome)*th
        if(wthome>0) {        
          xmhome=Dstep*sin(thx)
          ymhome=Dstep*cos(thy)
        }
      }
      Xend[ifish]<<-(X-Xs)*celllength
      Yend[ifish]<<-(Y-Ys)*celllength
      Dmove[ifish] <<- sqrt((X - Xs) ^ 2 + (Y - Ys) ^ 2)
    }
    Pexit=sum(ex)
    meanD=sum(Dmove) * steps
    distturn=swimspeed*1000/nday * steps
    
    Pexit = Pexit / nfish / 4    #probability of exiting the cell in one day
    meanD = celllength * meanD / input$nfish  #mean distance moved in one day km
    annmix=365*Pexit * steps
    mixecospace=swimspeed*365/(pi*celllength) * steps
    
    data.list = list("a" = distturn, "b" = Pexit, "c" = annmix)
    return(data.list)
  }
  
  prepare_plot = reactive({
    beginning <- Sys.time()
    data.out = do.calc()
    distturn = data.out$a
    Pexit = data.out$b
    annmix = data.out$c
    
    output$distturn <- renderText({paste("Distance between turns:  ",format(distturn,digits=1),"m")})
    output$Pexit <-    renderText({paste("Exit probability per cell face:  ",format(Pexit,digits=4),"per day")})
    if(Pexit>0.05) { output$Pexit_warning <- renderText({exit_txt})
    } else output$Pexit_warning <- renderText({""})
    output$diffus <-   renderText({paste("Diffusivity:  ",format(Pexit*input$celllength^2*input$steps,digits=1),"km^2/day")})
    output$annmix <-   renderText({paste("Mix rate for spatial models:  ",format(annmix,digits=3),"per year")})
    #output$swimspeed<- renderText({paste("Movement speed:  ",format(swimspeed*365,digits=1),"km/yr")})
    #output$mixecospace=renderText({paste("Ecospace mixing rate from annual speed:  ",format(mixecospace,digits=2))})
    #output$meanD=      rendeerText({paste("Mean movement dist:  ",format(meanD,digits=2),"km/day")})
    
    output$movedist <- renderText({paste("Ecospace movement dispersal:",format(annmix*input$celllength*pi,digits=1),"km/yr")})
    #output$meanDyr <-  renderText({paste("Ecospace movement dist from mean daily dist:",format(meanD*365,digits=1),"km/yr")})
    
    # store the settings and results for this run
    #c("time", "body length", "swim speed", "hours active", "turns","turn sd","home base","cell length","# fish","Ecospace dispersal")
    end <- Sys.time()
    output$time <-  renderText({paste("Run time:",format(end-beginning,digits=2))})
    # if updating this list, remember that it has to be done with names/dimensioning and in two places down at bottom
    dat = c(input$name,input$bl,input$bl.cv,input$bls,input$bls.cv,input$hours,input$hours.cv,input$turns,input$wthome,
            input$celllength,input$steps,input$nfish,Pexit*input$celllength^2*input$steps,annmix,annmix*input$celllength*pi,end-beginning)
    
    res <<- rbind(res,dat)
    print(end-beginning)    
  })  
  
  randomVals <- observeEvent(input$batch, {
    
    for(j in 1:20) {
      beginning <- Sys.time()
      data.out = do.calc()
      distturn = data.out$a
      Pexit = data.out$b
      annmix = data.out$c
      
      output$distturn <- renderText({paste("Distance between turns:  ",format(distturn,digits=1),"m")})
      output$Pexit <-    renderText({paste("Exit prob. per cell face:  ",format(Pexit,digits=4),"per day")})
      output$diffus <-   renderText({paste("Diffusivity (D):  ",format(Pexit*input$celllength^2*input$steps,digits=1),"km^2/day")})
      output$annmix <-   renderText({paste("Mix rate for spatial models:  ",format(annmix,digits=3),"per year")})
      #output$swimspeed<- renderText({paste("Movement speed:  ",format(swimspeed*365,digits=1),"km/yr")})
      #output$mixecospace=renderText({paste("Ecospace mixing rate from annual speed:  ",format(mixecospace,digits=2))})
      #output$meanD=      rendeerText({paste("Mean movement dist:  ",format(meanD,digits=2),"km/day")})
      
      output$movedist <- renderText({paste("Best Ecospace movement dist:",format(annmix*input$celllength*pi,digits=1),"km/yr")})
      #output$meanDyr <-  renderText({paste("Ecospace movement dist from mean daily dist:",format(meanD*365,digits=1),"km/yr")})
      
      end <- Sys.time()
      output$time <-  renderText({paste("Run time:",format(end-beginning,digits=2))})
      # store the settings and results for this run
      # if updating this list, remember that it has to be done with names/dimensioning and in two places down at bottom
      dat = c(input$name,input$bl,input$bl.cv,input$bls,input$bls.cv,input$hours,input$hours.cv,input$turns,input$wthome,
              input$celllength,input$steps,input$nfish,Pexit*input$celllength^2*input$steps,annmix,annmix*input$celllength*pi,end-beginning)
      res <<- rbind(res,dat)
      print(end-beginning)
    }
  })
  
  
  output$distPlot <- renderPlot({
    prepare_plot()
    plot(Xend,Yend,xlab="Distance, km",ylab="Distance, km",cex=0.5,col="darkblue")
    title("Distances moved over one day")
  })
  
  output$trajPlot <- renderPlot({
    prepare_plot()    # is reactive, so won't do anything
    plot(input$celllength*xone,input$celllength*yone,type="l",xlab="Distance, km",ylab="Distance, km")
    title("Sample fish trajectory")
  })
  
  output$histPlot <- renderPlot({
    prepare_plot()   # is reactive, so won't do anything
    points(0,0,pch=16,col="black")
    hist(Dmove*input$celllength,main="Distances moved per day",xlab="Distance, km",ylab="Proportion of fish",freq=FALSE)
  })
  
  
  output$downFile <- downloadHandler(
    filename = function() {
      paste0("Ecospace dispersal estimates", " ", Sys.time(), ".csv") 
    },
    content = function(file) {
      write.csv(res, file, row.names = FALSE)
    }
  )
  
  # output$downloadPlot <- downloadHandler(
  #   filename = function() { "output.pdf" },
  #   content = function(file) {
  #     pdf(file, paper = "default")
  #     output$histPlot()
  #     dev.off()
  #   }
  # )
}


# Run the application 
shinyApp(ui = ui, server = server)

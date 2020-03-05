library(shiny)
library(shinythemes)
library(DT)
library(survival)
library(SPREDA)
library(MASS)
library(actuar)

#Function for predicting quantiles
predictData <- function (p, location=0, scale=1,
                         shape=NULL, rate=NULL, dist) {
    if (dist=="weibull") {exp(qsev(p)*scale + location)}
    else if (dist=="sev") {qsev(p)*scale + location}
    else if (dist=="frechet") {exp(qlev(p)*scale + location)}
    else if (dist=="lev") {qlev(p)*scale + location}
    else if (dist=="lognormal") {exp(qnorm(p)*scale + location)}
    else if (dist=="gaussian") {qnorm(p)*scale + location}
    else if (dist=="loglogistic") {exp(qlogis(p)*scale + location)}
    else if (dist=="logistic") {qlogis(p)*scale + location}
    else if (dist=="exponential") {qexp(p)*scale}
    else if (dist=="gamma") {qgamma(p, shape=shape)*scale}
    else if (dist=="beta") {qbeta(p, shape1=shape[1],
                                  shape2=shape[2])*scale + location}
}



#Function for constructing probability plots for
# Datos completos (sin censura) y con censura a derecha
ProbPlot <- function (data, event=rep(1, length(data)),
                      weight=rep(1, length(data)), dist,
                      shape=NULL, simplify=TRUE){
    #-data un vector con los puntos observados.
    #-event es un indicador de estado (0 y 1) censura 0, sin censura 1
    # se omite para datos completos
    #-weight vector opcional contiene el peso de las observaciones
    #-dist distribucion asumida. Distribuciones "lognormal",
    # "gaussian", "exponential", "weibull", "sev", "frechet",
    # "lev", "logistic" ,"loglogistic", "gamma" y "beta"
    #-shape vector especificando los parametros de la distribucion
    #-simplify: if TRUE then do not return a matrix containing the data points
    # (=data), standardized quantiles (=standardQuantile), and corresponding
    # probabilities.
    
    #plotting positions
    dataSorted <- data[order(data)]
    eventSorted <- event[order(data)]
    datai <- unique(dataSorted[which(eventSorted==1)])
    
    cumSurvRaw <- survfit(Surv(data, event)~ 1, weights=weight)
    cumSurv <- unique(rep(cumSurvRaw$surv, cumSurvRaw$n.event))
    cumFail <- 1-cumSurv
    lagFail <- c(0, cumFail[-length(cumFail)])
    Prob <- .5*(cumFail+lagFail)
    
    #labels for ticks along y-axis
    tick.probs <- c(.001,.005,.01,.02,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.99,.999)
    
    #implement distributions
    if (dist=="lognormal" | dist=="gaussian") {
        yi <- qnorm(Prob) #plotting positions
        tick.pos <- qnorm(tick.probs) #positions for ticks along y-axis
        ylimr <- qnorm(range(tick.probs)) #range along y-axis
    }
    else if (dist=="exponential") {
        yi <- qexp(Prob) #plotting positions
        tick.pos <- qexp(tick.probs) #positions for ticks along y-axis
        ylimr <- qexp(range(tick.probs)) #range along y-axis
    }
    else if (dist=="weibull" | dist=="sev") {
        yi <- qsev(Prob) #plotting positions
        tick.pos <- qsev(tick.probs) #positions for ticks along y-axis
        ylimr <- qsev(range(tick.probs)) #range along y-axis
    }
    else if (dist=="frechet" | dist=="lev") {
        yi <- qlev(Prob) #plotting positions
        tick.pos <- qlev(tick.probs) #positions for ticks along y-axis
        ylimr <- qlev(range(tick.probs)) #range along y-axis
    }
    else if (dist=="loglogistic" | dist=="logistic") {
        yi <- qlogis(Prob) #plotting positions
        tick.pos <- qlogis(tick.probs) #positions for ticks along y-axis
        ylimr <- qlogis(range(tick.probs)) #range along y-axis
    }
    else if (dist=="gamma") {
        yi <- qgamma(Prob, shape=shape) #plotting positions
        tick.pos <- qgamma(tick.probs, shape=shape) #positions for ticks along y-axis
        ylimr <- qgamma(range(tick.probs), shape=shape) #range along y-axis
    }
    else if (dist=="beta") {
        yi <- qbeta(Prob, shape1=shape[1], shape2=shape[2]) #plotting positions
        tick.pos <- qbeta(tick.probs, shape1=shape[1], shape2=shape[2]) #positions for ticks along y-axis
        ylimr <- qbeta(range(tick.probs), shape1=shape[1], shape2=shape[2]) #range along y-axis
    }
    
    #determine range along x-axis
    rangeData <- range(data)
    
    #construct plot
    #distributions that require a log transform of the data scale
    if (dist=="weibull" | dist=="lognormal" | dist=="frechet" |
        dist=="loglogistic") {
        plot(0, type="n", log="x",
             xlim=c(rangeData[1], rangeData[2]),
             ylim=c(ylimr[1], ylimr[2]),
             xlab="Data", ylab="Probability", main=paste(dist, "distribution"),
             axes=FALSE, frame.plot=TRUE)}
    else {
        plot(0, type="n",
             xlim=c(rangeData[1], rangeData[2]),
             ylim=c(ylimr[1], ylimr[2]),
             xlab="Data", ylab="Probability", main=paste(dist, "distribution"),
             axes=FALSE, frame.plot=TRUE)}
    
    axis(2, at=tick.pos, labels=tick.probs)
    axis(1)
    points(datai, yi, col="red")
    #draw raster
    abline(h = tick.pos, lty=3, col="gray")
    abline(v = axTicks(1), lty=3, col="gray")
    
    #return plotting positions
    if (!simplify) cbind(data=datai, standardQuantile=yi, probability=Prob)
}
lineas <- function(x, cens=rep(1, length(x)),
                   weight=rep(1, length(x)), dist,
                   shape=NULL, simplify=TRUE){
    ps <- seq(.001, .999, length.out=200)
    #implement distributions
    if (dist=="lognormal") {
        mleLognormal <- survreg(Surv(x, cens) ~ 1,
                                weights=weight, dist="lognormal")
        predLognormal <- sapply(ps, predictData, location=coef(mleLognormal),
                                scale=mleLognormal$scale, dist="lognormal")
        lines(predLognormal, qnorm(ps), col="blue")
        
    }
    else if (dist=="gaussian"){
        mleNormal <- survreg(Surv(x, cens) ~ 1,
                             weights=weight, dist="gaussian")
        predNormal <- sapply(ps, predictData, location=coef(mleNormal),
                             scale=mleNormal$scale, dist="gaussian")
        lines(predNormal, qnorm(ps), col="blue")
    }
    else if (dist=="exponential") {
        mleExponential <- fitdistr(x, densfun="exponential")
        predExponential <- sapply(ps, predictData,
                                  scale=1/mleExponential$estimate,
                                  dist="exponential")
        lines(predExponential, qexp(ps), col="blue")
        
    }
    else if (dist=="weibull") {
        mleWeibull <- survreg(Surv(x, cens) ~ 1,
                              weights=weight, dist="weibull")
        predWeibull <- sapply(ps, predictData, location=coef(mleWeibull),
                              scale=mleWeibull$scale, dist="weibull")
        lines(predWeibull, qsev(ps), col="blue")
        
    }
    else if (dist=="sev") {
        mleGumbel <- survreg(Surv(x, cens) ~ 1,
                             weights=weight, dist="extreme")
        predGumbel <- sapply(ps, predictData, location=coef(mleGumbel),
                             scale=mleGumbel$scale, dist="sev")
        lines(predGumbel, qsev(ps), col="blue")
        
    }
    else if (dist=="frechet") {
        mleFrechet <- Lifedata.MLE(Surv(x, cens) ~ 1,
                                   weights=weight, dist="frechet")
        predFrechet <- sapply(ps, predictData, location=coef(mleFrechet)[1],
                              scale=coef(mleFrechet)[2], dist="frechet")
        lines(predFrechet, qlev(ps), col="blue")
        
    }
    else if (dist=="lev") {
        mleFrechet <- Lifedata.MLE(Surv(x, cens) ~ 1,
                                   weights=weight, dist="frechet")
        predFrechet <- sapply(ps, predictData, location=coef(mleFrechet)[1],
                              scale=coef(mleFrechet)[2], dist="lev")
        lines(predFrechet, qlev(ps), col="blue")
    }
    else if (dist=="loglogistic") {
        mleLoglogis <- survreg(Surv(x, cens) ~ 1,
                               weights=weight, dist="loglogistic")
        predLoglogis <- sapply(ps, predictData, location=coef(mleLoglogis),
                               scale=mleLoglogis$scale, dist="loglogistic")
        lines(predLoglogis, qlogis(ps), col="blue")
    }
    else if (dist=="logistic") {
        mleLogis <- survreg(Surv(x, cens) ~ 1,
                            weights=weight, dist="logistic")
        predLogis <- sapply(ps, predictData, location=coef(mleLogis),
                            scale=mleLogis$scale, dist="logistic")
        lines(predLogis, qlogis(ps), col="blue") 
    }
    
} 

# UI para aplicacion que realiza graficos de probabilidad
ui <- fluidPage(theme = shinytheme("united"),

    # Titulo de la aplicación

    titlePanel("Aplicación para realizar gráficos de probabilidad"),

    # Sidebar 
    sidebarLayout(
        sidebarPanel(
            radioButtons("formato", "Formato de los datos",
                         c("Completos"="comp", "Con cesura a derecha"="censura")),
            
            # Seleccionar archivo con los datos
            fileInput("file1", "Seleccionar datos (archivo CSV)",
                      accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv")),
            h6("Asegúrese de que su archivo CSV contenga una columna con los datos observados
            si los datos son completos, y dos columnas si tienen censura a derecha. La 
               primera columna contendrá los datos observados y la segunda columna
               ceros y unos, donde 0 representa censura y 1 representa falla", align = "left"),
            tags$hr(),
            h6("Si usted no ha seleccionado ningún archivo podrá utilizar la
               aplicación con datos simulados de una exponencial con rate=1 . Las fallas
               se simularan de una binomial con n=1 y p=0.8", align = "left"),
            tags$hr(),

            #Seleccionar la distribución
            selectInput("distribucion", "Distribución",
                        c(Normal = "gaussian",
                          Lognormal = "lognormal",
                          Exponencial = "exponential",
                          Weibull = "weibull",
                          SEV = "sev",
                          LEV = "lev",
                          Frechet = "frechet",
                          Logistica = "logistic",
                          Loglogistica = "loglogistic"
                          # Gamma = "gamma",
                          # Beta = "beta"
                          ),selected = "gaussian")

            ),
        # Datos complentos o con censura a derecha

        
        mainPanel(
            tabsetPanel(
             # Muestra el gráfico de probabilidad
             tabPanel("Grafico Probabilidad", plotOutput("grafico_probabilidad")),
             # Muestra los datos, cuantiles estandarizados con sus probabilidades
             tabPanel("Cuantiles y probabilidades",DTOutput("posiciones"))
             
            ),
            h4("Datos", align = "left"),
            DTOutput("contents")
        )
        
        
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    # Para cargar archivo
    filedata <- reactive({
        infile <- input$file1
        if (is.null(infile)){
            set.seed(1234)
            x <- rexp(30)
            set.seed(1234)
            cens <- rbinom(30,size=1,0.8)
            return(data.frame(x=sort(x),
                             cens=sort(cens,decreasing = T)))      
        }
        read.csv(infile$datapath)
    })
    # Para mostrar los datos cargados
    output$contents = renderDT({
        df <- filedata()
        return(df)
    })

    # Para obtener las posiciones
    output$posiciones = renderDT({
        df <- filedata()
        # Verifica el formato de los datos para utilizar ProbPlot
        if(input$formato == "comp"){ProbPlot(df[,1],dist = input$distribucion,simplify = FALSE)}
        else {ProbPlot(df[,1],event = df[,2],
             dist = input$distribucion,simplify = FALSE)}
    })
    
    # Realiza el grafico de probabilidad
    output$grafico_probabilidad <- renderPlot({
        df <- filedata()

        # Verifica el formato de los datos para utilizar ProbPlot
        if(input$formato == "comp"){ProbPlot(df[,1],dist = input$distribucion,simplify = FALSE)
                                    lineas(df[,1],dist = input$distribucion)}
        else{ProbPlot(df[,1],event = df[,2],
                      dist = input$distribucion,simplify = FALSE)
                      lineas(df[,1],cens = df[,2],dist = input$distribucion)}

    })
    

}

# Correr la aplicacion
shinyApp(ui = ui, server = server)

shinyUI(pageWithSidebar(
    
  # Application title
  headerPanel("SAM - Significance Analysis of Microarrays"),
  
  # Control panel;
  sidebarPanel(
    fileInput(inputId = "iFile", label = "", accept="application/vnd.ms-excel"),
    uiOutput(outputId = "ui"),
    actionButton("goButton", "Run"),

    tags$hr(),
    
    
    conditionalPanel(condition = "input.goButton != 0",
    
    radioButtons("deltaChoice", "Delta", c("Delta Slider" = "Delta Slider", "Manually Enter Data" = "Manually Enter Data")),
    conditionalPanel(condition = "input.deltaChoice == 'Delta Slider'",
      sliderInput("deltaSlider", label = "Delta value", min = 0, max = 10, value = 0, step = 0.01)
    ),
    
    conditionalPanel(condition = "input.deltaChoice == 'Manually Enter Data'",
      numericInput("deltaInput", label = "Delta value", min = 0, max = 10, value = 0, step = 0.01)
    ),
    
    conditionalPanel(condition = "(input.assayType == 'array' && (input.responseType_array == 'Two class unpaired' || input.responseType_array == 'Two class paired' || input.responseType_array == 'Two class unpaired timecourse' || input.responseType_array == 'Two class paired timecourse')) || (input.assayType == 'seq' && (input.responseType_seq == 'Two class unpaired' || input.responseType_seq == 'Two class paired'))",
      numericInput("min.foldchange", label = "Minimum fold change", min = 1, value = 0, step = 0.01) 
    ),
    
    conditionalPanel(condition = "(input.assayType == 'array' && (input.responseType_array == 'Two class unpaired' || input.responseType_array == 'Two class paired' || input.responseType_array == 'One class' || input.responseType_array == 'Survival')) || (input.assayType == 'seq' && (input.responseType_seq == 'Two class unpaired' || input.responseType_seq == 'Two class paired' || input.responseType_seq == 'Survival'))", 
      numericInput("dif", label = "Hypothesized mean difference in expression", value = 0),
      textInput("sampleSizeFactors", label = "Same size factors - four comma separated values", value = "1,2,3,5")
    ),
    
    radioButtons("localFDR", "Output local FDRs", c("No" = FALSE, "Yes" = TRUE)),
    
    downloadButton("downloadData","Save Result as XLSX File"),
    tags$hr()
    
    ),
  
    radioButtons("assayType", "Data type", c("Array" = "array", "Sequencing" = "seq")),
    
    conditionalPanel(condition = "input.assayType == 'seq'",
      selectInput("responseType_seq","Response Type", c('Quantitative','Two class unpaired', 'Survival', 'Multiclass', 'Two class paired'))
    ),
      
    conditionalPanel(condition = "input.assayType == 'array'",
      selectInput("responseType_array","Response Type", c('Quantitative','Two class unpaired', 'Survival', 'Multiclass', 'One class', 'Two class paired', 'Two class unpaired timecourse', 'One class timecourse', 'Two class paired timecourse','Pattern discovery')),       
      
      conditionalPanel(condition = paste("(input.responseType_array == 'Two class unpaired' ||", "input.responseType_array == 'Two class unpaired timecourse')"),
        radioButtons("testStatistic", "Test Statistic", c("T-statistic"="standard", "Wilcoxon"="wilcoxon"))                 
      ),
      
      conditionalPanel(condition = paste("(input.responseType_array == 'Two class unpaired timecourse' ||", "input.responseType_array == 'One class timecourse' ||", "input.responseType_array == 'Two class paired timecourse')"),
        radioButtons("timeSummaryType", "Time Summary Type", c("Slope"='slope', "Signed Area"='signed.area'))               
      ),
      
      radioButtons("centerArrays", "Median center the arrays", c("No" = FALSE, "Yes" = TRUE)),
      conditionalPanel(condition = paste("(input.responseType_array == 'Two class unpaired' ||", "input.responseType_array == 'Two class paired')"),
        radioButtons("dataLogged", "Are data in log scale (base 2)", c("No" = FALSE, "Yes" = TRUE))        
      ),
        
      conditionalPanel(condition = "input.responseType_array == 'Quantitative'",
        radioButtons("regressionMethod", "Regression method", c("Standard"='standard', "Ranks"='ranks'))
      ),
      
      radioButtons("s0", "Estimate of s0 factor for denominator", c("Automatic"="Automatic", "Use fixed percentile (eg 50)"="useFixed")),
      conditionalPanel(condition = "input.s0 == 'useFixed'",
        numericInput("s0.perc", "(-1 implies s0 is zero)", step=1, value = 50)  
      ),
      numericInput("numberOfNeighbors", "K-Nearest Neighbors Imputer: Number of Neighbors", value= 10, step=1)
    ),
 
    numericInput("nperms", "Number of Permutations", value=100, min=25, max=5000, step=1),
    uiOutput("rseed"),
    actionButton("randomButton", "Generate Random Seed")
       
  ),

  mainPanel(
#    tags$style(type="text/css",
#               ".shiny-output-error { visibility: hidden; }",
#               ".shiny-output-error:before { visibility: hidden; }"
#    ),
    tabsetPanel(id='SAM',
      tabPanel("SAM Plot", h3(textOutput("samPlotText")), plotOutput("samrPlot")), 
      tabPanel("Delta Table", h3(textOutput("deltaTableText")), tableOutput("deltaTable")),
      tabPanel("Significant Genes", h3(textOutput("positiveGenesText")), dataTableOutput("siggenes.table.up"), h3(textOutput("negativeGenesText")), dataTableOutput("siggenes.table.lo")),
      tabPanel("Sample Size", h3(textOutput("sampleSizePlotText")), plotOutput("samplePlot"),  h3(textOutput("sampleTableText1")), tableOutput("sampleTable1"), h3(textOutput("sampleTableText2")), tableOutput("sampleTable2"), h3(textOutput("sampleTableText3")), tableOutput("sampleTable3"), h3(textOutput("sampleTableText4")), tableOutput("sampleTable4"))
      
    )
    

  )
))
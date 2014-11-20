shinyUI(fluidPage(
    
  # Application title
  titlePanel("SAM - Significance Analysis of Microarrays"),
  
  # Control panel;
  fluidRow(
    column(3, 
    wellPanel(
      
    fileInput(inputId = "iFile", label = "", accept="application/vnd.ms-excel"),
    uiOutput(outputId = "ui"),
    conditionalPanel(condition = "input.analysisType == 'Standard'",
      actionButton("goButton", "Run")
    ),
    conditionalPanel(condition = "input.analysisType == 'Gene sets'",
      actionButton("goButton2", "Run")
    ),
       
    conditionalPanel(condition = "input.goButton != 0 && input.analysisType == 'Standard'",
    
    tags$hr(),                    
    radioButtons("deltaChoice", "Delta", c("Delta Slider" = "Delta Slider", "Manually Enter Delta" = "Manually Enter Delta")),
    conditionalPanel(condition = "input.deltaChoice == 'Delta Slider'",
      sliderInput("deltaSlider", label = "Delta value", min = 0, max = 10, value = 0, step = 0.01)
    ),
    
    conditionalPanel(condition = "input.deltaChoice == 'Manually Enter Delta'",
      numericInput("deltaInput", label = "Delta value", min = 0, max = 10, value = 0, step = 0.01)
    ),
    
    conditionalPanel(condition = "(input.assayType == 'array' && (input.responseType_array == 'Two class unpaired' || input.responseType_array == 'Two class paired' )) || (input.assayType == 'seq' && (input.responseType_seq == 'Two class unpaired' || input.responseType_seq == 'Two class paired'))",
      numericInput("min.foldchange", label = "Minimum fold change", min = 1, value = 0, step = 0.01) 
    ),
    
    conditionalPanel(condition = "(input.assayType == 'array' && (input.responseType_array == 'Two class unpaired' || input.responseType_array == 'Two class paired' || input.responseType_array == 'One class' || input.responseType_array == 'Survival')) || (input.assayType == 'seq' && (input.responseType_seq == 'Two class unpaired' || input.responseType_seq == 'Two class paired' || input.responseType_seq == 'Survival'))", 
      numericInput("dif", label = "Hypothesized mean difference in expression", value = 0),
      textInput("sampleSizeFactors", label = "Same size factors - four comma separated values", value = "1,2,3,5")
    ),
    
    radioButtons("localFDR", "Output local FDRs", c("No" = FALSE, "Yes" = TRUE)),
    
    textInput("dir", "Paste the filepath to save the output", value = getwd()),
    textInput("fname", "Type the file name you would like to save as", value = "result"),
    actionButton("saveButton", "Save")   
 
    ),
    
    
    conditionalPanel(condition = "input.goButton2 != 0 && input.analysisType == 'Gene sets'",
      
      tags$hr(),                    
      radioButtons("fdrChoice", "FDR", c("FDR Slider" = "FDR Slider", "Manually Enter FDR" = "Manually Enter FDR")),
      conditionalPanel(condition = "input.fdrChoice == 'FDR Slider'",
        sliderInput("fdrSlider", label = "FDR value", min = 0, max =1, value = 0.5, step = 0.01)
      ),
      
      conditionalPanel(condition = "input.fdrChoice == 'Manually Enter FDR'",
        numericInput("fdrInput", label = "FDR value", min = 0, max = 1, value = 0.5, step = 0.01)
      ),
      textInput("dir2", "Paste the filepath to save the output", value = getwd()),
      textInput("fname2", "Type the file name you would like to save as", value = "result"),
      actionButton("saveButton2", "Save")   
      
    )
    ),
    
    wellPanel(
    radioButtons("assayType", "Data type", c("Array" = "array", "Sequencing" = "seq")),
    
    conditionalPanel(condition = "input.assayType == 'seq'",
      selectInput("responseType_seq","Response Type", c('Quantitative','Two class unpaired', 'Survival', 'Multiclass', 'Two class paired'))
    ),    
    
    conditionalPanel(condition = "input.assayType == 'array'",
      selectInput("responseType_array","Response Type", c('Quantitative','Two class unpaired', 'Survival', 'Multiclass', 'One class', 'Two class paired', 'Two class unpaired timecourse', 'One class timecourse', 'Two class paired timecourse','Pattern discovery')),       
      
      conditionalPanel(condition = paste("(input.responseType_array == 'Quantitative' ||", "input.responseType_array == 'Two class unpaired' ||", "input.responseType_array == 'Survival' ||", "input.responseType_array == 'Multiclass' ||", "input.responseType_array == 'Two class paired' ||", "input.responseType_array == 'Two class unpaired')"),
        radioButtons("analysisType", "Analysis Type", c("Standard (genes)" = "Standard", "Gene sets" = "Gene sets"))
      ),
      
      conditionalPanel(condition = paste("input.analysisType == 'Gene sets' && ", "(input.responseType_array == 'Quantitative' ||", "input.responseType_array == 'Two class unpaired' ||", "input.responseType_array == 'Survival' ||", "input.responseType_array == 'Multiclass' ||", "input.responseType_array == 'Two class paired' ||", "input.responseType_array == 'Two class unpaired')"),
        fileInput(inputId = "gmtFile", label = "Browse for Gene Sets .gmt file"),
        numericInput("minGeneSet", label = "Minimum gene set size", min = 1, value = 15, step = 1), 
        numericInput("maxGeneSet", label = "Maximum gene set size", min = 1, value = 500, step = 1)  
      ),
      
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
    )
    
    ),
  
  
  column(9, 
#    tags$style(type="text/css",
#               ".shiny-output-error { visibility: hidden; }",
#               ".shiny-output-error:before { visibility: hidden; }"
#    ),

    conditionalPanel(condition = "input.analysisType == 'Gene sets'",
      tabsetPanel(id = 'gene set',
        tabPanel("GSA Plot", h3(textOutput("GSAPlotText")), plotOutput("geneSetPlot")),    
        tabPanel("Significant Gene Set", h3(textOutput("geneSetPositive")), tableOutput("positiveGeneSet"), h3(textOutput("geneSetNegative")), tableOutput("negativeGeneSet"), h3(textOutput("GeneScoreText")), numericInput("geneset.number", "geneset number", value= 1, step=1), tableOutput("testGeneSet")),
        tabPanel("All Gene Set",  h3(textOutput("geneSetFullPositive")), tableOutput("positiveFullGeneSet"), h3(textOutput("geneSetFullNegative")), tableOutput("negativeFullGeneSet")),
        tabPanel("Gene Set Collection Info", h3(textOutput("geneSetInfoText")), tableOutput("geneSetTable"), h3(textOutput("geneSetInfoText2")), tableOutput("geneSetQuantile"), h3(textOutput("geneSetInfoText3")), tableOutput("geneSetTableGenes")),
        tabPanel("Current Settings", h3(textOutput("inputParametersText2")), tableOutput("inputParameters2"), h3(textOutput("computedValuesText2")), tableOutput("computedValues2"))
        )      
    ),
    
    conditionalPanel(condition = "input.analysisType == 'Standard'", 
      tabsetPanel(id='SAM',
        tabPanel("SAM Plot", h3(textOutput("samPlotText")), plotOutput("samrPlot")), 
        tabPanel("Delta Table", h3(textOutput("deltaTableText")), tableOutput("deltaTable"), h3(textOutput("missRateText")), tableOutput("missRate")),
        tabPanel("Significant Genes", h3(textOutput("positiveGenesText")), dataTableOutput("siggenes.table.up"), h3(textOutput("negativeGenesText")), dataTableOutput("siggenes.table.lo")),
        tabPanel("All Genes", h3(textOutput("allPositiveGenesText")), dataTableOutput("Allgenes.table.up"), h3(textOutput("allNegativeGenesText")), dataTableOutput("Allgenes.table.lo")),
        tabPanel("Sample Size", h3(textOutput("sampleSizePlotText")), plotOutput("samplePlot"),  h3(textOutput("sampleTableText1")), tableOutput("sampleTable1"), h3(textOutput("sampleTableText2")), tableOutput("sampleTable2"), h3(textOutput("sampleTableText3")), tableOutput("sampleTable3"), h3(textOutput("sampleTableText4")), tableOutput("sampleTable4")),
        tabPanel("Current Settings", h3(textOutput("inputParametersText")), tableOutput("inputParameters"),  h3(textOutput("computedValuesText")), tableOutput("computedValues"), h3(textOutput("eigengeneText")), tableOutput("eigengene") )
      )
    )
    
  )

)
))
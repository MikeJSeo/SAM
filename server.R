#file:///C:/Users/mike/Desktop/MAIN/Examples/

library(XLConnect)
library(samr)

shinyServer(function(input, output) {  
  
  
  findDelta = reactive({
    
    if(input$deltaChoice == "Delta Slider")
      input$deltaSlider
    else if (input$deltaChoice == "Manually Enter Data")
      input$deltaInput
    
  })
  
  survival = function(firstrow){
    
    y = vector(length = length(firstrow), mode ="integer")
    censoring.status = vector(length = length(firstrow), mode ="integer")        
    for(i in seq.int(length(y))) {
      split = strsplit(as.character(firstrow[,i]), ",")  
      y[i] = as.character(data.frame(split)[1,])
      y[i] = gsub(" ","", y[i])
      y[i] = substring(y[i], 2)
      
      censoring.status[i] = as.character(data.frame(split)[2,])
      censoring.status[i] = gsub(" ","", censoring.status[i])
      censoring.status[i] = substring(censoring.status[i], 1, 1)
    }
    
    y = as.numeric(y)
    censoring.status = as.numeric(censoring.status)
    
    list(y= y, censoring.status = censoring.status)
  }
  
  output$negativeGenesText <- renderText({
    
    if(!is.null(getSiggenesTable())){
      siggenes.table = getSiggenesTable()
      paste("Negative Genes (", siggenes.table$ngenes.lo, ")", sep ="")
    }
  })
  
  output$positiveGenesText <- renderText({
    
    if(!is.null(getSiggenesTable())){
      siggenes.table = getSiggenesTable()
      paste("Positive Genes (", siggenes.table$ngenes.up, ")", sep = "")
    }
  })
  
  output$deltaTableText <- renderText({
    if(!is.null(getResult())){
      "Delta Table"
    }
  })
  
  output$samPlotText <- renderText({
    if(!is.null(getResult())){
      "SAM Plot"
    }
  })
  
  output$rseed <- renderUI({
  
    if(input$randomButton == 0){
      value = 1234567
    }else{
      value = round(runif(1, 0, 10^6))
    }
    
    numericInput("random.seed", "Random Seed", value= value)

  })
  
  
  getSiggenesTable = reactive({
    
    if(!is.null(getResult())){
      result = getResult()
      delta = findDelta()
      samr.obj = result$samr.obj
      data = result$data
      delta.table = result$delta.table
      min.foldchange = input$min.foldchange
    
      siggenes.table = samr.compute.siggenes.table(samr.obj,delta, data, delta.table, min.foldchange = min.foldchange)
      siggenes.table
    }
  })
  
  
  getResult <- reactive({
    
    input$goButton
    
    objFile = chooseFile()
    if(!is.null(objFile)){
      
      Sheet = input$sheet
      if (!is.null(Sheet)){
        wb = loadWorkbook(objFile$path)
        dat = readWorksheet(wb, Sheet, header = FALSE)
        x = dat[-1, c(-1,-2)]
        
        x = as.matrix(x)
        class(x) = "numeric"
        
#        if (sum(is.na(x)) > 0) {
#          require(impute)
#          x = impute.knn(x, k = input$numberOfNeighbors)
#          if (!is.matrix(x)) {
#            x = x$data
#          }
#        }

        geneid = dat[-1,2]
        genenames = dat[-1,1]
           
        isolate({
          
        firstrow = as.vector(dat[1,c(-1,-2)])
        censoring.status = NULL
        eigengene.number = NULL
        
        if( (input$responseType_seq == "Survival" && input$assayType == "seq")
           || (input$responseType_array == "Survival" && input$assayType == "array")){
          survival = survival(firstrow)
          y = survival$y
          censoring.status = survival$censoring.status
        }
        else if( (input$responseType_array == "Pattern discovery" && input$assayType == "array") ){
          y = NULL
          eigengene.number = firstrow[1]
          last = nchar(eigengene.number)
          eigengene.number = as.numeric(substr(eigengene.number, last, last))
        }
        
        else{
          y = firstrow
        } 
        
        
        data =list(x=x,y=y, genenames=genenames, geneid=geneid, logged2= as.logical(input$dataLogged), censoring.status = censoring.status, eigengene.number = eigengene.number)        
        resp.type = if(input$assayType == "seq"){input$responseType_seq}else{input$responseType_array} 
        s0.perc = if(is.na(input$s0.perc) || input$s0 == "Automatic"){NULL}else{input$s0.perc}
        center.arrays = as.logical(input$centerArrays)
        
        samr.obj = samr(data, resp.type = resp.type, assay.type = input$assayType, s0.perc = NULL, nperms = input$nperms, center.arrays = center.arrays, testStatistic = input$testStatistic, time.summary.type = input$timeSummaryType,  regression.method = input$regressionMethod, random.seed = input$random.seed)           
      
        delta.table = samr.compute.delta.table(samr.obj, min.foldchange = input$min.foldchange)
        
        list (data = data, samr.obj = samr.obj, delta.table = delta.table)  
        })
      }
    }  
        
   
  })
  
  

  chooseFile = reactive({
    
    input$goButton
    
    inFile = isolate(input$iFile)
    if (!is.null(inFile)) {
      # Determine document format;
      ptn = "\\.[[:alnum:]]{1,5}$"
      suf = tolower(regmatches(inFile$name, regexpr(ptn, inFile$name)))
 
      # Options for Excel documents;
      if (suf %in% c('.xls', '.xlsx')) {
        wb = loadWorkbook(inFile$datapath)
        sheets = getSheets(wb)
        output$ui <- renderUI({
          list(
            selectInput(inputId = "sheet", label = "Select a sheet:", choices = sheets),
            textInput(inputId = 'arg', label = 'Additional Arguments:', value = ' ')
          )
        })
        return(list(path = inFile$datapath, suf = suf))
      } 
      
      # Options for txt documents;
      if (suf %in% c('.txt', '.csv')) {
        output$ui <- renderUI({
          list(
            checkboxInput(inputId = 'header', label = 'First line as header', value = TRUE),
            textInput(inputId = 'sep', label = 'Separator', value = " "),
            textInput(inputId = 'quote', label = 'Quote', value = '\"'),
            textInput(inputId = 'arg', label = 'Additional Arguments:', value = ' '),
            tags$hr()
          )
        })
        return(list(path = inFile$datapath, suf = suf))
      }
    } else {return(NULL)}
  })
  

  
  output$samrPlot <- renderPlot({

      result = getResult()
      samr.obj = result$samr.obj
      delta = findDelta()
      min.foldchange = input$min.foldchange
      if(!is.null(result))
        samr.plot(samr.obj, delta, min.foldchange = min.foldchange)
  })
  
  output$deltaTable <- renderDataTable({
    
    
    result = getResult()
    delta.table = result$delta.table
    delta.table
  
  })
  
  output$siggenes.table.up <- renderDataTable({
    
    siggenes.table = getSiggenesTable()
    #write.csv(siggenes.table$genes.up, file = "siggenestable.csv")
    result = getResult()
    delta.table = result$delta.table
#  writeWorksheetToFile("result.xlsx", data = list(i1 = as.data.frame(delta.table), i2 = as.data.frame(siggenes.table$genes.up), i3 = as.data.frame(siggenes.table$genes.lo)),
 #                        sheet = c("FirstSheet", "SecondSheet", "FirstSheet"))
    if(!is.null(siggenes.table$genes.up)){
      siggenes.table$genes.up
    }
  #  else if(!is.null(siggenes.table$genes.up[,-3])){
    #  siggenes.table$genes.up[,-3]
  #  }
  
  })
  
  output$siggenes.table.lo <- renderDataTable({
    
    siggenes.table = getSiggenesTable()
    
    if(!is.null(siggenes.table$genes.lo)){
      
      siggenes.table$genes.lo
    #  siggenes.table$genes.lo[,-3]
    }
  })
  
  
  output$contents <- renderTable({
    objFile <- chooseFile()
    if (!is.null(objFile)) {
      suf <- objFile$suf
      # For Excel documents;
      if (suf %in% c('.xls', '.xlsx')) {
        Sheet <- input$sheet
        if (!is.null(Sheet)){
          
          if (input$arg %in% c(' ', '')) {
            wb <- loadWorkbook(objFile$path)
            dat <- readWorksheet(wb, Sheet, header = FALSE)
            return(dat)
          } else {
            wb <- loadWorkbook(objFile$path)
            expr <- paste('readWorksheet(wb, Sheet,', input$arg, ')', sep = '')
            print(expr)
            dat <- eval(parse(text = expr))
            return(dat)
          }
          
        } else {return(NULL)}
      }
      # For .txt and .csv documents;
      if (suf %in% c('.txt', '.csv')) {
        if (is.null(input$header)) {
          dat <- read.table(objFile$path)
          return(dat)
        } else {
          if (input$arg %in% c(' ', '')) {
            dat <- read.table(objFile$path, header=input$header, sep=input$sep, quote=input$quote)
            return(dat)
          } else {
            expr.1 <- paste('"', gsub('\\', '/', objFile$path, fixed = TRUE), '"', sep = '')
            expr.2 <- paste(expr.1, 
                            paste('header =', input$header), 
                            paste('sep =', paste("'", input$sep, "'", sep = '')), 
                            paste('quote =', paste("'", input$quote, "'", sep = '')), input$arg,  sep = ', ')
            print(expr.2)
            expr <- paste('read.table(', expr.2, ')', sep = '')
            print(expr)
            dat <- eval(parse(text = expr))
            return(dat)
          }
        }
      }
      
    } else {return(NULL)}
    
  })



  output$downloadData <- downloadHandler(
    filename = function() { "data.xlsx" },
    content = function(file) {
      siggenes.table = getSiggenesTable()
      result = getResult()
      delta.table = result$delta.table
      
      fname = paste(file, "xlsx", sep = ".")
#      wb = loadWorkbook(fname, create = TRUE)
#      createSheet(wb, name = c("Sheet1"))
#      writeWorksheet(wb, c(1:3), sheet = "Sheet1")
#      saveWorkbook(wb)
#      file.rename(fname, file)
      writeWorksheetToFile(fname, data = list(i1 = delta.table, i2 = siggenes.table$genes.up, i3 = siggenes.table$genes.lo), sheet = c("Delta Table", "Positive Genes", "Negative Genes"))
  #    saveWorkbook(wb)
      file.rename(fname, file)
  

      }
    )

  
})

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
  
  
  output$sampleSizePlotText <- renderText({
    if(!is.null(getSampleSize())){
      "Sample Size Plot"
    }
  })
  
  
  output$sampleTableText1 = renderText({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(getSampleSize())){
      paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[1] * samr.assess.samplesize.obj$n, sep = "") 
    }
    
  }) 
  
  output$sampleTableText2 = renderText({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(getSampleSize())){
      paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[2] * samr.assess.samplesize.obj$n, sep = "") 
    }
    
  }) 
  
  output$sampleTableText3 = renderText({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(getSampleSize())){
      paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[3] * samr.assess.samplesize.obj$n, sep = "") 
    }
    
  }) 
  
  output$sampleTableText4 = renderText({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(getSampleSize())){
      paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[4] * samr.assess.samplesize.obj$n, sep = "") 
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
  
  getSampleSize = reactive({
    
    
    ar = isolate(input$responseType_array)
    seq = isolate(input$responseType_seq)
    
    if(!is.null(getResult()) && ((input$assayType == "array" && (ar == "Two class unpaired" || ar == "Two class paired" || ar == "One class" || ar == "Survival")) || (input$assayType == "seq" && (seq == "Two class unpaired" || seq == "Two class paired" || seq == "Survival")))){
      
      result = getResult()
      samr.obj = result$samr.obj
      data = result$data
      dif = input$dif
      samplesize.factors = input$sampleSizeFactors
      samplesize.factors = as.numeric(unlist(strsplit(samplesize.factors, ",")))
      samr.assess.samplesize.obj =  samr.assess.samplesize(samr.obj, data, dif, samplesize.factors)
      samr.assess.samplesize.obj
    }
    
    
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
    
    objFile = isolate(input$iFile)
    
    if(!is.null(objFile)){
      wb = loadWorkbook(objFile$datapath)
      sheets = getSheets(wb)
      dat = readWorksheet(wb, sheets, header = FALSE)
      x = dat[-1, c(-1,-2)]
        
      x = as.matrix(x)
      class(x) = "numeric"

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
  })
  
  output$sampleTable1 = renderTable({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(getSampleSize())){
      samr.assess.samplesize.obj$results[,,1] 
    }

  })
  
  output$sampleTable2 = renderTable({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(getSampleSize())){
      samr.assess.samplesize.obj$results[,,2] 
    }    
    
  })
  
  output$sampleTable3 = renderTable({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(getSampleSize())){
      samr.assess.samplesize.obj$results[,,3]  
    }  
  })
  
  output$sampleTable4 = renderTable({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(getSampleSize())){
      samr.assess.samplesize.obj$results[,,4]  
    }
  })
  
  output$samplePlot = renderPlot({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(getSampleSize())){
      samr.assess.samplesize.plot(samr.assess.samplesize.obj)  
    }
    
  }) 
  
  output$samrPlot = renderPlot({

      result = getResult()
      samr.obj = result$samr.obj
      delta = findDelta()
      min.foldchange = input$min.foldchange
      if(!is.null(result))
        samr.plot(samr.obj, delta, min.foldchange = min.foldchange)
  })
  
  output$deltaTable <- renderTable({
    
    result = getResult()
    delta.table = result$delta.table
    delta.table
  
  })
  
  output$siggenes.table.up <- renderDataTable({
    
    siggenes.table = getSiggenesTable()
    result = getResult()
    delta.table = result$delta.table

    if(!is.null(siggenes.table$genes.up)){
      siggenes.table$genes.up
    }
  
  })
  
  output$siggenes.table.lo <- renderDataTable({
    
    siggenes.table = getSiggenesTable()
    
    if(!is.null(siggenes.table$genes.lo)){     
      siggenes.table$genes.lo
    }
  })
  
  
  output$contents <- renderTable({
    objFile <- chooseFile()
    if (!is.null(objFile)) {
      suf <- objFile$suf
  
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
    } else {return(NULL)}
    
  })



  output$downloadData <- downloadHandler(
    filename = function() { "result.xlsx" },
    content = function(file) {
      siggenes.table = getSiggenesTable()
      result = getResult()
      delta.table = result$delta.table
      samr.assess.samplesize.obj =  getSampleSize()

    #  samr.assess.samplesize.obj$results[,,4]  
      
      fname = paste(file, "xlsx", sep = ".")
      data = list()
      dataname = c()
    
      if(!is.null(delta.table)){
        data$delta.table = delta.table
        dataname = c(dataname, "delta table")
      }
      if(!is.null(siggenes.table$genes.up)){
        data$siggenes.table.genes.up = as.data.frame(siggenes.table$genes.up)
        dataname = c(dataname, "Postivie Genes")
      }
      if(!is.null(siggenes.table$genes.lo)){
        data$siggenes.table.genes.lo = as.data.frame(siggenes.table$genes.lo)
        dataname = c(dataname, "Negative Genes")
      }
      if(!is.null(samr.assess.samplesize.obj)){
        data$sampleTable1 = data$samr.assess.samplesize.obj$results[,,1]
        data$sampleTable2 = data$samr.assess.samplesize.obj$results[,,2]
        data$sampleTable3 = data$samr.assess.samplesize.obj$results[,,3]
        data$sampleTable4 = data$samr.assess.samplesize.obj$results[,,4]
        
        dataSample1 = paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[1] * samr.assess.samplesize.obj$n, sep = "")
        dataSample2 = paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[2] * samr.assess.samplesize.obj$n, sep = "")
        dataSample3 = paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[3] * samr.assess.samplesize.obj$n, sep = "")
        dataSample4 = paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[4] * samr.assess.samplesize.obj$n, sep = "")
        dataname = c(dataname, c(dataSample1, dataSample2, dataSample3, dataSample4))
      }
      
      writeWorksheetToFile(fname, data = data, sheet = dataname)
      file.rename(fname, file)

      }
    )

})

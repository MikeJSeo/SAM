#file:///C:/Users/mike/Desktop/MAIN/Examples/

library(XLConnect)
library(samr)
library(GSA)

source("GSA.listsets.revised.R")
source("GSA.correlate.revised.R")
source("GSA.plot.revised.R")

shinyServer(function(input, output) {  
  
  output$geneSetInfo = renderPrint({
    
    GSA = getGSA()
    if(!is.null(GSA)){      
      geneset.obj = GSA$geneset.obj
      genenames = GSA$genenames
      GSA.correlate.revised(geneset.obj, genenames)
    }
    
    
  })
  
  output$geneSetInfoText <- renderText({
    if(!is.null(getGSAList())){
      "Information for gene set collection"
    }
  })
  
  output$GSAPlotText <- renderText({
    if(!is.null(getGSAList())){
      "GSA Plot"
    }
  })
  
  output$GeneScoreText <- renderText({
    if(!is.null(getGSAList())){
      "Gene Score"
    }
  })
  
  output$geneSetPositive = renderText({
        
    GSA.list = getGSAList()
    if(!is.null(GSA.list$nsets.pos)){
      paste("Positive Gene Sets (", GSA.list$nsets.pos, ")", sep = "")      
    }
  })
  
  output$geneSetNegative = renderText({
    
    GSA.list = getGSAList()
    if(!is.null(GSA.list$nsets.neg)){
      paste("Negative Gene Sets (", GSA.list$nsets.neg, ")", sep = "")      
    }
  })
  
  
  output$geneSetFullPositive = renderText({
    
    GSAFullList = getGSAFullList()
    if(!is.null(GSAFullList)){
      paste("Positive Gene Sets (", GSAFullList$nsets.pos, ")", sep = "")      
    }
  })
  
  output$geneSetFullNegative = renderText({
    
    GSAFullList = getGSAFullList()
    if(!is.null(GSAFullList)){
      paste("Negative Gene Sets (", GSAFullList$nsets.neg, ")", sep = "")      
    }
  })
  
  
  output$geneSetPlot = renderPlot({
    
    GSA = getGSA()
    
    if(!is.null(GSA)){
      GSA.obj = GSA$GSA.obj
      GSA.plot.revised(GSA.obj, FDRcut = findFDR(), fac = 0)
    }
  })
  
  output$positiveGeneSet = renderTable({
    
    GSA.list = getGSAList()
    
    if(!is.na(GSA.list$positive[1])){
      GSA.list$positive
    }    
  })
  
  output$negativeGeneSet = renderTable({
    
    GSA.list = getGSAList()
    
    if(!is.na(GSA.list$negative[1])){
      GSA.list$negative      
    }    
  })
  
  output$positiveFullGeneSet = renderTable({
    
    GSAFullList = getGSAFullList()
    
    if(!is.null(GSAFullList)){
      GSAFullList$positive
    }    
  })
  
  output$negativeFullGeneSet = renderTable({
    
    GSAFullList = getGSAFullList()
    
    if(!is.null(GSAFullList)){
      GSAFullList$negative
    }    
  })
  
  
  output$testGeneSet = renderTable({
    
    GSA = getGSA()
    genesets = GSA$genesets
    GSA.obj = GSA$GSA.obj
    genenames = GSA$genenames
    
    if(!is.null(GSA)){
       GSA.genescores(input$geneset.number, genesets, GSA.obj, genenames)
    }    
  })
  
  getGSAList = reactive({
    GSA = getGSA()
    
    if(!is.null(GSA)){
      GSA.obj = GSA$GSA.obj
      geneset.names = GSA$geneset.names
      GSA.list = GSA.listsets.revised(GSA.obj, geneset.names = geneset.names, FDRcut = findFDR())
      GSA.list
    }
  })
  
  getGSAFullList = reactive({
    GSA = getGSA()
    
    if(!is.null(GSA)){
      GSA.obj = GSA$GSA.obj
      geneset.names = GSA$geneset.names
      GSA.list = GSA.listsets.revised(GSA.obj, geneset.names = geneset.names, FDRcut = 1)
      GSA.list
    }
    
  })
  
  getGSA = reactive({
    
    input$goButton
    
    isolate({
    
    data = getData()
    gmtFile = input$gmtFile
      
    if(!is.null(data) && !is.null(gmtFile)){
      
      x = data$x
      y = data$y
      data = getData()
      genenames = data$genenames
      
      geneset.obj = GSA.read.gmt(gmtFile$datapath)
      genesets = geneset.obj$genesets
      geneset.names = geneset.obj$geneset.names
      
      GSA.obj = GSA(x,y, genenames=genenames, genesets=genesets, method = "maxmean", resp.type="Two class unpaired", nperms=100, minsize = input$minGeneSet, maxsize = input$maxGeneSet)
      list(GSA.obj  = GSA.obj, geneset.names = geneset.names, genenames = genenames, genesets = genesets, geneset.obj = geneset.obj)
    }
    })
  })
  
  findFDR = reactive({
    
    if(input$fdrChoice == "FDR Slider")
      input$fdrSlider
    else if (input$fdrChoice == "Manually Enter FDR")
      input$fdrInput
  })


  findDelta = reactive({
    
    if(input$deltaChoice == "Delta Slider")
      input$deltaSlider
    else if (input$deltaChoice == "Manually Enter Delta")
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
  
  output$allPositiveGenesText <- renderText({
    
    allgenes = getAllgenesTable()
    
    if(!is.null(allgenes)){
      paste("Positive Genes (", allgenes$ngenes.up, ")", sep = "")
    }
  })
  
  
  
  output$allNegativeGenesText = renderText({
    
    allgenes = getAllgenesTable()
    
    if(!is.null(allgenes)){
      paste("Negative Genes (", allgenes$ngenes.lo, ")", sep = "")
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
    
      siggenes.table = samr.compute.siggenes.table(samr.obj,delta, data, delta.table, min.foldchange = min.foldchange, compute.localfdr= input$localFDR)
      siggenes.table
    }
  })
  
  
  getAllgenesTable = reactive({
    
    if(!is.null(getResult())){
      result = getResult()
      delta = findDelta()
      samr.obj = result$samr.obj
      data = result$data
      delta.table = result$delta.table
      min.foldchange = input$min.foldchange
      
      Allgenes = samr.compute.siggenes.table(samr.obj,delta, data, delta.table, min.foldchange = min.foldchange, compute.localfdr= input$localFDR, all.genes= TRUE)
      Allgenes
    }
  })
  
  
  getData = reactive({
    
    objFile = input$iFile
    if(!is.null(objFile)){
      wb = loadWorkbook(objFile$datapath)
      sheets = getSheets(wb)
      dat = readWorksheet(wb, sheets, header = FALSE)
      x = dat[-1, c(-1,-2)]
      
      x = as.matrix(x)
      class(x) = "numeric"
      
      geneid = dat[-1,1]
      genenames = dat[-1,2]
        
      firstrow = as.vector(dat[1,c(-1,-2)])
      censoring.status = NULL
      eigengene.number = NULL
        
      if( (input$responseType_seq == "Survival" && input$assayType == "seq") || (input$responseType_array == "Survival" && input$assayType == "array")){
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
      data
    }
    
  })
  
  getResult = reactive({
    
    input$goButton
    
    isolate({
    data = getData()

    if(!is.null(data)){
      resp.type = if(input$assayType == "seq"){input$responseType_seq}else{input$responseType_array} 
      s0.perc = if(is.na(input$s0.perc) || input$s0 == "Automatic"){NULL}else{input$s0.perc}
      center.arrays = as.logical(input$centerArrays)
        
      samr.obj = samr(data, resp.type = resp.type, assay.type = input$assayType, s0.perc = NULL, nperms = input$nperms, center.arrays = center.arrays, testStatistic = input$testStatistic, time.summary.type = input$timeSummaryType,  regression.method = input$regressionMethod, random.seed = input$random.seed)           
      delta.table = samr.compute.delta.table(samr.obj, min.foldchange = input$min.foldchange)
        
      list (data = data, samr.obj = samr.obj, delta.table = delta.table)  
    }
    
    })
  })

  
  output$sampleTable1 = renderTable({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(samr.assess.samplesize.obj)){
      samr.assess.samplesize.obj$results[,,1] 
    }

  })
  
  output$sampleTable2 = renderTable({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(samr.assess.samplesize.obj)){
      samr.assess.samplesize.obj$results[,,2] 
    }    
    
  })
  
  output$sampleTable3 = renderTable({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(samr.assess.samplesize.obj)){
      samr.assess.samplesize.obj$results[,,3]  
    }  
  })
  
  output$sampleTable4 = renderTable({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(samr.assess.samplesize.obj)){
      samr.assess.samplesize.obj$results[,,4]  
    }
  })
  
  output$samplePlot = renderPlot({
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(samr.assess.samplesize.obj)){
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
  
  
  
  output$deltaTable = renderTable({
    
    result = getResult()
    delta.table = result$delta.table
    delta.table
  
  })
  
  output$Allgenes.table.up = renderDataTable({
    Allgenes = getAllgenesTable()  
    
    if(!is.null(Allgenes$genes.up)){
      Allgenes$genes.up
    }
  })
  
  output$Allgenes.table.lo = renderDataTable({
    Allgenes = getAllgenesTable()  
    
    if(!is.null(Allgenes$genes.lo)){
      Allgenes$genes.lo
    }
  })
  
  
  output$siggenes.table.up <- renderDataTable({
    
    siggenes.table = getSiggenesTable()

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
  
  
  output$allGenes = renderDataTable({
    
    result = getResult()
    
    if(!is.null(result)){     
      data = result$data
      x = cbind(data$geneid, data$genenames, data$x)
      x
    }
  })

  output$downloadData <- downloadHandler(
    filename = function() { "result.xlsx" },
    content = function(file) {
      siggenes.table = getSiggenesTable()
      result = getResult()
      delta.table = result$delta.table
      samr.assess.samplesize.obj =  getSampleSize()
      Allgenes = getAllgenesTable()
      
      fname = paste(file, "xlsx", sep = ".")
      wb = loadWorkbook(fname, create = TRUE)
    
      samr.obj = result$samr.obj
      delta = findDelta()
      min.foldchange = input$min.foldchange
      
      png(file = "SAMPlot.png")
      samr.plot(samr.obj, delta, min.foldchange = min.foldchange)
      dev.off()
      
      if(!is.null(GSA)){
        createSheet(wb, name = "SAMPlot")
        createName(wb, name = "SAMPlot", formula = "SAMPlot!$B$2")
        addImage(wb, filename = "SAMPlot.png", name = "SAMPlot", originalSize = TRUE) 
      }
            
      if(!is.null(delta.table)){        
        createSheet(wb, name = "Delta Table")
        writeWorksheet(wb, delta.table, sheet = "Delta Table")            
      }
      if(!is.null(siggenes.table$genes.up)){
        createSheet(wb, name = "Significant Positive Genes")
        writeWorksheet(wb, siggenes.table$genes.up, sheet = "Significant Positive Genes")
      }
      if(!is.null(siggenes.table$genes.lo)){
        createSheet(wb, name = "Significant Negative Genes")
        writeWorksheet(wb, siggenes.table$genes.lo, sheet = "Significant Negative Genes")
      }
      
      if(!is.null(Allgenes$genes.up)){
        createSheet(wb, name = "All Positive Genes")
        writeWorksheet(wb, Allgenes$genes.up, sheet = "All Positive Genes")
      }
      
      if(!is.null(Allgenes$genes.lo)){
        createSheet(wb, name = "All Negative Genes")
        writeWorksheet(wb, Allgenes$genes.lo, sheet = "All Negative Genes")
      }
      
      if(!is.null(samr.assess.samplesize.obj)){        
        dataSample1 = paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[1] * samr.assess.samplesize.obj$n, sep = "")
        dataSample2 = paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[2] * samr.assess.samplesize.obj$n, sep = "")
        dataSample3 = paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[3] * samr.assess.samplesize.obj$n, sep = "")
        dataSample4 = paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[4] * samr.assess.samplesize.obj$n, sep = "")
        
        createSheet(wb, name = dataSample1)
        createSheet(wb, name = dataSample2)
        createSheet(wb, name = dataSample3)
        createSheet(wb, name = dataSample4)
        
        writeWorksheet(wb, samr.assess.samplesize.obj$results[,,1], sheet = dataSample1)
        writeWorksheet(wb, samr.assess.samplesize.obj$results[,,2], sheet = dataSample2)
        writeWorksheet(wb, samr.assess.samplesize.obj$results[,,3], sheet = dataSample3)
        writeWorksheet(wb, samr.assess.samplesize.obj$results[,,4], sheet = dataSample4)
      }
      
      saveWorkbook(wb)
      file.rename(fname, file)

      }
    )
  
  #download button for gene set data
  output$downloadData2 <- downloadHandler(
    filename = function() { "result.xlsx" },
    content = function(file) {
      
      GSA = getGSA()
      GSA.obj = GSA$GSA.obj
      GSA.list = getGSAList()
      GSAFullList = getGSAFullList()
      
      fname = paste(file, "xlsx", sep = ".")
      wb = loadWorkbook(fname, create = TRUE)
      
      png(file = "GSAPlot.png")
      GSA.plot.revised(GSA.obj, FDRcut = findFDR(), fac = 0)
      dev.off()
      
      if(!is.null(GSA)){
        createSheet(wb, name = "GSAPlot")
        createName(wb, name = "GSAPlot", formula = "GSAPlot!$B$2")
        addImage(wb, filename = "GSAPlot.png", name = "GSAPlot", originalSize = TRUE) 
      }
      
      if(!is.na(GSA.list$positive[1])){
        createSheet(wb, name = "Significant Positive Gene Sets")
        writeWorksheet(wb, GSA.list$positive, sheet = "Significant Positive Gene Sets")            
      }
      
      if(!is.na(GSA.list$negative[1])){
        createSheet(wb, name = "Significant Negative Gene Sets")
        writeWorksheet(wb, GSA.list$negative, sheet = "Significant Negative Gene Sets")      
      }
      
      if(!is.null(GSAFullList)){
        createSheet(wb, name = "Full Positive Gene Sets")
        writeWorksheet(wb, GSAFullList$positive, sheet = "Full Positive Gene Sets")
      }
      
      if(!is.null(GSAFullList)){
        createSheet(wb, name = "Full Negative Gene Sets")
        writeWorksheet(wb, GSAFullList$negative, sheet = "Full Negative Gene Sets")
      }
      
      saveWorkbook(wb)
      file.rename(fname, file)
      
    }
  )

})

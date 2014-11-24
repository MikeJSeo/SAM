#file:///C:/Users/mike/Desktop/MAIN/Examples/
options(shiny.maxRequestSize=10000*1024^2)

library(XLConnect)
library(samr)
library(GSA)

source("GSA.listsets.revised.R")
source("GSA.correlate.revised.R")
source("GSA.plot.revised.R")

shinyServer(function(input, output) {  
 
  
  getGeneSetTable = reactive({
    
    GSA = getGSA()
    if(!is.null(GSA)){      
      geneset.obj = GSA$geneset.obj
      genenames = GSA$genenames
      GSA.correlate.list = GSA.correlate.revised(geneset.obj, genenames)
      GSA.correlate.list$result
    }
  })

  output$geneSetTable = renderTable({
    geneSetTable = getGeneSetTable()
    if(!is.null(geneSetTable)){
      geneSetTable
    }
  })

  getGeneSetQuantile = reactive({
    GSA = getGSA()
    if(!is.null(GSA)){      
      geneset.obj = GSA$geneset.obj
      genenames = GSA$genenames
      GSA.correlate.list = GSA.correlate.revised(geneset.obj, genenames)
      GSA.correlate.list$QuantileCoverage
    }
    
  })
  
  output$geneSetQuantile = renderTable({
    geneSetQuantile = getGeneSetQuantile()
    if(!is.null(geneSetQuantile)){
      geneSetQuantile
    }
    
  })
  
  getGeneSetTableGenes = reactive({
    GSA = getGSA()
    if(!is.null(GSA)){      
      geneset.obj = GSA$geneset.obj
      genenames = GSA$genenames
      GSA.correlate.list = GSA.correlate.revised(geneset.obj, genenames)
      GSA.correlate.list$tableGenes
    }
    
  })
  
  output$geneSetTableGenes = renderTable({
    geneSetTableGenes = getGeneSetTableGenes()
    if(!is.null(geneSetTableGenes)){
      geneSetTableGenes
    }
 
  })
  
  output$geneSetInfoText <- renderText({
    if(!is.null(getGSAList())){
      "Information for gene set collection"
    }
  })
  
  output$geneSetInfoText2 <- renderText({
    if(!is.null(getGSAList())){
      "Quantiles of fraction coverage of gene-sets"
    }
  })
  
  output$geneSetInfoText3 <- renderText({
    if(!is.null(getGSAList())){
      "Table of number of genes in genesets"
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
      GSA.plot.revised(GSA.obj, FDRcut = 1, fac = 0)
    }
  })
  
  output$positiveGeneSet = renderTable({
    
    GSA.list = getGSAList()
    
    if(!is.null(GSA.list$positive)){
      if(!is.na(GSA.list$positive[1])){
        GSA.list$positive  
      }
    }    
  })
  
  output$negativeGeneSet = renderTable({
    
    GSA.list = getGSAList()
    if(!is.null(GSA.list$negative)){
      if(!is.na(GSA.list$negative[1])){
        GSA.list$negative  
      }  
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
    
    input$goButton2
    
    isolate({
    
    data = getData()
    gmtFile = input$gmtFile
      
    if(!is.null(data) && !is.null(gmtFile)){
      
      x = data$x
      y = data$y
      
      class(y) = "numeric"
      
      genenames = data$genenames
      censoring.status = data$censoring.status
      
      geneset.obj = GSA.read.gmt(gmtFile$datapath)
      genesets = geneset.obj$genesets
      geneset.names = geneset.obj$geneset.names
      
      s0.perc = if(is.na(input$s0.perc) || input$s0 == "Automatic"){NULL}else{input$s0.perc}   
      
      GSA.obj = GSA(x, y, genenames=genenames, genesets=genesets, method = "maxmean", resp.type= input$responseType_array, minsize = input$minGeneSet, maxsize = input$maxGeneSet, random.seed = input$random.seed, knn.neighbors = input$numberOfNeighbors, s0.perc = s0.perc, nperms = input$nperms, censoring.status = censoring.status)
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
  
  output$missRateText = renderText({
    if(!is.null(getResult())){
      paste("Estimated Miss rates for Delta = ", findDelta())
    }
    
  })
  
  output$eigengeneText = renderText({
    if(!is.null(getResult()) && input$responseType_array == "Pattern discovery"){
      "Eigengene"
    }
    
  })
  
  getEigengene = reactive({
    if(!is.null(getResult()) && input$responseType_array == "Pattern discovery"){
      result = getResult()
      samr.obj = result$samr.obj
      samr.obj$eigengene
    }
  })
  
  output$eigengene = renderTable({
    eigengene = getEigengene()
    if(!is.null(eigengene)){
      eigengene
    }
  })
  
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
  
  output$inputParametersText <- renderText({
    if(!is.null(getResult())){
      "Input Parameters"
    }
  })
  output$computedValuesText <- renderText({
    if(!is.null(getResult())){
      "Computed Values"
    }
  })
  
  output$inputParametersText2 <- renderText({
    if(!is.null(getGSA())){
      "Input Parameters"
    }
  })
  output$computedValuesText2 <- renderText({
    if(!is.null(getGSA())){
      "Computed Values"
    }
  })
  
  output$sampleSizePlotText2 <- renderText({
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
      y = data$y
      dif = input$dif
      if(!substring(y[1],2,6) %in% c("Block", "block")){
        samplesize.factors = input$sampleSizeFactors
        samplesize.factors = as.numeric(unlist(strsplit(samplesize.factors, ",")))
        samr.assess.samplesize.obj =  samr.assess.samplesize(samr.obj, data, dif, samplesize.factors)
        samr.assess.samplesize.obj
      }
    }
  })
  
  getMissRateTable = reactive({
    if(!is.null(getResult())){
      result = getResult()
      delta = findDelta()
      samr.obj = result$samr.obj
      delta.table = result$delta.table
      
      missrate.table = samr.missrate(samr.obj, delta, delta.table)
      missrate.table
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
    
      ptn <- "\\.[[:alnum:]]{1,4}$"
      Suf <- c(".txt", ".csv", ".xls", ".xlsx")
      suf <- tolower(regmatches(objFile, regexpr(ptn, objFile)))
      
      if(suf == ".xls" || suf == ".xlsx"){
        wb = loadWorkbook(objFile$datapath)      
        #read just first sheet of the file
        dat = readWorksheet(wb, 1, header = FALSE)
      }
      if(suf == ".csv"){
        dat = read.csv(objFile$datapath, header = FALSE)
      }
      
      geneid = dat[-1,1]
      genenames = dat[-1,2]
      
      imputedX = NULL
      x = dat[-1, c(-1,-2)]
      x = as.matrix(x)
      class(x) = "numeric"
      if(sum(is.na(x)) > 0){
        
        imputedummy = impute.knn(x, k = input$numberOfNeighbors)
        x = imputedummy$data
        imputedX = cbind(geneid, genenames, x)
        colnamesImputedX = as.vector(dat[1,])
        colnamesImputedX[1] = " "
        colnamesImputedX[2] = " "
        colnames(imputedX) = colnamesImputedX
      }
        
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

        
      data =list(x=x,y=y, genenames=genenames, geneid=geneid, logged2= as.logical(input$dataLogged), censoring.status = censoring.status, eigengene.number = eigengene.number, imputedX = imputedX)        
      data
    }
    
  })
  
  getImputedX = reactive({
    
    data = getData()
    if(!is.null(data)){
      imputedX = data$imputedX
      imputedX
    }
  })
  
  getResult = reactive({
    
    input$goButton
    input$min.foldchange
    
    isolate({
    data = getData()

    if(!is.null(data)){
      resp.type = if(input$assayType == "seq"){input$responseType_seq}else{input$responseType_array} 
      s0.perc = if(is.na(input$s0.perc) || input$s0 == "Automatic"){NULL}else{input$s0.perc}
      center.arrays = as.logical(input$centerArrays)
        
      samr.obj = samr(data, resp.type = resp.type, assay.type = input$assayType, s0.perc = s0.perc, nperms = input$nperms, center.arrays = center.arrays, testStatistic = input$testStatistic, time.summary.type = input$timeSummaryType,  regression.method = input$regressionMethod, random.seed = input$random.seed, knn.neighbors = input$numberOfNeighbors)           
      delta.table = samr.compute.delta.table(samr.obj, min.foldchange = input$min.foldchange)
        
      list (data = data, samr.obj = samr.obj, delta.table = delta.table)  
    }
    
    })
  })

  output$missRate = renderTable({
    
    if(!is.null(getMissRateTable())){
      getMissRateTable()      
    }
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
  
  capitalize = function(word){
    letters = strsplit(word,'')
    theletters = letters[[1]]
    theletters[1] = toupper(theletters[1])
    paste(theletters,collapse='')
  }
  
  output$computedValues = renderTable({
    if(!is.null(getComputedValues())){
      getComputedValues()      
    }
    
  })
  
  output$computedValues2 = renderTable({
    if(!is.null(getComputedValues2())){
      getComputedValues2()      
    }
    
  })
  
  getTailStrength = reactive({
    result = getResult()
    if(!is.null(result)){
      samr.obj = result$samr.obj
      samr.tail.strength(samr.obj)
    }    
  })
  
  getComputedValues2 = reactive({
    GSA = getGSA()
    if(!is.null(GSA)){
      GSA.obj = GSA$GSA.obj  
      computedValues = matrix(NA, nrow = 2, ncol = 1)
      colnames(computedValues) = "Values"
      rownames(computedValues) = c("Exchangibility factor s0", "s0 percentile")
        
      computedValues[1,1] = GSA.obj$s0  
      computedValues[2,1] = GSA.obj$s0.perc
      computedValues
    }
    
    
  })
  
  getComputedValues= reactive({
    
    result = getResult()
    if(!is.null(result)){
      
      computedValues = matrix(NA, nrow = 5, ncol = 1)
      
      colnames(computedValues) = "Values"
      rownames(computedValues) = c("Estimated proportion of non-null features (genes)", "s0 percentile", "False Discovery Rate", "Estimated tail strength", "Estimated standard error of tail strength")
      samr.obj = result$samr.obj
      tail.strength = getTailStrength()
      
      computedValues[1,1] = samr.obj$pi0  
      computedValues[2,1] = samr.obj$s0.perc
      delta.table = samr.compute.delta.table(samr.obj, min.foldchange= input$min.foldchange, dels= findDelta())[5]
      computedValues[3,1] = delta.table
      computedValues[4,1] = tail.strength$ts
      computedValues[5,1] = tail.strength$se.ts
      
      if(input$assayType == "array" && input$analysisType == "Standard" && input$responseType_array == "Multiclass"){
        multiclassAdd = matrix(NA, nrow = 2, ncol = 1)
        rownames(multiclassAdd) = c("2.5% null value for multiclass contrasts", "97.5% null value for multiclass contrasts")
        multiclassAdd[1,1] = samr.obj$stand.contrasts.95[1]
        multiclassAdd[2,1] = samr.obj$stand.contrasts.95[2]
        computedValues = rbind(computedValues, multiclassAdd)
      }
      
      if(input$assayType == "seq"){
        depth = samr.obj$depth
        seqAdd = matrix(NA, nrow = 1, ncol = 1)
        rownames(seqAdd) = "Quantiles of estimated sequencing depths (0%, 25%, 50%, 75%, 100%)"   
        seqAdd[1,1] = paste(quantile(depth, 0)[[1]], quantile(depth, 0.25)[[1]], quantile(depth, 0.5)[[1]], quantile(depth, 0.75)[[1]], quantile(depth, 1)[[1]],  sep = ", ")  
        computedValues = rbind(computedValues, seqAdd)
      }
      computedValues
      
      
    }
    
  })
  
  output$inputParameters = renderTable({
    
    if(!is.null(getInputParameters())){
      getInputParameters()
    }
    
  })
  
  output$inputParameters2 = renderTable({
    
    if(!is.null(getInputParameters2())){
      getInputParameters2()
    }
    
  })
  
  getInputParameters2 = reactive({
    
    input$goButton2
    findFDR()
        
    isolate({
      
      if(input$goButton2!= 0){
        current = matrix(NA, nrow = 10, ncol = 1)
        
        objFile = input$iFile
        gmtFile = input$gmtFile  
        s0.perc = if(is.na(input$s0.perc) || input$s0 == "Automatic"){"Automatic"}else{paste(input$s0.perc, " percentile")}
        
        current[1,1] = objFile$name
        current[2,1] = gmtFile$name
        current[3,1] = input$responseType_array
        current[4,1] = findFDR()
        current[5,1] = input$nperms
        current[6,1] = s0.perc
        current[7,1] = input$numberOfNeighbors
        current[8,1] = input$minGeneSet
        current[9,1] = input$maxGeneSet
        current[10,1] = input$random.seed
             
        rownames_current = c("File Name", "Gene set (gmt) file", "Data Type", "False discovery rate (in each tail)", "Number of permutations", "Input percentile for exchangeability factor s0", "Number of neighbors for KNN", "Minimum gene set size", "Maximum gene set size", "Seed for Random number generator")  
        
        rownames(current) = rownames_current
        colnames(current) = "value"
        current
      }
        
    })
  })
  
  getInputParameters = reactive({
    
    input$goButton
    findDelta()
    input$min.foldchange
    
    isolate({
    
    if(input$goButton!= 0){
        
      
    if(input$assayType == "array"){
      current = matrix(NA, nrow = 14, ncol = 1)
      
      objFile = input$iFile
      s0.perc = if(is.na(input$s0.perc) || input$s0 == "Automatic"){"Automatic"}else{paste(input$s0.perc, " percentile")}
      
      current[1,1] = objFile$name
      current[2,1] = input$responseType_array
      current[3,1] = capitalize(input$assayType)
      current[4,1] = input$centerArrays
      current[5,1] = findDelta()
      current[6,1] = input$min.foldchange
      current[7,1] = if(input$testStatistic == "standard"){"T-statistic"} else{"Wilcoxon"}
      current[8,1] = capitalize(input$regressionMethod)
      current[9,1] = input$dataLogged 
      current[10,1] = input$nperms
      current[11,1] = s0.perc
      current[12,1] = input$numberOfNeighbors
      current[13,1] = input$random.seed
      current[14,1] = capitalize(input$timeSummaryType)
      
      rownames_current = c("File Name", "Data Type", "Array or Seq data?", "Arrays centered?", "Delta", "Minimum fold change", "Test statistic", "Regression method", "Are data are log scale?", "Number of permutations", "Input percentile for exchangeability factor s0", "Number of neighbors for KNN", "Seed for Random number generator", "Time summary type")  
      
      if(input$responseType_array == "Quantitative"){
        current = matrix(current[c(-6, -7, -9, -14),], ncol = 1)
        rownames_current = rownames_current[c(-6,-7, -9, -14)]
      }
      else if(input$responseType_array == "Two class unpaired"){
        current = matrix(current[c(-8, -14),], ncol = 1)
        rownames_current = rownames_current[c(-8, -14)]
      }
      else if(input$responseType_array == "Survival"){
        current = matrix(current[c(-6, -7,-8,-9,-14),], ncol = 1)
        rownames_current = rownames_current[c(-6, -7,-8,-9, -14)]
      }
      else if(input$responseType_array == "Multiclass" || input$responseType_array == "One class" || input$responseType_array == "Pattern discovery"){
        current = matrix(current[c(-6,-7,-8,-9,-14),], ncol = 1)
        rownames_current = rownames_current[c(-6,-7,-8,-9, -14)]
      }
      else if(input$responseType_array == "Two class paired"){
        current = matrix(current[c(-7,-8,-14),], ncol = 1)
        rownames_current = rownames_current[c(-7,-8,-14)]
      }
      else if(input$responseType_array == "Two class unpaired timecourse"){
        current = matrix(current[c(-6, -8, -9),], ncol = 1)
        rownames_current = rownames_current[c(-6, -8, -9)]
      }
      else if(input$responseType_array == "Two class paired timecourse"){
        current = matrix(current[c(-6, -7, -8, -9),], ncol = 1)
        rownames_current = rownames_current[c(-6, -7, -8, -9)]
      }
      else if(input$responseType_array == "One class timecourse"){
        current = matrix(current[c(-6, -7, -8, -9),], ncol = 1)
        rownames_current = rownames_current[c(-6, -7, -8, -9)]
      }
      
      rownames(current) = rownames_current
      colnames(current) = "value"
    }
    
    if(input$assayType == "seq"){
      
      current = matrix(NA, nrow = 9, ncol = 1)
      
      objFile = input$iFile
      
      current[1,1] = objFile$name
      current[2,1] = input$responseType_seq
      current[3,1] = capitalize(input$assayType)
      current[4,1] = input$centerArrays
      current[5,1] = findDelta()
      current[6,1] = input$min.foldchange
      current[7,1] = input$nperms
      current[8,1] = input$numberOfNeighbors
      current[9,1] = input$random.seed
      
      rownames_current = c("File Name", "Data Type", "Array or Seq data?", "Arrays centered?", "Delta", "Minimum fold change", "Number of permutations", "Number of neighbors for KNN", "Seed for Random number generator")  

      if(input$responseType_seq == "Quantitative" || input$responseType_seq == "Survival" || input$responseType_seq == "Multiclass"){
        current = matrix(current[c(-6),], ncol = 1)
        rownames_current = rownames_current[c(-6)]
      }
      
      rownames(current) = rownames_current
      colnames(current) = "value"
    }
    
    current
    }
    })
  })
  
  savethis = observe({
  
    if(input$saveButton != 0){
      
    isolate({  
    
    dir = input$dir
    file = input$fname
    
    siggenes.table = getSiggenesTable()
    missrate.table = getMissRateTable()
    result = getResult()
    delta.table = result$delta.table
    samr.assess.samplesize.obj =  getSampleSize()
    Allgenes = getAllgenesTable()
    
    samr.obj = result$samr.obj
    delta = findDelta()
    min.foldchange = input$min.foldchange
    
    inputParameters = getInputParameters()
    computedValues = getComputedValues()
    eigengene = getEigengene()
    imputedX = getImputedX()
    
    fname = paste(file, "xlsx", sep = ".")
    
    if(file.exists(file.path(dir, fname))){
      file.remove(file.path(dir, fname))      
    }
    wb = loadWorkbook(file.path(dir,fname), create = TRUE)
    
    if(!is.null(samr.obj)){
      png(file = "SAMPlot.png")
      samr.plot(samr.obj, delta, min.foldchange = min.foldchange)
      dev.off()      
    }
    
    samr.assess.samplesize.obj =  getSampleSize()
    if(!is.null(samr.assess.samplesize.obj)){
      png(file = "samplesizePlot.png")
      samr.assess.samplesize.plot(samr.assess.samplesize.obj)  
      dev.off()  
    }
    
    if(!is.null(samr.obj)){      
      createSheet(wb, name = "SAMPlot")
      createName(wb, name = "SAMPlot", formula = "SAMPlot!$B$2")
      addImage(wb, filename = "SAMPlot.png", name = "SAMPlot", originalSize = TRUE) 
      if (file.exists("SAMPlot.png")) file.remove("SAMPlot.png")
    }
    
    if(!is.null(delta.table)){        
      createSheet(wb, name = "Delta Table")
      writeWorksheet(wb, delta.table, sheet = "Delta Table")            
    }
    
    if(!is.null(missrate.table)){
      createSheet(wb, name = "Miss Rate Table")
      writeWorksheet(wb, missrate.table, sheet = "Miss Rate Table")      
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
      createSheet(wb, name = "sampleSizePlot")
      createName(wb, name = "sampleSizePlot", formula = "sampleSizePlot!$B$2")
      addImage(wb, filename = "samplesizePlot.png", name = "sampleSizePlot", originalSize = TRUE) 
      if (file.exists("samplesizePlot.png")) file.remove("samplesizePlot.png")
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
    
    if(!is.null(inputParameters)){
      createSheet(wb, name = "Input Parameters")
      newInputParameters = cbind(rownames(inputParameters), inputParameters)
      colnames(newInputParameters) = c("Questions","values")
      writeWorksheet(wb, newInputParameters, sheet = "Input Parameters")
    }
    
    if(!is.null(eigengene)){
      createSheet(wb, name = "Eigengene")
      writeWorksheet(wb, eigengene, sheet = "Eigengene")
    }
    
    if(!is.null(imputedX)){
      createSheet(wb, name = "Imputed Data")
      writeWorksheet(wb, imputedX, sheet = "Imputed Data")
    }
    
    if(!is.null(computedValues)){
      createSheet(wb, name = "Computed Values")
      newComputedValues = cbind(rownames(computedValues), computedValues)
      colnames(newComputedValues) = c("Questions","values")
      
      writeWorksheet(wb, newComputedValues, sheet = "Computed Values")
      
    }
    
    if(!is.null(result)){
      saveWorkbook(wb, file.path(dir, fname))
    # shell(file.path(dir, fname))
    }
    
    
    })
    }
      
  })

  
  savethis2 = observe({
    
    if(input$saveButton2 != 0){
      
      isolate({  
        
        dir = input$dir2
        file = input$fname2
        
        GSA = getGSA()
        GSA.obj = GSA$GSA.obj
        GSA.list = getGSAList()
        GSAFullList = getGSAFullList()
        
        geneSetTable = getGeneSetTable()
        geneSetQuantile = getGeneSetQuantile()
        geneSetTableGenes = getGeneSetTableGenes()
        
        inputParameters = getInputParameters2()
        computedValues = getComputedValues2()
        
        fname = paste(file, "xlsx", sep = ".")
        
        if(file.exists(file.path(dir, fname))){
          file.remove(file.path(dir, fname))      
        }
        wb = loadWorkbook(file.path(dir,fname), create = TRUE)
        
        if(!is.null(GSA.obj)){
          png(file = "GSAPlot.png")
          GSA.plot.revised(GSA.obj, FDRcut = 1, fac = 0)
          dev.off()    
        }
        
        if(!is.null(GSA)){
          createSheet(wb, name = "GSAPlot")
          createName(wb, name = "GSAPlot", formula = "GSAPlot!$B$2")
          addImage(wb, filename = "GSAPlot.png", name = "GSAPlot", originalSize = TRUE) 
          if (file.exists("GSAPlot.png")) file.remove("GSAPlot.png")
        }
        
        if(!is.null(GSA.list$positive[1])){
          createSheet(wb, name = "Significant Positive Gene Sets")
          writeWorksheet(wb, GSA.list$positive, sheet = "Significant Positive Gene Sets")            
        }

        if(!is.null(GSA.list$negative[1])){
          createSheet(wb, name = "Significant Negative Gene Sets")
          writeWorksheet(wb, GSA.list$negative, sheet = "Significant Negative Gene Sets")      
        }
        
        if(!is.null(GSAFullList$positive)){
          createSheet(wb, name = "Full Positive Gene Sets")
          writeWorksheet(wb, GSAFullList$positive, sheet = "Full Positive Gene Sets")
        }
        
        if(!is.null(GSAFullList$negative)){
          createSheet(wb, name = "Full Negative Gene Sets")
          writeWorksheet(wb, GSAFullList$negative, sheet = "Full Negative Gene Sets")
        }
        
        if(!is.null(geneSetTable)){
          createSheet(wb, name = "Gene set collection")
          newGeneSetTable = cbind(rownames(geneSetTable), geneSetTable)
          colnames(newGeneSetTable) = c("Questions","values")
          writeWorksheet(wb, newGeneSetTable, sheet = "Gene set collection")          
        }
        if(!is.null(geneSetQuantile)){
          createSheet(wb, name = "Quantiles of fraction coverage")
          writeWorksheet(wb, geneSetQuantile, sheet = "Quantiles of fraction coverage")          
        }
        if(!is.null(geneSetTableGenes)){
          createSheet(wb, name = "Number of genes")
          writeWorksheet(wb, geneSetTableGenes, sheet = "Number of genes")                    
        }
        
        if(!is.null(inputParameters)){
          createSheet(wb, name = "Input Parameters")
          newInputParameters = cbind(rownames(inputParameters), inputParameters)
          colnames(newInputParameters) = c("Questions","values")
          writeWorksheet(wb, newInputParameters, sheet = "Input Parameters")
        }
        
        if(!is.null(computedValues)){
          createSheet(wb, name = "Computed Values")
          newComputedValues = cbind(rownames(computedValues), computedValues)
          colnames(newComputedValues) = c("Questions","values")
          
          writeWorksheet(wb, newComputedValues, sheet = "Computed Values")
          
        }
        
        if(!is.null(GSA)){
          saveWorkbook(wb, file.path(dir, fname))
          # shell(file.path(dir, fname))
        }
        
      })
    }
    
  })
  
  

})

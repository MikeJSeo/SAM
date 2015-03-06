options(shiny.maxRequestSize=10000*1024^2)

library(openxlsx)
library(samr)
library(GSA)

source("GSA.listsets.revised.R")
source("GSA.correlate.revised.R")
source("GSA.plot.revised.R")

shinyServer(function(input, output) {  
 
  
  ##########Read uploaded data!
  
  getData = reactive({
    
    objFile = input$iFile
    if(!is.null(objFile)){
   
      dat = read.xlsx(objFile$datapath, 1, colNames = FALSE)   
      
      geneid = dat[-1,1]
      genenames = dat[-1,2]
      
      imputedX = NULL
      x = as.matrix(dat[-1, c(-1,-2)])
      class(x) = "numeric"
      
      originalX = NULL
      originalX = cbind(geneid, genenames, x)
      colnamesOriginalX = as.vector(dat[1,])
      colnamesOriginalX[1] = " "
      colnamesOriginalX[2] = " "
      colnames(originalX) = colnamesOriginalX
      
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
      
      data =list(x=x,y=y, genenames=genenames, geneid=geneid, logged2= as.logical(input$dataLogged), censoring.status = censoring.status, eigengene.number = eigengene.number, imputedX = imputedX, originalX = originalX)        
      data
    }
  })
  
  #########Modify changes for survival datasets
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
  
  ##########GSA ##############################
  
  ######get results from GSA
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
  
  getGSACorrelateList = reactive({
    GSA = getGSA()
    if(!is.null(GSA)){
      geneset.obj = GSA$geneset.obj
      genenames = GSA$genenames
      GSA.correlate.list = GSA.correlate.revised(geneset.obj, genenames)
      GSA.correlate.list
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
  
  findFDR = reactive({
    if(input$fdrChoice == "FDR Slider")
      input$fdrSlider
    else if (input$fdrChoice == "Manually Enter FDR")
      input$fdrInput
  })
  
  #########Output GSA tables/plots
  output$geneSetTable = renderTable({
    GSA.correlate.list = getGSACorrelateList()
    if(!is.null(GSA.correlate.list)){
      GSA.correlate.list$result
    }
  })
  
  output$geneSetQuantile = renderTable({
    GSA.correlate.list = getGSACorrelateList()
    if(!is.null(GSA.correlate.list)){
      GSA.correlate.list$QuantileCoverage
    }
  })
  
  output$geneSetTableGenes = renderTable({
    GSA.correlate.list = getGSACorrelateList()
    if(!is.null(GSA.correlate.list)){
      GSA.correlate.list$tableGenes
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
    if(!is.null(GSA)){
      genesets = GSA$genesets
      GSA.obj = GSA$GSA.obj
      genenames = GSA$genenames
      GSA.genescores(input$geneset.number, genesets, GSA.obj, genenames)
    }    
  })
  
  output$computedValues2 = renderTable({
    if(!is.null(getComputedValues2())){
      getComputedValues2()      
    }
  })
  
  output$inputParameters2 = renderTable({
    inputParameters2 = getInputParameters2()
    if(!is.null(inputParameters2)){
      inputParameters2
    }
  })
  
  ######output texts in the GSA interface
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
  
  #################################SAM
  
  ############get results from SAM
  getSamrObj = reactive({

    input$goButton
    
    isolate({
      data = getData()
      
      if(!is.null(data)){
        resp.type = if(input$assayType == "seq"){input$responseType_seq}else{input$responseType_array} 
        s0.perc = if(is.na(input$s0.perc) || input$s0 == "Automatic"){NULL}else{input$s0.perc}
        center.arrays = as.logical(input$centerArrays)
        
        samr.obj = samr(data, resp.type = resp.type, assay.type = input$assayType, s0.perc = s0.perc, nperms = input$nperms, center.arrays = center.arrays, testStatistic = input$testStatistic, time.summary.type = input$timeSummaryType,  regression.method = input$regressionMethod, random.seed = input$random.seed, knn.neighbors = input$numberOfNeighbors)           
        samr.obj
      }
      
    })
  })
  
  getDeltaTable = reactive({  
    samr.obj = getSamrObj()
    if(!is.null(samr.obj)){
      delta.table = samr.compute.delta.table(samr.obj, min.foldchange = input$min.foldchange)
      delta.table
    }
    
  })
  
  getEigengene = reactive({
    samr.obj = getSamrObj()
    if(!is.null(samr.obj) && input$responseType_array == "Pattern discovery"){      
      eigengene = samr.obj$eigengene
      eigengene
    }
  })
  
  getSampleSize = reactive({
    
    ar = isolate(input$responseType_array)
    seq = isolate(input$responseType_seq)
    samr.obj = getSamrObj() 
    
    if(!is.null(samr.obj) && ((input$assayType == "array" && (ar == "Two class unpaired" || ar == "Two class paired" || ar == "One class" || ar == "Survival")) || (input$assayType == "seq" && (seq == "Two class unpaired" || seq == "Two class paired" || seq == "Survival")))){
      
      data = getData()   
      y = data$y
      dif = input$dif
      #quick fix not allowing blocks for sample size
      if(!substring(y[1],2,6) %in% c("Block", "block")){  
        samplesize.factors = input$sampleSizeFactors
        samplesize.factors = as.numeric(unlist(strsplit(samplesize.factors, ",")))
        samr.assess.samplesize.obj =  samr.assess.samplesize(samr.obj, data, dif, samplesize.factors)
        samr.assess.samplesize.obj
      }
    }
  })
  
  getMissRateTable = reactive({
    
    samr.obj = getSamrObj()
    if(!is.null(samr.obj)){
      
      delta = findDelta()
      delta.table = getDeltaTable()
      
      missrate.table = samr.missrate(samr.obj, delta, delta.table)
      missrate.table
    }
    
  })
  
  getSiggenesTable = reactive({
    
    samr.obj = getSamrObj()
    if(!is.null(samr.obj)){
      delta = findDelta()
      data = getData()
      delta.table = getDeltaTable()
      min.foldchange = input$min.foldchange
      
      siggenes.table = samr.compute.siggenes.table(samr.obj,delta, data, delta.table, min.foldchange = min.foldchange, compute.localfdr= input$localFDR)
      siggenes.table
    }
  })
  
  
  getAllgenesTable = reactive({
    
    samr.obj = getSamrObj()
    if(!is.null(samr.obj)){
      delta = findDelta()
      data = getData()
      delta.table = getDeltaTable()
      min.foldchange = input$min.foldchange
      
      Allgenes = samr.compute.siggenes.table(samr.obj,delta, data, delta.table, min.foldchange = min.foldchange, compute.localfdr= input$localFDR, all.genes= TRUE)
      Allgenes
    }
  })
  
  getImputedX = reactive({
    
    data = getData()
    imputedX = data$imputedX
    if(!is.null(imputedX)){
      imputedX
    }
  })
  
  getOriginalX = reactive({
    
    data = getData()
    originalX = data$originalX
    if(!is.null(originalX)){
      originalX
    }
  })
  
  getComputedValues= reactive({
    
    samr.obj = getSamrObj()
    if(!is.null(samr.obj)){
      
      computedValues = matrix(NA, nrow = 5, ncol = 1)
      
      colnames(computedValues) = "Values"
      rownames(computedValues) = c("Estimated proportion of non-null features (genes)", "s0 percentile", "False Discovery Rate", "Estimated tail strength", "Estimated standard error of tail strength")
      tail.strength = getTailStrength()
      
      computedValues[1,1] = samr.obj$pi0  
      computedValues[2,1] = samr.obj$s0.perc
      delta.table = samr.compute.delta.table(samr.obj, min.foldchange= input$min.foldchange, dels= findDelta())[5]
      computedValues[3,1] = if(is.na(delta.table)){0.00}else{delta.table}
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
  
  getTailStrength = reactive({
    samr.obj = getSamrObj()
    if(!is.null(samr.obj)){
      samr.tail.strength(samr.obj)
    }    
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
  
  findDelta = reactive({
    
    if(input$deltaChoice == "Delta Slider")
      input$deltaSlider
    else if (input$deltaChoice == "Manually Enter Delta")
      input$deltaInput
    
  })
  
  output$rseed <- renderUI({
    
    if(input$randomButton == 0){
      value = 1234567
    }else{
      value = round(runif(1, 0, 10^6))
    }
    numericInput("random.seed", "Random Seed", value= value)
  })
  
  output$slider <- renderUI({
    
    delta.table = getDeltaTable()
    if(!is.null(delta.table)){
      sliderInput("deltaSlider", label = "Delta value", min = 0, max = as.numeric(delta.table[50,1]), value = 0, step = 0.00001)  
    }
    
  })
  

  
  
  capitalize = function(word){
    letters = strsplit(word,'')
    theletters = letters[[1]]
    theletters[1] = toupper(theletters[1])
    paste(theletters,collapse='')
  }
  
  ############Write texts for SAM interface
  
  output$missRateText = renderText({
    if(!is.null(getSamrObj())){
      paste("Estimated Miss rates for Delta = ", findDelta())
    }
    
  })
  
  output$eigengeneText = renderText({
    if(!is.null(getSamrObj()) && input$responseType_array == "Pattern discovery"){
      "Eigengene"
    }
  })
  
  
  output$eigengene = renderTable({
    eigengene = getEigengene()
    if(!is.null(eigengene)){
      eigengene
    }
  })
  
  output$negativeGenesText <- renderText({
    
    siggenes.table = getSiggenesTable()
    if(!is.null(siggenes.table)){
      paste("Negative Genes (", siggenes.table$ngenes.lo, ")", sep ="")
    }
  })
  
  output$positiveGenesText <- renderText({
    
    siggenes.table = getSiggenesTable()
    if(!is.null(siggenes.table)){
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
    if(!is.null(getSamrObj())){
      "Delta Table"
    }
  })
  
  output$samPlotText <- renderText({
    if(!is.null(getSamrObj())){
      "SAM Plot"
    }
  })
  
  output$inputParametersText <- renderText({
    if(!is.null(getSamrObj())){
      "Input Parameters"
    }
  })
  
  output$computedValuesText <- renderText({
    if(!is.null(getSamrObj())){
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
  

  
  ########Output tables and plots for SAM
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
    samr.obj = getSamrObj()
    delta = findDelta()
    min.foldchange = input$min.foldchange
    if(!is.null(samr.obj)){
      samr.plot(samr.obj, delta, min.foldchange = min.foldchange)
    }
  })
  
  output$deltaTable = renderTable({
    delta.table = getDeltaTable()
    if(!is.null(delta.table)){
      delta.table
    }
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
  
  output$computedValues = renderTable({
    if(!is.null(getComputedValues())){
      getComputedValues()      
    }
  })
  
  output$inputParameters = renderTable({
    inputParameters = getInputParameters()
    if(!is.null(inputParameters)){
      inputParameters
    }
  })
  
  
  ##########Saving results into xlsx format
  
  #######SAM results save
  savethis = observe({
  
    if(input$saveButton != 0){
      
    isolate({  
    
    dir = input$dir
    file = input$fname
    
    siggenes.table = getSiggenesTable()
    missrate.table = getMissRateTable()
    delta.table = getDeltaTable()
    samr.assess.samplesize.obj =  getSampleSize()
    Allgenes = getAllgenesTable()
    
    samr.obj = getSamrObj()
    delta = findDelta()
    min.foldchange = input$min.foldchange
    
    inputParameters = getInputParameters()
    computedValues = getComputedValues()
    eigengene = getEigengene()
    imputedX = getImputedX()
    originalX = getOriginalX()
    
    titleStyle = createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center", fgFill = "#4F81BD")
    
    wb = createWorkbook()  
      
    
    
    if(!is.null(samr.obj)){
      png(file = "SAMPlot.png")
      samr.plot(samr.obj, delta, min.foldchange = min.foldchange)
      dev.off()      
    }
    
    if(!is.null(samr.assess.samplesize.obj)){
      png(file = "samplesizePlot.png")
      samr.assess.samplesize.plot(samr.assess.samplesize.obj)  
      dev.off()  
    }
    
    if(!is.null(originalX)){
      addWorksheet(wb, sheetName = "Original Data")
      writeData(wb, originalX, sheet = "Original Data")
    }
    
    if(!is.null(imputedX)){
      addWorksheet(wb, sheetName = "Imputed Data")
      writeData(wb, imputedX, sheet = "Imputed Data")
    }
    
    if(!is.null(samr.obj)){      
      addWorksheet(wb, sheetName = "SAM Plot")
      insertImage(wb, sheet = "SAM Plot", file = "SAMPlot.png", width = 5, height = 5)
    }
    
    if(!is.null(delta.table)){ 
      addWorksheet(wb, sheetName = "Delta Table")
      writeData(wb, sheet = "Delta Table", x = "Delta Table")      
      writeData(wb, sheet = "Delta Table", x = delta.table, startRow = 2)  
      
      addStyle(wb, sheet = "Delta Table", titleStyle, rows = 1, cols = 1)
      setColWidths(wb, "Delta Table", cols= 1:8, widths = 18)      
    }
    
    if(!is.null(missrate.table)){
      missraterow = nrow(delta.table) + 4
      writeData(wb, sheet = "Delta Table", x = "Miss Rates", startRow = missraterow)
      writeData(wb, sheet = "Delta Table", x = missrate.table, startRow = missraterow + 1)
      addStyle(wb, sheet = "Delta Table", titleStyle, rows = missraterow, cols = 1)
    }
    
    if((siggenes.table$ngenes.up != 0) || (siggenes.table$ngenes.lo != 0)){
      addWorksheet(wb, sheetName = "Significant Genes")  
      setColWidths(wb, "Significant Genes", cols= 1:11, widths = 18)
    }
    
    significantgenerow = 1
    if(!is.null(siggenes.table$genes.up)){
      writeData(wb, sheet = "Significant Genes", x = "Significant Pos", startRow = significantgenerow)
      addStyle(wb, sheet = "Significant Genes", titleStyle, rows = significantgenerow, cols = 1)
      writeData(wb, sheet = "Significant Genes", x = siggenes.table$genes.up, startRow = significantgenerow + 1)
      significantgenerow = siggenes.table$ngenes.up + 4
      
    }    
    
    if(!is.null(siggenes.table$genes.lo)){
      writeData(wb, sheet = "Significant Genes", x = "Significant Neg", startRow = significantgenerow)
      writeData(wb, sheet = "Significant Genes", x = siggenes.table$genes.lo, startRow = significantgenerow + 1)
      addStyle(wb, sheet = "Significant Genes", titleStyle, rows = significantgenerow, cols = 1)   
    }    
    
    
    if(!is.null(Allgenes)){
      addWorksheet(wb, sheetName = "All Genes")  
      setColWidths(wb, "All Genes", cols= 1:11, widths = 18)  
    }  
    
    allgenerow = 1
    if(!is.null(Allgenes$genes.up)){
      writeData(wb, sheet = "All Genes", x = "All Pos", startRow = allgenerow)
      addStyle(wb, sheet = "All Genes", titleStyle, rows = allgenerow, cols = 1)
      writeData(wb, sheet = "All Genes", x = Allgenes$genes.up, startRow = allgenerow + 1)
      allgenerow = Allgenes$ngenes.up + 4
    }    
    
    if(!is.null(Allgenes$genes.lo)){
      writeData(wb, sheet = "All Genes", x = "All Neg", startRow = allgenerow)
      addStyle(wb, sheet = "All Genes", titleStyle, rows = allgenerow, cols = 1)
      writeData(wb, sheet = "All Genes", x = Allgenes$genes.lo, startRow = allgenerow + 1)
    }
    
    if(!is.null(samr.assess.samplesize.obj)){
      addWorksheet(wb, sheetName = "Sample Size") 
      setColWidths(wb, "Sample Size", cols= 1:8, widths = 20)
      insertImage(wb, sheet = "Sample Size", file = "samplesizePlot.png", width = 5, height = 5, startRow = 1)
      
      dataSample1 = paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[1] * samr.assess.samplesize.obj$n, sep = "")
      writeData(wb, sheet = "Sample Size", x = dataSample1, startRow = 26)
      addStyle(wb, sheet = "Sample Size", titleStyle, rows = 26, cols = 1)
      writeData(wb, sheet = "Sample Size", x = samr.assess.samplesize.obj$results[,,1], startRow = 27)
      
      dataSample2 = paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[2] * samr.assess.samplesize.obj$n, sep = "")
      writeData(wb, sheet = "Sample Size", x = dataSample2, startRow = 39)
      addStyle(wb, sheet = "Sample Size", titleStyle, rows = 39, cols = 1)
      writeData(wb, sheet = "Sample Size", x = samr.assess.samplesize.obj$results[,,2], startRow = 40)
      
      dataSample3 = paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[3] * samr.assess.samplesize.obj$n, sep = "")
      writeData(wb, sheet = "Sample Size", x = dataSample3, startRow = 52)
      addStyle(wb, sheet = "Sample Size", titleStyle, rows = 52, cols = 1)
      writeData(wb, sheet = "Sample Size", x = samr.assess.samplesize.obj$results[,,3], startRow = 53)
      
      dataSample4 = paste("Sample Size = ", samr.assess.samplesize.obj$samplesize.factor[4] * samr.assess.samplesize.obj$n, sep = "")
      writeData(wb, sheet = "Sample Size", x = dataSample4, startRow = 65)
      addStyle(wb, sheet = "Sample Size", titleStyle, rows = 65, cols = 1)
      writeData(wb, sheet = "Sample Size", x = samr.assess.samplesize.obj$results[,,4], startRow = 66)
    }
    
    if(!is.null(inputParameters)){
      addWorksheet(wb, sheetName = "Current Settings")
      setColWidths(wb, "Current Settings", cols= 1, widths = 46)
      setColWidths(wb, "Current Settings", cols= 2, widths = 18)
      newInputParameters = cbind(rownames(inputParameters), inputParameters)
      colnames(newInputParameters) = c("Questions","Values")
      writeData(wb, sheet = "Current Settings", x = "Input Parameters")
      addStyle(wb, sheet = "Current Settings", titleStyle, rows = 1, cols = 1)
      writeData(wb, sheet = "Current Settings", x = newInputParameters, startRow = 2)
      
      currentsettingsrow = nrow(newInputParameters) + 4
      writeData(wb, sheet = "Current Settings", x = "Computed Values", startRow = currentsettingsrow)
      addStyle(wb, sheet = "Current Settings", titleStyle, rows = currentsettingsrow, cols = 1)
      
      newComputedValues = cbind(rownames(computedValues), computedValues)
      colnames(newComputedValues) = c("Questions","Values")
      writeData(wb, sheet = "Current Settings", x = newComputedValues, startRow = currentsettingsrow + 1)
      
      
      if(!is.null(eigengene)){
        currentsettingsrow = currentsettingsrow + nrow(newComputedValues)
        writeData(wb, sheet = "Current Settings", x = "Eigengene", startRow = currentsettingsrow + 3)
        addStyle(wb, sheet = "Current Settings", titleStyle, rows = currentsettingsrow + 3, cols = 1)
        writeData(wb, sheet = "Current Settings", x = eigengene, startRow = currentsettingsrow + 4)
      }
    }
  
    if(!is.null(samr.obj)){
      fname = paste(file, "xlsx", sep = ".")
      saveWorkbook(wb, file.path(dir, fname), overwrite = TRUE)
      
      if(Sys.info()[['sysname']] == "Windows"){
        shell(file.path(dir, fname))  
      } else if(Sys.info()[['sysname']] == "Darwin"){
        system(paste("open", fname))
      }
          
    }
    
    
    })
    }
      
  })

  ####### Save .xlsx for GSA methods
  savethis2 = observe({
    
    if(input$saveButton2 != 0){
      
      isolate({  
        
        dir = input$dir2
        file = input$fname2
        
        GSA = getGSA()
        GSA.obj = GSA$GSA.obj
        GSA.list = getGSAList()
        GSAFullList = getGSAFullList()
        
        GSA.correlate.list = getGSACorrelateList()
        geneSetTable = GSA.correlate.list$result
        geneSetQuantile = GSA.correlate.list$QuantileCoverage
        geneSetTableGenes = GSA.correlate.list$tableGenes
               
        inputParameters = getInputParameters2()
        computedValues = getComputedValues2()
        
        originalX = getOriginalX()
        imputedX = getImputedX()
        
        titleStyle = createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center", fgFill = "#4F81BD")
        
        wb = createWorkbook()  
        
        if(!is.null(GSA.obj)){
          png(file = "GSAPlot.png")
          GSA.plot.revised(GSA.obj, FDRcut = 1, fac = 0)
          dev.off()    
        }
        
        if(!is.null(originalX)){
          addWorksheet(wb, sheetName = "Original Data")
          writeData(wb, originalX, sheet = "Original Data")
        }
        
        if(!is.null(imputedX)){
          addWorksheet(wb, sheetName = "Imputed Data")
          writeData(wb, imputedX, sheet = "Imputed Data")
        }
        
        if(!is.null(GSA)){
          addWorksheet(wb, sheetName = "GSA Plot")
          insertImage(wb, sheet = "GSA Plot", file = "GSAPlot.png", width = 5, height = 5)
        }
        
        if(!is.null(GSA.list$positive[1]) || !is.null(GSA.list$negative[1])){
          addWorksheet(wb, sheetName = "Significant Gene Sets")
          setColWidths(wb, "Significant Gene Sets", cols= 1:5, widths = 18)
        }
        
        significantgenerow = 1
        if(!(GSA.list$nsets.pos == 0)){
          writeData(wb, sheet = "Significant Gene Sets", x = "Significant Pos", startRow = significantgenerow)
          addStyle(wb, sheet = "Significant Gene Sets", titleStyle, rows = significantgenerow, cols = 1)
          writeData(wb, sheet = "Significant Gene Sets", x = GSA.list$positive, startRow = significantgenerow + 1)
          significantgenerow = GSA.list$nsets.pos + 4          
        }

        if(!(GSA.list$nsets.neg == 0)){
          writeData(wb, sheet = "Significant Gene Sets", x = "Significant Neg", startRow = significantgenerow)
          writeData(wb, sheet = "Significant Gene Sets", x = GSA.list$negative, startRow = significantgenerow + 1)
          addStyle(wb, sheet = "Significant Gene Sets", titleStyle, rows = significantgenerow, cols = 1)   
        }    
        
        
        if(!is.null(GSAFullList$positive) || !is.null(GSAFullList$negative)){
          addWorksheet(wb, sheetName = "Full Gene Sets")
          setColWidths(wb, "Full Gene Sets", cols= 1:5, widths = 18)
        }
        
        allgenerow = 1
        if(!is.null(GSAFullList$positive)){
          writeData(wb, sheet = "Full Gene Sets", x = "All Positive", startRow = allgenerow)
          addStyle(wb, sheet = "Full Gene Sets", titleStyle, rows = allgenerow, cols = 1)
          writeData(wb, sheet = "Full Gene Sets", x = GSAFullList$positive, startRow = allgenerow + 1)
          allgenerow = GSAFullList$nsets.pos + 4
        }
        
        if(!is.null(GSAFullList$negative)){
          writeData(wb, sheet = "Full Gene Sets", x = "All Negative", startRow = allgenerow)
          addStyle(wb, sheet = "Full Gene Sets", titleStyle, rows = allgenerow, cols = 1)
          writeData(wb, sheet = "Full Gene Sets", x = GSAFullList$negative, startRow = allgenerow + 1)     
        }
        
        if(!is.null(geneSetTable)){
          addWorksheet(wb, sheetName = "Gene Set Collection")
          setColWidths(wb, "Gene Set Collection", cols= 1:2, widths = 40)
          
          genesetinforow = 1
          newGeneSetTable = cbind(rownames(geneSetTable), geneSetTable)
          colnames(newGeneSetTable) = c("Questions","Values")
          writeData(wb, sheet = "Gene Set Collection", x = "Information for gene set collection", startRow = genesetinforow)
          addStyle(wb, sheet = "Gene Set Collection", titleStyle, rows = genesetinforow, cols = 1)
          writeData(wb, sheet = "Gene Set Collection", x = newGeneSetTable, startRow = genesetinforow + 1)
          genesetinforow = nrow(newGeneSetTable) + 4
          
          newGeneSetQuantile = cbind(colnames(geneSetQuantile), geneSetQuantile[1,])
          colnames(newGeneSetQuantile) = c("Quantiles","Values")
          writeData(wb, sheet = "Gene Set Collection", x = "Quantiles of fraction coverage", startRow = genesetinforow)
          addStyle(wb, sheet = "Gene Set Collection", titleStyle, rows = genesetinforow, cols = 1)
          writeData(wb, sheet = "Gene Set Collection", x = newGeneSetQuantile, startRow = genesetinforow + 1)
          genesetinforow = genesetinforow + nrow(newGeneSetQuantile) + 3
          
          newGeneSetTableGenes = cbind(rownames(geneSetTableGenes), geneSetTableGenes)
          colnames(newGeneSetTableGenes) = c("Number of genes in a set","count")
          writeData(wb, sheet = "Gene Set Collection", x = "Number of genes", startRow = genesetinforow)
          addStyle(wb, sheet = "Gene Set Collection", titleStyle, rows = genesetinforow, cols = 1)
          writeData(wb, sheet = "Gene Set Collection", x = newGeneSetTableGenes, startRow = genesetinforow + 1)
        }

        if(!is.null(inputParameters)){
          addWorksheet(wb, sheetName = "Current Settings")
          setColWidths(wb, "Current Settings", cols= 1:2, widths = 40)
          
          newInputParameters = cbind(rownames(inputParameters), inputParameters)
          colnames(newInputParameters) = c("Questions","Values")
          
          currentsettingsrow = 1
          writeData(wb, sheet = "Current Settings", x = "Input Parameters", startRow = currentsettingsrow)
          addStyle(wb, sheet = "Current Settings", titleStyle, rows = currentsettingsrow, cols = 1)
          writeData(wb, sheet = "Current Settings", x = newInputParameters, startRow = currentsettingsrow + 1)
          currentsettingsrow = nrow(newInputParameters) + 4
          
          writeData(wb, sheet = "Current Settings", x = "Computed Values", startRow = currentsettingsrow)
          addStyle(wb, sheet = "Current Settings", titleStyle, rows = currentsettingsrow, cols = 1)
          newComputedValues = cbind(rownames(computedValues), computedValues)
          colnames(newComputedValues) = c("Questions","Values")
          writeData(wb, sheet = "Current Settings", x = newComputedValues, startRow = currentsettingsrow + 1)
        }

        if(!is.null(GSA)){
          fname = paste(file, "xlsx", sep = ".")
          saveWorkbook(wb, file.path(dir, fname), overwrite = TRUE)
          
          if(Sys.info()[['sysname']] == "Windows"){
            shell(file.path(dir, fname))  
          } else if(Sys.info()[['sysname']] == "Darwin"){
            system(paste("open", fname))
          }
          
        }
        
      })
    }
    
  })
  
  

})


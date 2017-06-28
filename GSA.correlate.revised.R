GSA.correlate.revised =function(GSA.genesets.obj, genenames){
  nsets=length(GSA.genesets.obj$genesets)
  ngenes=unlist(lapply(GSA.genesets.obj$genesets,length))
  allgenes=unlist(GSA.genesets.obj$genesets)
  sets.in.exp=match(unique(allgenes),genenames)
  exp.in.sets=match(genenames,allgenes)
  
  result = matrix(nrow = 6, ncol = 1)
  colnames(result) = "value"
  rownames(result) = c("Number of gene-sets", "Total number of genes in gene-set collection", "Total number of unique genes in gene-set collection", "Total number of genes in  genenames list", "Total number of unique genes in genenames list", "Number of unique genes in both collections")
  result[1,1] = nsets
  result[2,1] = sum(ngenes)
  result[3,1] = length(unique(allgenes))
  result[4,1] = length(genenames)
  result[5,1] = length(unique(genenames))
  result[6,1] = sum(!is.na(sets.in.exp))
  
  nn=rep(NA,nsets)
  for(i in 1:nsets){
    nn[i]=sum(!is.na(match(GSA.genesets.obj$genesets[[i]],genenames)))
  }
  
  QuantileCoverage =  quantile(nn/ngenes, seq(0,1,by=.1))
  res = matrix(as.numeric(NA), 1, 11)
  res[1,] = QuantileCoverage
  dimnames(res) = list("value", names(QuantileCoverage))
  QuantileCoverage = res
  
  tableGenes = table(ngenes, dnn = "Table of number of genes in genesets")
  tableGenesN = length(names(tableGenes))
  res2 = matrix(as.numeric(NA), tableGenesN, 1)
  dimnames(res2) = list(names(tableGenes), "count")
  for(i in 1:tableGenesN){
    res2[i,1] = tableGenes[[i]]
  }
  tableGenes = res2
  #  cat(c("Quantiles of fraction coverage of gene-sets"),fill=T)
  #  print(quantile(nn/ngenes, seq(0,1,by=.1)),digits=4, fill = T)
  #  cat("", fill = T)
  #  table(ngenes, dnn = "Table of number of genes in genesets")
  list(result = result, QuantileCoverage = QuantileCoverage, tableGenes = tableGenes)
}

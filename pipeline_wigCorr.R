#present output from wigCorr
#transform into matrix
#plot heatmap
#output files for replicates that don't correlate best with their counterpart

source('scripts/heatmap.3-split.R')

wigCorr.asMatrix = function(fileName){
  wigCorr=as.matrix(read.table(file=fileName))
  wigCorr = sub(x = wigCorr, pattern = "__treat_pileup", replacement = "")
  colA = wigCorr[,1]
  colB = wigCorr[,2]
  
  simpleNames = function(names){
    n = strsplit(names, split = "[/._]")
    n = matrix(unlist(n), nrow = length(n[[1]]),byrow = F)
    return(paste(n[1,], n[2,], n[3,]))
  }
  
  colA = simpleNames(colA)
  
  colB = simpleNames(colB)
  
  names = unique(c(colA,colB))
  o = order(names)
  names = names[o]
  len = length(names)
  
  final = matrix(1, nrow = len, ncol = len)
  rownames(final) = names
  colnames(final) = names
  
  for(i in 1:length(colA)){
    a = colA[i]
    b = colB[i]
    val = as.numeric(wigCorr[i,3])
    final[a,b] = val
    final[b,a] = val
  }
  return(final)
}

plotWigCorr = function(fileName){
  corrMatrix = wigCorr.asMatrix(fileName)
  heatmap.2(corrMatrix, Colv = NA, Rowv = NA, scale = 'none', revC = T, trace = "none", margins = c(10,12), cexCol = .7, cexRow = .6)
  layout(1)
  text(x = .5,y=.9, fileName)
  return(corrMatrix)
}

tp = 'wigCorr.txt'
pdf('wigCorr plots.pdf')
  output = plotWigCorr(tp)
  
  outFile = sub(tp, pattern = '.txt', replacement = '.csv')
  write.table(format(output, digits = 3), outFile, sep = ',' )
dev.off()

wc = wigCorr.asMatrix(tp)
toPlot = ifelse(wc == wc, 0, 0)
for(i in 1:ncol(wc)){
  o = order(wc[,i], decreasing = T)
  n = colnames(wc)[i]
  n = strsplit(n, ' ')[[1]]
  n = paste(n[1], n[2])
  #print(colnames(wc)[o])
  cnames = colnames(wc)[o]
  x = grepl(pattern = n,x = cnames)
  toPlot[x,i] = 1
}

heatmap.2(toPlot, Colv = F, Rowv = F, margins = c(8,8), trace = 'n')
wc_order = function(i){
  return(wc[,i,drop = F][order(wc[,i], decreasing = T),,drop = F])
}
bad = c(4,8, 15, 16)
for(i in bad){
  tmp = wc_order(i)
  tmp = cbind(rownames(tmp), format(round(tmp,3), width = 3))
  write.csv(tmp, paste('wigCorr.', colnames(wc)[i],'.csv',sep =""), quote = F, row.names = F)
}

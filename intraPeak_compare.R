#compare intra rep called peaks to peaks called vs input.
scoreBy = 8 #column to score by, 8 is pval 9 is fdr?

realPeaks = list()
for(f in dir('replicate_peaks/', full.names = T)){
  tmp = read.table(f)
  realPeaks[[length(realPeaks) + 1]] = tmp[,scoreBy]
}
n = dir('replicate_peaks/', full.names = T)
n = strsplit(n, '[/_]')
n = matrix(unlist(n), nrow = length(n[[1]]))
n = paste(n[3,], n[4,], n[5,])
names(realPeaks) = n


falsePeaks = list()
for(f in dir('intra_narrowPeaks/', full.names = T)){
  tmp = read.table(f)
  falsePeaks[[length(falsePeaks) + 1]] = tmp[,scoreBy]
}
n = dir('intra_narrowPeaks/', full.names = T)
n = strsplit(n, '[/_]')
n = matrix(unlist(n), nrow = length(n[[1]]))
n = paste(n[3,], n[4,], n[5,])
names(falsePeaks) = n

plot(c(0,100), c(0,10000), type = 'n')
toCheck = c(0,5,20)
dataR = matrix(0, nrow = length(realPeaks), ncol = length(toCheck))
rownames(dataR) = names(realPeaks)
for(rep in names(realPeaks)){
  rP = realPeaks[[rep]]
  
  yrP = sapply(toCheck, function(x)return(sum(rP > x)))
  dataR[rep,] = yrP
  
  #  lines(x, yfP, type = 'l')
#   print(rep)
#   print(nR / nF)
#   dR = density(rP)
#   dF = density(fP)
#   plot(dR)
#   lines(dF, col = 'red')
#   boxplot(list(rP, fP))
}

dataF = matrix(0, nrow = length(falsePeaks), ncol = length(toCheck))
rownames(dataF) = names(falsePeaks)
for(rep in names(falsePeaks)){
  fP = falsePeaks[[rep]]
  yfP = sapply(toCheck, function(x)return(sum(fP > x)))
  dataF[rep,] = yfP
  #lines(x, yfP, type = 'l')
  #   print(rep)
  #   print(nR / nF)
  #   dR = density(rP)
  #   dF = density(fP)
  #   plot(dR)
  #   lines(dF, col = 'red')
  #   boxplot(list(rP, fP))
}

write.csv(dataR, 'realpeaks.csv')
write.csv(dataF, 'falsepeaks.csv')
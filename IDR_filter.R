if(!exists('wd')){
  wd = ""
}
indexDict = read.table("ref/hg38.index2symbol.1kb_ext_promoters.txt", row.names = 1, stringsAsFactors = F)
lines = c("MCF10A", "MCF7", "MDA231")
mods = c("H3K4AC", "H3K4ME3")
prefixes = list()
for (l in lines) {
  for (m in mods) {
    fname = paste(wd, "/IDR_results/", l, "_", m, sep = "")
    if(wd == ""){
      fname = paste("IDR_results/", l, "_", m, sep = "")
    }
    prefixes[[length(prefixes) + 1]] = fname
  }
}
prefixes = unlist(prefixes)

#extract chr position
keep = rownames(indexDict) != "55637"
indexDict = indexDict[keep,]
indexes = rownames(indexDict)
str = indexDict[,4]
str = matrix(unlist(strsplit(str, split = ':')), nrow = 2)
chroms = str[1,]

str = str[2,]
str = matrix(unlist(strsplit(str, '-')), nrow = 2)
starts = as.numeric(str[1,])
ends = as.numeric(str[2,])
chroms = as.factor(chroms)

ranges = cbind(starts, ends)
rownames(ranges) = indexes
names(chroms) = indexes

isPassed = matrix(F, nrow = length(indexes), ncol = length(prefixes))
rownames(isPassed) = indexes
tmp = strsplit(prefixes, '/')
n = matrix(unlist(tmp), nrow = length(tmp[[1]]))[length(tmp[[1]]),]
colnames(isPassed) = n

chrRanges = list()
for(chr in unique(chroms)){
  keep = chroms == chr
  chrRanges[[length(chrRanges) + 1]] = ranges[keep,]
}
names(chrRanges) = unique(chroms)

checkOverlaps = function(chr, start, end){
  ranges = chrRanges[[chr]]
  keep = !(start > ranges[,2] | end < ranges[,1])
  if(sum(keep) > 0){
    return(rownames(ranges)[keep])
  }
}


targetIDR = .01
print(paste('filtering peaks for IDR cutoff', targetIDR, '...'))
for(pref in prefixes){
  tmp = read.table(paste(pref, '.aboveIDR.txt', sep = ''))
  thresh = tmp[,3]
  num = tmp[,4]
  
  keep = thresh == targetIDR
  targetNum = as.numeric(num[keep])
  
  print(paste(pref,targetNum))
}
  
  idrpref = sub('IDR_results/', 'pooled_narrowPeaks/', pref)
  narrowPeaks = read.table(paste(idrpref, '.narrowPeak', sep = ''), stringsAsFactors = F)
  pvals = narrowPeaks[,8]
  o = order(pvals, decreasing = T)
  narrowPeaks = narrowPeaks[o,]
  narrowPeaks = narrowPeaks[1:targetNum,]
  
  o = order(narrowPeaks[,2])#filter by position then by chrom
  narrowPeaks = narrowPeaks[o,]
  o = order(narrowPeaks[,1])
  narrowPeaks = narrowPeaks[o,]
  nam = basename(pref)
  for(i in 1:nrow(narrowPeaks)){
    hitIndex = checkOverlaps(narrowPeaks[i,1], narrowPeaks[i,2], narrowPeaks[i,3])
    if(is.null(hitIndex))
      next
    for(ind in hitIndex){
      isPassed[ind, nam] = T
    }
  }
}

idrPassed_byIndex = isPassed[,c(1,3,5,2,4,6)]
save(idrPassed_byIndex, file = 'idrPassed_byIndex.save')

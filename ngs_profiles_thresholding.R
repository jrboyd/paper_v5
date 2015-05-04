lines = c('MCF10A', 'MCF7', 'MDA231')
mods = c('H3K4AC', 'H3K4ME3')
other_mod = function(m){
  if(m == mods[1])
    return(mods[2])
  else if(m == mods[2])
    return(mods[1])
  return('unrecognized mod')
}
load('mycounts_data.save')
load('dictionaries.save')

primeColors = list(
  c(27,158,119),
  c(217,95,2),
  c(117,112,176),
  c(102,166,30),
  c(230,171,2),
  c(231,41,138)
)
names(primeColors) = colnames(markData_4me3_4ac)

secondColors = list(
  c(102,194,165),
  c(252,141,98),
  c(141,160,203),
  c(166,216,84),
  c(255,217,47),
  c(231,138,195)
)
names(secondColors) = colnames(markData_4me3_4ac)
secondColors = primeColors

applyWindow = function(dat, win = 10){
  if(win < 2)
    return(dat)
  out = matrix(0, nrow = nrow(dat), ncol = ncol(dat)/win)
  for(i in 1:ncol(out)){
    start = (i-1)*win+1
    end = i*win
    #out[,i] = apply(dat[,start:end],1,median)
    out[,i] = rowMeans(dat[,start:end])
  }
  return(out)
}

layout(matrix(1:6, ncol = 2, byrow = T))
par(mai = rep(0.3, 4))
for(l in lines){
  ac_dir = dir('ngsplot_data', full.names = T, pattern = paste(l,'.+','AC',sep = ""))
  me_dir = dir('ngsplot_data', full.names = T, pattern = paste(l,'.+','ME3',sep = ""))
  ac_prof = paste(ac_dir, '/hm1.txt', sep = "")
  me_prof = paste(me_dir, '/hm1.txt', sep = "")
  
  tmp = read.table(ac_prof, stringsAsFactors = F)
  ensgs = tmp[2:nrow(tmp),1]
  
  strand = tmp[2:nrow(tmp),4]
  ac_dat = as.matrix(tmp[2:nrow(tmp), 5:ncol(tmp)])
  rownames(ac_dat) = ensgs
  
  tmp = read.table(me_prof, stringsAsFactors = F)
  me_dat = as.matrix(tmp[2:nrow(tmp), 5:ncol(tmp)])
  rownames(me_dat) = ensgs
  
  sortBy = markData_4me3_4ac[,paste(l, 'H3K4AC')]
  o = order(sortBy, decreasing = T)
  top1k = rownames(markData_4me3_4ac)[o][1:1000]
  top1k = ensg2cut[top1k]
  top1k = intersect(top1k, rownames(ac_dat))
  plot(colMeans(ac_dat[top1k,]), type = 'l', col = 'blue', xlim = c(0,100), ylim = c(0,6))
  lines(colMeans(me_dat[top1k,]), col = 'green')
  
  sortBy = markData_4me3_4ac[,paste(l, 'H3K4ME3')]
  o = order(sortBy, decreasing = T)
  top1k = rownames(markData_4me3_4ac)[o][1:1000]
  top1k = ensg2cut[top1k]
  top1k = intersect(top1k, rownames(ac_dat))
  plot(colMeans(ac_dat[top1k,]), type = 'l', col = 'blue', xlim = c(0,100), ylim = c(0,6))
  lines(colMeans(me_dat[top1k,]), col = 'green')
  
}
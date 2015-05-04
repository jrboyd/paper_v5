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
pdf('ngs_custom.pdf')
for(l in lines){
  ac_dir = dir('ngsplot_data', full.names = T, pattern = paste(l,'.+','AC',sep = ""))
  me_dir = dir('ngsplot_data', full.names = T, pattern = paste(l,'.+','ME3',sep = ""))
  ac_prof = paste(ac_dir, '/hm1.txt', sep = "")
  me_prof = paste(me_dir, '/hm1.txt', sep = "")
  
  tmp = read.table(ac_prof, stringsAsFactors = F)
  ensgs = tmp[2:nrow(tmp),1]
  strand = tmp[2:nrow(tmp),4]
  ac_dat = as.matrix(tmp[2:nrow(tmp), 5:ncol(tmp)])
  
  
  tmp = read.table(me_prof, stringsAsFactors = F)
  me_dat = as.matrix(tmp[2:nrow(tmp), 5:ncol(tmp)])
  
  
  dat = cbind(ac_dat, me_dat)
  rownames(dat) = ensgs
  keep = intersect(ensgs, ensg2cut[rownames(markData_4me3_4ac)])
  dat = dat[keep,]
  ac_val = rowSums(dat[,1:(ncol(dat)/2)])
  me_val = rowSums(dat[,(ncol(dat)/2+1):ncol(dat)])
  o = order(ac_val,decreasing = T)
  cr = colorRamp(colors = c('blue', 'black', 'yellow'))
  m = max(abs(dat))
  toPlot = dat[o,]
  cap = 5
  toPlot = ifelse(toPlot > cap, cap, toPlot)
  toPlot = ifelse(toPlot < -cap, -cap, toPlot)
  heatmap.2(toPlot, Colv = NA, Rowv = NA,dendrogram = 'none', labRow = "", labCol = "", scale = 'n', col = rgb(cr(0:100/100)/255), symm = T, trace = 'n', density.info = 'none', key.xlab = 'log2 FE', main = paste(l,'\nH3K4ac                H3K4me3'))
}
dev.off()


#do some spiffy line plots
source('functions_movingAverage.R')
mat_lay = matrix(1:6, nrow = 3, ncol = 2, byrow = T)
mat_lay = cbind((1:3), mat_lay+3)#line labels
mat_lay = rbind((c(-1,1,2)), mat_lay+2)#mark labels
mat_lay = rbind(c(0,1,1), mat_lay+1)#mark labels

nf = layout(mat_lay, widths = c(.2,1,1), heights = c(.5,.2,1,1,1))
par(mai = c(.4,.5,.1,.1))
textPlot = function(txt){
  plot(c(0,1), c(0,1), type = 'n', axes = F)
  text(.5,.5, txt)
}
forePlot = function(x, lm){
  c = (primeColors[[lm]])/255
  c = rgb(c[1], c[2], c[3])
  lines(x, col = c, lwd = 3)
}
toptops = list()
tops = list()
bgPlot = function(x, lm){
  x = movingAverage(x, n = 3, centered = T)
  c = (secondColors[[lm]])/255
  c = rgb(c[1], c[2], c[3])
  plot(x, ylim = c(0,6), type = 'l', lwd = 1, col = c, xlab = '', ylab = '')
  perc_95 = length(x)/25*5
  perc_97 = length(x)/25*3
  perc_99 = length(x)/25*1
  perc_98 = length(x)/25*2
  perc_kept = perc_98
  lines(c(perc_kept, perc_kept),c(0,6), lty = 3)
  toptop = names(x)[1:perc_kept]
  toptops[[length(toptops) + 1]] <<- ensg2sym[toptop]
  top = names(x)
  tops[[length(tops) + 1]] <<- ensg2sym[top]
  #print(ensg2sym[top])
}
par(mai = rep(0,4))
textPlot(paste('Highly methylated promoters are acetylated.\nHighly acetylated promoters are not necessarily methylated.'))
for(m in mods){
  textPlot(paste(m, 'over 75th percentile sorted by', m,'\n',other_mod(m),'in background'))
}
for(l in lines){
  textPlot(l)
}
par(mai = c(.4,.5,.1,.1))
for(l in lines){
  ac_dat = markData_4me3_4ac[, paste(l,'H3K4AC')]
  me_dat = markData_4me3_4ac[, paste(l,'H3K4ME3')]
  
  q_thresh = .75
  ac_thresh = quantile(ac_dat, q_thresh)
  
  keep = ac_dat > ac_thresh
  ac_dat = ac_dat[keep]
  me_dat = me_dat[keep]
  
  o = order(ac_dat, decreasing = T)
  bgPlot(me_dat[o], paste(l, mods[2]))
  forePlot(ac_dat[o], paste(l, mods[1]))
  
  
  ac_dat = markData_4me3_4ac[, paste(l,'H3K4AC')]
  me_dat = markData_4me3_4ac[, paste(l,'H3K4ME3')]
  
  me_thresh = quantile(me_dat, q_thresh)
  keep = me_dat > me_thresh
  ac_dat = ac_dat[keep]
  me_dat = me_dat[keep]
  
  o = order(me_dat, decreasing = T)
  bgPlot(ac_dat[o], paste(l, mods[1]))
  forePlot(me_dat[o], paste(l, mods[2]))
}

thresh98 = apply(markData_4me3_4ac, 2, function(x)quantile(x, .98))
thresh75 = apply(markData_4me3_4ac, 2, function(x)quantile(x, .75))
pass98 = ifelse(markData_4me3_4ac, F, F)
for(i in 1:ncol(pass98)){
  pass98[,i] = markData_4me3_4ac[,i] >= thresh98[i]
}
pass75 = ifelse(markData_4me3_4ac, F, F)
for(i in 1:ncol(pass98)){
  pass75[,i] = markData_4me3_4ac[,i] >= thresh75[i]
}

library(limma)
library(venneuler)
pdf('percentile venns.pdf')
layout(1:2)
vennDiagram(pass98[,1:3], names = lines, )
title('98th percentile, H3K4ac')
vennDiagram(pass75[,1:3], names = lines)
title('75th percentile, H3K4ac')

perc_vennDiagram = function(x){
  tmp = vennCounts(x)
  tmp[1,4] = 0
  tmp[,4] = round(tmp[,4] / sum(tmp[,4]), digits = 2)
  tmp[1,4] = ""
  vennDiagram(tmp, names = lines)
}

perc_vennDiagram(pass98[,1:3])
title('98th percentile, H3K4ac')
perc_vennDiagram(pass75[,1:3])
title('75th percentile, H3K4ac')


vennDiagram(pass98[,4:6], names = lines)
title('98th percentile, H3K4me3')
vennDiagram(pass75[,4:6], names = lines)
title('75th percentile, H3K4me3')

perc_vennDiagram(pass98[,4:6])
title('98th percentile, H3K4me3')
perc_vennDiagram(pass75[,4:6])
title('75th percentile, H3K4me3')

plot(venneuler(pass98[,4:6]))
plot(venneuler(pass75[,4:6]))
dev.off()

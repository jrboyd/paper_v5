load('mycounts_data_bu.save')
ac_thresh = sapply(1:3, function(x)return(quantile(markData_4me3_4ac[,x], .75)))
me_thresh = sapply(4:6, function(x)return(quantile(markData_4me3_4ac[,x], .75)))

ac_keep = apply(markData_4me3_4ac[,1:3],1,function(x)return(all(x>ac_thresh)))
me_keep = apply(markData_4me3_4ac[,4:6],1,function(x)return(all(x>me_thresh)))

source('scripts/heatmap.3-split.R')

pdf('over 75th heatmaps.pdf')
me_dist = function(x){
  return(dist(x[,4:6]))
}
#filter by ac, cluster on me
heatmap.3(markData_4me3_4ac[ac_keep,], distfun = me_dist, trace = 'n', main = paste(sum(ac_keep), ' with ac over 75th percentile in all 3 lines\nclustered by me3'))


ac_dist = function(x){
  return(dist(x[,1:3]))
}
#filter by me, cluster on ac
res = heatmap.3(markData_4me3_4ac[me_keep,], distfun = ac_dist, trace = 'n', main = paste(sum(me_keep), ' with me3 over 75th percentile in all 3 lines\nclustered by ac'))


both_keep = ac_keep & me_keep
heatmap.3(markData_4me3_4ac[both_keep,], trace = 'n', main = paste(sum(both_keep), ' with both over 75th percentile in all 3 lines'))
dev.off()

for(i in 1:3){
  keep = markData_4me3_4ac[,i] > ac_thresh[i]
  o = order(markData_4me3_4ac[,i+3])
  #heatmap.2(markData_4me3_4ac[o,c(i,i+3)], Rowv = F, Colv = F, trace = 'n')
  plot(markData_4me3_4ac[o,i+3])
  lines(markData_4me3_4ac[o,i], col = rgb(0,0,0,.15))
  title('sorted by', colnames(markData_4me3_4ac)[i+3], 'with ')
}

for(i in 1:3){
  keep = markData_4me3_4ac[,i+3] > me_thresh[i]
  o = order(markData_4me3_4ac[,i])
  #heatmap.2(markData_4me3_4ac[o,c(i,i+3)], Rowv = F, Colv = F, trace = 'n')
  plot(markData_4me3_4ac[o,i])
  lines(markData_4me3_4ac[o,i+3], col = rgb(0,0,0,.15))
  title(colnames(markData_4me3_4ac)[i])
}
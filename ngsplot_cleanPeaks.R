downers = c(1:15,49:50,87:100)
uppers1 = 40:45
uppers2 = 55:60
colors = rep('black',ncol(d))
colors[downers] = 'red'
colors[uppers1] = 'green'
colors[uppers2] = 'green'
boxplot(d, col = colors)

highPeaks = d[,uppers2]
lowPeaks = d[,uppers1]
floor = d[,downers]
highPeakTest = (rowMeans(highPeaks) - (rowMeans(lowPeaks))) > .1
floorTest = (rowMeans(lowPeaks) - apply(d,1,function(x)quantile(x,.8))) > 0

ceilingTest = apply(highPeaks,1,min) > 1
highPeakTest = apply(highPeaks,1,max) > apply(lowPeaks,1,max)
floorTest = apply(lowPeaks,1,max) > apply(floor,1,max)
keep = highPeakTest & floorTest & ceilingTest
tmp = heatmap(d[keep,], Colv = NA, , labRow = "", labCol = "", scale = 'n')#, col = c('black', 'green'))


doMostChanged = function(base, test, markName, threshold = 4) {
    plotUnderOver = F
    if(doText)
    for (i in 1:3) {
        # write lists of changing genes
        p1 = base[i]
        p2 = test[i]
        dat = markData_4me3_4ac[, p2] - markData_4me3_4ac[, p1]
        names(dat) = ensg2sym[rownames(markData_4me3_4ac)]
        
        keep = dat > log2(threshold)
        print(paste(sum(keep), markName, "going UP  - ", colnames(markData_4me3_4ac)[p1], " to ", colnames(markData_4me3_4ac)[p2]))
        write.table(names(dat)[keep], file = paste("up - ", colnames(markData_4me3_4ac)[p1], " to ", colnames(markData_4me3_4ac)[p2], 
            ".txt", sep = ""), quote = F, row.names = F, col.names = F)
        
        keep = dat < -log2(threshold)
        print(paste(sum(keep), markName, "going DOWN  - ", colnames(markData_4me3_4ac)[p1], " to ", colnames(markData_4me3_4ac)[p2]))
        write.table(names(dat)[keep], file = paste("down - ", colnames(markData_4me3_4ac)[p1], " to ", colnames(markData_4me3_4ac)[p2], 
            ".txt", sep = ""), quote = F, row.names = F, col.names = F)
        
    }
    if(doPlots) pdf(paste("most_changed", markName, "heatmaps.pdf", sep = "_"))
    for (i in 1:3) {
        # plot heatmap of up/down
        p1 = base[i]
        p2 = test[i]
        dat = markData_4me3_4ac[, p2] - markData_4me3_4ac[, p1]
        keep = dat > log2(threshold) | dat < -log2(threshold)
        plotTitle = paste(
          sum(keep),'promoters\n',
          'meeting', threshold, 'fold change\nbetween', colnames(markData_4me3_4ac)[p1], 'and', colnames(markData_4me3_4ac)[p2])
        # print(sum(keep))
        options(warn = -1)
        myclust = function(x){
          attr(x, 'Size') = nrow(x)
          hclust(dist(x/rowSums(x)))
        }
        res = heatmap.3(markData_4me3_4ac[keep, ], hclustfun = hclust, trace = "n", classCount = 4, main = plotTitle)
        options(warn = 0)
        res_dat = res[[3]]
        rownames(res_dat) = ensg2sym[rownames(res_dat)]
        res[[3]] = res_dat
        plot.HeatmapLists(res)
        
    }
    if(doPlots) dev.off()
    if(plotUnderOver){
      
      if(doPlots) pdf(paste("most_changed", markName, "scatterplots.pdf", sep = "_"))
      for (i in 1:3) {
        # plot over under threshold
        p1 = base[i]
        p2 = test[i]
        dat = markData_4me3_4ac[, p2] - markData_4me3_4ac[, p1]
        
        plot(markData_4me3_4ac[, p1], markData_4me3_4ac[, p2], xlim = c(0, 6), ylim = c(0, 6), pch = 16, col = rgb(0, 0, 0, 0.1), 
             cex = 0.7, xlab = colnames(markData_4me3_4ac)[p1], ylab = colnames(markData_4me3_4ac)[p2])
        lines(c(0, 6), c(0 + log2(threshold), 6 + log2(threshold)))
        lines(c(0, 6), c(0 - log2(threshold), 6 - log2(threshold)))
        
      }
      if(doPlots) dev.off()
    }
    
    
    print(paste("wrote most changed lists for", markName))
}




base = c(1, 1, 2)
test = c(2, 3, 3)
doMostChanged(base, test, "K4ac")

base = base + 3
test = test + 3
doMostChanged(base, test, "K4me3")

 

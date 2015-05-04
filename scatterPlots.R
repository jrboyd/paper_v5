d = markData_4me3_4ac
fname = "scatterplots.pdf"
pdf(fname)
for (i in 1:3) {
    xy = d[, c(i + 3, i)]
    
    color_1slope = "#377eb8"
    color_trendLine = "#e41a1c"
    
    plot(xy, pch = 16, col = rgb(0, 0, 0, 0.1), cex = 0.7, type = "n", xlim = c(0, max(markData_4me3_4ac)), ylim = c(0, max(markData_4me3_4ac)))
    
    
    
    points(xy, pch = 16, col = rgb(0, 0, 0, 0.1), cex = 0.7)
    
    lines(c(-1, 10), c(-1, 10), lty = 2, col = color_1slope, lwd = 5)
    
    fit <- glm(xy[, 2] ~ xy[, 1])
    co <- coef(fit)
    abline(fit, col = color_trendLine, lwd = 2, )
    
    
    
    
    r2 = cor(xy)[1, 2]
    r2 = round(r2, digits = 2)
    r2 = format(r2, nsmall = 2)
    legend("topleft", legend = c("Slope = 1", "Linear Regression"), fill = c(color_1slope, color_trendLine))
    legend("bottomright", legend = bquote(R^2 == .(rs), list(rs = r2)))
}
dev.off()
print(paste("wrote scatter plots to", fname)) 

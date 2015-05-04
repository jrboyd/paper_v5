d = markData_4me3_4ac

fname = "scatterplots_overUnder.pdf"
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
    
    b = fit$coefficients[1]
    m = fit$coefficients[2]
    
    lines(c(0, 6), c(0, 6) * m + b + 1)
    lines(c(0, 6), c(0, 6) * m + b - 1)
    
    total = nrow(xy)
    isOver = xy[, 2] > ((xy[, 1] * m) + b + 1)
    isUnder = xy[, 2] < ((xy[, 1] * m) + b - 1)
    nOver = sum(isOver)
    nUnder = sum(isUnder)
    
    pOver = nOver/total
    pUnder = nUnder/total
    
    txtOver = paste(nOver, "/", total, "=", round(pOver, 4))
    txtUnder = paste(nUnder, "/", total, "=", round(pUnder, 4))
    
    points(xy[isOver, ], col = "red")
    points(xy[isUnder, ], col = "green")
    
    
    r2 = cor(xy)[1, 2]
    r2 = round(r2, digits = 2)
    r2 = format(r2, nsmall = 2)
    # legend('topleft', legend = c('Slope = 1', 'Linear Regression'), fill = c(color_1slope, color_trendLine))
    legend("topleft", legend = c(txtOver, txtUnder), fill = c("red", "green"))
    legend("bottomright", legend = bquote(R^2 == .(rs), list(rs = r2)))
}
dev.off()
print(paste("wrote over under plots to", fname)) 

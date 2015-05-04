plot_IDR_results = function(prefixes, output.file.prefix = "IDR_results") {
  npair <- length(prefixes)  # number of curves to plot on the same figure
  df.txt <- 10
  ntemp <- as.numeric(npair)
  saved.file.prefix <- list()  # identifier of filenames that contain the em and URI results
  
  uri.list <- list()
  uri.list.match <- list()
  ez.list <- list()
  legend.txt <- c()
  em.output.list <- list()
  uri.output.list <- list()
  
  for (i in 1:npair) {
    saved.file.prefix[i] <- prefixes[i]
    
    load(paste(saved.file.prefix[i], ".uri.sav", sep = ""))
    load(paste(saved.file.prefix[i], ".em.sav", sep = ""))
    
    uri.output.list[[i]] <- uri.output
    em.output.list[[i]] <- em.output
    
    ez.list[[i]] <- get.ez.tt.all(em.output, uri.output.list[[i]]$data12.enrich$merge1, uri.output.list[[i]]$data12.enrich$merge2)  # reverse =T for error rate
    
    # URI for all peaks
    uri.list[[i]] <- uri.output$uri.n
    # URI for matched peaks
    uri.match <- get.uri.matched(em.output$data.pruned, df = df.txt)
    uri.list.match[[i]] <- uri.match$uri.n
    
    file.name <- unlist(strsplit(as.character(saved.file.prefix[i]), "/"))
    
    legend.txt[i] <- paste(file.name[length(file.name)])
    
  }
  
  plot.uri.file <- paste(output.file.prefix, "-plot.pdf", sep = "")
  
  library(RColorBrewer)
  
  colors = RColorBrewer::brewer.pal(n = npair, name = "Dark2")
  ############# plot and report output plot correspondence curve for each pair, plot number of selected peaks vs IDR plot all into 1 file
  pdf(plot.uri.file, width = 14)
  par(mfcol = c(2, 3), mar = c(5, 6, 4, 2) + 0.1)
  plot.uri.group(uri.list, NULL, file.name = NULL, title.txt = "all peaks", col.txt = colors)
  plot.uri.group(uri.list.match, NULL, file.name = NULL, title.txt = "matched peaks", col.txt = colors)
  plot.ez.group(ez.list, plot.dir = NULL, file.name = NULL, y.lim = c(0, 0.6), col.txt = colors)
  plot(0, 1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", xaxt = "n", yaxt = "n")  # legends
  legend("center", legend.txt, cex = 2.5, fill = colors, box.col = "white")
  
  dev.off()
} 

# change the scale of uri from based on t (percentage) to n (number of peaks or basepairs) this is for plotting multiple
# pairwise URI's on the same plot
scale.t2n <- function(uri) {
    
    ntotal <- uri$ntotal
    tv <- uri$tv * uri$ntotal
    uri.uri <- uri$uri * uri$ntotal
    jump.left <- uri$jump.left
    uri.spl <- uri$uri.spl
    uri.spl$x <- uri$uri.spl$x * uri$ntotal
    uri.spl$y <- uri$uri.spl$y * uri$ntotal
    
    t.binned <- uri$t.binned * uri$ntotal
    uri.slope <- uri$uri.slope
    uri.der <- uri$uri.der
    uri.der$x <- uri$uri.der$x * uri$ntotal
    uri.der$y <- uri$uri.der$y
    
    uri.n <- list(tv = tv, uri = uri.uri, t.binned = t.binned, uri.slope = uri.slope, uri.spl = uri.spl, uri.der = uri.der, ntotal = ntotal, 
        jump.left = jump.left)
    return(uri.n)
}

############### plot correspondence profile

# plot multiple comparison wrt one template uri.list contains the total number of peaks plot.missing=F: not plot the missing
# points on the right
plot.uri.group <- function(uri.n.list, plot.dir, file.name = NULL, xlab.txt = "num of significant peaks", ylab.txt = "num of peaks in common", 
    col.start = 0, col.txt = NULL, plot.missing = F, title.txt = NULL) {
    
    if (is.null(col.txt)) 
        col.txt <- c("black", "red", "purple", "green", "blue", "cyan", "magenta", "orange", "grey")
    
    n <- length(uri.n.list)
    
    ntotal <- c()
    for (i in 1:n) ntotal[i] <- uri.n.list[[i]]$ntotal
    
    jump.left <- c()
    jump.left.der <- c()
    
    ncommon <- c()
    for (i in 1:n) {
        # jump.left[i] <- which.max(uri.n.list[[i]]$uri[-1]-uri.n.list[[i]]$uri[-length(uri.n.list[[i]]$uri)]) if(jump.left[i] < 6)
        # jump.left[i] <- length(uri.n.list[[i]]$uri)
        
        ## reversed.index <- seq(length(uri.n.list[[i]]$tv[,1]), 1, by=-1) nequal <- sum(uri.n.list[[i]]$uri[reversed.index]==
        ## uri.n.list[[i]]$tv[reversed.index,1]) temp <- which(uri.n.list[[i]]$uri[reversed.index]==
        ## uri.n.list[[i]]$tv[reversed.index,1])[nequal] jump.left[i] <- length(uri.n.list[[i]]$tv[,1])-temp print(uri.n.list[[i]]$uri)
        ## print(uri.n.list[[i]]$tv[,1]) jump.left[i] <- uri.n.list[[i]]$jump.left
        
        # jump.left.der[i] <- sum(uri.n.list[[i]]$t.binned < uri.n.list[[i]]$uri.der$x[length(uri.n.list[[i]]$uri.der$x)])
        
        jump.left[i] <- uri.n.list[[i]]$jump.left
        jump.left.der[i] <- jump.left[i]
        ncommon[i] <- uri.n.list[[i]]$tv[jump.left[i], 1]
    }
    
    
    if (plot.missing) {
        max.peak <- max(ntotal)
    } else {
        max.peak <- max(ncommon) * 1.05
    }
    
    if (!is.null(file.name)) {
        postscript(paste(plot.dir, "uri.", file.name, sep = ""))
        par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
    }
    
    plot(uri.n.list[[1]]$tv[, 1], uri.n.list[[1]]$uri, type = "n", xlab = xlab.txt, ylab = ylab.txt, xlim = c(0, max.peak), ylim = c(0, 
        max.peak), cex.lab = 1)
    
    for (i in 1:n) {
        
        if (plot.missing) {
            points(uri.n.list[[i]]$tv[, 1], uri.n.list[[i]]$uri, col = col.txt[i + col.start], cex = 0.5)
        } else {
            points(uri.n.list[[i]]$tv[1:jump.left[i], 1], uri.n.list[[i]]$uri[1:jump.left[i]], col = col.txt[i + col.start], cex = 0.5)
        }
        lines(uri.n.list[[i]]$uri.spl, col = col.txt[i + col.start], lwd = 4)
    }
    abline(coef = c(0, 1), lty = 3)
    # legend(0, max.peak, legend=legend.txt, col=col.txt[(col.start+1):length(col.txt)], lty=1, lwd=3, cex=2)
    if (!is.null(title)) 
        title(title.txt)
    
    if (!is.null(file.name)) {
        dev.off()
    }
    
    if (!is.null(file.name)) {
        postscript(paste(plot.dir, "duri.", file.name, sep = ""))
        par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
    }
    plot(uri.n.list[[1]]$t.binned, uri.n.list[[1]]$uri.slope, type = "n", xlab = xlab.txt, ylab = "slope", xlim = c(0, max.peak), 
        ylim = c(0, 1.5), cex.lab = 1)
    
    for (i in 1:n) {
        # if(plot.missing){ points(uri.n.list[[i]]$t.binned, uri.n.list[[i]]$uri.slope, col=col.txt[i+col.start], cex=0.5) } else {
        # points(uri.n.list[[i]]$t.binned[1:jump.left.der[i]], uri.n.list[[i]]$uri.slope[1:jump.left.der[i]],
        # col=col.txt[i+col.start], cex=0.5) }
        lines(uri.n.list[[i]]$uri.der, col = col.txt[i + col.start], lwd = 4)
    }
    abline(h = 1, lty = 3)
    # legend(0.5*max.peak, 1.5, legend=legend.txt, col=col.txt[(col.start+1):length(col.txt)], lty=1, lwd=3, cex=2)
    
    if (!is.null(title)) 
        title(title.txt)
    
    if (!is.null(file.name)) {
        dev.off()
    }
    
}

# compute uri for matched sample
get.uri.matched <- function(data12, df = 10) {
    
    tt <- seq(0.01, 1, by = 0.01)
    vv <- tt
    uri <- get.uri.2d(data12$sample1$sig.value, data12$sample2$sig.value, tt, vv, spline.df = df)
    
    # change scale from t to n
    uri.n <- scale.t2n(uri)
    
    return(list(uri = uri, uri.n = uri.n))
    
}

get.uri.2d <- function(x1, x2, tt, vv, spline.df = NULL) {
    
    o <- order(x1, x2, decreasing = T)
    
    # sort x2 by the order of x1
    x2.ordered <- x2[o]
    
    tv <- cbind(tt, vv)
    ntotal <- length(x1)  # number of peaks
    
    uri <- apply(tv, 1, comp.uri, x = x2.ordered)
    
    # compute the derivative of URI vs t using small bins
    uri.binned <- uri[seq(1, length(uri), by = 4)]
    tt.binned <- tt[seq(1, length(uri), by = 4)]
    uri.slope <- (uri.binned[2:(length(uri.binned))] - uri.binned[1:(length(uri.binned) - 1)])/(tt.binned[2:(length(uri.binned))] - 
        tt.binned[1:(length(tt.binned) - 1)])
    
    # smooth uri using spline first find where the jump is and don't fit the jump this is the index on the left jump.left.old <-
    # which.max(uri[-1]-uri[-length(uri)])
    short.list.length <- min(sum(x1 > 0)/length(x1), sum(x2 > 0)/length(x2))
    
    if (short.list.length < max(tt)) {
        jump.left <- which(tt > short.list.length)[1] - 1
    } else {
        jump.left <- which.max(tt)
    }
    
    # reversed.index <- seq(length(tt), 1, by=-1) nequal <- sum(uri[reversed.index]== tt[reversed.index]) temp <-
    # which(uri[reversed.index]== tt[reversed.index])[nequal] jump.left <- length(tt)-temp
    
    if (jump.left < 6) {
        jump.left <- length(tt)
    }
    
    
    if (is.null(spline.df)) 
        uri.spl <- smooth.spline(tt[1:jump.left], uri[1:jump.left], df = 6.4) else {
        uri.spl <- smooth.spline(tt[1:jump.left], uri[1:jump.left], df = spline.df)
    }
    # predict the first derivative
    uri.der <- predict(uri.spl, tt[1:jump.left], deriv = 1)
    
    invisible(list(tv = tv, uri = uri, uri.slope = uri.slope, t.binned = tt.binned[2:length(uri.binned)], uri.spl = uri.spl, uri.der = uri.der, 
        jump.left = jump.left, ntotal = ntotal))
}
# compute upper rank intersection for one t tv: the upper percentile x is sorted by the order of paired variable
comp.uri <- function(tv, x) {
    n <- length(x)
    qt <- quantile(x, prob = 1 - tv[1])  # tv[1] is t
    # sum(x[1:ceiling(n*tv[2])] >= qt)/n/tv[2]- tv[1]*tv[2] #tv[2] is v
    sum(x[1:ceiling(n * tv[2])] >= qt)/n
    
}

plot.ez.group <- function(ez.list, plot.dir, file.name = NULL, y.lim = NULL, xlab.txt = "num of significant peaks", ylab.txt = "IDR", 
    col.txt = NULL, title.txt = NULL) {
    
    if (is.null(col.txt)) 
        col.txt <- c("black", "red", "purple", "green", "blue", "cyan", "magenta", "orange", "grey")
    
    n.entry <- length(ez.list)
    x <- rep(NA, n.entry)
    y.max <- rep(NA, n.entry)
    
    for (i in 1:n.entry) {
        x[i] <- max(ez.list[[i]]$n)
        
        y.max[i] <- max(ez.list[[i]]$IDR)
        
    }
    
    if (is.null(y.lim)) 
        y.lim <- c(0, max(y.max))
    
    if (!is.null(file.name)) {
        postscript(paste(plot.dir, "ez.", file.name, sep = ""))
        par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
    }
    
    
    
    plot(c(0, max(x)), y.lim, ylim = y.lim, type = "n", xlab = xlab.txt, ylab = ylab.txt, lwd = 5, cex = 5, cex.axis = 2, cex.lab = 2)
    
    q <- seq(0.01, 0.99, by = 0.01)
    
    for (i in 1:length(ez.list)) {
        
        n.plot <- round(quantile(ez.list[[i]]$n, prob = q))
        IDR.plot <- ez.list[[i]]$IDR[n.plot]
        lines(n.plot, IDR.plot, col = col.txt[i], cex = 2, lwd = 5)
    }
    
    
    # legend(0, y.lim[2], legend=legend.txt, col=col.txt[1:length(col.txt)], lty=1, lwd=5, cex=2)
    
    if (!is.null(title)) 
        title(title.txt)
    
    if (!is.null(file.name)) {
        dev.off()
    }
    
}

get.ez.tt.all <- function(em.fit, all.data1, all.data2, idr.level = c(0.01, 0.05, 0.1)) {
    
    u <- em.fit$data.pruned$sample1$sig.value
    v <- em.fit$data.pruned$sample2$sig.value
    # u <- em.fit$data.pruned$sample1 v <- em.fit$data.pruned$sample2
    
    e.z <- 1 - em.fit$em.fit$e.z  # this is the error prob
    
    o <- order(e.z)
    e.z.ordered <- e.z[o]
    n.select <- c(1:length(e.z))
    IDR <- cumsum(e.z.ordered)/n.select
    
    u.o <- u[o]
    v.o <- v[o]
    
    n.level <- length(idr.level)
    # sig.value1 <- rep(NA, n.level) sig.value2 <- rep(NA, n.level)
    ez.cutoff <- rep(NA, n.level)
    n.selected <- rep(NA, n.level)
    npeak.rep1 <- rep(NA, n.level)
    npeak.rep2 <- rep(NA, n.level)
    
    for (i in 1:length(idr.level)) {
        
        # find which uri.ez is closet to fdr.level
        index <- which.min(abs(IDR - idr.level[i]))
        # sig.value1[i] <- min(u.o[1:index]) sig.value2[i] <- min(v.o[1:index])
        ez.cutoff[i] <- e.z.ordered[index]  # fixed on 02/20/10
        n.selected[i] <- sum(e.z <= ez.cutoff[i])
        # npeak.rep1[i] <- sum(all.data1['sig.value'] >= sig.value1[i]) npeak.rep2[i] <- sum(all.data2['sig.value'] >= sig.value2[i])
    }
    
    # output the cutoff of posterior probability, number of selected overlapped peaks
    map.uv <- cbind(ez.cutoff, n.selected)
    
    return(list(n = n.select, IDR = IDR, idr.level = idr.level, map.uv = map.uv))
} 

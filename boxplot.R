if(!exists('markData_4me3_4ac')) load('mycounts_data.save')
d = markData_4me3_4ac
keep = apply(d, 1, max) > 1
d= d[keep,]
fname = "boxplot.pdf"
pdf(fname, width = 12, height = 6)
boxplot(d, ylab = 'log2 Fold Enrichment', notch = T)




#http://stats.stackexchange.com/questions/81864/hypothesis-test-for-difference-in-medians-among-more-than-two-samples
#http://digitalcommons.wayne.edu/cgi/viewcontent.cgi?article=1031&context=jmasm
#http://dornsife.usc.edu/assets/sites/239/docs/update_info.pdf
if(T){
  nq = 14
  high_qs = (nq/2):nq/nq #top 50
  low_qs = 0:(nq/2)/nq#lower 5
  highexp_qs = 1-(1/2^(1:nq)) #log slices at tippy top
  lowexp_qs = (1/2^(1:nq))
  qs = c(low_qs, high_qs)
  
  
  sliceQ = function(x){
    return(quantile(x, qs))
  }
  dq = apply(d, 2, sliceQ)
  sliceLines = function(x){
    lines(1:3, x[1:3])
    lines(4:6, x[4:6])
  }
  plot(c(0,7), c(0, max(d)), type = 'n')
  apply(dq, 1, sliceLines)
  
  plotQS = function(){
    dq = apply(d, 2, sliceQ)
    plot(c(0,7), c(0, max(d)), type = 'n')
    apply(dq, 1, sliceLines)
  }
  qs = c(low_qs, high_qs)
  plotQS()
  qs = c(lowexp_qs, highexp_qs)
  plotQS()
  
  #filter by q
  qval = .9
  sliceQ = function(x){
    return(quantile(x, qval))
  }
  dq = apply(d, 2, sliceQ)
  overQ = function(x){
    return(any(x > dq))
  }
  keep = apply(d, 1, overQ)
  boxplot(d[keep,], notch = T)
  
  #test filtered
  median.test <- function(x, y){
    z <- c(x, y)
    g <- rep(1:2, c(length(x), length(y)))
    m <- median(z)
    fisher.test(z < m, g)$p.value
  }
  a = c(1,1,2)
  b = c(2,3,3)
  for(i in 1:3){
    t = median.test(d[,a[i]], d[,b[i]])
    print(t)
  }
  a = a + 3
  b = b + 3
  for(i in 1:3){
    t = median.test(d[,a[i]], d[,b[i]])
    print(t)
  }
  
  a = c(1,1,2)
  b = c(2,3,3)
  for(i in 1:3){
    t = median.test(d[keep,a[i]], d[keep,b[i]])
    print(t)
  }
  a = a + 3
  b = b + 3
  for(i in 1:3){
    t = median.test(d[keep,a[i]], d[keep,b[i]])
    print(t)
  }
  
  qqplot(d[keep,1],d[keep,2])
  lines(c(0,10),c(0,10))
  
  
  library(snpar)
  # one-sample test
  x <- c(14.22, 15.83, 17.74, 19.88, 20.42, 21.96, 22.33, 22.79, 23.56, 24.45)
  ## normal approximation test
  quant.test(x, q = 14)
  ## exact quantile test 
  quant.test(x, q = 19, exact = TRUE)
  
  # two-sample test
  y <- c(5.54, 5.52, 5.00, 4.89, 4.95, 4.85, 4.80, 4.78, 4.82, 4.85, 4.72, 4.48, 
         4.39, 4.36, 4.30, 4.26, 4.25, 4.22)
  group <- as.numeric(gl(2,9))
  ## independent two-sample test
  quant.test(y, group, exact = TRUE)
  ## paired two-sample test
  quant.test(y,group, paired = TRUE)
  
  #library(WRS)
}
dev.off()
print(paste("wrote boxplots to", fname)) 
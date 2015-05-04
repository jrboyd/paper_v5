source('counts_to_FE.R')
load('dictionaries.save')
# keep = !(duplicated(index2ensg[rownames(data)]))
# data = data[keep,]
# keep = !is.na(index2ensg[rownames(data)])
# data = data[keep,]
b = agg_data
rownames(b) = index2ensg[rownames(b)]
b = ensg.expand(b)
colnames(b) = colnames(agg_data)




tmp = read.table(dir('counts_mine/',full.names = T)[1])
my_counts = tmp[,2, drop = F]
rownames(my_counts) = tmp[,1]
i = 0
for(f in dir('counts_mine/', full.names = T)){
  i = i + 1
  if(i == 1)
    next
  my_counts = cbind(my_counts, read.table(f)[,2])
  
}
n = dir('counts_mine/')
n = strsplit(n, '[_\\.]')
n = matrix(unlist(n), nrow = length(n[[1]]))
n = paste(n[1,], n[2,], n[3,], sep = ' ')
colnames(my_counts) = n

common = intersect(rownames(my_counts), rownames(b))

pdf('new vs old.pdf', height = 40)
layout(matrix(1:ncol(my_counts), ncol = 2, byrow = T))
for(c in 1:ncol(b)){
  n = colnames(b)[c]
  smoothScatter(log2(b[common,c]+1), log2(my_counts[common,paste(n,'R1')]+1), xlab = paste(n, 'old'), ylab = paste(n,'R1', 'new'), nrpoints = 0)
  smoothScatter(log2(b[common,c]+1), log2(my_counts[common,paste(n,'R2')]+1), xlab = paste(n, 'old'), ylab = paste(n,'R2 new'), nrpoints = 0)
}
dev.off()

b = agg_data
rownames(b) = index2ensg[rownames(b)]
keep = !grepl(',', rownames(b))
b = b[keep,]
common = intersect(rownames(my_counts), rownames(b))
pdf('new vs oldsimple.pdf', height = 40)
layout(matrix(1:ncol(my_counts), ncol = 2, byrow = T))
for(c in 1:ncol(b)){
  n = colnames(b)[c]
  smoothScatter(log2(b[common,c]+1), log2(my_counts[common,paste(n,'R1')]+1), xlab = paste(n, 'old'), ylab = paste(n,'R1', 'new'), nrpoints = 0)
  smoothScatter(log2(b[common,c]+1), log2(my_counts[common,paste(n,'R2')]+1), xlab = paste(n, 'old'), ylab = paste(n,'R2 new'), nrpoints = 0)
}
dev.off()

pdf('new vs oldnotsure.pdf', height = 40)
layout(matrix(1:ncol(my_counts), ncol = 2, byrow = T))
c = 1
n = colnames(b)[c]
x = log2(b[common,c]+1)
y = log2(my_counts[common,paste(n,'R2')]+1)
keep = y > x + 1
for(c in 1:ncol(b)){
  n = colnames(b)[c]
  smoothScatter(log2(b[common,c]+1)[keep], log2(my_counts[common,paste(n,'R1')]+1)[keep], xlab = paste(n, 'old'), ylab = paste(n,'R1', 'new'), nrpoints = 0, xlim = c(0,10), ylim = c(0,10))
  smoothScatter(log2(b[common,c]+1)[keep], log2(my_counts[common,paste(n,'R2')]+1)[keep], xlab = paste(n, 'old'), ylab = paste(n,'R2 new'), nrpoints = 0, xlim = c(0,10), ylim = c(0,10))
}
dev.off()

source('scripts/PCA.R')
n = strsplit(colnames(my_counts), ' ')
n = matrix(unlist(n), nrow = 3)
plotPCA(log2(my_counts+1), conditionLabels = n[2,], replicateLabels = n[3,], secondaryConditionLabels = n[1,], n = 3)

new_d = log2(my_counts[common,paste(n,'R1')]+1)[keep]
old_d = log2(b[common,c]+1)[keep]
keep = old_d < 1
old_zero_ensg = names(old_d)[keep]
old_zero_ensgAttribs = ensgAttribs_V21[old_zero_ensg,]

keep = old_d > 1
old_1plus = old_d[keep]
plot(old_d[keep], new_d[keep], xlim = c(0,10), ylim = c(0,10))
old_1plus_ensg = names(old_d)[keep]
old_1plus_ensgAttribs = ensgAttribs_V21[old_1plus_ensg,]

getTypeCounts = function(dat, column, min_count = -1){#get counts of Attrib table dat for column
  
  uniq = unique(ensgAttribs_V21[,column])
  output = rep(0, length(uniq))
  for(i in 1:length(uniq)){
    t = uniq[i]
    n_match = sum(dat[,column] == t)
    #print(paste(t,n_match, sep = '     '))
    output[i] = n_match
  }
  names(output) = uniq
  keep = output > min_count
  output = output[keep]
  writeClipboard(paste(names(output),output, sep = '     '))
  return(output)
}
bg = getTypeCounts(ensgAttribs_V21, 2)
zeros = getTypeCounts(old_zero_ensgAttribs, 2)
pluses = getTypeCounts(old_1plus_ensgAttribs, 2)
load('halfSecondary_ensg.save')
halfSec_ensgAttribs = ensgAttribs_V21[intersect(halfSecondary_ensg,rownames(ensgAttribs_V21)),]
halfsec = getTypeCounts(halfSec_ensgAttribs, 2)
grpCounts = rep(0,length(unique(type2groups)))
names(grpCounts) = unique(type2groups)
atrb = halfSec_ensgAttribs
toDo = list(bg, zeros+pluses, halfsec)
layout(matrix(1:3,nrow = 1))
for(tc in toDo){
for(n in names(grpCounts)){
  keep = type2groups[names(tc)] == n
  keep = keep[!is.na(keep)]
  grpCounts[n] = sum(tc[keep])
}
pie(grpCounts)
}
layout(1)
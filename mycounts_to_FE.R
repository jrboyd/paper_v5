indexDict = read.table("ref/hg38.index2symbol.1kb_ext_promoters.txt", row.names = 1, stringsAsFactors = F)
ref2sym = as.matrix(indexDict[, 1])
rownames(ref2sym) = rownames(indexDict)

ref2lengths = indexDict[, 2]
names(ref2lengths) = rownames(indexDict)

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

l_row = 1
m_row = 2
r_row = 3

indi_header = strsplit(colnames(my_counts), split = ' ')
indi_header = matrix(unlist(indi_header), nrow = 3)

lines = unique(indi_header[l_row, ])
mods = unique(indi_header[m_row, ])
reps = unique(indi_header[r_row, ])

data = my_counts

print("aggregating counts...")
agg_data = matrix(0, nrow = nrow(data), ncol = length(lines) * length(mods))
colnames(agg_data) = 1:ncol(agg_data)
rownames(agg_data) = rownames(data)
i = 1
for (l in lines) {
    for (m in mods) {
        keep = (indi_header[l_row, ] == l & indi_header[m_row, ] == m)
        agg_data[, i] = rowSums(data[, keep, drop = F])
        colnames(agg_data)[i] = paste(l, m)
        i = i + 1
    }
}



# calc cpm
print("normalizing per million reads...")
mod_data = agg_data
cpm_factors = 1e+06/colSums(mod_data)
mod_data = t(t(mod_data) * cpm_factors)

#norm for feature length
featureLengths = 2000
print("normalizing for kb feature length...")
mod_data = mod_data/featureLengths * 1000



print("applying log2 transform...")
rpkm_data = mod_data
fe_data = matrix(0, nrow = nrow(rpkm_data), ncol = length(lines) * (length(mods) - 1))
keep = !grepl("input", colnames(rpkm_data))
colnames(fe_data) = colnames(rpkm_data)[keep]
rownames(fe_data) = rownames(rpkm_data)

doFloor = F
if(doFloor)
  print("flooring data below 4 RPKM...")
print("calculating FE over input...")
for (i in 1:ncol(fe_data)) {
    key = colnames(fe_data)[i]
    l = strsplit(key, " ")[[1]][1]
    key_input = paste(l, "input")
    # print(paste(key, 'and', key_input))
    treat_data = rpkm_data[, key]
    if(doFloor){
      treat_data = ifelse(treat_data - 2 < 0, yes = 0, no = treat_data)
    }
    fe_data[, i] = treat_data / rpkm_data[, key_input]
}

print("flooring data at 1 FE, log2(1) = 0...")
fe_data = ifelse(fe_data < 0, 0, fe_data)
print("removing undetected entries...")
keep = apply(fe_data, 1, max) > 0
fe_data = fe_data[keep, ]
print("sorting columns by histone mark (4ac 4me3) then by line (mcf10a, mcf7, mda231)...")
markData_4me3_4ac = fe_data[, c(1, 3, 5, 2, 4, 6)]
#save(markData_4me3_4ac, file = "mycounts_data.save")
print("markData_4me3_4ac is set NOT saved!") 

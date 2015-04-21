library(ggplot2)
source('jrb_R_scripts/read_countDir.R')
source('jrb_R_scripts/norm_counts.R')
source('jrb_R_scripts/jrb_matrix.R')
source('jrb_R_scripts/heatmap.3-split.R')
load('ref//enst_dicts.save')
ensgs = enst_dict$ensg_id
chrms = enst_dict$chrms
keep = !duplicated(ensgs)

ensg2chrm = chrms[keep]
names(ensg2chrm) = ensgs[keep]

savePlot = T

exData = countDir.example()
exData = norm.read_depth(exData)
exData = norm.remove_unaligned(exData)

chrSizes = read.table('jrb_R_scripts/hg38.chrom.sizes', row.names = 1)
exData = norm.feature_length(exData, chrSizes)


tmp = jrb.split_colnames(exData)
cell_lines = tmp[1,]
histone_mods = tmp[2,]
colnames(exData) = paste(cell_lines, histone_mods)


exData = norm.fe(exData)
exData = log2(exData)
tmp = jrb.split_colnames(exData, split = ' ')
cell_lines = tmp[1,]
histone_mods = tmp[2,]

chrm_avg = numeric(length = 0)
for(hm in unique(plot_data$histone_mods)){
  for(cl in unique(plot_data$cell_lines)){
    keep = (plot_data$histone_mods == hm) & (plot_data$cell_lines == cl)
    chrm_avg = c(chrm_avg, log2(mean(2^plot_data[keep,1])))
    names(chrm_avg)[length(chrm_avg)] = paste(cl, hm)
  }
}
tmp = matrix(unlist(strsplit(names(chrm_avg), split = ' ')), nrow = 2)
cl = tmp[1,]
hm = tmp[2,]
chrm_averages = data.frame(values = chrm_avg, cell_lines = cl, histone_mods = hm, chrms = 'all', is_diff = 'mean', is_input = 'no')

matrix2data.frame = function(mat){
  chrms = factor(sub('chr', '', rep(colnames(mat), nrow(mat))), levels = sub('chr', '', colnames(mat)))
  
  i = unlist(lapply(1:nrow(mat), function(x)rep(x, ncol(mat))))#appends each row into single vector
  
  is_diff = rep('other', length(i))
  is_diff[chrms == '16'] = 'chr16'
  is_diff[chrms == '17'] = 'chr17'
  is_diff[chrms == '19'] = 'chr19'
  is_diff[chrms == '20'] = 'chr20'
  is_diff[chrms == '22'] = 'chr22'
  
  cl = jrb.split_colnames(t(mat), split = ' ')[1,]
  cl = cl[i]
  
  hm = jrb.split_colnames(t(mat), split = ' ')[2,]
  hm = hm[i]
  
  plot_data = data.frame(values = t(mat)[1:length(mat)], chrms = chrms, cell_lines = cl, is_input = ifelse(histone_mods[i]=='input', 'input','histone mark'),  histone_mods = hm, is_diff = is_diff)
  plot_data = rbind(plot_data, chrm_averages)
  return(plot_data)
}




tmp = t(exData)
o = mixedorder(rownames(exData))
tmp = tmp[,o]
tmp = tmp[,colnames(tmp) != 'chrY']
k4me3_dat = tmp[histone_mods == 'H3K4ME3',]
k4ac_dat = tmp[histone_mods == 'H3K4AC',]
tmp = tmp[histone_mods == 'H3K4AC' | histone_mods == 'H3K4ME3',]
plot_data = matrix2data.frame(tmp)



p <- ggplot(data = plot_data, mapping = aes(x = cell_lines, y = values, color = is_diff))
p = p + scale_colour_manual(values=c(RColorBrewer::brewer.pal(length(unique(plot_data$is_diff))-1, 'Set1'), '#595959'))

if(savePlot) pdf('k4 reads per chromosome.pdf', height = 3, width = 6, useDingbats = F)
print(p + geom_point(position = position_jitter(width = .04, height = 0)) + facet_grid(. ~ histone_mods, scales = "free"))
if(savePlot) dev.off()

plot_bychr_whighlight = function(plot_data){
  o = order(plot_data$is_diff, decreasing = T)
  
  plot_data = plot_data[o,]
  
  shapes = c(16,1,16,16,16,16,16)
  names(shapes) = unique(plot_data$is_diff)
  
  sizes = rep(2.5, length(unique(plot_data$is_diff)))
  names(sizes) = unique(plot_data$is_diff)
  sizes['other'] = 2
  
  colors = c(RColorBrewer::brewer.pal(length(unique(plot_data$is_diff))-1, 'Set1'), '#AAAAAA')
  colors = colors[c(1:(length(colors)-2), length(colors), length(colors)-1)]
  
  p <- ggplot(data = plot_data, mapping = aes(x = cell_lines, y = values, color = is_diff, shape = is_diff, size = is_diff)) + theme_bw()
  p = p + 
    scale_colour_manual(values = colors) + 
    scale_shape_manual(values = shapes) +
    scale_size_manual(values = sizes)
  p = p + layer(geom = "point", position = position_jitter(width = .12, height = 0))+ facet_grid(. ~ histone_mods, scales = "free")
#   stat_sum_single <- function(fun, geom="point", ...) {
#     stat_summary(fun.y=fun, fun.ymin = fun, fun.ymax = fun, colour="red", geom=geom, size = 3)
#   }
#   p = p + stat_sum_single(mean)
  print(p)
}

if(savePlot) pdf('k4 fe per chromosome.pdf', height = 3, width = 6, useDingbats = F)

p <- ggplot(data = plot_data, mapping = aes(x = chrms, y = values, color = cell_lines))+ theme_bw()
p = p + scale_colour_manual(values=c(RColorBrewer::brewer.pal(length(unique(plot_data$is_diff))-1, 'Set1'), '#595959'))
print(p + geom_point(position = position_jitter(width = .04, height = 0)) + facet_grid(. ~ histone_mods, scales = "free"))

plot_bychr_whighlight(plot_data)

matrix2data.frame = function(mat){
  chrms = factor(sub('chr', '', rep(colnames(mat), nrow(mat))), levels = sub('chr', '', colnames(mat)))
  
  i = unlist(lapply(1:nrow(mat), function(x)rep(x, ncol(mat))))#appends each row into single vector
  
  is_diff = rep('other', length(i))
  is_diff[chrms == '3'] = 'chr3'
  is_diff[chrms == '4'] = 'chr4'
  is_diff[chrms == '13'] = 'chr13'
  is_diff[chrms == '18'] = 'chr18'
  is_diff[chrms == 'X'] = 'chrX'

  cl = jrb.split_colnames(t(mat), split = ' ')[1,]
  cl = cl[i]
  
  hm = jrb.split_colnames(t(mat), split = ' ')[2,]
  hm = hm[i]
  
  plot_data = data.frame(values = t(mat)[1:length(mat)], chrms = chrms, cell_lines = cl, is_input = ifelse(histone_mods[i]=='input', 'input','histone mark'),  histone_mods = hm, is_diff = is_diff)
  plot_data = rbind(plot_data, chrm_averages)
  return(plot_data)
}

tmp = t(exData)
o = mixedorder(rownames(exData))
tmp = tmp[,o]
tmp = tmp[,colnames(tmp) != 'chrY']
k4me3_dat = tmp[histone_mods == 'H3K4ME3',]
k4ac_dat = tmp[histone_mods == 'H3K4AC',]
tmp = tmp[histone_mods == 'H3K4AC' | histone_mods == 'H3K4ME3',]
plot_data = matrix2data.frame(tmp)
plot_bychr_whighlight(plot_data)
if(savePlot) dev.off()

library(ggplot2)
source('jrb_R_scripts/read_countDir.R')
source('jrb_R_scripts/norm_counts.R')
source('jrb_R_scripts/jrb_matrix.R')
source('jrb_R_scripts/heatmap.3-split.R')
load('data_raw/my_fe_corrected.save')
load('ref//enst_dicts.save')
ensgs = enst_dict$ensg_id
chrms = enst_dict$chrms
keep = !duplicated(ensgs)

ensg2chrm = chrms[keep]
names(ensg2chrm) = ensgs[keep]

savePlot = T

chrSizes = read.table('jrb_R_scripts/hg38.chrom.sizes', row.names = 1)





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
  return(plot_data)
}

header = jrb.split_colnames(my_fe, split = ' ')
cell_lines = header[1,]
histone_mods = header[2,]

fe_per_chrm = matrix(0, nrow = ncol(my_fe), ncol = nrow(chrSizes))
rownames(fe_per_chrm) = colnames(my_fe)
colnames(fe_per_chrm) = rownames(chrSizes)

for(chr in rownames(chrSizes)){
  keep =  ensg2chrm[rownames(my_fe)] == chr
  chrm_fe = my_fe[keep,]
  dat = apply(chrm_fe, 2, quantile)[4,]
  fe_per_chrm[,chr] = dat
}


o = mixedorder(colnames(fe_per_chrm))
fe_per_chrm = fe_per_chrm[,o]
fe_per_chrm = fe_per_chrm[,colnames(fe_per_chrm) != 'chrY']
plot_data = matrix2data.frame(fe_per_chrm)

if(savePlot) pdf('k4 promoter fe per chromosome.pdf', height = 3, width = 6, useDingbats = F)

p <- ggplot(data = plot_data, mapping = aes(x = chrms, y = values, color = cell_lines))+ theme_bw()
p = p + scale_colour_manual(values=c(RColorBrewer::brewer.pal(length(unique(plot_data$is_diff))-1, 'Set1'), '#595959'))
print(p + geom_point(position = position_jitter(width = .04, height = 0)) + facet_grid(. ~ histone_mods, scales = "free"))


#if(savePlot) dev.off()

plot_bychr_whighlight = function(plot_data){
  o = order(plot_data$is_diff, decreasing = T)
  
  plot_data = plot_data[o,]
  
  shapes = c(1,16,6,3,15,8)
  names(shapes) = unique(plot_data$is_diff)
  
  sizes = rep(2.5, length(unique(plot_data$is_diff)))
  names(sizes) = unique(plot_data$is_diff)
  sizes['other'] = 2
  
  p <- ggplot(data = plot_data, mapping = aes(x = cell_lines, y = values, color = is_diff, shape = is_diff, size = is_diff)) + theme_bw()
  p = p + 
    scale_colour_manual(values = c(RColorBrewer::brewer.pal(length(unique(plot_data$is_diff))-1, 'Set1'), '#AAAAAA')) + 
    scale_shape_manual(values = shapes) +
    scale_size_manual(values = sizes)
  print(p + layer(geom = "point", position = position_jitter(width = .12, height = 0))+ facet_grid(. ~ histone_mods, scales = "free"))
}

#if(savePlot) pdf('k4 promoter fe per chromosome.pdf', height = 3, width = 6, useDingbats = F)

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
  return(plot_data)
}

plot_data = matrix2data.frame(fe_per_chrm)
plot_bychr_whighlight(plot_data)
if(savePlot) dev.off()

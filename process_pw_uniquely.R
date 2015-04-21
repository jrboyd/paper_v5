library(limma)
source('jrb_R_scripts//jrb_matrix.R')
load('data_raw/my_fe_corrected.save')
load('ref/dictionaries.save')
cell_lines = jrb.split_colnames(my_fe, split = ' ')[1,]
histone_marks = jrb.split_colnames(my_fe, split = ' ')[2,]

sym2dat = function(sym){
  ensg = names(ensg2sym)[ensg2sym == sym]
  print(my_fe[ensg,])
  print(mycounts[ensg,])
  return(my_fe[ensg,])
}

get_uniqSets = function(a = 1, b = 2, min_foldChange = 3, me3_factor = .8){
  minDiff = log2(min_foldChange)
  
  isMarked = matrix(F, nrow = nrow(my_fe), ncol = 0)
  isMarked = cbind(isMarked, my_fe[,a] > (my_fe[,b] + minDiff))
  isMarked = cbind(isMarked, my_fe[,b] > (my_fe[,a] + minDiff))
  
  minDiff = minDiff + log2(me3_factor)
  isMarked = cbind(isMarked, my_fe[,a+3] > (my_fe[,b+3] + minDiff))
  isMarked = cbind(isMarked, my_fe[,b+3] > (my_fe[,a+3] + minDiff))
  
  colnames(isMarked) = colnames(my_fe)[c(a,b,a+3,b+3)]
  
  layout(matrix(1:4, nrow = 2))
  vennDiagram(isMarked[,1:2], cex = .7)
  vennDiagram(isMarked[,3:4], cex = .7)
  vennDiagram(isMarked[,c(1,3)], cex = .7)
  vennDiagram(isMarked[,c(2,4)], cex = .7)
  
  uniquely_K4_membership = matrix(F, nrow = nrow(isMarked), ncol = 6)
  uniquely_K4_membership[,1] = isMarked[,1] & !isMarked[,3]
  uniquely_K4_membership[,2] = isMarked[,1] & isMarked[,3]
  uniquely_K4_membership[,3] = !isMarked[,1] & isMarked[,3]
  uniquely_K4_membership[,4] = isMarked[,2] & !isMarked[,4]
  uniquely_K4_membership[,5] = isMarked[,2] & isMarked[,4]
  uniquely_K4_membership[,6] = !isMarked[,2] & isMarked[,4]
  
  rownames(uniquely_K4_membership) = ensg2sym[rownames(isMarked)]
  cnames = colnames(my_fe)
  colnames(uniquely_K4_membership) = c(cnames[a], paste(cell_lines[a], 'both'), cnames[a+3], cnames[b], paste(cell_lines[b], 'both'), cnames[b+3])
  keep = rowSums(uniquely_K4_membership) > 0
  uniquely_K4_membership = uniquely_K4_membership[keep,]
  return(uniquely_K4_membership)
}

uniq_10a_v_7 = get_uniqSets(1, 2)
uniq_10a_v_231 = get_uniqSets(1, 3)
uniq_7_v_231 = get_uniqSets(2, 3)

save(uniq_10a_v_7, uniq_10a_v_231, uniq_7_v_231, file = 'data_intermediate/pw_uniquely.save')

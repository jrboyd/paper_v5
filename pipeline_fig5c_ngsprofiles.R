#outputs all ngs profiles for gsea lists in gsea_passing.save.
#separately plots groups in uniq_passing.save, highlighting whether mark trend hold for full gsea list.

##dependencies
source('functions_fig5c_ngsprofiles.R')

a = ngsplot_passing(res_10a_v_231, 'MCF10A', 'MDA231')
b = ngsplot_passing(res_10a_v_7, 'MCF10A', 'MCF7')
c = ngsplot_passing(res_7_v_231, 'MCF7', 'MDA231')
# 
# abc = union(union(a, b), c)
# 
# tmp = sapply(ensg2sym, function(x)return(any(x == abc)))
# tmp = names(ensg2sym)[tmp]
# toIPA = my_fe[tmp,]
# 
# gsea_list = gsea_lists$C2$GOZGIT_ESR1_TARGETS_DN
# tmp = sapply(ensg2sym, function(x)return(any(x == gsea_list)))
# tmp = names(ensg2sym)[tmp]
# tmp = my_fe[intersect(tmp, rownames(my_fe)),]
# write.table(tmp[,2] - tmp[,1], file = 'tmp2.txt', quote = F, row.names = T, col.names = T, sep = '\t')
# 
# tmp = sapply(ensg2sym, function(x)return(any(x == sym)))
# tmp = names(ensg2sym)[tmp]
# my_fe[tmp,]
# 
# uniqs = rownames(res_10a_v_7[[2]])[res_10a_v_7[[2]][,5]]
# inters = intersect(uniqs, gsea_list)
# tmp = sapply(ensg2sym, function(x)return(any(x == inters)))
# tmp = names(ensg2sym)[tmp]
# toIPA = my_fe[tmp,]
# write.table(rownames(toIPA), file = 'tmp5.txt', quote = F, row.names = F, col.names = F, sep = '\t')

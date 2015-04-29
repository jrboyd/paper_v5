load('data_raw/my_fe_corrected.save')
load('ref//dictionaries.save')
load('ref//gsea_dataset.save')
toOut = c(
  'CHARAFE_BREAST_CANCER_LUMINAL_VS_MESENCHYMAL_UP',
  'CHARAFE_BREAST_CANCER_LUMINAL_VS_MESENCHYMAL_DN',
  'CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_UP',
  'CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN',
  'CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_UP',
  'CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_DN',
  'GOZGIT_ESR1_TARGETS_DN',
  'GOZGIT_ESR1_TARGETS_UP')

for(o in toOut){
  l = rownames(gsea_master)[gsea_master[,o]]
  write(l, file = paste(o, '.txt', sep =''))
}


#my_fe = ifelse(my_fe < 0, 0, my_fe)
dat = as.data.frame(my_fe)

dat$gene_symbols = ensg2sym[rownames(dat)]
dat$ensg_ids = rownames(dat)
key = rep(T, ncol(dat))
ac_key = key; ac_key[4:6] = F
me_key = key; me_key[1:3] = F
ac_dat = dat[,ac_key]
me_dat = dat[,me_key]
write.table(ac_dat, file = 'ac_dat.txt',row.names = F, sep = '\t', quote = F)
write.table(me_dat, file = 'me_dat.txt',row.names = F, sep = '\t', quote = F)
fc_7_v_10a_dat = cbind(dat[,2] - dat[,1], dat[,7:8]); colnames(fc_ac_dat)[1] = 'MCF7 / MCF10A H3K4ac'
fc_231_v_10a_dat = cbind(dat[,3] - dat[,1], dat[,7:8]); colnames(fc_ac_dat)[1] = 'MDA231 / MCF10A H3K4ac'
write.table(fc_7_v_10a_dat, file = 'fc_7_v_10a_ac_dat.txt',row.names = F, sep = '\t', quote = F)
write.table(fc_231_v_10a_dat, file = 'fc_231_v_10a_ac_dat.txt',row.names = F, sep = '\t', quote = F)

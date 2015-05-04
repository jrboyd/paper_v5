if (exists("wd")) setwd(wd)
remove(list = ls())

filter4protein_coding = T
filter4idr_passed = F
forceRecount = F

source('load_dictionaries.R')
source("scripts/heatmap.3-split.R")

if(file.exists('mycounts_data.save') && !forceRecount){
  load('mycounts_data.save') #load previously processed data
}else{
  print('processing data...')
  source('mycounts_to_FE.R')  #process the data
}

keep = apply(markData_4me3_4ac, 1, max) > 1 #require one detected mark
markData_4me3_4ac = markData_4me3_4ac[keep, ]


if(filter4idr_passed){
  load('idrPassed_byIndex.save')
  idrPassed_byIndex = idrPassed_byIndex[rownames(markData_4me3_4ac),]
  markData_4me3_4ac = ifelse(idrPassed_byIndex, markData_4me3_4ac, 0)
}
#rownames(markData_4me3_4ac) = ensg.cut(rownames(markData_4me3_4ac)) #cutting no longer necessary
if(filter4protein_coding){
  keep = ensg2type[rownames(markData_4me3_4ac)] == 'protein_coding'
  tmp = rownames(markData_4me3_4ac)[is.na(keep)]
  if(length(tmp) > 0){
    warning(paste('problem with ensg2type,', length(tmp), ' ensg NOT recognized'))
    a = !is.na(keep)
    keep = keep[a]
  }
  markData_4me3_4ac = markData_4me3_4ac[keep,]
}
dirName = "paper_figs_2-26-15"
if (!file.exists(dirName)) dir.create(dirName)
wd = getwd()
setwd(dirName)
sourcewd = function(src) {
    source(paste(wd, src, sep = "/"))
}

doPlots = T
doText = T

#sourcewd('functions_batch-consistency-plot.R')
#sourcewd("batch-consistency-plot.r")
sourcewd('boxplot.R')
sourcewd("vennDiagrams.R")
sourcewd("scatterPlots.R")
sourcewd("scatterPlots_overUnder.R")
sourcewd("most_changed.R")



setwd(wd) 

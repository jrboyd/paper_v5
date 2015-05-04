source('')
do_oldway = F

ensgTable = read.table('gencode.v21.annotation_genes.gff3', stringsAsFactors = F)
head(ensgTable)
tmp = (ensgTable[,9])
tmp = strsplit(tmp, ';')
subsplit = function(x){
  x = strsplit(x, '=')
  x = matrix(unlist(x), nrow = 2)
  keys = x[1,]
  vals = x[2,]
  names(vals) = keys
  return(vals)
}
tmp = sapply(tmp, subsplit)
cnames = names(tmp[[1]])
getentry = function(x){
  gene_id = x[1]
  gene_type = x[4]
  gene_status = x[5]
  gene_name = x[6]
  return(c(gene_id, gene_type, gene_status, gene_name))
  
}
ensgAttribs = t(sapply(tmp, getentry))
rownames(ensgAttribs) = ensgAttribs[,1]
ensgAttribs_V21 = ensgAttribs
ensgAttribs_V21[,1] = ensg.cut(ensgAttribs_V21[,1])
save(ensgAttribs_V21, file = 'ensgAttribs_V21.save')

if(do_oldway){
  load('mycounts_data.save')
  ensgs = rownames(markData_4me3_4ac)
  ensgs = ensg.cut(ensgs)
  library(biocLite)
  source("http://bioconductor.org/biocLite.R")
  library(RCurl)
  library(biomaRt)
  # define biomart object
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  # query biomart
  results <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"), filters = "ensembl_gene_id", 
                   values = ensgs, mart = mart)
  results 
}
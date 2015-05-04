#set how ensg ids are convertd to gene names and other metadata - biotype etc.
#if dictionaries need to be updated use update_ensg.R

load('ensgAttribs_V21.save') #this is based directly on enseml V21 gene annotation
ensg2type = ensgAttribs_V21[,2]
names(ensg2type) = rownames(ensgAttribs_V21)

load('type2groups.save') #consolidates biotypes to broad categories

ensg2sym = ensgAttribs_V21[,4]
names(ensg2sym) = rownames(ensgAttribs_V21)

cutensg = rownames(ensgAttribs_V21)
cutensg = strsplit(cutensg, '\\.')
cutensg = matrix(unlist(cutensg), nrow = 2)[1,]
ensg2cut = cutensg
names(ensg2cut) = rownames(ensgAttribs_V21)
save(ensg2cut, ensg2sym, ensg2type, type2groups, file = 'dictionaries.save')
